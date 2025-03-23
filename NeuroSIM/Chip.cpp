#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "MaxPooling.h"
#include "Sigmoid.h"
#include "BitShifter.h"
#include "AdderTree.h"
#include "Buffer.h"
#include "HTree.h"
#include "ProcessingUnit.h"
#include "Tile.h"
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "Chip.h"
#include "XYBus.h"			  

using namespace std;

extern Param *param;
double globalBusWidth = 0;
int    numBufferCore  = 0;				  

/*** Circuit Modules ***/
Buffer *globalBuffer;
HTree *GhTree;
XYBus *GBus;							  
AdderTree *Gaccumulation;
Sigmoid *Gsigmoid;
BitShifter *GreLu;
MaxPooling *maxPool;

void ChipDesignInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell){
	globalBuffer  = new Buffer(inputParameter, tech, cell)    ;
	GhTree        = new HTree(inputParameter, tech, cell)     ;
	GBus          = new XYBus(inputParameter, tech, cell)     ;												
	Gaccumulation = new AdderTree(inputParameter, tech, cell) ;
	Gsigmoid      = new Sigmoid(inputParameter, tech, cell)   ;
	GreLu         = new BitShifter(inputParameter, tech, cell);
	maxPool       = new MaxPooling(inputParameter, tech, cell);
}

vector<double> Mapping(bool findsubArrayDup, bool findpeDup, bool findtileDup, bool findNumTile, bool find_no_tileDup_NumTile, bool findUtilization, bool findSpeedUp, 
								const vector<vector<double> > &netStructure, double *desiredNumTileCM, double *desiredTileSizeCM, double *desiredPESizeCM, int *numTileRow, 
								int *numTileCol, int raw_processing, int layer_by_layer, int full_pipeline, int partial_pipeline) {
	
	int numRowPerSynapse, numColPerSynapse;

	numRowPerSynapse = param->numRowPerSynapse;
	numColPerSynapse = param->numColPerSynapse;
	
	vector<double> tileDup    ;
	vector<double> peDup      ;
	vector<double> subArrayDup;

	vector<double> numTileEachLayer           ;
	vector<double> no_tileDup_numTileEachLayer;
	vector<double> utilizationEachLayer       ;
	vector<double> speedUpEachLayer           ;
	
	*desiredNumTileCM  = param->numTileRow*param->numTileCol                                    ; // get numTileTotal
	*desiredTileSizeCM = param->numPERowPerTile*param->numSubArrayRowPerPE*param->numRowSubArray;
	*desiredPESizeCM   = param->numSubArrayRowPerPE*param->numRowSubArray                       ;
	*numTileRow        = param->numTileRow                                                      ; // get the number row of Tile
	*numTileCol        = param->numTileCol                                                      ; // get the number column of Tile

	//max(duplication time of each layer): (netStructure[i][0]-netStructure[i][3]+1)/netStructure[i][7]*(netStructure[i][1]-netStructure[i][4]+1)/netStructure[i][7];
	/*** SubArray Duplication (consider duplication inside one PE) ***/
	subArrayDup = SubArrayDup((*desiredPESizeCM), netStructure, numRowPerSynapse, numColPerSynapse);
	//cout << 3 << endl;
	/*** PE Duplication (consider duplication inside one Tile) ***/ 
	peDup = PeDup(subArrayDup, (*desiredPESizeCM), (*desiredTileSizeCM), netStructure, numRowPerSynapse, numColPerSynapse);
	//cout << 4 << endl;
	/*** PE Duplication (consider duplication inside one Tile) ***/ 
	tileDup = TileDup(*desiredTileSizeCM, subArrayDup, peDup, netStructure, numRowPerSynapse, numColPerSynapse, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
	//cout << 5 << endl;
	
	numTileEachLayer            = OverallEachLayer(false, false, true , tileDup, peDup, subArrayDup, (*desiredTileSizeCM), netStructure, numRowPerSynapse, numColPerSynapse);
	speedUpEachLayer            = OverallEachLayer(false, true , false, tileDup, peDup, subArrayDup, (*desiredTileSizeCM), netStructure, numRowPerSynapse, numColPerSynapse);
	utilizationEachLayer        = OverallEachLayer(true , false, false, tileDup, peDup, subArrayDup, (*desiredTileSizeCM), netStructure, numRowPerSynapse, numColPerSynapse);
	no_tileDup_numTileEachLayer = OverallEachLayer(false, false, false, tileDup, peDup, subArrayDup, (*desiredTileSizeCM), netStructure, numRowPerSynapse, numColPerSynapse);

	if (findsubArrayDup) {
		return subArrayDup;
	}
	else if (findpeDup) {
		return peDup;
	}
	else if (findtileDup) {
		return tileDup;
	}
	else if (findNumTile) {
		return numTileEachLayer;
	} 
	else if (find_no_tileDup_NumTile) {
		return no_tileDup_numTileEachLayer;
	}
	else if (findUtilization) {
		return utilizationEachLayer;
	} 
	else if (findSpeedUp) {
		return speedUpEachLayer;
	}
}

void ChipInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell, const vector<vector<double> > &netStructure, const vector<double> &numTileEachLayer,
					double desiredNumTileCM, double desiredTileSizeCM, double desiredPESizeCM, int numTileRow, int numTileCol, int full_pipeline, const vector<double> &tileDup) { 
	
	int numRowPerSynapse;

	numRowPerSynapse = param->numRowPerSynapse;

	/*** Initialize Tile ***/
	//ceil((double)(desiredTileSizeCM)/(double)(desiredPESizeCM)): the number of PE in one Tile
	TileInitialize(inputParameter, tech, cell, ceil((double)(desiredTileSizeCM)/(double)(desiredPESizeCM)), desiredPESizeCM); 

	// find max layer and define the global buffer: enough to hold the max layer inputs
	double maxLayerInput = 0;
	// find max # tiles needed to be added at the same time
	double maxTileAdded = 0;

	for (int i=0; i<netStructure.size(); i++) {
		//netStructure[i][0]: represent IFM Length, netStructure[i][1]: represent IFM Width, netStructure[i][2]: represent IFM Channel Depth
		double input = netStructure[i][0]*netStructure[i][1]*netStructure[i][2]; 
		if (!full_pipeline) {
			if (input > maxLayerInput) {
				maxLayerInput = input;
			}
			globalBusWidth += (desiredTileSizeCM)+(desiredTileSizeCM)/param->numColMuxed; // numColMuxed: How many columns share 1 ADC (for eNVM and FeFET) or parallel SRAM
		}
		else {
			maxLayerInput += netStructure[i][0]*netStructure[i][1]*netStructure[i][2]/2;
			globalBusWidth += ((desiredTileSizeCM)+(desiredTileSizeCM)/param->numColMuxed)*numTileEachLayer[i];
		}

		if (tileDup[i]*ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/desiredTileSizeCM) > maxTileAdded) {
			maxTileAdded = tileDup[i]*ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/desiredTileSizeCM);
		}
	}
	// have to limit the global bus width --> cannot grow dramatically with num of tile
	while (globalBusWidth > param->maxGlobalBusWidth) {
		globalBusWidth /= 2;
	}
	
	// define bufferSize for inference operation
	int bufferSize;

	// numBitInput: precision of input activation
	bufferSize = param->numBitInput*maxLayerInput;	
	// calculate that for the maxLayerInput, it needs the number of Buffer core
	numBufferCore = ceil((double) bufferSize/((double) param->globalBufferCoreSizeRow*(double) param->globalBufferCoreSizeCol));
	
	globalBuffer->Initialize((param->globalBufferCoreSizeRow*param->globalBufferCoreSizeCol), param->globalBufferCoreSizeCol, 1, 
							param->unitLengthWireResistance, param->clkFreq, param->globalBufferType);

	maxPool->Initialize(param->numBitInput, 2*2, (desiredTileSizeCM/(double) param->numColMuxed), param->clkFreq);
	
	// globalBusType->true: represent using H-tree
	if (param->globalBusType) {
		GhTree->Initialize((numTileRow), (numTileCol), param->globalBusDelayTolerance, globalBusWidth, param->clkFreq);
	} 
	else {
		GBus->Initialize((numTileRow), (numTileCol), param->globalBusDelayTolerance, globalBusWidth, param->clkFreq);
	}
	
	int numBitSubarrayOutput = 0;

	if (param->parallelRead) {		
		numBitSubarrayOutput = log2((double)param->levelOutput)+ceil(log2(ceil(param->numRowSubArray/param->numRowParallel)))+
									param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	} 
	else{
		numBitSubarrayOutput = ceil(log2((double)param->numRowSubArray))+param->cellBit+param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	}

	int maxAddFromSubArray;

	maxAddFromSubArray = (int) ceil((double)(desiredPESizeCM)/(double)param->numRowSubArray); // from subArray to ProcessingUnit
	maxAddFromSubArray *= (int) ceil((double)(desiredTileSizeCM)/(double)(desiredPESizeCM)); // from ProcessingUnit to Tile
	
	if (full_pipeline) {
		maxAddFromSubArray *= (netStructure.size()+1);
	}
	
	Gaccumulation->Initialize((int) maxTileAdded, numBitSubarrayOutput+ceil(log2((double)maxAddFromSubArray)), 
							ceil((double)(desiredTileSizeCM)/(double)param->numColMuxed), param->clkFreq);

	if (param->reLu) {
		GreLu->Initialize(ceil((double)(desiredTileSizeCM)/(double)param->numColMuxed), param->numBitInput, param->clkFreq);
	} 
	else {
		Gsigmoid->Initialize(param->Activationtype, param->numBitInput, numBitSubarrayOutput+ceil(log2((double) maxAddFromSubArray))+ceil(log2((double) maxTileAdded)), 
							ceil((double) (desiredTileSizeCM)/(double) param->numColMuxed), param->clkFreq);
	}
}

vector<double> ChipCalculateArea(InputParameter& inputParameter, Technology& tech, MemCell& cell, double desiredNumTileCM, double desiredTileSizeCM, 
						double desiredPESizeCM, int numTileRow, double *CMTileheight, double *CMTilewidth) {
	
	vector<double> areaResults;
	
	double area      = 0;
	double areaIC    = 0;
	double areaADC   = 0;
	double areaAccum = 0;
	double areaOther = 0;
	double areaArray = 0;
	
	double CMheight = 0;
	double CMwidth  = 0;
	
	*CMTileheight = 0;
	*CMTilewidth  = 0;
	
	vector<double> areaCMTile;
	// (double) desiredTileSizeCM/(double) desiredPESizeCM: represent the number of PE inside one Tile
	areaCMTile = TileCalculateArea(pow(ceil((double) desiredTileSizeCM/(double) desiredPESizeCM), 2), desiredPESizeCM, &CMheight, &CMwidth);
	
	double CMTileArea      = areaCMTile[0];
	double CMTileAreaIC    = areaCMTile[1]; // hTreeCM->area (interconnect)
	double CMTileAreaADC   = areaCMTile[2];
	double CMTileAreaAccum = areaCMTile[3];
	double CMTileAreaOther = areaCMTile[4];
	double CMTileAreaArray = areaCMTile[5];

	area      += CMTileArea*desiredNumTileCM     ;
	areaIC    += CMTileAreaIC*desiredNumTileCM   ;
	areaADC   += CMTileAreaADC*desiredNumTileCM  ;
	areaAccum += CMTileAreaAccum*desiredNumTileCM;
	areaOther += CMTileAreaOther*desiredNumTileCM;
	areaArray += CMTileAreaArray*desiredNumTileCM;
	
	*CMTileheight = CMheight; //CMheight = sqrt(area) CMwidth = area/(*height) area-> the area of one Tile
	*CMTilewidth  = CMwidth ;
	
	// global buffer is made up by multiple cores
	globalBuffer->CalculateArea(numTileRow*CMheight, NULL, NONE);

	double globalBufferArea = globalBuffer->area*numBufferCore;
	double globalBufferHeight = numTileRow*CMheight;
	double globalBufferWidth = globalBufferArea/globalBufferHeight;														

	// globalBusType->true: represent using H-tree
	if (param->globalBusType) {
		GhTree->CalculateArea(CMheight, CMwidth, param->treeFoldedRatio);

		area   += GhTree->area;
		areaIC += GhTree->area;
	} 
	else {													
		GBus->CalculateArea(CMheight, CMwidth, param->treeFoldedRatio);

		area   += GBus->area;
		areaIC += GBus->area;
	}

	maxPool->CalculateUnitArea(NONE);
	maxPool->CalculateArea(globalBufferWidth);
	Gaccumulation->CalculateArea(NULL, globalBufferHeight/3, NONE);
	
	double areaGreLu = 0;
	double areaGsigmoid = 0;
	
	if (param->reLu) {
		GreLu->CalculateArea(NULL, globalBufferWidth/3, NONE);
		area      += GreLu->area;
		areaGreLu += GreLu->area;
	} 
	else {
		Gsigmoid->CalculateUnitArea(NONE);
		Gsigmoid->CalculateArea(NULL, globalBufferWidth/3, NONE);
		area         += Gsigmoid->area;
		areaGsigmoid += Gsigmoid->area;
	}
	
	area += globalBufferArea + maxPool->area + Gaccumulation->area;
	
	areaResults.push_back(area);
	areaResults.push_back(areaIC);
	areaResults.push_back(areaADC);
	areaResults.push_back(areaAccum + Gaccumulation->area);
	areaResults.push_back(areaOther + globalBufferArea + maxPool->area + areaGreLu + areaGsigmoid);
	areaResults.push_back(areaArray);
	
	return areaResults;
}

double ChipCalculatePerformance(InputParameter& inputParameter, Technology& tech, MemCell& cell, int layerNumber, const string &newweightfile, const string &oldweightfile, 
							const string &inputfile, bool followedByMaxPool, const vector<vector<double> > &netStructure, const vector<double> &numTileEachLayer, 
							const vector<double> &speedUpEachLayer, const vector<vector<double> > &tileLocaEachLayer, double desiredTileSizeCM, double desiredPESizeCM, double CMTileheight, 
							double CMTilewidth, double *readLatency, double *readDynamicEnergy, double *writeLatency, double *writeDynamicEnergy, double *leakage, double *leakageSRAMInUse, 
							double *bufferLatency, double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy, double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, 
							double *coreEnergyADC, double *coreEnergyAccum, double *coreEnergyOther, bool CalculateclkFreq, double *clkPeriod, int *partition_col, int *partition_row, 
							const vector<double> &tileDup, const vector<double> &peDup, const vector<double> &subArrayDup, double *act_poolLatency, double *computationLatency) {
	
	int numRowPerSynapse, numColPerSynapse;

	numRowPerSynapse = param->numRowPerSynapse;
	numColPerSynapse = param->numColPerSynapse;

	// only get performance of single layer
	int l = layerNumber;

	// get weight matrix file Size
	int weightMatrixRow = netStructure[l][2]*netStructure[l][3]*netStructure[l][4]*numRowPerSynapse;
	int weightMatrixCol = netStructure[l][5]*numColPerSynapse;

	// load in whole file 
	vector<vector<double> > inputVector;
	inputVector = LoadInInputData(inputfile); 

	vector<vector<double> > newMemory;
	newMemory = LoadInWeightData(newweightfile, numRowPerSynapse, numColPerSynapse, param->maxConductance, param->minConductance);

	vector<vector<double> > oldMemory;
	oldMemory = LoadInWeightData(oldweightfile, numRowPerSynapse, numColPerSynapse, param->maxConductance, param->minConductance);
	
	*readLatency         = 0;
	*readDynamicEnergy   = 0;
	*writeLatency        = 0;
	*writeDynamicEnergy  = 0;
	*leakage             = 0;
	*leakageSRAMInUse    = 0;
	*bufferLatency       = 0;
	*bufferDynamicEnergy = 0;
	*icLatency           = 0;
	*icDynamicEnergy     = 0;
	*coreEnergyADC       = 0;
	*coreEnergyAccum     = 0;
	*coreEnergyOther     = 0;
	*coreLatencyADC      = 0;
	*coreLatencyAccum    = 0;
	*coreLatencyOther    = 0;

	*partition_row = 0;
	*partition_col = 0;

	*act_poolLatency    = 0;
	*computationLatency = 0;
	
	double tileLeakage          = 0;
	double tileLeakageSRAMInUse = 0;

	double tileReadLatency         = 0;
	double tileReadDynamicEnergy   = 0;
	double tileWriteLatency        = 0;
	double tileWriteDynamicEnergy  = 0;
	double tilebufferLatency       = 0;
	double tilebufferDynamicEnergy = 0;
	double tileicLatency           = 0;
	double tileicDynamicEnergy     = 0;
	double tileLatencyADC          = 0;
	double tileLatencyAccum        = 0;
	double tileLatencyOther        = 0;
	double tileEnergyADC           = 0;
	double tileEnergyAccum         = 0;
	double tileEnergyOther         = 0;

	double tile_computationLatency       = 0;
	// double computation_latency_per_layer = 0;

	if(!CalculateclkFreq) {	
		globalBuffer->clkFreq  = param->clkFreq; 
		GhTree->clkFreq        = param->clkFreq; 						  
		Gaccumulation->clkFreq = param->clkFreq; 
		Gsigmoid->clkFreq      = param->clkFreq; 
		GreLu->clkFreq         = param->clkFreq; 
		maxPool->clkFreq       = param->clkFreq; 
	}
	// netStructure[l][7]: represent the kernel stride of each layer 
	// int numInVector = (netStructure[l][0]-netStructure[l][3]+1)/netStructure[l][7]*(netStructure[l][1]-netStructure[l][4]+1)/netStructure[l][7];
	int numInVector  = inputVector[0].size()/param->numBitInput;
	int totalNumTile = param->numTileRow*param->numTileCol;
	
	double totalnumTileRow = param->numTileRow;
	double totalnumTileCol = param->numTileCol;																			 

	if (speedUpEachLayer[l] > 1) {
		if ((tileDup[l] == totalnumTileRow*totalnumTileCol)) {
			// duplication in Tile or PE or subArray --> tell each Tile to take the whole assigned weight  --> "fully" duplication
			// assign weight and input to specific tile
			vector<vector<double> > tileMemory;
			tileMemory = CopyArray(newMemory, 0, 0, weightMatrixRow, weightMatrixCol);

			vector<vector<double> > tileMemoryOld;
			tileMemoryOld = CopyArray(oldMemory, 0, 0, weightMatrixRow, weightMatrixCol);
			
			vector<vector<double> > tileInput;
			tileInput = CopyInput(inputVector, 0, numInVector*param->numBitInput, weightMatrixRow);
	
			TileCalculatePerformance(tileMemory, tileMemoryOld, tileInput, ceil((double)desiredTileSizeCM/(double)desiredPESizeCM), desiredPESizeCM, ceil(speedUpEachLayer[l]/(totalnumTileRow*totalnumTileCol)), 
									peDup[l], weightMatrixRow, weightMatrixCol, numInVector*param->numBitInput, tech, 
									cell, &tileReadLatency, &tileReadDynamicEnergy, &tileWriteLatency, &tileWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &tilebufferLatency, &tilebufferDynamicEnergy, 
									&tileicLatency, &tileicDynamicEnergy, &tileLatencyADC, &tileLatencyAccum, &tileLatencyOther, &tileEnergyADC, &tileEnergyAccum, &tileEnergyOther, 
									CalculateclkFreq, clkPeriod, &tile_computationLatency);
			// further speed up in Tile Level
			*readLatency          = tileReadLatency/(totalnumTileRow*totalnumTileCol)  ;
			*readDynamicEnergy    = tileReadDynamicEnergy                              ;
			*bufferLatency        = tilebufferLatency/(totalnumTileRow*totalnumTileCol);
			*bufferDynamicEnergy  = tilebufferDynamicEnergy                            ;
			*icLatency            = tileicLatency/(totalnumTileRow*totalnumTileCol)    ;
			*icDynamicEnergy      = tileicDynamicEnergy                                ;
			
			*coreLatencyADC   = tileLatencyADC/(totalnumTileRow*totalnumTileCol)  ;
			*coreLatencyAccum = tileLatencyAccum/(totalnumTileRow*totalnumTileCol);
			*coreLatencyOther = tileLatencyOther/(totalnumTileRow*totalnumTileCol);
			
			*coreEnergyADC   = tileEnergyADC  ;
			*coreEnergyAccum = tileEnergyAccum;
			*coreEnergyOther = tileEnergyOther;

			*computationLatency = tile_computationLatency/(totalnumTileRow*totalnumTileCol);

			//cout << "*readLatency: " << *readLatency << endl;

			*writeLatency = tileWriteLatency                                              ;
			*writeDynamicEnergy = tileWriteDynamicEnergy*(totalnumTileRow*totalnumTileCol);
			//cout << "layer " << l+1 << " computation latency: " << computation_latency_per_layer << endl;
		}
		else {
			// # duplication is smaller then # Tile, means only a group of Tile take the assigned weight  --> not "fully" duplication
			// also need to redefine a few data-grab start-point
			//cout <<"layer(in 2):"<< l <<endl;
			for (int i=0; i<ceil((double) weightMatrixRow/(double) desiredTileSizeCM); i++) { // # of tiles in row
				for (int j=0; j<ceil((double) weightMatrixCol/(double) desiredTileSizeCM); j++) { // # of tiles in Column
					if ((i*desiredTileSizeCM < weightMatrixRow) && (j*desiredTileSizeCM < weightMatrixCol)) {
						int numRowMatrix = min(desiredTileSizeCM, weightMatrixRow-i*desiredTileSizeCM);
						int numColMatrix = min(desiredTileSizeCM, weightMatrixCol-j*desiredTileSizeCM);
						
						// assign weight and input to specific tile
						vector<vector<double> > tileMemory;
						tileMemory = CopyArray(newMemory, i*desiredTileSizeCM, j*desiredTileSizeCM, numRowMatrix, numColMatrix);

						vector<vector<double> > tileMemoryOld;
						tileMemoryOld = CopyArray(oldMemory, i*desiredTileSizeCM, j*desiredTileSizeCM, numRowMatrix, numColMatrix);
						
						vector<vector<double> > tileInput;
						tileInput = CopyInput(inputVector, i*desiredTileSizeCM, numInVector*param->numBitInput, numRowMatrix);

						TileCalculatePerformance(tileMemory, tileMemoryOld, tileInput, ceil((double)desiredTileSizeCM/(double)desiredPESizeCM), desiredPESizeCM,
												ceil(speedUpEachLayer[l]/tileDup[l]), peDup[l], numRowMatrix, numColMatrix, 
												numInVector*param->numBitInput, tech, cell, &tileReadLatency, &tileReadDynamicEnergy, &tileWriteLatency, &tileWriteDynamicEnergy, &tileLeakage, 
												&tileLeakageSRAMInUse, &tilebufferLatency, &tilebufferDynamicEnergy, &tileicLatency, &tileicDynamicEnergy, &tileLatencyADC, &tileLatencyAccum, 
												&tileLatencyOther, &tileEnergyADC, &tileEnergyAccum, &tileEnergyOther, CalculateclkFreq, clkPeriod, &tile_computationLatency);

						*readLatency          = MAX(tileReadLatency, (*readLatency))    ;
						*readDynamicEnergy   += tileReadDynamicEnergy                   ;
						*bufferLatency        = MAX(tilebufferLatency, (*bufferLatency));
						*bufferDynamicEnergy += tilebufferDynamicEnergy                 ;
						*icLatency            = MAX(tileicLatency, (*icLatency))        ;
						*icDynamicEnergy     += tileicDynamicEnergy                     ;
						
						*coreLatencyADC   = MAX(tileLatencyADC, (*coreLatencyADC))    ;
						*coreLatencyAccum = MAX(tileLatencyAccum, (*coreLatencyAccum));
						*coreLatencyOther = MAX(tileLatencyOther, (*coreLatencyOther));
						
						*coreEnergyADC   += tileEnergyADC  ;
						*coreEnergyAccum += tileEnergyAccum;
						*coreEnergyOther += tileEnergyOther;

						*computationLatency = MAX(tile_computationLatency, *computationLatency);

						*writeLatency       += tileWriteLatency      ;
						*writeDynamicEnergy += tileWriteDynamicEnergy;
					}
				}
			}
			// further speed up in Tile Level	
			*readLatency      /= (tileDup[l]);   
			*coreLatencyADC   /= (tileDup[l]);
			*coreLatencyAccum /= (tileDup[l]);
			*coreLatencyOther /= (tileDup[l]);
			*bufferLatency    /= (tileDup[l]);
			*icLatency        /= (tileDup[l]);

			*computationLatency = *computationLatency/(tileDup[l]);

			//*writeLatency       *= tileDup[l];
			*writeDynamicEnergy = (*writeDynamicEnergy)*tileDup[l];

			// whether go through accumulation?
			if (ceil((double) weightMatrixRow/(double) desiredTileSizeCM) > 1) {  
				Gaccumulation->CalculateLatency(ceil(weightMatrixCol/numColPerSynapse*(numInVector/(double) Gaccumulation->numAdderTree)), tileDup[l]*ceil((double) weightMatrixRow/(double) desiredTileSizeCM), 0);
				Gaccumulation->CalculatePower(ceil(weightMatrixCol/numColPerSynapse*(numInVector/(double) Gaccumulation->numAdderTree)), tileDup[l]*ceil((double) weightMatrixRow/(double) desiredTileSizeCM));
				
				*readDynamicEnergy += Gaccumulation->readDynamicEnergy;
				*readLatency += Gaccumulation->readLatency;
				*coreLatencyAccum += Gaccumulation->readLatency;
				*coreEnergyAccum += Gaccumulation->readDynamicEnergy;
			}
			//cout << "*readLatency: " << *readLatency << endl;
		}
	}
	else {
		if (numTileEachLayer[l] <= totalnumTileRow*totalnumTileCol) {
			// no duplication --> tell Tile to further partition the weight and grab data (redefine a few data-grab start-point)
			for (int i=0; i<totalnumTileRow; i++) { // # of tiles in row
				for (int j=0; j<totalnumTileCol; j++) { // # of tiles in Column
					if ((i*desiredTileSizeCM < weightMatrixRow) && (j*desiredTileSizeCM < weightMatrixCol)){
						
						int numRowMatrix = min(desiredTileSizeCM, weightMatrixRow-i*desiredTileSizeCM);
						int numColMatrix = min(desiredTileSizeCM, weightMatrixCol-j*desiredTileSizeCM);
						// assign weight and input to specific tile
						vector<vector<double> > tileMemory;
						tileMemory = CopyArray(newMemory, i*desiredTileSizeCM, j*desiredTileSizeCM, numRowMatrix, numColMatrix);

						vector<vector<double> > tileMemoryOld;
						tileMemoryOld = CopyArray(oldMemory, i*desiredTileSizeCM, j*desiredTileSizeCM, numRowMatrix, numColMatrix);
						
						vector<vector<double> > tileInput;
						tileInput = CopyInput(inputVector, i*desiredTileSizeCM, numInVector*param->numBitInput, numRowMatrix);

						TileCalculatePerformance(tileMemory, tileMemoryOld, tileInput, ceil((double)desiredTileSizeCM/(double)desiredPESizeCM), desiredPESizeCM, 1, 
												1, numRowMatrix, numColMatrix, numInVector*param->numBitInput, tech, cell, &tileReadLatency, &tileReadDynamicEnergy, &tileWriteLatency, 
												&tileWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &tilebufferLatency, &tilebufferDynamicEnergy, &tileicLatency, &tileicDynamicEnergy, 
												&tileLatencyADC, &tileLatencyAccum, &tileLatencyOther, &tileEnergyADC, &tileEnergyAccum, &tileEnergyOther, 
												CalculateclkFreq, clkPeriod, &tile_computationLatency);
						
						*readLatency          = MAX(tileReadLatency, (*readLatency))    ;
						*readDynamicEnergy   += tileReadDynamicEnergy                   ;
						*bufferLatency        = MAX(tilebufferLatency, (*bufferLatency));
						*bufferDynamicEnergy += tilebufferDynamicEnergy                 ;
						*icLatency            = MAX(tileicLatency, (*icLatency))        ;
						*icDynamicEnergy     += tileicDynamicEnergy                     ;
						
						*coreLatencyADC   = MAX(tileLatencyADC, (*coreLatencyADC))    ;
						*coreLatencyAccum = MAX(tileLatencyAccum, (*coreLatencyAccum));
						*coreLatencyOther = MAX(tileLatencyOther, (*coreLatencyOther));
						
						*coreEnergyADC   += tileEnergyADC  ;
						*coreEnergyAccum += tileEnergyAccum;
						*coreEnergyOther += tileEnergyOther;

						*computationLatency = MAX(tile_computationLatency, (*computationLatency));

						*writeLatency       += tileWriteLatency      ;
						*writeDynamicEnergy += tileWriteDynamicEnergy;
					}
				}
			}
			// whether go through accumulation?
			if (ceil((double) weightMatrixRow/(double) desiredTileSizeCM) > 1) {  
				Gaccumulation->CalculateLatency(ceil(weightMatrixCol/numColPerSynapse  *(numInVector/(double) Gaccumulation->numAdderTree)), ceil((double) weightMatrixRow/(double) desiredTileSizeCM), 0);
				Gaccumulation->CalculatePower(ceil(weightMatrixCol/numColPerSynapse  *(numInVector/(double) Gaccumulation->numAdderTree)), ceil((double) weightMatrixRow/(double) desiredTileSizeCM));
				
				*readDynamicEnergy += Gaccumulation->readDynamicEnergy;
				*readLatency += Gaccumulation->readLatency;
				*coreLatencyAccum += Gaccumulation->readLatency;
				*coreEnergyAccum += Gaccumulation->readDynamicEnergy;
			}
			//cout << "*readLatency: " << *readLatency << endl;
		}
		else {
			// this part needs to be revised later, because different parts latency need to be accumulated. (Done)
			// the defined accelerator cannot store all the weights of this layer. (weight partition)
			// cout << "OMG" <<endl;
			/*int weightMatrixCol_reduce_record = 0;
			int weightMatrixRow_reduce_record = 0;
			if (weightMatrixCol > desiredTileSizeCM*param->numTileCol*param->numTileRow) {
				while (weightMatrixCol > desiredTileSizeCM*param->numTileCol*param->numTileRow) {
					weightMatrixCol = weightMatrixCol/2;
					weightMatrixCol_reduce_record = weightMatrixCol_reduce_record+1;
				}
			}
			if (weightMatrixRow > desiredTileSizeCM*param->numTileRow*param->numTileCol) {
				while (weightMatrixRow > desiredTileSizeCM*param->numTileRow*param->numTileCol) {
					weightMatrixRow = weightMatrixRow/2;
					weightMatrixRow_reduce_record = weightMatrixRow_reduce_record+1;
				}
			}

			*partition_col = pow(2, weightMatrixCol_reduce_record);
			*partition_row = pow(2, weightMatrixRow_reduce_record);*/

			int flag = ceil(numTileEachLayer[l]/(totalnumTileRow*totalnumTileCol));
			
			flag = ceil(flag/2);

			int col = 0;
			int row = 0;

			for (int i=0; i<flag; i++) {
				if (weightMatrixCol > weightMatrixRow) {
					weightMatrixCol = weightMatrixCol/2;
					col = col+1;
				}
				else {
					weightMatrixRow = weightMatrixRow/2;
					row = row+1;
				}
			}

			*partition_col = pow(2, col);
			*partition_row = pow(2, row);

			int partition_col_copy = *partition_col;
			int partition_row_copy = *partition_row;
			int i, j, k, z;

			double t_readLatency      = 0;   
			double t_bufferLatency    = 0;   
			double t_icLatency        = 0;   
			double t_coreLatencyADC   = 0;   
			double t_coreLatencyAccum = 0;   
			double t_coreLatencyOther = 0;   

			for (i=0; i<partition_row_copy; i++) {
				for (j=0; j<partition_col_copy; j++) {
					for (k=0; k<totalnumTileRow; k++) {
						for (z=0; z<totalnumTileCol; z++) {
							if ((k*desiredTileSizeCM < weightMatrixRow) && (z*desiredTileSizeCM < weightMatrixCol)){
								int numRowMatrix = min(desiredTileSizeCM, weightMatrixRow-k*desiredTileSizeCM);
								int numColMatrix = min(desiredTileSizeCM, weightMatrixCol-z*desiredTileSizeCM);
								// assign weight and input to specific tile
								vector<vector<double> > tileMemory;
								tileMemory = CopyArray(newMemory, k*desiredTileSizeCM+i*weightMatrixRow, z*desiredTileSizeCM+j*weightMatrixCol, numRowMatrix, numColMatrix);

								vector<vector<double> > tileMemoryOld;
								tileMemoryOld = CopyArray(oldMemory, k*desiredTileSizeCM+i*weightMatrixRow, z*desiredTileSizeCM+j*weightMatrixCol, numRowMatrix, numColMatrix);
								
								vector<vector<double> > tileInput;
								tileInput = CopyInput(inputVector, k*desiredTileSizeCM+i*weightMatrixRow, numInVector*param->numBitInput, numRowMatrix);

								TileCalculatePerformance(tileMemory, tileMemoryOld, tileInput, ceil((double)desiredTileSizeCM/(double)desiredPESizeCM), desiredPESizeCM, 1, 
														1, numRowMatrix, numColMatrix, numInVector*param->numBitInput, tech, cell, &tileReadLatency, &tileReadDynamicEnergy, &tileWriteLatency, 
														&tileWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &tilebufferLatency, &tilebufferDynamicEnergy, &tileicLatency, &tileicDynamicEnergy, 
														&tileLatencyADC, &tileLatencyAccum, &tileLatencyOther, &tileEnergyADC, &tileEnergyAccum, &tileEnergyOther, 
														CalculateclkFreq, clkPeriod, &tile_computationLatency);
								
								t_readLatency         = MAX(tileReadLatency, (t_readLatency))    ;
								*readDynamicEnergy   += tileReadDynamicEnergy                    ;
								t_bufferLatency       = MAX(tilebufferLatency, (t_bufferLatency));
								*bufferDynamicEnergy += tilebufferDynamicEnergy                  ;
								t_icLatency           = MAX(tileicLatency, (t_icLatency))        ;
								*icDynamicEnergy     += tileicDynamicEnergy                      ;
								
								t_coreLatencyADC   = MAX(tileLatencyADC, (t_coreLatencyADC))    ;
								t_coreLatencyAccum = MAX(tileLatencyAccum, (t_coreLatencyAccum));
								t_coreLatencyOther = MAX(tileLatencyOther, (t_coreLatencyOther));
								
								*coreEnergyADC   += tileEnergyADC  ;
								*coreEnergyAccum += tileEnergyAccum;
								*coreEnergyOther += tileEnergyOther;

								*computationLatency = MAX(tile_computationLatency, (*computationLatency));

								*writeLatency       += tileWriteLatency      ;
								*writeDynamicEnergy += tileWriteDynamicEnergy;
							}
						}
					}
					if (k >= 1) {
						Gaccumulation->CalculateLatency(ceil(weightMatrixCol/numColPerSynapse*(numInVector/(double) Gaccumulation->numAdderTree)), k+1, 0);
						Gaccumulation->CalculatePower(ceil(weightMatrixCol/numColPerSynapse*(numInVector/(double) Gaccumulation->numAdderTree)), k+1);
						
						*readDynamicEnergy += Gaccumulation->readDynamicEnergy;
						*readLatency += Gaccumulation->readLatency;
						*coreLatencyAccum += Gaccumulation->readLatency;
						*coreEnergyAccum += Gaccumulation->readDynamicEnergy;
					}
					*readLatency      = *readLatency + t_readLatency          ;
					*bufferLatency    = *bufferLatency + t_bufferLatency      ;
					*icLatency        = *icLatency + t_icLatency              ;
					*coreLatencyADC   = *coreLatencyADC + t_coreLatencyADC    ;
					*coreLatencyAccum = *coreLatencyAccum + t_coreLatencyAccum;
					*coreLatencyOther = *coreLatencyOther + t_coreLatencyOther;
				}
				if (partition_col_copy > 1)  {
					Gaccumulation->CalculateLatency(ceil(weightMatrixCol/numColPerSynapse  *(numInVector/(double) Gaccumulation->numAdderTree)), partition_col_copy, 0);
					Gaccumulation->CalculatePower(ceil(weightMatrixCol/numColPerSynapse  *(numInVector/(double) Gaccumulation->numAdderTree)), partition_col_copy);
					
					*readDynamicEnergy += Gaccumulation->readDynamicEnergy;
					*readLatency += Gaccumulation->readLatency;
					*coreLatencyAccum += Gaccumulation->readLatency;
					*coreEnergyAccum += Gaccumulation->readDynamicEnergy;
				}
			}
		}
	}

	if(!CalculateclkFreq){
		if (param->reLu) {
			GreLu->CalculateLatency(ceil(numInVector*netStructure[l][5]/(double) GreLu->numUnit));
			GreLu->CalculatePower(ceil(numInVector*netStructure[l][5]/(double) GreLu->numUnit));
			
			*readLatency += GreLu->readLatency;
			*coreLatencyOther += GreLu->readLatency;
			*readDynamicEnergy += GreLu->readDynamicEnergy;
			*coreEnergyOther += GreLu->readDynamicEnergy;

			*act_poolLatency = *act_poolLatency + GreLu->readLatency;
		}
		else {
			Gsigmoid->CalculateLatency(ceil(numInVector*netStructure[l][5]/Gsigmoid->numEntry));
			Gsigmoid->CalculatePower(ceil(numInVector*netStructure[l][5]/Gsigmoid->numEntry));
			
			*readLatency += Gsigmoid->readLatency;
			*coreLatencyOther += Gsigmoid->readLatency;
			*readDynamicEnergy += Gsigmoid->readDynamicEnergy;
			*coreEnergyOther += Gsigmoid->readDynamicEnergy;

			*act_poolLatency = *act_poolLatency + Gsigmoid->readLatency;
		}
		
		// if this layer is followed by Max Pool
		if (followedByMaxPool) {
			maxPool->CalculateLatency(1e20, 0, ceil((double) (numInVector * weightMatrixCol/numColPerSynapse/(double) maxPool->window)/maxPool->numMaxPooling));
			maxPool->CalculatePower(ceil((double) (numInVector * weightMatrixCol/numColPerSynapse/(double) maxPool->window)/maxPool->numMaxPooling));

			*readLatency += maxPool->readLatency;
			*coreLatencyOther += maxPool->readLatency; 
			*readDynamicEnergy += maxPool->readDynamicEnergy;
			*coreEnergyOther += maxPool->readDynamicEnergy;

			*act_poolLatency = *act_poolLatency + maxPool->readLatency;
		}							  
		
		double numBitToLoadOut = weightMatrixRow*param->numBitInput*numInVector;
		double numBitToLoadIn = ceil(weightMatrixCol/param->numColPerSynapse)*param->numBitInput*numInVector/(netStructure[l][6]? 4:1);
		
		if (param->globalBusType) {
			double fraction = netStructure[l][2]*netStructure[l][3]*netStructure[l][4] * (netStructure[l][0]-netStructure[l][3]+1)/netStructure[l][7]*(netStructure[l][1]-netStructure[l][4]+1)/netStructure[l][7];
		
			GhTree->CalculateLatency(0, 0, tileLocaEachLayer[0][l], tileLocaEachLayer[1][l], CMTileheight, CMTilewidth, ceil((numBitToLoadOut+numBitToLoadIn)/ceil(GhTree->busWidth*(numTileEachLayer[l]/totalNumTile))));
			GhTree->CalculatePower(0, 0, tileLocaEachLayer[0][l], tileLocaEachLayer[1][l], CMTileheight, CMTilewidth, ceil(GhTree->busWidth*(numTileEachLayer[l]/totalNumTile)), 
						ceil((numBitToLoadOut*param->inputtoggle)/ceil(GhTree->busWidth*(numTileEachLayer[l]/totalNumTile) )) + ceil((numBitToLoadIn*param->outputtoggle)/ceil(GhTree->busWidth*(numTileEachLayer[l]/totalNumTile) ) ) );
		} 

		double chip_bufferclk = max(GhTree->critical_latency, 1.0); 
		
		globalBuffer->CalculateLatency(globalBuffer->interface_width, numBitToLoadOut/globalBuffer->interface_width,
								globalBuffer->interface_width, numBitToLoadIn/globalBuffer->interface_width);
		
		globalBuffer->CalculatePower(globalBuffer->interface_width, numBitToLoadOut/globalBuffer->interface_width,
								globalBuffer->interface_width, numBitToLoadIn/globalBuffer->interface_width);
		
		// since multi-core buffer has improve the parallelism
		globalBuffer->readLatency /= MIN(numBufferCore, ceil(globalBusWidth/globalBuffer->interface_width));
		globalBuffer->writeLatency /= MIN(numBufferCore, ceil(globalBusWidth/globalBuffer->interface_width));			
	} 
	if(!CalculateclkFreq){
		*bufferLatency += globalBuffer->readLatency + globalBuffer->writeLatency;
		*icLatency += GhTree->readLatency;
		*readLatency += globalBuffer->readLatency + globalBuffer->writeLatency + GhTree->readLatency;
		*coreLatencyOther += globalBuffer->readLatency + globalBuffer->writeLatency + GhTree->readLatency;
		*bufferDynamicEnergy += globalBuffer->readDynamicEnergy + globalBuffer->writeDynamicEnergy;
		*icDynamicEnergy += GhTree->readDynamicEnergy;
		*readDynamicEnergy += globalBuffer->readDynamicEnergy + globalBuffer->writeDynamicEnergy + GhTree->readDynamicEnergy;
		*coreEnergyOther += globalBuffer->readDynamicEnergy + globalBuffer->writeDynamicEnergy + GhTree->readDynamicEnergy;
		*leakage = tileLeakage;
		*leakageSRAMInUse = tileLeakageSRAMInUse;
	}
	return 0;
}

vector<double> TileDup(double tileSize, const vector<double> &subArrayDup, const vector<double> &peDup, const vector<vector<double> > &netStructure, 
									int numRowPerSynapse, int numColPerSynapse, int raw_processing, int layer_by_layer, int full_pipeline, int partial_pipeline) {
	
	vector<double> tileDup;

	if (raw_processing) {
		for (int i=0; i<netStructure.size(); i++) {
			tileDup.push_back(1);
		}
	}
	else if (layer_by_layer) {
		for (int i=0; i<netStructure.size(); i++) {
			int actualDup = 0;

			int numTileRow = param->numTileRow; //number row of Tile in the accelerator
			int numTileCol = param->numTileCol; //number column of Tile in the accelerator

			int maxDup = netStructure[i][0]*netStructure[i][1];

			if (peDup[i] > 1) {
				actualDup = numTileRow*numTileCol;
			}
			// netStructure[i][2]: represent IFM Channel Depth, netStructure[i][3]: represent Kernel Length, netStructure[i][4]: represent Kernel Width, netStructure[i][5]: represent Kernel Depth
			else if ((netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*numRowPerSynapse <= tileSize*param->numTileRow) && (netStructure[i][5]*numColPerSynapse <= tileSize*param->numTileCol)) {
				int tileForOneMatrixRow = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/(double) tileSize); // this layer needs the number row of PE
				int tileForOneMatrixCol = ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) tileSize); // this layer needs the number column of Tile

				actualDup = floor((numTileRow*numTileCol)/(tileForOneMatrixRow*tileForOneMatrixCol))==0? 1:floor((numTileRow*numTileCol)/(tileForOneMatrixRow*tileForOneMatrixCol));
			} 
			else {
				actualDup = 1;
			}

			int flag = actualDup;

			if ((netStructure[i][0] == 1) && (netStructure[i][1] == 1)) {
				actualDup = 1;
				tileDup.push_back(actualDup);
			}
			else if ((1 <= actualDup*subArrayDup[i]*peDup[i]) && (actualDup*subArrayDup[i]*peDup[i] <= maxDup)) {
				tileDup.push_back(actualDup);
			}
			else {
				for(int j=1; j<=(flag-1); j++) {
					actualDup = actualDup-1;
					if ((1 <= actualDup*subArrayDup[i]*peDup[i]) && (actualDup*subArrayDup[i]*peDup[i] <= maxDup)) {
						tileDup.push_back(actualDup);
						break;
					}
				}
			}
		}
	}
	else if (full_pipeline) {
		int before_dup_used_tile_num = 0;
		int min_tile_num_required_layer = 10000000000;

		vector<int> tile_num_each_layer;

		int numTileRow = param->numTileRow; //number row of Tile in the accelerator
		int numTileCol = param->numTileCol; //number column of Tile in the accelerator

		// tileDup initial
		for (int i=0; i<netStructure.size(); i++) {
			tileDup.push_back(1);
		}

		for (int i=0; i<netStructure.size(); i++) {
			int tileForOneMatrixRow = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/(double) tileSize);
			int tileForOneMatrixCol = ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) tileSize);

			if (min_tile_num_required_layer > (tileForOneMatrixRow*tileForOneMatrixCol)) {
				min_tile_num_required_layer = tileForOneMatrixRow*tileForOneMatrixCol;
			}

			before_dup_used_tile_num = before_dup_used_tile_num + tileForOneMatrixRow*tileForOneMatrixCol;

			tile_num_each_layer.push_back(tileForOneMatrixRow*tileForOneMatrixCol);
		}

		int extra_tile_num = numTileCol*numTileRow - before_dup_used_tile_num;

		vector<int> time_predictor;

		vector<int> record_cannot_accelerate_layer;

		//initial inference time predictor
		for (int i=0; i<netStructure.size(); i++) {
			int maxDup = netStructure[i][0]*netStructure[i][1];
			
			time_predictor.push_back(maxDup/(subArrayDup[i]*peDup[i]));
		}

		int flag = 1;

		while (extra_tile_num >= min_tile_num_required_layer) {
			// update predictor
			for (int i=0; i<netStructure.size(); i++) {
				int maxDup = netStructure[i][0]*netStructure[i][1];
				
				time_predictor[i] = maxDup/(subArrayDup[i]*peDup[i]*tileDup[i]);
			}

			int bottleneck_layer = 0;

			for (int i=0; i<netStructure.size(); i++) {
				if (flag == 1) {
					if (time_predictor[bottleneck_layer] < time_predictor[i]) {
						bottleneck_layer = i;
					}
				}
				else if (flag == 0) {
					for (int j=0; j<record_cannot_accelerate_layer.size(); j++) {
						if ((time_predictor[bottleneck_layer] < time_predictor[i]) && (i != record_cannot_accelerate_layer[j])) {
							bottleneck_layer = i;
						}
					}
				}
			}

			int maxDup = netStructure[bottleneck_layer][0]*netStructure[bottleneck_layer][1];

			if (tile_num_each_layer[bottleneck_layer] <= extra_tile_num) {
				flag = 1;
				extra_tile_num = extra_tile_num - tile_num_each_layer[bottleneck_layer];

				tileDup[bottleneck_layer] = tileDup[bottleneck_layer] +1;

				int temp = tileDup[bottleneck_layer];

				if ((netStructure[bottleneck_layer][0] == 1) && (netStructure[bottleneck_layer][1] == 1)) {
					tileDup[bottleneck_layer] = 1;
					break;
				}
				else if (tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer] > maxDup) {
					for (int j=0; j<=(temp-1); j++) {
						tileDup[bottleneck_layer] = tileDup[bottleneck_layer]-1;
						if ((1<=tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer]) && (tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer]<= maxDup)) {
							break;
						}
					}
				}
			}
			else {
				flag = 0;
				record_cannot_accelerate_layer.push_back(bottleneck_layer);
			}
		}
		if (!param->pipeline_dynamic_dup) {
			for (int i=0; i<netStructure.size(); i++) {
				tileDup[i]=1;
			}
		}
	}
	else if (partial_pipeline) {
		vector<int> tile_num_each_layer;

		for (int i=0; i<netStructure.size(); i++) {
			int tileForOneMatrixRow = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/(double) tileSize);
			int tileForOneMatrixCol = ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) tileSize);

			tile_num_each_layer.push_back(tileForOneMatrixRow*tileForOneMatrixCol);
		}

		vector<vector<int>> record;
		vector<int> recordRow;

		for (int i=0; i<netStructure.size(); i++) {
			if (tile_num_each_layer[i] <= param->numTileRow*param->numTileCol) {
				double p = tile_num_each_layer[i];
                recordRow.push_back(i);
                
                for(int j=i+1; j<netStructure.size(); j++) {
                    p = p + tile_num_each_layer[j];
                    if (p <= (param->numTileRow*param->numTileCol)) {
                        recordRow.push_back(j);
                    }
                    else {
                        break;
                    }
                }
                record.push_back(recordRow);
                i = recordRow.size()-1;
			}
			else {
                recordRow.push_back(i);
                record.push_back(recordRow);
                i = recordRow.size()-1;
            }
		}
		
		int temp1 = 0;
		int temp2 = 0;
		int temp3 = 0;
		int temp4 = 0;

		//cout << 6 << endl;

		// tileDup initial
		for (int i=0; i<netStructure.size(); i++) {
			tileDup.push_back(1);
		}

		vector<int> time_predictor;

		//initial inference time predictor
		for (int i=0; i<netStructure.size(); i++) {
			int maxDup = netStructure[i][0]*netStructure[i][1];
			
			time_predictor.push_back(maxDup/(subArrayDup[i]*peDup[i]));
		}

		//cout << 7 << endl;

		//cout << "record.size(): " << record.size() << endl;

		for (int i=0; i<record.size(); i++) {
			int before_dup_used_tile_num = 0;
			int min_tile_num_required_layer = 10000000000;

			int numTileRow = param->numTileRow; //number row of Tile in the accelerator
			int numTileCol = param->numTileCol; //number column of Tile in the accelerator

			for (int j=temp1; j<record[i].size(); j++) {
				int tileForOneMatrixRow = ceil((double) netStructure[j][2]*(double) netStructure[j][3]*(double) netStructure[j][4]*(double) numRowPerSynapse/(double) tileSize);
				int tileForOneMatrixCol = ceil((double) netStructure[j][5]*(double) numColPerSynapse/(double) tileSize);

				if (min_tile_num_required_layer > (tileForOneMatrixRow*tileForOneMatrixCol)) {
					min_tile_num_required_layer = tileForOneMatrixRow*tileForOneMatrixCol;
				}

				before_dup_used_tile_num = before_dup_used_tile_num + tileForOneMatrixRow*tileForOneMatrixCol;

				temp1 = j + 1;
			}

			//cout << 9 << endl;

			int extra_tile_num = numTileCol*numTileRow - before_dup_used_tile_num;

			vector<int> record_cannot_accelerate_layer;

			int flag = 1;
			int tag = 1;
			//cout << 10 << endl;
			while ((extra_tile_num >= min_tile_num_required_layer) && (tag == 1)) {
				// update predictor
				for (int j=temp3; j<record[i].size(); j++) {
					int maxDup = netStructure[j][0]*netStructure[j][1];
					
					time_predictor[j] = ceil(maxDup/(subArrayDup[j]*peDup[j]*tileDup[j]));
					//time_predictor[j] = maxDup/(subArrayDup[j]*peDup[j]*tileDup[j]);
				}
				int bottleneck_layer = temp1-1;
				for (int j=temp4; j<record[i].size(); j++) {
					if (flag == 1) {
						if (time_predictor[bottleneck_layer] < time_predictor[j]) {
							bottleneck_layer = j;
						}
					}
					else if (flag == 0) {
						for (int k=0; k<record_cannot_accelerate_layer.size(); k++) {
							if ((time_predictor[bottleneck_layer] < time_predictor[j]) && (j != record_cannot_accelerate_layer[k])) {
								bottleneck_layer = j;
								tag = 1;
							}
							else {
								tag = 0;
							}
						}
					}
				}
				//cout << 11 << endl;
				//cout << "bottleneck_layer: " << bottleneck_layer << endl;
				//cout << "extra_tile_num: " << extra_tile_num << endl;
				int maxDup = netStructure[bottleneck_layer][0]*netStructure[bottleneck_layer][1];

				//cout << "tile_num_each_layer[bottleneck_layer]: " << tile_num_each_layer[bottleneck_layer] << endl;

				if (tile_num_each_layer[bottleneck_layer] <= extra_tile_num) {
					flag = 1;
					extra_tile_num = extra_tile_num - tile_num_each_layer[bottleneck_layer];
					//cout << "tile_num_each_layer[bottleneck_layer]: " << tile_num_each_layer[bottleneck_layer] << endl;
					tileDup[bottleneck_layer] = tileDup[bottleneck_layer] +1;

					int temp = tileDup[bottleneck_layer];

					if ((netStructure[bottleneck_layer][0] == 1) && (netStructure[bottleneck_layer][1] == 1)) {
						tileDup[bottleneck_layer] = 1;
						flag = 0;
						break;
					}
					else if (tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer] > maxDup) {
						flag = 0;
						for (int j=0; j<=(temp-1); j++) {
							tileDup[bottleneck_layer] = tileDup[bottleneck_layer]-1;
							if ((1<=tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer]) && (tileDup[bottleneck_layer]*subArrayDup[bottleneck_layer]*peDup[bottleneck_layer]<= maxDup)) {
								break;
							}
						}
					}
				}
				else {
					flag = 0;
					record_cannot_accelerate_layer.push_back(bottleneck_layer);
				}
			}
			//cout << 11 << endl;
			temp3 = temp1;
			temp4 = temp1;
		}

		//cout << 8 << endl;

		if (!param->pipeline_dynamic_dup) {
			for (int i=0; i<netStructure.size(); i++) {
				tileDup[i]=1;
			}
		}
	}
	return tileDup;
}

vector<double> PeDup(const vector<double> &subArrayDup, double peSize, double desiredTileSize, 
								const vector<vector<double> > &netStructure, int numRowPerSynapse, int numColPerSynapse) {

	vector<double> peDup;

	for (int i=0; i<netStructure.size(); i++) {
		int actualDup = 0;

		int maxDup = netStructure[i][0]*netStructure[i][1];

		int numPERow = ceil((double) desiredTileSize/(double) peSize); //number row of PE in one Tile
		int numPECol = ceil((double) desiredTileSize/(double) peSize); //number column of PE in one Tile

		if (subArrayDup[i] > 1) {
			actualDup = numPERow*numPECol;
		}
		// netStructure[i][2]: represent IFM Channel Depth, netStructure[i][3]: represent Kernel Length, netStructure[i][4]: represent Kernel Width, netStructure[i][5]: represent Kernel Depth
		else if ((netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*numRowPerSynapse <= desiredTileSize) || (netStructure[i][5]*numColPerSynapse <= desiredTileSize)) {
			int peForOneMatrixRow = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/(double) peSize); // this layer needs the number row of PE
			int peForOneMatrixCol = ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) peSize); // this layer needs the number column of PE

			actualDup = floor((numPERow*numPECol)/(peForOneMatrixRow*peForOneMatrixCol))==0? 1:floor((numPERow*numPECol)/(peForOneMatrixRow*peForOneMatrixCol));
		}
		else {
			actualDup = 1;
		}

		int flag = actualDup;

		if ((netStructure[i][0] == 1) && (netStructure[i][1] == 1)) {
			actualDup = 1;
			peDup.push_back(actualDup);
		}
		else if ((1 <= actualDup*subArrayDup[i]) && (actualDup*subArrayDup[i] <= maxDup)) {
			peDup.push_back(actualDup);
		}
		else {
			for(int j=1; j<=(flag-1); j++) {
				actualDup = actualDup-1;
				if ((1 <= actualDup*subArrayDup[i]) && (actualDup*subArrayDup[i] <= maxDup)) {
					peDup.push_back(actualDup);
					break;
				}
			}
		}
	}
	return peDup;
}

vector<double> SubArrayDup(double desiredPESizeCM, const vector<vector<double> > &netStructure, int numRowPerSynapse, int numColPerSynapse) {
	vector<double> subArrayDup;
	
	for (int i=0; i<netStructure.size(); i++) {
		int actualDup = 0;

		//int maxDup = (netStructure[i][0]-netStructure[i][3]+1)/netStructure[i][7]*(netStructure[i][1]-netStructure[i][4]+1)/netStructure[i][7];
		int maxDup = netStructure[i][0]*netStructure[i][1];
		
		if ((netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*numRowPerSynapse <= desiredPESizeCM) || (netStructure[i][5]*numColPerSynapse <= desiredPESizeCM)) {
			int arrayForOneMatrixRow = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/(double) param->numRowSubArray); // this layer needs the number row of subarray
			int arrayForOneMatrixCol = ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) param->numColSubArray); // this layer needs the number column of subarray
			int numSubArrayRow = ceil((double) desiredPESizeCM/(double) param->numRowSubArray); // number row of Subarray in one PE
			int numSubArrayCol = ceil((double) desiredPESizeCM/(double) param->numColSubArray); // number column of Subarray in one PE

			actualDup = floor((numSubArrayRow*numSubArrayCol)/(arrayForOneMatrixRow*arrayForOneMatrixCol))==0? 1:floor((numSubArrayRow*numSubArrayCol)/(arrayForOneMatrixRow*arrayForOneMatrixCol)); // the same as PEDesign
		}
		else {
			actualDup = 1;
		}

		int flag = actualDup;

		if ((netStructure[i][0] == 1) && (netStructure[i][1] == 1)) {
			actualDup = 1;
			subArrayDup.push_back(actualDup);
		}
		else if ((1 <= actualDup) && (actualDup <= maxDup)) {
			subArrayDup.push_back(actualDup);
		}
		else {
			for(int j=1; j<=(flag-1); j++) {
				actualDup = actualDup-1;
				if ((1 <= actualDup) && (actualDup <= maxDup)) {
					subArrayDup.push_back(actualDup);
					break;
				}
			}
		}
	}
	return subArrayDup;
}

vector<double> OverallEachLayer(bool utilization, bool speedUp, bool numTile, const vector<double> &tileDup, const vector<double> &peDup, 
										const vector<double> &subArrayDup, double desiredTileSizeCM, const vector<vector<double> > &netStructure, 
										int numRowPerSynapse, int numColPerSynapse) {
	vector<double> numTileEachLayer;
	vector<double> no_tileDup_numTileEachLayer;
	vector<double> utilizationEachLayer;	
	vector<double> speedUpEachLayer;
	
	for (int i=0; i<netStructure.size(); i++) {

		double numtileEachLayer, utilizationEach;
		double no_tileDup_numtileEachLayer;
		// netStructure[i][2]: represent IFM Channel Depth, netStructure[i][3]: represent Kernel Length, netStructure[i][4]: represent Kernel Width, netStructure[i][5]: represent Kernel Depth
		numtileEachLayer = tileDup[i]*ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/desiredTileSizeCM)*
							ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) desiredTileSizeCM);
		no_tileDup_numtileEachLayer = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) numRowPerSynapse/desiredTileSizeCM)*
										ceil((double) netStructure[i][5]*(double) numColPerSynapse/(double) desiredTileSizeCM);
		
		utilizationEach = (tileDup[i]*peDup[i]*subArrayDup[i]*netStructure[i][2]*netStructure[i][3]*netStructure[i][4]
									*numRowPerSynapse*netStructure[i][5]*numColPerSynapse)/(param->numTileRow*param->numTileCol*desiredTileSizeCM*desiredTileSizeCM);
		
		numTileEachLayer.push_back(numtileEachLayer);
		no_tileDup_numTileEachLayer.push_back(no_tileDup_numtileEachLayer);
		utilizationEachLayer.push_back(utilizationEach);
		speedUpEachLayer.push_back(tileDup[i]*peDup[i]*subArrayDup[i]);
	}

	if (utilization) {
		return utilizationEachLayer;
	} 
	else if (speedUp) {
		return speedUpEachLayer;
	} 
	else if (numTile) {
		return numTileEachLayer;
	}
	else {
		return no_tileDup_numTileEachLayer;
	}
}

vector<vector<double> > LoadInWeightData(const string &weightfile, int numRowPerSynapse, int numColPerSynapse, double maxConductance, double minConductance) {
	
	ifstream fileone(weightfile.c_str());                           
	string lineone;
	string valone;
	
	int ROW = 0;
	int COL = 0;
	
	if (!fileone.good()) {                                       
		cerr << "Error: the fileone cannot be opened!" << endl;
		exit(1);
	}else{
		while (getline(fileone, lineone, '\n')) {                   
			ROW++;                                             
		}
		fileone.clear();
		fileone.seekg(0, ios::beg);                               
		if (getline(fileone, lineone, '\n')) {                      
			istringstream iss (lineone);                         
			while (getline(iss, valone, ',')) {                   
				COL++;
			}
		}	
	}
	fileone.clear();
	fileone.seekg(0, ios::beg);                   
	
	
	double NormalizedMin = 0;
	double NormalizedMax = pow(2, param->synapseBit);
	
	double RealMax = param->algoWeightMax;
	double RealMin = param->algoWeightMin;
	
	vector<vector<double> > weight;            
	// load the data into a weight matrix ...
	for (int row=0; row<ROW; row++) {	
		vector<double> weightrow;
		vector<double> weightrowb;
		getline(fileone, lineone, '\n');              
		istringstream iss;
		iss.str(lineone);
		for (int col=0; col<COL; col++) {       
			while(getline(iss, valone, ',')){	
				istringstream fs;
				fs.str(valone);
				double f=0;
				fs >> f;	
				//normalize weight to integer
				double newdata = ((NormalizedMax-NormalizedMin)/(RealMax-RealMin)*(f-RealMax)+NormalizedMax);
				if (newdata >= 0) {
					newdata += 0.5;
				}else {
					newdata -= 0.5;
				}
				//cout << "newdata: " << newdata << endl;
				// map and expend the weight in memory array
				int cellrange = pow(2, param->cellBit);
				vector<double> synapsevector(numColPerSynapse);       
				int value = newdata; 
				
				if (param->BNNparallelMode) {
					if (value == 1) {
						weightrow.push_back(maxConductance);
						weightrow.push_back(minConductance);
					} else {
						weightrow.push_back(minConductance);
						weightrow.push_back(maxConductance);
					}
				} else if (param->XNORparallelMode || param->XNORsequentialMode) {
					if (value == 1) {
						weightrow.push_back(maxConductance);
						weightrowb.push_back(minConductance);
					} else {
						weightrow.push_back(minConductance);
						weightrowb.push_back(maxConductance);
					}
				} else {
					int remainder;   
					for (int z=0; z<numColPerSynapse; z++) {   
						remainder = ceil((double)(value%cellrange));
						value = ceil((double)(value/cellrange));
						synapsevector.insert(synapsevector.begin(), remainder);
					}
					for (int u=0; u<numColPerSynapse; u++) {
						double cellvalue = synapsevector[u];
						double conductance = cellvalue/(cellrange-1) * (maxConductance-minConductance) + minConductance;
						//cout << "conductance: " << conductance << endl; 
						weightrow.push_back(conductance);
					}
				}
			}
		}
		if (param->XNORparallelMode || param->XNORsequentialMode) {
			weight.push_back(weightrow);
			weightrow.clear();
			weight.push_back(weightrowb);
			weightrowb.clear();
		} else {
			weight.push_back(weightrow);
			weightrow.clear();
		}
	}
	fileone.close();
	
	return weight;
	weight.clear();
}

vector<vector<double> > CopyArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol) {
	
	vector<vector<double> > copy;
	for (int i=0; i<numRow; i++) {
		vector<double> copyRow;
		for (int j=0; j<numCol; j++) {
			copyRow.push_back(orginal[positionRow+i][positionCol+j]);
		}
		copy.push_back(copyRow);
		copyRow.clear();
	}
	
	return copy;
	copy.clear();
} 

vector<vector<double> > ReshapeArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol, int numPE, int weightMatrixRow) {
	
	vector<vector<double> > copy;

	for (int k=0; k<numPE; k++) {
		for (int i=0; i<numRow; i++) {
			vector<double> copyRow;
			for (int j=0; j<numCol; j++) {
				copyRow.push_back(orginal[positionRow+k*weightMatrixRow+i][positionCol+j]);
			}
			copy.push_back(copyRow);
			copyRow.clear();
		}
	}
	
	return copy;
	copy.clear();
}

vector<vector<double> > LoadInInputData(const string &inputfile) {
	
	ifstream infile(inputfile.c_str());     
	string inputline;
	string inputval;
	
	int ROWin=0, COLin=0;      
	if (!infile.good()) {       
		cerr << "Error: the input file cannot be opened!" << endl;
		exit(1);
	}else{
		while (getline(infile, inputline, '\n')) {      
			ROWin++;                               
		}
		infile.clear();
		infile.seekg(0, ios::beg);    
		if (getline(infile, inputline, '\n')) {        
			istringstream iss (inputline);      
			while (getline(iss, inputval, ',')) {       
				COLin++;
			}
		}	
	}
	infile.clear();
	infile.seekg(0, ios::beg);          

	vector<vector<double> > inputvector;              
	// load the data into inputvector ...
	for (int row=0; row<ROWin; row++) {	
		vector<double> inputvectorrow;
		vector<double> inputvectorrowb;
		getline(infile, inputline, '\n');             
		istringstream iss;
		iss.str(inputline);
		for (int col=0; col<COLin; col++) {
			while(getline(iss, inputval, ',')){	
				istringstream fs;
				fs.str(inputval);
				double f=0;
				fs >> f;
				
				if (param->BNNparallelMode) {
					if (f == 1) {
						inputvectorrow.push_back(1);
					} else {
						inputvectorrow.push_back(0);
					}
				} else if (param->XNORparallelMode || param->XNORsequentialMode) {
					if (f == 1) {
						inputvectorrow.push_back(1);
						inputvectorrowb.push_back(0);
					} else {
						inputvectorrow.push_back(0);
						inputvectorrowb.push_back(1);
					}
				} else {
					inputvectorrow.push_back(f);
				}
			}
		}
		if (param->XNORparallelMode || param->XNORsequentialMode) {
			inputvector.push_back(inputvectorrow);
			inputvectorrow.clear();
			inputvector.push_back(inputvectorrowb);
			inputvectorrowb.clear();
		} else {
			inputvector.push_back(inputvectorrow);
			inputvectorrow.clear();
		}
	}
	// close the input file ...
	infile.close();
	
	return inputvector;
	inputvector.clear();
}

vector<vector<double> > CopyInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow) {
	
	vector<vector<double> > copy;
	for (int i=0; i<numRow; i++) {
		vector<double> copyRow;
		for (int j=0; j<numInputVector; j++) {
			copyRow.push_back(orginal[positionRow+i][j]);
		}
		copy.push_back(copyRow);
		copyRow.clear();
	}
	
	return copy;
	copy.clear();
	
} 

vector<vector<double> > ReshapeInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow, int numPE, int weightMatrixRow) {
	
	vector<vector<double> > copy;

	for (int k=0; k<numPE; k++) {
		for (int i=0; i<numRow; i++) {
			vector<double> copyRow;
			for (int j=0; j<numInputVector; j++) {
				copyRow.push_back(orginal[positionRow+k*weightMatrixRow+i][j]);
			}
			copy.push_back(copyRow);
			copyRow.clear();
		}
	}
	
	return copy;
	copy.clear();
} 











