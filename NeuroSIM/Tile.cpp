#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "Sigmoid.h"
#include "BitShifter.h"
#include "AdderTree.h"
#include "Buffer.h"
#include "HTree.h"
#include "ProcessingUnit.h"
#include "SubArray.h"
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "Tile.h"

using namespace std;

extern Param *param;
// Anni update
int numInBufferCoreCM = 0;
int numOutBufferCoreCM = 0;	
int numInBufferCoreNM = 0;
int numOutBufferCoreNM = 0;									 

SubArray *subArrayInPE;
Buffer *inputBufferCM;
Buffer *outputBufferCM;
HTree *hTreeCM;
AdderTree *accumulationCM;
Sigmoid *sigmoidCM;
BitShifter *reLuCM;
HTree *hTreeNM;
AdderTree *accumulationNM;
Sigmoid *sigmoidNM;
BitShifter *reLuNM;

void TileInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell, double _numPECM, double _peSizeCM){
	
	subArrayInPE = new SubArray(inputParameter, tech, cell);
	hTreeNM = new HTree(inputParameter, tech, cell);
	accumulationNM = new AdderTree(inputParameter, tech, cell);
	inputBufferCM = new Buffer(inputParameter, tech, cell);
	outputBufferCM = new Buffer(inputParameter, tech, cell);
	hTreeCM = new HTree(inputParameter, tech, cell);
	accumulationCM = new AdderTree(inputParameter, tech, cell);
	// Anni update
	reLuNM = new BitShifter(inputParameter, tech, cell);
	reLuCM = new BitShifter(inputParameter, tech, cell);
	sigmoidNM = new Sigmoid(inputParameter, tech, cell);
	sigmoidCM = new Sigmoid(inputParameter, tech, cell);

	
	/*** Parameters ***/
	double numPECM, peSizeCM, numSubArrayCM;
	int numRowPerSynapse, numColPerSynapse;
	
	numPECM = _numPECM;
	peSizeCM = _peSizeCM;
	numRowPerSynapse = param->numRowPerSynapse;
	numColPerSynapse = param->numColPerSynapse;
	
	/*** Initialize ProcessingUnit ***/
	numSubArrayCM = ceil((double)peSizeCM/(double)param->numRowSubArray)*ceil((double)peSizeCM/(double)param->numColSubArray);

	ProcessingUnitInitialize(subArrayInPE, inputParameter, tech, cell, ceil(sqrt(numSubArrayCM)), ceil(sqrt(numSubArrayCM)));

	// Anni update: numBitPEOutput
	int numBitSubarrayOutput, numBitPEOutputCM, numBitPEOutputNM;
	if (param->parallelRead) {		
		numBitSubarrayOutput = log2((double)param->levelOutput)+ceil(log2(ceil(param->numRowSubArray/param->numRowParallel)))+param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	} else{
		numBitSubarrayOutput = ceil(log2((double)param->numRowSubArray))+param->cellBit+param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	}

	// Anni update: numBitPEOutputCM, numOutBufferCoreCM, numInBufferCoreCM
	numBitPEOutputCM = numBitSubarrayOutput + ceil(sqrt(numSubArrayCM));	

    accumulationCM->Initialize(numPECM, numBitPEOutputCM, ceil((double)numPECM*(double)peSizeCM/(double)param->numColMuxed), param->clkFreq);
	numOutBufferCoreCM = ceil(((numBitPEOutputCM+ceil(log2((double)numPECM)))*numPECM*param->numColSubArray/param->numColMuxed)/(param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol));

	if (((numBitPEOutputCM+ceil(log2((double)numPECM)))*numPECM*param->numColSubArray/param->numColMuxed) < (param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol)) {
		outputBufferCM->Initialize((numBitPEOutputCM+ceil(log2((double)numPECM)))*numPECM*param->numColSubArray/param->numColMuxed, (numBitPEOutputCM+ceil(log2((double)numPECM)))*numPECM, 1, param->unitLengthWireResistance, param->clkFreq, param->peBufferType);
	} 
	else {
		outputBufferCM->Initialize((param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol), param->tileBufferCoreSizeCol, 1, param->unitLengthWireResistance, param->clkFreq, param->peBufferType);
	}	

	numInBufferCoreCM = ceil((numPECM*param->numBitInput*param->numRowSubArray)/(param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol));
	
    // 230920 update
	if ((numPECM*param->numBitInput*param->numRowSubArray) < (param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol)) {
		inputBufferCM->Initialize(numPECM*param->numBitInput*param->numRowSubArray, numPECM*param->numRowSubArray, 1, param->unitLengthWireResistance, param->clkFreq, param->peBufferType);
	} 
	else {
		inputBufferCM->Initialize((param->tileBufferCoreSizeRow*param->tileBufferCoreSizeCol), param->tileBufferCoreSizeCol, 1, param->unitLengthWireResistance, param->clkFreq, param->peBufferType);
	}   

	 hTreeCM->Initialize(numPECM, numPECM, param->localBusDelayTolerance, numPECM*(double)peSizeCM, param->clkFreq);

}

vector<double> TileCalculateArea(double numPE, double peSize, double *height, double *width) {
	double area = 0;
	double PEheight, PEwidth, PEbufferArea;
	*height = 0;
	*width = 0;
	vector<double> areaResults;
	vector<double> peAreaResults;
	double areareLu = 0;
	double areasigmoid = 0;
	int numSubArray = ceil((double) peSize/(double) param->numRowSubArray)*ceil((double) peSize/(double) param->numColSubArray);
	
	peAreaResults = ProcessingUnitCalculateArea(subArrayInPE, ceil((double)sqrt((double)numSubArray)), ceil((double)sqrt((double)numSubArray)), &PEheight, &PEwidth, &PEbufferArea);
	
	double PEarea = peAreaResults[0];
	double PEareaADC = peAreaResults[1];
	double PEareaAccum = peAreaResults[2];
	double PEareaOther = peAreaResults[3];
	double PEareaArray = peAreaResults[4];

	accumulationCM->CalculateArea(NULL, ceil(sqrt((double)numPE))*PEwidth, NONE);
	inputBufferCM->CalculateArea(ceil(sqrt((double)numPE))*PEheight, NULL, NONE);
	outputBufferCM->CalculateArea(NULL, ceil(sqrt((double)numPE))*PEwidth, NONE);
	
	inputBufferCM->area *= numInBufferCoreCM;
	outputBufferCM->area *= numOutBufferCoreCM;												  
	hTreeCM->CalculateArea(PEheight, PEwidth, 16);

	area += PEarea*numPE + accumulationCM->area + inputBufferCM->area + outputBufferCM->area + hTreeCM->area; // the area of one Tile
	
	*height = sqrt(area);
	*width = area/(*height);
	
	areaResults.push_back(area);
	areaResults.push_back(hTreeCM->area);
	areaResults.push_back(PEareaADC*numPE);
	areaResults.push_back(PEareaAccum*numPE + accumulationCM->area);
	areaResults.push_back(PEareaOther*numPE + inputBufferCM->area + outputBufferCM->area + areareLu + areasigmoid);
	areaResults.push_back(PEareaArray*numPE);
	
	return areaResults;
}

void TileCalculatePerformance(const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, const vector<vector<double> > &inputVector, double numPE, 
							double peSize, int speedUp, int PespeedUp, int weightMatrixRow, int weightMatrixCol, int numInVector, Technology& tech, MemCell& cell,
							double *readLatency, double *readDynamicEnergy, double *writeLatency, double *writeDynamicEnergy, double *leakage, double *leakageSRAMInUse, double *bufferLatency, 
							double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy, double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, double *coreEnergyADC, 
							double *coreEnergyAccum, double *coreEnergyOther, bool CalculateclkFreq, double *clkPeriod, double *computationLatency) {

	/*** sweep PE ***/
	int numRowPerSynapse, numColPerSynapse;

	numRowPerSynapse = param->numRowPerSynapse;
	numColPerSynapse = param->numColPerSynapse;

	double PEreadLatency, PEreadDynamicEnergy, PEwriteLatency, PEwriteDynamicEnergy, PEleakage, PEleakageSRAMInUse, PEbufferLatency, PEbufferDynamicEnergy, PEicLatency, PEicDynamicEnergy;
	double peLatencyADC, peLatencyAccum, peLatencyOther, peEnergyADC, peEnergyAccum, peEnergyOther, peComputationLatency;

	int numSubArrayRow = ceil((double)peSize/(double)param->numRowSubArray);
	int numSubArrayCol = ceil((double)peSize/(double)param->numColSubArray);
	
	*readLatency = 0;
	*readDynamicEnergy = 0;
	*writeLatency = 0;
	*writeDynamicEnergy = 0;
	*leakage = 0;
	*leakageSRAMInUse = 0;
	*bufferLatency = 0;
	*bufferDynamicEnergy = 0;
	*icLatency = 0;
	*icDynamicEnergy = 0;
	*coreEnergyADC = 0;
	*coreEnergyAccum = 0;
	*coreEnergyOther = 0;
	*coreLatencyADC = 0;
	*coreLatencyAccum = 0;
	*coreLatencyOther = 0;

	*computationLatency = 0;

	// Anni update: update Clock frequency
	if(!CalculateclkFreq) {	
		inputBufferCM->clkFreq = param->clkFreq; 
		outputBufferCM->clkFreq = param->clkFreq; 
		hTreeCM->clkFreq = param->clkFreq; 
		accumulationCM->clkFreq = param->clkFreq; 
		reLuCM->clkFreq = param->clkFreq; 
		sigmoidCM->clkFreq = param->clkFreq; 
		
		hTreeNM->clkFreq = param->clkFreq; 
		accumulationNM->clkFreq = param->clkFreq; 	
		reLuNM->clkFreq = param->clkFreq; 
		sigmoidNM->clkFreq = param->clkFreq; 	
	}

	if (speedUp > 1) {
		if (PespeedUp == numPE*numPE) {
			// duplication in PE or subArray --> tell each PE to take the whole assigned weight  --> "fully" duplication
			// assign weight and input to specific pe
			vector<vector<double> > pEMemory;
			pEMemory = CopyPEArray(newMemory, 0, 0, weightMatrixRow, weightMatrixCol);

			vector<vector<double> > pEMemoryOld;
			pEMemoryOld = CopyPEArray(oldMemory, 0, 0, weightMatrixRow, weightMatrixCol);

			vector<vector<double> > pEInput;
			pEInput = CopyPEInput(inputVector, 0, numInVector, weightMatrixRow);
		
			ProcessingUnitCalculatePerformance(subArrayInPE, tech, cell, pEMemory, pEMemoryOld, pEInput, ceil((double)speedUp/(double)PespeedUp), 
										numSubArrayRow, numSubArrayCol, weightMatrixRow, weightMatrixCol, numInVector,
										&PEreadLatency, &PEreadDynamicEnergy, &PEwriteLatency, &PEwriteDynamicEnergy, &PEleakage, &PEleakageSRAMInUse,
										&PEbufferLatency, &PEbufferDynamicEnergy, &PEicLatency, &PEicDynamicEnergy,
										&peLatencyADC, &peLatencyAccum, &peLatencyOther, &peEnergyADC, &peEnergyAccum, &peEnergyOther, CalculateclkFreq, clkPeriod, &peComputationLatency);
			*readLatency       = PEreadLatency/(numPE*numPE); // further speed up in PE level
			*readDynamicEnergy = PEreadDynamicEnergy        ; // since subArray.cpp takes all input vectors, no need to *numPE here

			*bufferLatency       = PEbufferLatency/(numPE*numPE); //#cycles of PE-level buffers (DFF)
			*bufferDynamicEnergy = PEbufferDynamicEnergy        ;
			*icLatency           = PEicLatency/(numPE*numPE)    ; //s
			*icDynamicEnergy     = PEicDynamicEnergy            ;
			
			*coreLatencyADC   = peLatencyADC/(numPE*numPE)  ; //#sensing cycles
			*coreLatencyAccum = peLatencyAccum/(numPE*numPE); //#cycles
			*coreLatencyOther = peLatencyOther/(numPE*numPE);
			
			*coreEnergyADC   = peEnergyADC  ;
			*coreEnergyAccum = peEnergyAccum;
			*coreEnergyOther = peEnergyOther;

			*computationLatency = peComputationLatency/(numPE*numPE);

			*writeLatency       = PEwriteLatency                    ;
			*writeDynamicEnergy = PEwriteDynamicEnergy*(numPE*numPE);

			//cout << "*readLatency: " << *readLatency << endl;
			// no accumulation access
		} 
		else {
			// # duplication is smaller then # PE, means only a group of PE take the assigned weight  --> not "fully" duplication
			// also need to redefine a few data-grab start-point
			for (int i=0; i<ceil((double)weightMatrixRow/(double)peSize); i++) {
				for (int j=0; j<ceil((double)weightMatrixCol/(double)peSize); j++) {
					if ( (i*peSize < weightMatrixRow) && (j*peSize < weightMatrixCol) ) {
						int numRowMatrix = min(peSize, (double) weightMatrixRow-i*peSize);
						int numColMatrix = min(peSize, (double) weightMatrixCol-j*peSize);
				
						// assign weight and input to specific tile
						vector<vector<double> > pEMemory;
						pEMemory = CopyPEArray(newMemory, i*peSize, j*peSize, numRowMatrix, numColMatrix);

						vector<vector<double> > pEMemoryOld;
						pEMemoryOld = CopyPEArray(oldMemory, i*peSize, j*peSize, numRowMatrix, numColMatrix);

						vector<vector<double> > pEInput;
						pEInput = CopyPEInput(inputVector, i*peSize, numInVector, numRowMatrix);
						
						ProcessingUnitCalculatePerformance(subArrayInPE, tech, cell, pEMemory, pEMemoryOld, pEInput, ceil((double)speedUp/(double)PespeedUp), 
											numSubArrayRow, numSubArrayCol, numRowMatrix, numColMatrix, numInVector,
											&PEreadLatency, &PEreadDynamicEnergy, &PEwriteLatency, &PEwriteDynamicEnergy, &PEleakage, &PEleakageSRAMInUse,
											&PEbufferLatency, &PEbufferDynamicEnergy, &PEicLatency, &PEicDynamicEnergy,
											&peLatencyADC, &peLatencyAccum, &peLatencyOther, &peEnergyADC, &peEnergyAccum, &peEnergyOther, CalculateclkFreq, clkPeriod, &peComputationLatency);

						*readLatency          = MAX(PEreadLatency, (*readLatency))    ;
						*readDynamicEnergy   += PEreadDynamicEnergy                   ;
						*bufferLatency        = MAX(PEbufferLatency, (*bufferLatency));
						*bufferDynamicEnergy += PEbufferDynamicEnergy                 ;
						*icLatency            = MAX(PEicLatency,(*icLatency))         ;
						*icDynamicEnergy     += PEicDynamicEnergy                     ;
						
						*coreLatencyADC   = MAX(peLatencyADC, (*coreLatencyADC))    ;
						*coreLatencyAccum = MAX(peLatencyAccum, (*coreLatencyAccum));
						*coreLatencyOther = MAX(peLatencyOther, (*coreLatencyOther));
						
						*coreEnergyADC   += peEnergyADC  ;
						*coreEnergyAccum += peEnergyAccum;
						*coreEnergyOther += peEnergyOther;

						*computationLatency = MAX(peComputationLatency, (*computationLatency));

						*writeLatency       += PEwriteLatency      ;
						*writeDynamicEnergy += PEwriteDynamicEnergy;
					}
				}
			}
			*readLatency      /= (PespeedUp);   // further speedup in PE level
			*coreLatencyADC   /= (PespeedUp);
			*coreLatencyAccum /= (PespeedUp);
			*coreLatencyOther /= (PespeedUp);
			*bufferLatency    /= (PespeedUp);
			*icLatency        /= (PespeedUp);

			*computationLatency = *computationLatency/(PespeedUp);

			//*writeLatency       *= PespeedUp;
			*writeDynamicEnergy = (*writeDynamicEnergy)*PespeedUp;

			//cout << "*readLatency: " << *readLatency << endl;

			// whether go through accumulation?
			if (ceil((double)weightMatrixRow/(double)peSize) > 1) {
				accumulationCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil((double)param->numColMuxed/(double)param->numColPerSynapse), PespeedUp*ceil((double)weightMatrixRow/(double)peSize), 0);
				accumulationCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil((double)param->numColMuxed/(double)param->numColPerSynapse), PespeedUp*ceil((double)weightMatrixRow/(double)peSize));

				*readLatency       += accumulationCM->readLatency      ; 
				*coreLatencyAccum  += accumulationCM->readLatency      ; 
				*readDynamicEnergy += accumulationCM->readDynamicEnergy;					
				*coreEnergyAccum   += accumulationCM->readDynamicEnergy;
			}
		}
	} 
	else {
		// no duplication --> tell PE to further partition the weight and grab data (redefine a few data-grab start-point)
		for (int i=0; i<numPE; i++) {
			for (int j=0; j<numPE; j++) {
				// each cycle assign to different PE
				if ( (i*peSize < weightMatrixRow) && (j*peSize < weightMatrixCol) ) {
					// assign weight and input to specific pe
					int numRowMatrix = min(peSize, (double) weightMatrixRow-i*peSize);
					int numColMatrix = min(peSize, (double) weightMatrixCol-j*peSize);

					vector<vector<double> > pEMemory;
					pEMemory = CopyPEArray(newMemory, i*peSize, j*peSize, numRowMatrix, numColMatrix);

					vector<vector<double> > pEMemoryOld;
					pEMemoryOld = CopyPEArray(oldMemory, i*peSize, j*peSize, numRowMatrix, numColMatrix);

					vector<vector<double> > pEInput;
					pEInput = CopyPEInput(inputVector, i*peSize, numInVector, numRowMatrix);
					
					ProcessingUnitCalculatePerformance(subArrayInPE, tech, cell, pEMemory, pEMemoryOld, pEInput, 1, numSubArrayRow, numSubArrayCol, numRowMatrix,
											numColMatrix, numInVector, &PEreadLatency, &PEreadDynamicEnergy, &PEwriteLatency, &PEwriteDynamicEnergy, &PEleakage, &PEleakageSRAMInUse,
											&PEbufferLatency, &PEbufferDynamicEnergy, &PEicLatency, &PEicDynamicEnergy, 
											&peLatencyADC, &peLatencyAccum, &peLatencyOther, &peEnergyADC, &peEnergyAccum, &peEnergyOther, CalculateclkFreq, clkPeriod, &peComputationLatency);
					
					*readLatency        = max(PEreadLatency, (*readLatency));
					*readDynamicEnergy += PEreadDynamicEnergy               ;
					
					*bufferLatency        = max(PEbufferLatency, (*bufferLatency));
					*bufferDynamicEnergy += PEbufferDynamicEnergy                 ;
					*icLatency            = max(PEicLatency,(*icLatency))         ;
					*icDynamicEnergy     += PEicDynamicEnergy                     ;
					
					*coreLatencyADC   = MAX(peLatencyADC, (*coreLatencyADC))    ;
					*coreLatencyAccum = MAX(peLatencyAccum, (*coreLatencyAccum));
					*coreLatencyOther = MAX(peLatencyOther, (*coreLatencyOther));
					
					*coreEnergyADC   += peEnergyADC  ;
					*coreEnergyAccum += peEnergyAccum;
					*coreEnergyOther += peEnergyOther;

					*computationLatency = MAX(peComputationLatency, (*computationLatency));

					*writeLatency       += PEwriteLatency      ;
					*writeDynamicEnergy += PEwriteDynamicEnergy;
				}
			}
		}

		accumulationCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil((double)param->numColMuxed/(double)param->numColPerSynapse), numPE, 0);
		accumulationCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil((double)param->numColMuxed/(double)param->numColPerSynapse), numPE);
		
		*readLatency       += accumulationCM->readLatency      ;
		*coreLatencyAccum  += accumulationCM->readLatency      ;	
		*readDynamicEnergy += accumulationCM->readDynamicEnergy;			
		*coreEnergyAccum   += accumulationCM->readDynamicEnergy;

		//cout << "*readLatency: " << *readLatency << endl;
	}
	//cout << "*readLatency1: " << *readLatency << endl;
	if(!CalculateclkFreq){
		double numBitToLoadOut, numBitToLoadIn;
		
		numBitToLoadIn = MAX(ceil(weightMatrixCol/param->numColPerSynapse)*(accumulationCM->numStage+accumulationCM->numAdderBit)*numInVector/param->numBitInput, 0);

		outputBufferCM->CalculateLatency(outputBufferCM->interface_width, numBitToLoadIn/outputBufferCM->interface_width/2.0, outputBufferCM->interface_width, numBitToLoadIn/outputBufferCM->interface_width/2.0);
		outputBufferCM->CalculatePower(outputBufferCM->interface_width, numBitToLoadIn/outputBufferCM->interface_width/2.0, outputBufferCM->interface_width, numBitToLoadIn/outputBufferCM->interface_width/2.0);
		//considering buffer activation: no matter speedup or not, the total number of data transferred is fixed
		numBitToLoadOut = MAX(weightMatrixRow*numInVector, 0);

		inputBufferCM->CalculateLatency(inputBufferCM->interface_width, numBitToLoadOut/inputBufferCM->interface_width/2.0, inputBufferCM->interface_width, numBitToLoadOut/inputBufferCM->interface_width/2.0);
		inputBufferCM->CalculatePower(inputBufferCM->interface_width, numBitToLoadOut/inputBufferCM->interface_width/2.0, inputBufferCM->interface_width, numBitToLoadOut/inputBufferCM->interface_width/2.0);
		// since multi-core buffer has improve the parallelism
		// Anni update: numInBufferCoreCM, numOutBufferCoreCM
		//cout<<"numBitToLoadOut: "<<numBitToLoadOut<<"	numBitToLoadIn: "<<numBitToLoadIn<<endl;
		//cout<<"inputBufferCM->readLatency before: "<<inputBufferCM->readLatency<<"	outputBufferCM->readLatency: "<<outputBufferCM->readLatency<<endl;
		inputBufferCM->readLatency /= MIN(numInBufferCoreCM, ceil(hTreeCM->busWidth/inputBufferCM->interface_width));
		inputBufferCM->writeLatency /= MIN(numInBufferCoreCM, ceil(hTreeCM->busWidth/inputBufferCM->interface_width));
		outputBufferCM->readLatency /= MIN(numOutBufferCoreCM, ceil(hTreeCM->busWidth/outputBufferCM->interface_width));
		outputBufferCM->writeLatency /= MIN(numOutBufferCoreCM, ceil(hTreeCM->busWidth/outputBufferCM->interface_width));
		//cout<<"numInBufferCoreCM: "<<numInBufferCoreCM<<"	ceil(hTreeCM->busWidth/inputBufferCM->interface_width): "<<ceil(hTreeCM->busWidth/inputBufferCM->interface_width)<<endl;
		//cout<<"numOutBufferCoreCM: "<<numOutBufferCoreCM<<"	ceil(hTreeCM->busWidth/outputBufferCM->interface_width): "<<ceil(hTreeCM->busWidth/outputBufferCM->interface_width)<<endl;	
		//cout<<"inputBufferCM->interface_width: "<<inputBufferCM->interface_width<<"	outputBufferCM->interface_width: "<<outputBufferCM->interface_width<<endl;
		//cout<<"inputBufferCM->readLatency after: "<<inputBufferCM->readLatency<<"	outputBufferCM->readLatency: "<<outputBufferCM->readLatency<<endl;

		//cout << "inputBufferCM->readLatency + inputBufferCM->writeLatency: " << inputBufferCM->readLatency + inputBufferCM->writeLatency << endl;
		//cout << "outputBufferCM->readLatency + outputBufferCM->writeLatency: " << outputBufferCM->readLatency + outputBufferCM->writeLatency << endl;
		*readLatency += (inputBufferCM->readLatency + inputBufferCM->writeLatency);
		*readLatency += (outputBufferCM->readLatency + outputBufferCM->writeLatency);

		*readDynamicEnergy += inputBufferCM->readDynamicEnergy + inputBufferCM->writeDynamicEnergy;
		*readDynamicEnergy += outputBufferCM->readDynamicEnergy + outputBufferCM->writeDynamicEnergy;
		// used to define travel distance
		double PEheight, PEwidth, PEbufferArea;
		int numSubArray = ceil((double) peSize/(double) param->numRowSubArray)*ceil((double) peSize/(double) param->numColSubArray);
		vector<double> PEarea;
		
		PEarea = ProcessingUnitCalculateArea(subArrayInPE, ceil((double)sqrt((double)numSubArray)), ceil((double)sqrt((double)numSubArray)), &PEheight, &PEwidth, &PEbufferArea);
		hTreeCM->CalculateLatency(NULL, NULL, NULL, NULL, PEheight, PEwidth, (numBitToLoadOut+numBitToLoadIn)/hTreeCM->busWidth);
		hTreeCM->CalculatePower(NULL, NULL, NULL, NULL, PEheight, PEwidth, hTreeCM->busWidth, (numBitToLoadOut * param->inputtoggle +numBitToLoadIn * param->outputtoggle )/hTreeCM->busWidth);	  
		//cout << "hTreeCM->readLatency: " << hTreeCM->readLatency << endl;
		*readLatency += hTreeCM->readLatency;
		*bufferLatency += (inputBufferCM->readLatency + outputBufferCM->readLatency + inputBufferCM->writeLatency + outputBufferCM->writeLatency);
		*icLatency += hTreeCM->readLatency;
		*coreLatencyOther += (inputBufferCM->readLatency + inputBufferCM->writeLatency + outputBufferCM->readLatency + outputBufferCM->writeLatency + hTreeCM->readLatency);

		*readDynamicEnergy += hTreeCM->readDynamicEnergy;
		*bufferDynamicEnergy += inputBufferCM->readDynamicEnergy + outputBufferCM->readDynamicEnergy + inputBufferCM->writeDynamicEnergy + outputBufferCM->writeDynamicEnergy;
		*icDynamicEnergy += hTreeCM->readDynamicEnergy;
		*coreEnergyOther += inputBufferCM->readDynamicEnergy + inputBufferCM->writeDynamicEnergy + outputBufferCM->readDynamicEnergy + outputBufferCM->writeDynamicEnergy + hTreeCM->readDynamicEnergy;
		
		*leakage = PEleakage*numPE*numPE + accumulationCM->leakage + inputBufferCM->leakage + outputBufferCM->leakage + hTreeCM->leakage;
		*leakageSRAMInUse = PEleakageSRAMInUse*numPE*numPE;
	} 
	//cout << "*readLatency2: " << *readLatency << endl;
}

vector<vector<double> > CopyPEArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol) {
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

vector<vector<double> > CopyPEInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow) {
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

