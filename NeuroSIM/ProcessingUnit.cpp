#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "Bus.h"
#include "SubArray.h"
#include "constant.h"
#include "formula.h"
#include "ProcessingUnit.h"
#include "Param.h"
#include "AdderTree.h"
#include "Bus.h"
#include "DFF.h"

using namespace std;

extern Param *param;

AdderTree *adderTreeNM;
Bus *busInputNM;
Bus *busOutputNM;
DFF *bufferInputNM;
DFF *bufferOutputNM;

AdderTree *adderTreeCM;
Bus *busInputCM;
Bus *busOutputCM;
DFF *bufferInputCM;
DFF *bufferOutputCM;

void ProcessingUnitInitialize(SubArray *& subArray, InputParameter& inputParameter, Technology& tech, MemCell& cell, int _numSubArrayRowCM, int _numSubArrayColCM) {

	/*** circuit level parameters ***/
	switch(param->memcelltype) {
		case 3:     cell.memCellType = Type::FeFET; break;
		case 2:	    cell.memCellType = Type::RRAM; break;
		case 1:	    cell.memCellType = Type::SRAM; break;
		case -1:	break;
		default:	exit(-1);
	}
	switch(param->accesstype) {
		case 4:	    cell.accessType = none_access;  break;
		case 3:	    cell.accessType = diode_access; break;
		case 2:	    cell.accessType = BJT_access;   break;
		case 1:	    cell.accessType = CMOS_access;  break;
		case -1:	break;
		default:	exit(-1);
	}				
					
	switch(param->transistortype) {
		case 3:	    inputParameter.transistorType = TFET;          break;
		case 2:	    inputParameter.transistorType = FET_2D;        break;
		case 1:	    inputParameter.transistorType = conventional;  break;
		case -1:	break;
		default:	exit(-1);
	}
	
	switch(param->deviceroadmap) {
		case 2:	    inputParameter.deviceRoadmap = LSTP;  break;
		case 1:	    inputParameter.deviceRoadmap = HP;    break;
		case -1:	break;
		default:	exit(-1);
	}
	
	subArray = new SubArray(inputParameter, tech, cell);
	adderTreeNM = new AdderTree(inputParameter, tech, cell);
	busInputNM = new Bus(inputParameter, tech, cell);
	busOutputNM = new Bus(inputParameter, tech, cell);
	bufferInputNM = new DFF(inputParameter, tech, cell);
	bufferOutputNM = new DFF(inputParameter, tech, cell);
	adderTreeCM = new AdderTree(inputParameter, tech, cell);
	busInputCM = new Bus(inputParameter, tech, cell);
	busOutputCM = new Bus(inputParameter, tech, cell);
	bufferInputCM = new DFF(inputParameter, tech, cell);
	bufferOutputCM = new DFF(inputParameter, tech, cell);
		
	/* Create SubArray object and link the required global objects (not initialization) */
	inputParameter.temperature = param->temp;   // Temperature (K)
	inputParameter.processNode = param->technode;    // Technology node
	tech.Initialize(inputParameter.processNode, inputParameter.deviceRoadmap, inputParameter.transistorType);
	
	cell.resistanceOn = param->resistanceOn;	                                // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
	cell.resistanceOff = param->resistanceOff;	                                // Roff resistance at Vr in the reported measurement dat (need to recalculate below if considering the nonlinearity)
	cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;            // Average resistance (for energy estimation)
	cell.readVoltage = param->readVoltage;	                                    // On-chip read voltage for memory cell
	cell.readPulseWidth = param->readPulseWidth;
	cell.accessVoltage = param->accessVoltage;                                       // Gate voltage for the transistor in 1T1R
	cell.resistanceAccess = param->resistanceAccess;
	cell.featureSize = param->featuresize; 
	//cell.writeVoltage = param->writeVoltage;

	cell.maxNumLevelLTP = param->maxNumLevelLTP;	                            // Maximum number of conductance states during LTP or weight increase
	cell.maxNumLevelLTD = param->maxNumLevelLTD;	                            // Maximum number of conductance states during LTD or weight decrease
	double writeVoltageLTP = param->writeVoltage;
	double writeVoltageLTD = param->writeVoltage;
	cell.writeVoltage = sqrt(writeVoltageLTP * writeVoltageLTP + writeVoltageLTD * writeVoltageLTD);    // Use an average value of write voltage 
	double writePulseWidthLTP = param->writePulseWidth;
	double writePulseWidthLTD = param->writePulseWidth;
	cell.writePulseWidth = (writePulseWidthLTP + writePulseWidthLTD) / 2;

	if (cell.memCellType == Type::SRAM) {   // SRAM
		cell.heightInFeatureSize = param->heightInFeatureSizeSRAM;                   // Cell height in feature size
		cell.widthInFeatureSize = param->widthInFeatureSizeSRAM;                     // Cell width in feature size
		cell.widthSRAMCellNMOS = param->widthSRAMCellNMOS;
		cell.widthSRAMCellPMOS = param->widthSRAMCellPMOS;
		cell.widthAccessCMOS = param->widthAccessCMOS;
		cell.minSenseVoltage = param->minSenseVoltage;
	} else {
		cell.heightInFeatureSize = (cell.accessType==CMOS_access)? param->heightInFeatureSize1T1R : param->heightInFeatureSizeCrossbar;         // Cell height in feature size
		cell.widthInFeatureSize = (cell.accessType==CMOS_access)? param->widthInFeatureSize1T1R : param->widthInFeatureSizeCrossbar;            // Cell width in feature size
	} 

	subArray->XNORparallelMode = param->XNORparallelMode;               
	subArray->XNORsequentialMode = param->XNORsequentialMode;             
	subArray->BNNparallelMode = param->BNNparallelMode;                
	subArray->BNNsequentialMode = param->BNNsequentialMode;              
	subArray->conventionalParallel = param->conventionalParallel;                  
	subArray->conventionalSequential = param->conventionalSequential;                 
	subArray->numRow = param->numRowSubArray;
	subArray->numCol = param->numRowSubArray;
	subArray->levelOutput = param->levelOutput;
	subArray->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
    subArray->clkFreq = param->clkFreq;                       // Clock frequency
	subArray->relaxArrayCellHeight = param->relaxArrayCellHeight;
	subArray->relaxArrayCellWidth = param->relaxArrayCellWidth;
	subArray->numReadPulse = param->numBitInput;
	subArray->avgWeightBit = param->cellBit;
	subArray->numCellPerSynapse = param->numColPerSynapse;
	subArray->SARADC = param->SARADC;
	subArray->currentMode = param->currentMode;
	subArray->validated = param->validated;
	subArray->spikingMode = NONSPIKING;
	// Anni update
	subArray->numRowParallel = param->numRowParallel;
	subArray->numAdd = ceil(param->numRowSubArray/param->numRowParallel);
	
	int numRow = param->numRowSubArray;
	int numCol = param->numColSubArray;
	
	if (subArray->numColMuxed > numCol) {                      // Set the upperbound of numColMuxed
		subArray->numColMuxed = numCol;
	}

	subArray->numReadCellPerOperationFPGA = numCol;	           // Not relevant for IMEC
	subArray->numWriteCellPerOperationFPGA = numCol;	       // Not relevant for IMEC
	subArray->numReadCellPerOperationMemory = numCol;          // Define # of SRAM read cells in memory mode because SRAM does not have S/A sharing (not relevant for IMEC)
	subArray->numWriteCellPerOperationMemory = numCol/8;       // # of write cells per operation in SRAM memory or the memory mode of multifunctional memory (not relevant for IMEC)
	subArray->numReadCellPerOperationNeuro = numCol;           // # of SRAM read cells in neuromorphic mode
	subArray->numWriteCellPerOperationNeuro = numCol;	       // For SRAM or analog RRAM in neuro mode
    subArray->maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);
	
	int numSubArrayRowCM = _numSubArrayRowCM;
	int numSubArrayColCM = _numSubArrayColCM;
	/*** initialize modules ***/
	subArray->Initialize(numRow, numCol, param->unitLengthWireResistance);        // initialize subArray
	subArray->CalculateArea();
	// Anni update: numBitSubarrayOutput
	int numBitSubarrayOutput = 0;
	if (param->parallelRead) {		
		numBitSubarrayOutput = log2((double)param->levelOutput)+ceil(log2(ceil(param->numRowSubArray/param->numRowParallel)))+param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	} 
	else{
		numBitSubarrayOutput = ceil(log2((double)param->numRowSubArray))+param->cellBit+param->numBitInput+(param->numColPerSynapse-1)*param->cellBit+1;
	}
	
	// Anni update: numBitSubarrayOutput	
	adderTreeCM->Initialize(numSubArrayRowCM, numBitSubarrayOutput, ceil((double)numSubArrayColCM*(double)numCol/(double)param->numColMuxed), param->clkFreq);

	busInputCM->Initialize(HORIZONTAL, numSubArrayRowCM, numSubArrayColCM, 0, numRow, subArray->height, subArray->width, param->clkFreq);
	busOutputCM->Initialize(VERTICAL, numSubArrayRowCM, numSubArrayColCM, 0, numCol, subArray->height, subArray->width, param->clkFreq);
	bufferOutputCM->Initialize((numCol/param->numColMuxed)*(numBitSubarrayOutput+adderTreeCM->numStage), param->clkFreq);   
	bufferInputCM->Initialize(param->numBitInput*numRow, param->clkFreq);	
}

vector<double> ProcessingUnitCalculateArea(SubArray *subArray, int numSubArrayRow, int numSubArrayCol, double *height, double *width, double *bufferArea) {
	vector<double> areaResults;
	*height = 0;
	*width = 0;
	*bufferArea = 0;
	double area = 0;
	
	subArray->CalculateArea();
	adderTreeCM->CalculateArea(NULL, subArray->width, NONE);
	bufferInputCM->CalculateArea(numSubArrayRow*subArray->height, NULL, NONE);
	bufferOutputCM->CalculateArea(NULL, numSubArrayCol*subArray->width, NONE);

	busInputCM->CalculateArea(1, true); 
	busOutputCM->CalculateArea(1, true);

	area += subArray->area  * (numSubArrayRow*numSubArrayCol) + adderTreeCM->area + bufferInputCM->area + bufferOutputCM->area;
	
	*height = sqrt(area);
	*width = area/(*height);
	
	areaResults.push_back(area);
	areaResults.push_back(subArray->areaADC*(numSubArrayRow*numSubArrayCol));
	areaResults.push_back(subArray->areaAccum*(numSubArrayRow*numSubArrayCol)+adderTreeCM->area);
	areaResults.push_back(subArray->areaOther*(numSubArrayRow*numSubArrayCol)+ bufferInputCM->area + bufferOutputCM->area);
	areaResults.push_back(subArray->areaArray*(numSubArrayRow*numSubArrayCol));
	
	return areaResults;
}

double ProcessingUnitCalculatePerformance(SubArray *subArray, Technology& tech, MemCell& cell, const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, 
											const vector<vector<double> > &inputVector,
											int arrayDup, int numSubArrayRow, int numSubArrayCol, int weightMatrixRow,
											int weightMatrixCol, int numInVector, double *readLatency, double *readDynamicEnergy, double *writeLatency, double *writeDynamicEnergy, 
											double *leakage, double *leakageSRAMInUse, double *bufferLatency, double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy,
											double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, double *coreEnergyADC, double *coreEnergyAccum, double *coreEnergyOther, 
											bool CalculateclkFreq, double *clkPeriod, double *computationLatency) {
	
	*readLatency        = 0;
	*readDynamicEnergy  = 0;
	*writeLatency       = 0;
	*writeDynamicEnergy = 0;
	*leakage            = 0;
	
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



	if(!CalculateclkFreq) {	
		subArray->clkFreq = param->clkFreq; 
		adderTreeNM->clkFreq = param->clkFreq; 
		busInputNM->clkFreq = param->clkFreq; 
		busOutputNM->clkFreq = param->clkFreq; 
		bufferInputNM->clkFreq = param->clkFreq; 
		bufferOutputNM->clkFreq = param->clkFreq; 
		adderTreeCM->clkFreq = param->clkFreq; 
		busInputCM->clkFreq = param->clkFreq; 
		busOutputCM->clkFreq = param->clkFreq; 
		bufferInputCM->clkFreq = param->clkFreq; 
		bufferOutputCM->clkFreq = param->clkFreq; 
	}

	double subArrayReadLatency, subArrayReadDynamicEnergy, subArrayLeakage, subArrayLeakageSRAMInUse, subArrayLatencyADC, subArrayLatencyAccum, subArrayLatencyOther;

	if (arrayDup > 1) {
		// weight matrix is duplicated among subArray
		if (arrayDup < numSubArrayRow*numSubArrayCol) {
			// a couple of subArrays are mapped by the matrix and none full duplication
			// need to redefine the data-grab start-point
			for (int i=0; i<ceil((double) weightMatrixRow/(double) param->numRowSubArray); i++) {
				for (int j=0; j<ceil((double) weightMatrixCol/(double) param->numColSubArray); j++) {
					int numRowMatrix = min(param->numRowSubArray, weightMatrixRow-i*param->numRowSubArray);
					int numColMatrix = min(param->numColSubArray, weightMatrixCol-j*param->numColSubArray);
					
					if ((i*param->numRowSubArray < weightMatrixRow) && (j*param->numColSubArray < weightMatrixCol)) {
						// assign weight and input to specific subArray
						vector<vector<double> > subArrayMemory;
						subArrayMemory = CopySubArray(newMemory, i*param->numRowSubArray, j*param->numColSubArray, numRowMatrix, numColMatrix);

						vector<vector<double> > subArrayMemoryOld;
						subArrayMemoryOld = CopySubArray(oldMemory, i*param->numRowSubArray, j*param->numColSubArray, numRowMatrix, numColMatrix);

						vector<vector<double> > subArrayInput;
						subArrayInput = CopySubInput(inputVector, i*param->numRowSubArray, numInVector, numRowMatrix);
						
						subArrayReadLatency  = 0;
						subArrayLatencyADC   = 0;
						subArrayLatencyAccum = 0;
						subArrayLatencyOther = 0;

						double activityColWrite = 0;
						double activityRowWrite = 0;
						int numWritePulseAVG=0;
						int totalNumWritePulse = 0;
						double writeDynamicEnergyArray = 0;
						
						GetWriteEstimation(subArray, tech, cell, subArrayMemory, subArrayMemoryOld, 
							&activityColWrite, &activityRowWrite, &numWritePulseAVG, &totalNumWritePulse, &writeDynamicEnergyArray);
						
						subArray->activityColWrite = activityColWrite;
						subArray->activityRowWrite = activityRowWrite;
						subArray->numWritePulseAVG = numWritePulseAVG;
						subArray->totalNumWritePulse = totalNumWritePulse;
						subArray->writeDynamicEnergyArray = writeDynamicEnergyArray;

						for (int k=0; k<numInVector; k++) {                 // calculate single subArray through the total input vectors
							double activityRowRead = 0;
							vector<double> input; 
							input = GetInputVector(subArrayInput, k, &activityRowRead);
							subArray->activityRowRead = activityRowRead;
							
							int cellRange = pow(2, param->cellBit);
							if (param->parallelRead) {
								subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
							} else {
								subArray->levelOutput = cellRange;
							}
							
							vector<double> columnResistance;
							columnResistance = GetColumnResistance(input, subArrayMemory, cell, param->parallelRead, subArray->resCellAccess);
							
							subArray->CalculateLatency(1e20, columnResistance, CalculateclkFreq);
							if(CalculateclkFreq && (*clkPeriod < subArray->readLatency)){
								*clkPeriod = subArray->readLatency;					//clk freq is decided by the longest sensing latency
							}							
							
							if(!CalculateclkFreq){
								subArray->CalculatePower(columnResistance);
								*readDynamicEnergy += subArray->readDynamicEnergy;
								subArrayLeakage = subArray->leakage;
								// Anni update: 
								subArrayLeakageSRAMInUse = subArray->leakageSRAMInUse;

								subArrayLatencyADC += subArray->readLatencyADC;			//sensing cycle
								subArrayLatencyAccum += subArray->readLatencyAccum;		//#cycles
								subArrayReadLatency += subArray->readLatency;		//#cycles + sensing cycle
								subArrayLatencyOther += subArray->readLatencyOther;		
								
								*coreEnergyADC += subArray->readDynamicEnergyADC;
								*coreEnergyAccum += subArray->readDynamicEnergyAccum;
								*coreEnergyOther += subArray->readDynamicEnergyOther;
							}
						}
						// accumulate write latency as array need to be write sequentially (worst case)
						// limitation by on-chip buffer, write latency will be divided by numArrayWriteParallel (real case)
						*writeLatency += subArray->writeLatency;
						*writeDynamicEnergy += subArray->writeDynamicEnergy;

						adderTreeCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray), 0);
						adderTreeCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray));
						//cout << "*readLatency: " << *readLatency << endl;
						*readLatency = MAX(subArrayReadLatency + adderTreeCM->readLatency, (*readLatency));
						*readDynamicEnergy += adderTreeCM->readDynamicEnergy;
						*coreLatencyADC = MAX(subArrayLatencyADC, (*coreLatencyADC));
						*coreLatencyAccum = MAX(subArrayLatencyAccum + adderTreeCM->readLatency, (*coreLatencyAccum));
						*coreLatencyOther = MAX(subArrayLatencyOther, (*coreLatencyOther));
						*coreEnergyAccum += adderTreeCM->readDynamicEnergy;
					}
				}
			}
			*writeDynamicEnergy = *writeDynamicEnergy * arrayDup;

			//cout << "*readLatency: " << *readLatency << endl;

			// considering speedup, the latency of processing each layer is decreased
			*readLatency      = (*readLatency     )/(arrayDup);
			*coreLatencyADC   = (*coreLatencyADC  )/(arrayDup);
			*coreLatencyAccum = (*coreLatencyAccum)/(arrayDup);
			*coreLatencyOther = (*coreLatencyOther)/(arrayDup);

			*computationLatency = *readLatency;
		} 
		else {
			// assign weight and input to one subArray and full duplication
			vector<vector<double> > subArrayMemory;
			subArrayMemory = CopySubArray(newMemory, 0, 0, weightMatrixRow, weightMatrixCol);

			vector<vector<double> > subArrayMemoryOld;
			subArrayMemoryOld = CopySubArray(oldMemory, 0, 0, weightMatrixRow, weightMatrixCol);

			vector<vector<double> > subArrayInput;
			subArrayInput = CopySubInput(inputVector, 0, numInVector, weightMatrixRow);

			subArrayReadLatency  = 0;
			subArrayLatencyADC   = 0;
			subArrayLatencyAccum = 0;
			subArrayLatencyOther = 0;

			double activityColWrite = 0;
			double activityRowWrite = 0;
			int numWritePulseAVG=0;
			int totalNumWritePulse = 0;
			double writeDynamicEnergyArray = 0;
			
			GetWriteEstimation(subArray, tech, cell, subArrayMemory, subArrayMemoryOld, 
				&activityColWrite, &activityRowWrite, &numWritePulseAVG, &totalNumWritePulse, &writeDynamicEnergyArray);
			
			subArray->activityColWrite = activityColWrite;
			subArray->activityRowWrite = activityRowWrite;
			subArray->numWritePulseAVG = numWritePulseAVG;
			subArray->totalNumWritePulse = totalNumWritePulse;
			subArray->writeDynamicEnergyArray = writeDynamicEnergyArray;

			for (int k=0; k<numInVector; k++) {                 // calculate single subArray through the total input vectors
				double activityRowRead = 0;
				vector<double> input;
				input = GetInputVector(subArrayInput, k, &activityRowRead);
				subArray->activityRowRead = activityRowRead;
				int cellRange = pow(2, param->cellBit);
				
				if (param->parallelRead) {
					subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
				} 
				else {
					subArray->levelOutput = cellRange;
				}
				
				vector<double> columnResistance;
				columnResistance = GetColumnResistance(input, subArrayMemory, cell, param->parallelRead, subArray->resCellAccess);
				
				subArray->CalculateLatency(1e20, columnResistance, CalculateclkFreq);
				if(CalculateclkFreq && (*clkPeriod < subArray->readLatency)){
					*clkPeriod = subArray->readLatency;					//clk freq is decided by the longest sensing latency
				}
				
				if(!CalculateclkFreq){
					subArray->CalculatePower(columnResistance);
					*readDynamicEnergy += subArray->readDynamicEnergy;
					subArrayLeakage = subArray->leakage;
					// Anni update: 
					subArrayLeakageSRAMInUse = subArray->leakageSRAMInUse;
					
					subArrayLatencyADC += subArray->readLatencyADC;			//sensing cycle
					subArrayLatencyAccum += subArray->readLatencyAccum;		//#cycles
					subArrayReadLatency += subArray->readLatency;		//#cycles + sensing cycle
					subArrayLatencyOther += subArray->readLatencyOther;
					
					*coreEnergyADC += subArray->readDynamicEnergyADC;
					*coreEnergyAccum += subArray->readDynamicEnergyAccum;
					*coreEnergyOther += subArray->readDynamicEnergyOther;
				}
			}
			//cout << "subArrayReadLatency: " << subArrayReadLatency << endl;
			*writeLatency = subArray->writeLatency;
			*writeDynamicEnergy = subArray->writeDynamicEnergy*arrayDup;
			// do not pass adderTree 
			*readLatency      = subArrayReadLatency /(arrayDup);
			*coreLatencyADC   = subArrayLatencyADC  /(arrayDup);
			*coreLatencyAccum = subArrayLatencyAccum/(arrayDup);
			*coreLatencyOther = subArrayLatencyOther/(arrayDup);

			*computationLatency = *readLatency;
		}
	}
	else {
		// weight matrix is further partitioned inside PE (among subArray) --> no duplicated
		for (int i=0; i<numSubArrayRow/*ceil((double) weightMatrixRow/(double) param->numRowSubArray)*/; i++) {
			for (int j=0; j<numSubArrayCol/*ceil((double) weightMatrixCol/(double) param->numColSubArray)*/; j++) {
				if ((i*param->numRowSubArray < weightMatrixRow) && (j*param->numColSubArray < weightMatrixCol)) {
					int numRowMatrix = min(param->numRowSubArray, weightMatrixRow-i*param->numRowSubArray);
					int numColMatrix = min(param->numColSubArray, weightMatrixCol-j*param->numColSubArray);
					// assign weight and input to specific subArray
					vector<vector<double> > subArrayMemory;
					subArrayMemory = CopySubArray(newMemory, i*param->numRowSubArray, j*param->numColSubArray, numRowMatrix, numColMatrix);
					
					vector<vector<double> > subArrayMemoryOld;
					subArrayMemoryOld = CopySubArray(oldMemory, i*param->numRowSubArray, j*param->numColSubArray, numRowMatrix, numColMatrix);

					vector<vector<double> > subArrayInput;
					subArrayInput = CopySubInput(inputVector, i*param->numRowSubArray, numInVector, numRowMatrix);
					
					subArrayReadLatency = 0;
					subArrayLatencyADC = 0;
					subArrayLatencyAccum = 0;
					subArrayLatencyOther = 0;
					
					double activityColWrite = 0;
					double activityRowWrite = 0;
					int numWritePulseAVG=0;
					int totalNumWritePulse = 0;
					double writeDynamicEnergyArray = 0;
					
					GetWriteEstimation(subArray, tech, cell, subArrayMemory, subArrayMemoryOld, 
						&activityColWrite, &activityRowWrite, &numWritePulseAVG, &totalNumWritePulse, &writeDynamicEnergyArray);
					
					subArray->activityColWrite = activityColWrite;
					subArray->activityRowWrite = activityRowWrite;
					subArray->numWritePulseAVG = numWritePulseAVG;
					subArray->totalNumWritePulse = totalNumWritePulse;
					subArray->writeDynamicEnergyArray = writeDynamicEnergyArray;

					for (int k=0; k<numInVector; k++) {                 // calculate single subArray through the total input vectors
						double activityRowRead = 0;

						vector<double> input;
						input = GetInputVector(subArrayInput, k, &activityRowRead);

						subArray->activityRowRead = activityRowRead;
						
						int cellRange = pow(2, param->cellBit);
						if (param->parallelRead) {
							subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
						} else {
							subArray->levelOutput = cellRange;
						}
						
						vector<double> columnResistance;
						columnResistance = GetColumnResistance(input, subArrayMemory, cell, param->parallelRead, subArray->resCellAccess);

						
						subArray->CalculateLatency(1e20, columnResistance, CalculateclkFreq);
						
						if(CalculateclkFreq && (*clkPeriod < subArray->readLatency)){
							*clkPeriod = subArray->readLatency;					//clk freq is decided by the longest sensing latency
						}
						
						if(!CalculateclkFreq){
							subArray->CalculatePower(columnResistance);
							*readDynamicEnergy += subArray->readDynamicEnergy;
							subArrayLeakage = subArray->leakage;
							// Anni update: 
							subArrayLeakageSRAMInUse = subArray->leakageSRAMInUse;
							
							subArrayLatencyADC += subArray->readLatencyADC;			//sensing cycle
							subArrayLatencyAccum += subArray->readLatencyAccum;		//#cycles
							subArrayReadLatency += subArray->readLatency;		//#cycles + sensing cycle
							subArrayLatencyOther += subArray->readLatencyOther;
							
							*coreEnergyADC += subArray->readDynamicEnergyADC;
							*coreEnergyAccum += subArray->readDynamicEnergyAccum;
							*coreEnergyOther += subArray->readDynamicEnergyOther;
						}
					}
					*writeLatency       += subArray->writeLatency      ;
					*writeDynamicEnergy += subArray->writeDynamicEnergy;

					*readLatency      = MAX(subArrayReadLatency, (*readLatency))      ;
					*coreLatencyADC   = MAX(subArrayLatencyADC, (*coreLatencyADC))    ;
					*coreLatencyAccum = MAX(subArrayLatencyAccum, (*coreLatencyAccum));
					*coreLatencyOther = MAX(subArrayLatencyOther, (*coreLatencyOther));
				}
			}
		}
		adderTreeCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray), 0);
		adderTreeCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray));

		*readLatency += adderTreeCM->readLatency;
		*coreLatencyAccum += adderTreeCM->readLatency;
		*readDynamicEnergy += adderTreeCM->readDynamicEnergy;
		*coreEnergyAccum += adderTreeCM->readDynamicEnergy;
		//cout << "*readLatency: " << *readLatency << endl;
		*computationLatency = *readLatency;
	}

	if(!CalculateclkFreq){
		// considering buffer activation: no matter speedup or not, the total number of data transferred is fixed
		// input buffer: total num of data loaded in = weightMatrixRow*numInVector
		// output buffer: total num of data transferred = weightMatrixRow*numInVector/param->numBitInput (total num of IFM in the PE) *adderTree->numAdderTree*adderTree->numAdderBit (bit precision of OFMs) 
		
		bufferInputCM->CalculateLatency(0, weightMatrixRow/param->numRowPerSynapse*numInVector/(bufferInputCM->numDff));
		bufferOutputCM->CalculateLatency(0, weightMatrixCol/param->numColPerSynapse*adderTreeCM->numAdderBit*numInVector/param->numBitInput/(bufferOutputCM->numDff));
		bufferInputCM->CalculatePower(weightMatrixRow/param->numRowPerSynapse*numInVector/(bufferInputCM->numDff), bufferInputCM->numDff, false);
		bufferOutputCM->CalculatePower(weightMatrixCol/param->numColPerSynapse*adderTreeCM->numAdderBit*numInVector/param->numBitInput/(bufferOutputCM->numDff), bufferOutputCM->numDff, false);
		
		busInputCM->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth)); 
		busInputCM->CalculatePower(busInputCM->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth), true);
		
		// Anni update: adderTreeCM->numStage+adderTreeCM->numAdderBit
		busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*(adderTreeCM->numStage+adderTreeCM->numAdderBit)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
		busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*(adderTreeCM->numStage+adderTreeCM->numAdderBit)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth), false);

		// 230920 update
		
		*bufferLatency = bufferInputCM->readLatency + bufferOutputCM->readLatency;	//considered in ic
		
		*icLatency = busInputCM->readLatency + busOutputCM->readLatency;	

		*bufferDynamicEnergy += bufferInputCM->readDynamicEnergy + bufferOutputCM->readDynamicEnergy;
		*icDynamicEnergy += busInputCM->readDynamicEnergy + busOutputCM->readDynamicEnergy;

		// 1.4 update : leakage energy of IC
		*leakage = subArrayLeakage*numSubArrayRow*numSubArrayCol + adderTreeCM->leakage + bufferInputCM->leakage + bufferOutputCM->leakage + busInputCM->leakage + busOutputCM->leakage;
		// Anni update
		*leakageSRAMInUse = subArrayLeakageSRAMInUse*numSubArrayRow*numSubArrayCol;

		// 230920 update
		*readLatency += (*bufferLatency) + (*icLatency);
		*coreLatencyOther += (*bufferLatency) + (*icLatency);
		*readDynamicEnergy += (*bufferDynamicEnergy) + (*icDynamicEnergy);
		*coreEnergyOther += (*bufferDynamicEnergy) + (*icDynamicEnergy);		
	}
	return 0;
}

vector<vector<double> > CopySubArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol) {
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


vector<vector<double> > CopySubInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow) {
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


vector<double> GetInputVector(const vector<vector<double> > &input, int numInput, double *activityRowRead) {
	vector<double> copy;
	for (int i=0; i<input.size(); i++) {
		double x = input[i][numInput];
		copy.push_back(x);   
	}  
	double numofreadrow = 0;  // initialize readrowactivity parameters
	for (int i=0; i<input.size(); i++) {
		if (copy[i] != 0) {
			numofreadrow += 1;
		}else {
			numofreadrow += 0;
		}
	}
	double totalnumRow = input.size();
	*(activityRowRead) = numofreadrow/totalnumRow;
	return copy;
	copy.clear();
} 


vector<double> GetColumnResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess) {
	vector<double> resistance;
	vector<double> conductance;
	double columnG = 0; 
	
	for (int j=0; j<weight[0].size(); j++) {
		int activatedRow = 0;
		columnG = 0;
		for (int i=0; i<weight.size(); i++) {
			if (cell.memCellType == Type::RRAM) {	// eNVM
				double totalWireResistance;
				if (cell.accessType == CMOS_access) {
					totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol + cell.resistanceAccess;
				} else {
					totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol;
				}
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
			} else if (cell.memCellType == Type::FeFET) {
				double totalWireResistance;
				totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol;
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
				
			} else if (cell.memCellType == Type::SRAM) {	
				// SRAM: weight value do not affect sense energy --> read energy calculated in subArray.cpp (based on wireRes wireCap etc)
				double totalWireResistance = (double) (resCellAccess + param->wireResistanceCol);
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
			}
		}
		
		// Anni update TODO: SRAM col resistance?
		if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (!parallelRead) {  
				conductance.push_back((double) columnG/activatedRow);
			} else {
				// Anni update
				int numAdd = ceil(param->numRowSubArray / param->numRowParallel);
				conductance.push_back(columnG/(double) numAdd);
			}
		} else {
			conductance.push_back(columnG);
		}
	}
	// covert conductance to resistance
	for (int i=0; i<weight[0].size(); i++) {
		resistance.push_back((double) 1.0/conductance[i]);
	}
		
	return resistance;
}


double GetWriteEstimation(SubArray *subArray, Technology& tech, MemCell& cell, const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, 
							double *activityColWrite, double *activityRowWrite, int *numWritePulseAVG, int *totalNumWritePulse, double *writeDynamicEnergyArray) {

	int maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);
	double minDeltaConductance = (double) (param->maxConductance-param->minConductance)/maxNumWritePulse;     // define the min delta weight
	int totalNumSetWritePulse = 0;
	int totalNumResetWritePulse = 0;
	
	*activityColWrite = 0;
	*activityRowWrite = 0;
	*numWritePulseAVG = 0;
	*totalNumWritePulse = 0;
	*writeDynamicEnergyArray = 0;
	
	int numSelectedRowSet = 0;							// used to calculate activityRowWrite
	int numSelectedRowReset = 0;						// used to calculate activityRowWrite
	int numSelectedColSet = 0;							// used to calculate activityColWrite
	int numSelectedColReset = 0;						// used to calculate activityColWrite
	for (int i=0; i<newMemory.size(); i++) {    		// update weight row-by-row
		int numSet = 0;          						// num of columns need to be set
		int numReset = 0;        						// num of columns need to be reset
		int numSetWritePulse = 0;						// num of set pulse of each row
		int numResetWritePulse = 0;						// num of reset pulse of each row
		bool rowSelected = false;
		
		// default: oldMemory[i][j] = 0, every element in newMemory needs to rewirte. (the worst case)
		for (int j=0; j<newMemory[0].size(); j++) {   	// sweep column for a row
			if (param->memcelltype != 1) { // eNVM
				if (abs(newMemory[i][j]-0) >= minDeltaConductance) {
					rowSelected = true;
					if (newMemory[i][j] > 0) {  // LTP
						numSet += 1;
						int thisPulse = (int)ceil(abs(newMemory[i][j]-oldMemory[i][j])/minDeltaConductance);
						numSetWritePulse = MAX( numSetWritePulse, thisPulse );
						// energy in each cell
						*writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / (abs(1/newMemory[i][j] + 1/oldMemory[i][j])/2) * cell.writePulseWidth * thisPulse *((cell.memCellType == Type::FeFET)==true? 0:1);
					} 
					else {   // LTD
						numReset += 1;
						int thisPulse = (int)ceil(abs(newMemory[i][j]-oldMemory[i][j])/minDeltaConductance);
						numResetWritePulse = MAX( numResetWritePulse, thisPulse );
						// energy in each cell
						*writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / (abs(1/newMemory[i][j] + 1/oldMemory[i][j])/2) * cell.writePulseWidth * thisPulse *((cell.memCellType == Type::FeFET)==true? 0:1);
					}
					if (cell.memCellType == Type::FeFET) { //FeFET
						double newPr = (newMemory[i][j]/minDeltaConductance-maxNumWritePulse/2)*(param->polarization*2/maxNumWritePulse);
						double oldPr = (0/minDeltaConductance-maxNumWritePulse/2)*(param->polarization*2/maxNumWritePulse);
						// assume pr and conductance are linear mapped
						double deltaPr = abs(newPr+(param->polarization))+abs(oldPr+(param->polarization));  // uC/cm^2 (assume erase before program)
						*writeDynamicEnergyArray += deltaPr*0.01*cell.writeVoltage*(2*tech.featureSize*tech.featureSize);
					}
				} 
				else { // no update
					numSet += 0;
					numReset += 0;
				}
			}
			else {  // SRAM
				//cout << "newMemory[i][j]: " << newMemory[i][j] << endl;
				if (newMemory[i][j] != 0) {
					rowSelected = true;
					if (newMemory[i][j] > 0) {  // LTP
						numSet += 1;
						numSetWritePulse = 1;
					} 
					else {   // LTD
						numReset += 1;
						numResetWritePulse = 1;
					}
				} 
				else { // no update
					numSet += 0;
					numReset += 0;
				}
			}
		}

		if (rowSelected && (numSet>0)) {  			 // if set happens in this row
			numSelectedRowSet += 1;
		} 
		else if (rowSelected && (numReset>0)) { 	 // if reset happens in this row
			numSelectedRowReset += 1;
		} 
		else {
			numSelectedRowSet += 0;
			numSelectedRowReset += 0;
		}
		numSelectedColSet += numSet;
		numSelectedColReset += numReset;
		totalNumSetWritePulse += numSetWritePulse;
		totalNumResetWritePulse += numResetWritePulse;
	}
	
	// get average num of selected column for set and reset
	numSelectedColSet = numSelectedRowSet==0? 0:ceil(numSelectedColSet/numSelectedRowSet);
	numSelectedColReset = numSelectedRowReset==0? 0:ceil(numSelectedColReset/numSelectedRowReset);
		
	*totalNumWritePulse = totalNumResetWritePulse + totalNumSetWritePulse;
	*numWritePulseAVG = (*totalNumWritePulse)/(MAX(1, (numSelectedRowSet+numSelectedRowReset)/2.0));
	*activityColWrite = ((numSelectedColSet+numSelectedColReset)/2.0)/newMemory[0].size();
	*activityRowWrite = ((numSelectedRowSet+numSelectedRowReset)/2.0)/newMemory.size();
	//cout << "*activityRowWrite: " << *activityRowWrite << endl;	
	
	// calculate WL BL and SL energy
	if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		if (cell.accessType == CMOS_access) {
			if (cell.memCellType == Type::FeFET) {
				// SET
				*writeDynamicEnergyArray += subArray->capRow2 * tech.vdd * tech.vdd * totalNumSetWritePulse;	
				*writeDynamicEnergyArray += (subArray->capCol + param->gateCapFeFET * numSelectedRowSet) * cell.writeVoltage * cell.writeVoltage * numSelectedColSet * totalNumSetWritePulse;
				// RESET
				*writeDynamicEnergyArray += subArray->capRow2 * tech.vdd * tech.vdd * totalNumResetWritePulse;	
				*writeDynamicEnergyArray += (subArray->capCol + param->gateCapFeFET * numSelectedRowReset) * cell.writeVoltage * cell.writeVoltage * numSelectedColReset * totalNumResetWritePulse;
			} 
			else {
				// SET
				*writeDynamicEnergyArray += subArray->capRow2 * tech.vdd * tech.vdd * totalNumSetWritePulse;																                // Selected WL
				*writeDynamicEnergyArray += subArray->capCol * cell.writeVoltage * cell.writeVoltage * (newMemory[0].size()>=numSelectedColSet? (newMemory[0].size()-numSelectedColSet):(newMemory[0].size())) * totalNumSetWritePulse;	                    // Unselected SLs
				*writeDynamicEnergyArray += subArray->capRow1 * cell.writeVoltage * cell.writeVoltage * numSelectedColSet * totalNumSetWritePulse;											// Selected BL
				// RESET
				*writeDynamicEnergyArray += subArray->capRow2 * tech.vdd * tech.vdd * totalNumResetWritePulse;																				// Selected WL
				*writeDynamicEnergyArray += subArray->capCol * cell.writeVoltage * cell.writeVoltage * numSelectedColReset * totalNumResetWritePulse;										// Selected SLs
				*writeDynamicEnergyArray += subArray->capRow1 * cell.writeVoltage * cell.writeVoltage * (newMemory[0].size()>=numSelectedColReset? (newMemory[0].size()-numSelectedColReset):(newMemory[0].size())) * totalNumResetWritePulse;				// Unselected BL
			}
		} 
		else {
			// SET
			*writeDynamicEnergyArray += subArray->capRow1 * cell.writeVoltage * cell.writeVoltage * totalNumSetWritePulse;   																// Selected WL
			*writeDynamicEnergyArray += subArray->capRow1 * cell.writeVoltage/2 * cell.writeVoltage/2 * (newMemory.size()>=numSelectedRowSet? (newMemory.size()-numSelectedRowSet):(newMemory.size())) * (*numWritePulseAVG);  						// Unselected WLs
			*writeDynamicEnergyArray += subArray->capCol * cell.writeVoltage/2 * cell.writeVoltage/2 * (newMemory[0].size()>=numSelectedColSet? (newMemory[0].size()-numSelectedColSet):(newMemory[0].size())) * totalNumSetWritePulse; 					// Unselected BLs
			*writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 
										* cell.writePulseWidth * (newMemory[0].size()>=numSelectedColSet? (newMemory[0].size()-numSelectedColSet):(newMemory[0].size())) * totalNumSetWritePulse;    										                // Half-selected (unselected) cells on the selected row
			*writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 
										* cell.writePulseWidth * (newMemory.size()>=numSelectedRowSet? (newMemory.size()-numSelectedRowSet):(newMemory.size())) * totalNumSetWritePulse;  											                // Half-selected (unselected) cells on the selected columns
			// RESET
			*writeDynamicEnergyArray += subArray->capRow1 * cell.writeVoltage/2 * cell.writeVoltage/2 * (newMemory.size()>=numSelectedRowReset? (newMemory.size()-numSelectedRowReset):(newMemory.size())) * (*numWritePulseAVG);  					    // Unselected WLs
			*writeDynamicEnergyArray += subArray->capCol * cell.writeVoltage * cell.writeVoltage * totalNumResetWritePulse; 																	// Selected BLs
			*writeDynamicEnergyArray += subArray->capCol * cell.writeVoltage/2 * cell.writeVoltage/2 * (newMemory[0].size()>=numSelectedColReset? (newMemory[0].size()-numSelectedColReset):(newMemory[0].size())) * totalNumResetWritePulse; 					// Unselected BLs
			*writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 
										* cell.writePulseWidth * (newMemory[0].size()>=numSelectedColReset? (newMemory[0].size()-numSelectedColReset):(newMemory[0].size())) * totalNumResetWritePulse;    									                    // Half-selected (unselected) cells on the selected row
			*writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 
										* cell.writePulseWidth * (newMemory.size()>=numSelectedRowReset? (newMemory.size()-numSelectedRowReset):(newMemory.size())) * totalNumResetWritePulse;   										                // Half-selected (unselected) cells on the selected columns			
		}
	} 
	else {   // SRAM
		*writeDynamicEnergyArray = 0; // leave to subarray.cpp 
	}
} 






