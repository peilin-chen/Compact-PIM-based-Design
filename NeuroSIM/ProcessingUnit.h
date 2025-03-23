#ifndef PROCESSINGUNIT_H_
#define PROCESSINGUNIT_H_
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "SubArray.h"
 
/*** Functions ***/
//2023.12.12 delete _numSubArrayRowNM _numSubArrayColNM parameter of ProcessingUnitInitialize
void ProcessingUnitInitialize(SubArray *& subArray, InputParameter& inputParameter, Technology& tech, MemCell& cell, int _numSubArrayRowCM, int _numSubArrayColCM);
//2023.12.12 delete NMpe of ProcessingUnitCalculateArea
vector<double> ProcessingUnitCalculateArea(SubArray *subArray, int numSubArrayRow, int numSubArrayCol, double *height, double *width, double *bufferArea);
// Anni update: double *leakageSRAMInUse
//2023.12.12 delete NMpe of ProcessingUnitCalculatePerformance
double ProcessingUnitCalculatePerformance(SubArray *subArray, Technology& tech, MemCell& cell, const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, const vector<vector<double> > &inputVector, 
										int arrayDup, int numSubArrayRow, int numSubArrayCol, int weightMatrixRow, int weightMatrixCol, 
										int numInVector, double *readLatency, double *readDynamicEnergy, double *writeLatency, double *writeDynamicEnergy, double *leakage, double *leakageSRAMInUse,
										double *bufferLatency, double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy,
										double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, double *coreEnergyADC, double *coreEnergyAccum, double *coreEnergyOther, 
										bool CalculateclkFreq, double*clkPeriod, double *computationLatency);

vector<vector<double> > CopySubArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol);
vector<vector<double> > CopySubInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow);
vector<double> GetInputVector(const vector<vector<double> > &input, int numInput, double *activityRowRead);
vector<double> GetColumnResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess);
double GetWriteEstimation(SubArray *subArray, Technology& tech, MemCell& cell, const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, 
						double *activityColWrite, double *activityRowWrite, int *numWritePulseAVG, int *totalNumWritePulse, double *writeDynamicEnergyArray);

#endif /* PROCESSINGUNIT_H_ */