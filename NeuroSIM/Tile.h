#ifndef TILE_H_
#define TILE_H_
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
 
using namespace std;

/*** Functions ***/
//2023.12.12 delete _numPENM parameter of TileInitialize
//2023.12.12 delete _peSizeNM parameter of TileInitialize
void TileInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell, double _numPECM, double _peSizeCM);
//2023.12.12 delete NMTile of TileCalculateArea
vector<double> TileCalculateArea(double numPE, double peSize, double *height, double *width);
// Anni update: add double *leakageSRAMInUse
//2023.12.12 delete novelMap of TileCalculatePerformance
void TileCalculatePerformance(const vector<vector<double> > &newMemory, const vector<vector<double> > &oldMemory, const vector<vector<double> > &inputVector, 
			double numPE, double peSize, int speedUp, int PespeedUp, int weightMatrixRow, int weightMatrixCol, int numInVector, Technology& tech, MemCell& cell, double *readLatency, double *readDynamicEnergy, 
			double *writeLatency, double *writeDynamicEnergy, double *leakage, double *leakageSRAMInUse, double *bufferLatency, double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy,
			double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, double *coreEnergyADC, double *coreEnergyAccum, double *coreEnergyOther, 
			bool CalculateclkFreq, double*clkPeriod, double *computationLatency);
		
vector<vector<double> > CopyPEArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol);
vector<vector<double> > CopyPEInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow);
	

#endif /* TILE_H_ */