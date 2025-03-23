#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Precharger.h"

using namespace std;

Precharger::Precharger(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Precharger::Initialize(int _numCol, double _resLoad, double _activityColWrite, int _numReadCellPerOperationNeuro, int _numWriteCellPerOperationNeuro) {
	if (initialized)
		cout << "[Precharger] Warning: Already initialized!" << endl;

	numCol = _numCol;
	resLoad = _resLoad;
	activityColWrite = _activityColWrite;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;
	numWriteCellPerOperationNeuro = _numWriteCellPerOperationNeuro;

	double temp = 0;
	widthPMOSBitlineEqual = MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&temp, &widthPMOSBitlineEqual, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	widthPMOSBitlinePrecharger = 6 * tech.featureSize;
	EnlargeSize(&temp, &widthPMOSBitlinePrecharger, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	initialized = true;
}

void Precharger::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		double hBitlinePrecharger, wBitlinePrecharger;
		double hBitlineEqual, wBitlineEqual;
		CalculateGateArea(INV, 1, 0, widthPMOSBitlinePrecharger, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBitlinePrecharger, &wBitlinePrecharger);
		CalculateGateArea(INV, 1, 0, widthPMOSBitlineEqual, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBitlineEqual, &wBitlineEqual);
		area = 0;
		height = 0;
		width = 0;

		// 230920 update
		double hUnit = hBitlinePrecharger*2 + hBitlineEqual ;
		double wUnit = MAX(wBitlinePrecharger, wBitlineEqual);

		if (_newWidth && _option==NONE) {
			int numRowUnit;  // Number of rows of unit
			int numUnitPerRow;
			if (_newWidth < wUnit) {
				cout << "[Precharger] Error: Precharger width is even larger than the assigned width !" << endl;
			}
			numUnitPerRow = (int)(_newWidth/wUnit);
			if (numUnitPerRow > numCol) {
				numUnitPerRow = numCol;
			}
			numRowUnit = (int)ceil((double)numCol/numUnitPerRow);
			width = _newWidth;
			height = numRowUnit * hUnit;
		} else {
			width = numCol * wUnit;
			height = hUnit;
		}
		area = height * width;
		
		// Modify layout
		newHeight = _newHeight;
		newWidth = _newWidth;
		switch (_option) {
			case MAGIC:
				MagicLayout();
				break;
			case OVERRIDE:
				OverrideLayout();
				break;
			default:    // NONE
				break;
		}
		
		// Capacitance
		capOutputBitlinePrecharger = CalculateDrainCap(widthPMOSBitlinePrecharger, PMOS, hBitlinePrecharger, tech) + CalculateDrainCap(widthPMOSBitlineEqual, PMOS, hBitlineEqual, tech);
	}
}

void Precharger::CalculateLatency(double _rampInput, double _capLoad, double numRead, double numWrite){
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		writeLatency = 0;

		rampInput = _rampInput;
		capLoad = _capLoad;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp;
		double tau;

		resPullUp = CalculateOnResistance(widthPMOSBitlinePrecharger, PMOS, inputParameter.temperature, tech);
		tau = resPullUp * (capLoad + capOutputBitlinePrecharger) + resLoad * capLoad / 2;
		gm = CalculateTransconductance(widthPMOSBitlinePrecharger, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tau, beta, 1e20, &rampOutput);
		writeLatency = readLatency;

		readLatency *= numRead;
		writeLatency *= numWrite;
	}
}

void Precharger::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		// 230925 update
		/* Leakage power */
		leakage = 2*CalculateGateLeakage(INV, 1, 0, widthPMOSBitlinePrecharger, inputParameter.temperature, tech) * tech.vdd * numCol;

		/* Dynamic energy */
		// Read
		readDynamicEnergy = capLoad * tech.vdd * tech.vdd * MIN(numReadCellPerOperationNeuro, numCol);   // BL and BL_bar
		readDynamicEnergy *= numRead;
		// Write

		// 1.4 update : update precharge write energy
		writeDynamicEnergy = capLoad * tech.vdd * tech.vdd * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite);
		writeDynamicEnergy *= numWrite;
	}
}

void Precharger::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

