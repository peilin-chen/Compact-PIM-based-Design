#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "SRAMWriteDriver.h"

using namespace std;

SRAMWriteDriver::SRAMWriteDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void SRAMWriteDriver::Initialize(int _numCol, double _activityColWrite, int _numWriteCellPerOperationNeuro){
	if (initialized)
		cout << "[SRAMWriteDriver] Warning: Already initialized!" << endl;
	
	numCol = _numCol;
	activityColWrite = _activityColWrite;
	numWriteCellPerOperationNeuro = _numWriteCellPerOperationNeuro;
	
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	initialized = true;
}

void SRAMWriteDriver::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv;
		// INV
		// 1.4 update: (added explanation) Assume 3 inverters are used - 2 for driving BL, 1 for driving BL bar
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		area = 0;
		height = 0;
		width = 0;
		double hUnit = hInv * 3;
		double wUnit = wInv;

		if (_newWidth && _option==NONE) {
			int numRowUnit;  // Number of rows of unit
			int numUnitPerRow;
			if (wUnit > _newWidth) {
				cout << "[SRAMWriteDriver] Error: SRAMWriteDriver width is even larger than the assigned width !" << endl;
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
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
	}
}

void SRAMWriteDriver::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numWrite){
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		writeLatency = 0;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double rampInvOutput;
		
		// 1st stage INV (Pullup)
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capInvOutput + capInvInput);
		gm = CalculateTransconductance(widthInvP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		writeLatency += horowitz(tr, beta, rampInput, &rampInvOutput);
		
		// 2nd stage INV (Pulldown)
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech);
		tr = resPullDown * (capLoad + capInvOutput) + resLoad * capLoad / 2;
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		writeLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);

		writeLatency *= numWrite;
	}
}

void SRAMWriteDriver::CalculatePower(double numWrite) {
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		/* Leakage power */
		leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 3 * numCol;

		/* Write Dynamic energy */
		// Energy consumption due to charging and discharging of BL & BL_bar already accounted for in precharger.cpp
		// Assuming the write data is 0, the energy consumption is on the output of 1st INV and the input of 2nd INV
		writeDynamicEnergy = (capInvOutput + capInvInput) * tech.vdd * tech.vdd * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite);
		writeDynamicEnergy *= numWrite;
	}
}

void SRAMWriteDriver::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

