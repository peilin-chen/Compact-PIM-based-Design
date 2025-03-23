#include <cmath>
#include <iostream>
#include "constant.h"
#include "typedef.h"
#include "formula.h"
#include "Bus.h"
#include "Param.h"

using namespace std;

extern Param *param;

Bus::Bus(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Bus::Initialize(BusMode _mode, int _numRow, int _numCol, double _delaytolerance, double _busWidth, double _unitHeight, double _unitWidth, double _clkFreq){
	if (initialized)
		cout << "[Bus] Warning: Already initialized!" << endl;
	
	mode = _mode;
	numRow = _numRow;
	numCol = _numCol;     // num of Row and Col in tile/pe level
	
	unitHeight = _unitHeight;
	unitWidth = _unitWidth;
	
	clkFreq = _clkFreq;

	delaytolerance = _delaytolerance;
	busWidth = _busWidth;
	unitLengthWireResistance = param->unitLengthWireResistance;
	unitLengthWireCap = 0.2e-15/1e-6;   // 0.2 fF/mm
	
	// define min INV resistance and capacitance to calculate repeater size
	widthMinInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthMinInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	CalculateGateArea(INV, 1, widthMinInvN, widthMinInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hMinInv, &wMinInv);
	CalculateGateCapacitance(INV, 1, widthMinInvN, widthMinInvP, hMinInv, tech, &capMinInvInput, &capMinInvOutput);
	// 1.4 update: change the formula
	double resOnRep = (CalculateOnResistance(widthMinInvN, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(widthMinInvP, PMOS, inputParameter.temperature, tech))/2;
	
	// optimal repeater design to achieve highest speed
	repeaterSize = floor(sqrt(resOnRep*unitLengthWireCap/capMinInvInput/unitLengthWireResistance));
	minDist = sqrt(2*resOnRep*(capMinInvOutput+capMinInvInput)/(unitLengthWireResistance*unitLengthWireCap));
	CalculateGateArea(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hRep, &wRep);
	CalculateGateCapacitance(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, hRep, tech, &capRepInput, &capRepOutput);
	// 1.4 update: change the formula
	resOnRep = (CalculateOnResistance(MIN_NMOS_SIZE * tech.featureSize * repeaterSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, PMOS, inputParameter.temperature, tech))/2;
	double minUnitLengthDelay = 0.7*(resOnRep*(capRepInput+capRepOutput+unitLengthWireCap*minDist)+0.54*unitLengthWireResistance*minDist*unitLengthWireCap*minDist+unitLengthWireResistance*minDist*capRepInput)/minDist;
	double maxUnitLengthEnergy = (capRepInput+capRepOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist;
	
	if (delaytolerance) {   // tradeoff: increase delay to decrease energy
		double delay = 0;
		double energy = 100;
		while ((delay<minUnitLengthDelay*(1+delaytolerance)) && (repeaterSize >= 1)) {
			repeaterSize -= 1;
			minDist *= 0.9;
			CalculateGateArea(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hRep, &wRep);
			CalculateGateCapacitance(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, hRep, tech, &capRepInput, &capRepOutput);
			// 1.4 update: change the formula
			resOnRep = (CalculateOnResistance(MIN_NMOS_SIZE * tech.featureSize * repeaterSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, PMOS, inputParameter.temperature, tech))/2;
			delay = 0.7*(resOnRep*(capRepInput+capRepOutput+unitLengthWireCap*minDist)+0.54*unitLengthWireResistance*minDist*unitLengthWireCap*minDist+unitLengthWireResistance*minDist*capRepInput)/minDist;
			energy = (capRepInput+capRepOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist;
		}
	}
	
	widthInvN = MAX(1,repeaterSize) * MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = MAX(1,repeaterSize) * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	if (mode == HORIZONTAL) {
		wireLength = unitWidth*(numCol-1);
	} else {
		wireLength = unitHeight*(numRow-1);
	}
	
	initialized = true;
}

void Bus::CalculateArea(double foldedratio, bool overLap) {
	if (!initialized) {
		cout << "[Bus] Error: Require initialization first!" << endl;
	} else {
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);

		numRepeater = ceil(wireLength/minDist);
	
		if (numRepeater > 0) {
			wireWidth += busWidth*wInv/foldedratio;  
		} else {
			wireWidth += busWidth*param->wireWidth/foldedratio;
		}
		
		if (!overLap) {
			area = numRow*wireLength*wireWidth;
		} else {
			area = 0;
		}

		// Capacitance
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
	}
}

void Bus::CalculateLatency(double numRead){
	if (!initialized) {
		cout << "[Bus] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		double resOnRep = (CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech))/2;
		unitLatencyRep = 0.7*(resOnRep*(capInvInput+capInvOutput+unitLengthWireCap*minDist)+0.54*unitLengthWireResistance*minDist*unitLengthWireCap*minDist+unitLengthWireResistance*minDist*capInvInput)/minDist;
		unitLatencyWire = 0.7*unitLengthWireResistance*minDist*unitLengthWireCap*minDist/minDist;
		
		if (numRepeater > 0) {
			readLatency += wireLength*unitLatencyRep;
		} else {
			readLatency += wireLength*unitLatencyWire;
		}
		if (param->synchronous) {
			readLatency = ceil(readLatency*clkFreq);
		}
		readLatency *= numRead;	
	}
}

void Bus::CalculatePower(double numBitAccess, double numRead, bool input) {
	if (!initialized) {
		cout << "[Bus] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		unitLengthLeakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd / minDist;
		// 230920 update
		leakage = unitLengthLeakage * wireLength * numBitAccess;

		// 1.4 update -updated interconnect energy

		if (input == true) {
			unitLengthEnergyRep = (capInvInput+capInvOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist * 0.5 * param->inputtoggle;
			unitLengthEnergyWire = (unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist * 0.5 * param->inputtoggle;
		}

		else {
			unitLengthEnergyRep = (capInvInput+capInvOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist * 0.5 * param->outputtoggle;
			unitLengthEnergyWire = (unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist * 0.5 * param->outputtoggle;
		}


		if (numRepeater > 0) {
			readDynamicEnergy += wireLength*unitLengthEnergyRep;
		} else {
			readDynamicEnergy += wireLength*unitLengthEnergyWire;
		}
		readDynamicEnergy *= numBitAccess*numRead;
	}
}

void Bus::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void Bus::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}


