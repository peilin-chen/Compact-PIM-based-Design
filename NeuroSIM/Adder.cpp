#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "Adder.h"

using namespace std;

extern Param *param;

Adder::Adder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Adder::Initialize(int _numBit, int _numAdder, double _clkFreq){
	if (initialized)
		cout << "[Adder] Warning: Already initialized!" << endl;
	
	numBit = _numBit;
	numAdder = _numAdder;
	clkFreq = _clkFreq;
	
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	EnlargeSize(&widthNandN, &widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	initialized = true;
}

void Adder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		double hNand, wNand;
		
		// 1.4 update: GAA special layout
		if ((tech.featureSize <= 2e-9) && param->speciallayout) {CalculateGateArea(NAND, 2, MIN_NMOS_SIZE * tech.featureSize, MIN_NMOS_SIZE * tech.featureSize, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);}
		else {CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);}
		

		area = 0;
		height = 0;
		width = 0;
		if (_newHeight && _option==NONE) {   // Adder in multiple columns given the total height
			hAdder = hNand;
			wAdder = wNand * 9 * numBit;
			
			if (hAdder > _newHeight) {
				cout << "[Adder] Error: A single adder height is even larger than the assigned height ! " << endl;
			} else {
				height = _newHeight;
				width = wAdder * hAdder * numAdder / _newHeight;
			}
		} else if (_newWidth && _option==NONE) { // Adder in multiple rows given the total width
			hAdder = hNand * numBit;
			wAdder = wNand * 9;
			
			if (wAdder > _newWidth) {
				cout << "[Adder] Error: A single adder width is even larger than the assigned width ! " << endl;
			} else {
				width = _newWidth;
				height = wAdder * hAdder * numAdder / _newWidth;
            }
		} else {    // Assume one row of adder by default
			hAdder = hNand;
			wAdder = wNand * 9 * numBit;
			width = wAdder * numAdder;
			height = hAdder;
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
		
		// 1.4 update: GAA special layout
		if ((tech.featureSize == 2e-9) && param->speciallayout) { 
			CalculateGateCapacitance_GAA(NAND, 2, MIN_NMOS_SIZE * tech.featureSize, MIN_NMOS_SIZE * tech.featureSize, hNand, tech, &capNandInput, &capNandOutput, 1.0, 22.0/15.0, 8.0/15.0); }
	    else if ((tech.featureSize == 1e-9) && param->speciallayout) {
			CalculateGateCapacitance_GAA(NAND, 2, MIN_NMOS_SIZE * tech.featureSize, MIN_NMOS_SIZE * tech.featureSize, hNand, tech, &capNandInput, &capNandOutput, 1.0, 23.0/15.0, 7.0/15.0); }
		else {
			CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		}
	}
}

void Adder::CalculateLatency(double _rampInput, double _capLoad, double numRead){
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		rampInput = _rampInput;
		capLoad = _capLoad;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double readLatencyIntermediate = 0;
		double ramp[10];
		
		ramp[0] = rampInput;

		// 1.4 update: GAA special layout
		// Calibration data pattern is A=1111111..., B=1000000... and Cin=1
		// 1st
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullDown = resPullDown *3.0/2.0;}	
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[0], &ramp[1]);
		
		// 2nd
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullUp= resPullUp *3.0/2.0;}		
		tr = resPullUp * (capNandOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[1], &ramp[2]);
		
		// 3rd
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullDown = resPullDown *3.0/2.0;}		
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[2], &ramp[3]);

		// 4th
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullUp= resPullUp *3.0/2.0;}		
		tr = resPullUp * (capNandOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[3], &ramp[4]);
		
		if (numBit > 2) {
			readLatency += readLatencyIntermediate * (numBit - 2);
		}
		
		// 5th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullDown = resPullDown *3.0/2.0;}		
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[4], &ramp[5]);

		// 6th
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullUp= resPullUp *3.0/2.0;}		
		tr = resPullUp * (capNandOutput + capNandInput);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[5], &ramp[6]);
		
		// 7th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		if ((tech.featureSize <= 2e-9) && param->speciallayout)
		{resPullDown = resPullDown *3.0/2.0;}		
		tr = resPullDown * (capNandOutput + capLoad);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[6], &ramp[7]);
		
		// Anni update
		// if (param->synchronous) {
		// 	readLatency  = ceil(readLatency*clkFreq);	//#cycles
		// }
		readLatency *= numRead;
		rampOutput = ramp[7];
	}
}

void Adder::CalculatePower(double numRead, int numAdderPerOperation) {
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		// 1.4 update - updated
		/* Leakage power */
		if ((tech.featureSize == 2e-9) && param->speciallayout) { 
				/* Leakage power */
		leakage += CalculateGateLeakage(NAND, 2, MIN_NMOS_SIZE * tech.featureSize, MIN_NMOS_SIZE * tech.featureSize, inputParameter.temperature, tech) * tech.vdd * 9 * numBit * numAdder * 2/3; }
	    else if ((tech.featureSize == 1e-9) && param->speciallayout) {
				/* Leakage power */
		leakage += CalculateGateLeakage(NAND, 2, MIN_NMOS_SIZE * tech.featureSize, MIN_NMOS_SIZE * tech.featureSize, inputParameter.temperature, tech) * tech.vdd * 9 * numBit * numAdder * 2/3;}
		else {
				/* Leakage power */
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * 9 * numBit * numAdder;
		}
		
		/* Read Dynamic energy */
		// Calibration data pattern of critical path is A=1111111..., B=1000000... and Cin=1
		// Only count 0 to 1 transition for energy
		// First stage
		readDynamicEnergy += (capNandInput * 6) * tech.vdd * tech.vdd;    // Input of 1 and 2 and Cin
        readDynamicEnergy += (capNandOutput * 2) * tech.vdd * tech.vdd;  // Output of S[0] and 5
		// Second and later stages
		readDynamicEnergy += (capNandInput * 7) * tech.vdd * tech.vdd * (numBit-1);
		readDynamicEnergy += (capNandOutput * 3) * tech.vdd * tech.vdd * (numBit-1);
		
		// Hidden transition
		// First stage
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * 2;	// #2 and #3
		readDynamicEnergy += (capNandOutput + capNandInput * 2) * tech.vdd * tech.vdd;	// #4
		readDynamicEnergy += (capNandOutput + capNandInput * 3) * tech.vdd * tech.vdd;	// #5
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd;		// #6
		// Second and later stages
		readDynamicEnergy += (capNandOutput + capNandInput * 3) * tech.vdd * tech.vdd * (numBit-1);	// # 1
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * (numBit-1);		// # 3
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * 2 * (numBit-1);		// #6 and #7
	
		readDynamicEnergy *= MIN(numAdderPerOperation, numAdder) * numRead;
		
		if(param->validated){
				readDynamicEnergy *= param->delta; 	// switching activity of adder, delta = 0.15 by default
		}
	}
}

void Adder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void Adder::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}


