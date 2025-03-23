#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "MaxPooling.h"

using namespace std;

extern Param *param;

MaxPooling::MaxPooling(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), comparator(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void MaxPooling::Initialize(int _numBit, int _window, int _numMaxPooling, double _clkFreq) {    // able to assign multiple MPU to operate in parallel
	if (initialized)
		cout << "[MaxPooling] Warning: Already initialized!" << endl;
	
	numBit = _numBit;                 // # of comparing elements
	window = _window;                  // window size of max pool
	numMaxPooling = _numMaxPooling;   // # of Max Pooling Unit (MPU)
	clkFreq = _clkFreq;
	
	numComparator = 0;                // initialize the # of N-bit comparator in each MPU
	numStage = 0;                     // # of N-bit comparator stage in each MPU
	
	int n = window;
	int m = n%2;
	while (n != 0) {
		int add = n/2;
		numComparator += add;
		numStage += 1;
		n /= 2;
	}
	numComparator += m;
	numStage += m;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	// NAND
	widthNandN = 2*MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNandN, &widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	// NOR1
	widthNorN = 4*MIN_NMOS_SIZE * tech.featureSize;
	widthNorP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNorN, &widthNorP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	// NOR2 (numBit-1) inputs
	widthNorN2 = (numBit*2)*MIN_NMOS_SIZE * tech.featureSize;
	widthNorP2 = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNorN2, &widthNorP2, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	// 1-bit Comparator
	comparator.Initialize(1, 1);    // initialize single comparator

	initialized = true;
}

void MaxPooling::CalculateUnitArea(AreaModify _option) {
	if (!initialized) {
		cout << "[MaxPooling] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand, hNor, wNor, hNor2, wNor2;
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// NOR (2 inputs)
		CalculateGateArea(NOR, 2, widthNorN, widthNorP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNor, &wNor);		
		// NOR (numBit-1 inputs)
		CalculateGateArea(NOR, numBit, widthNorN2, widthNorP2, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNor2, &wNor2);		
		// 1-bit comparator
		comparator.CalculateUnitArea(_option);   // single 1-bit comparator area
		
		// for one MPU of size window and N-bit input
		areaUnit = ((comparator.areaUnit+(hInv*wInv)*5)*numBit);      // each N-bit comparator needs *numBit 1-bit comparator, each 1-bit comparator needs 4 TG and 1 INV as enable
		areaUnit += ((hInv*wInv)*6)+(hNand*wNand+hInv*wInv)+(hNor*wNor)+(hNor2*wNor2);     // each N-bit comparator needs 3*TG, 3*INV, 1*AND, 1*NOR and 1*OR
		areaUnit *= numComparator;       // each MPU need *(numComparator) N-bit omparators

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
		// NAND
		CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		// NOR1
		CalculateGateCapacitance(NOR, 2, widthNorN, widthNorP, hNor, tech, &capNorInput, &capNorOutput);
		// NOR2
		CalculateGateCapacitance(NOR, numBit, widthNorN2, widthNorP2, hNor2, tech, &capNor2Input, &capNor2Output);
		
	}
}

void MaxPooling::CalculateArea(double widthArray){
	if (!initialized) {
		cout << "[MaxPooling] Error: Require initialization first!" << endl;
	} else {
		area = 0;
		height = 0;
		width = 0;
		
		width= widthArray;
		area = areaUnit * numMaxPooling;      // able to assign multiple MPU to operate in parallel
		height = area/width;
	}
}


void MaxPooling::CalculateLatency(double _rampInput, double _capLoad, double numRead){
	if (!initialized) {
		cout << "[MaxPooling] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		readLatency = 0;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resNOR, resINV, resTG;
		double rampNOROutput, rampINVOutput, rampTGOutput;

		comparator.CalculateLatency(rampInput, capInvOutput*2, 1);    // 1-bit comparator read latency
		readLatency += comparator.readLatency * numBit/2;    // assume the comparator will go to the half way and stop
		
		// Gout pass NOR2, INV and TG
		// NOR2
		resNOR = CalculateOnResistance(widthNorN2, NMOS, inputParameter.temperature, tech) * 2;
		tr = resNOR * (capInvInput*2 + numBit*capInvOutput);
		gm = CalculateTransconductance(widthNorN2, NMOS, tech);
		beta = 1 / (resNOR * gm);
		readLatency += horowitz(tr, beta, 1e20, &rampNOROutput);
		// INV
		resINV = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resINV * (capInvInput + capNor2Output);
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resINV * gm);
		readLatency += horowitz(tr, beta, 1e20, &rampINVOutput);
		// TG
		resTG = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resTG * (capInvInput*2 + capLoad);
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resTG * gm);
		readLatency += horowitz(tr, beta, rampINVOutput, &rampTGOutput);
		
		if (param->synchronous) { 
			readLatency = ceil(readLatency*clkFreq);
		}
		
		readLatency *= numStage;
		readLatency *= numRead;
	}
}

void MaxPooling::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[MaxPooling] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		/* Leakage power */
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numComparator * numMaxPooling;
		comparator.CalculatePower(1,1);   // single 1-bit comparator read only once

		/* Read Dynamic energy */
		// AND
		readDynamicEnergy += (capNandInput + capInvOutput) * tech.vdd * tech.vdd;
		// NOR
		readDynamicEnergy += (capNorInput + capInvOutput) * tech.vdd * tech.vdd;
		// INV
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * 2;
		// NOR2
		readDynamicEnergy += (capNor2Input + capInvOutput) * tech.vdd * tech.vdd;
		
		if(param->validated){
			readDynamicEnergy *= param->epsilon; 	// switching activity of control circuits, epsilon = 0.05 by default
		}
		
		readDynamicEnergy += comparator.readDynamicEnergy * numBit/2;   // assume the comparator will go to the half way and stop

		readDynamicEnergy *= numComparator;    // need *numComparator N-bit comparator
		readDynamicEnergy *= numMaxPooling;
		readDynamicEnergy *= numRead;
	}
}

void MaxPooling::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}



