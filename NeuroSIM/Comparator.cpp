#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "Comparator.h"

using namespace std;

extern Param *param;

Comparator::Comparator(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Comparator::Initialize(int _numBit, int _numComparator) {
	if (initialized)
		cout << "[Comparator] Warning: Already initialized!" << endl;
	
	numBit = _numBit;
	numComparator = _numComparator;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	// NAND2
	widthNand2N = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNand2P = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNand2N, &widthNand2P, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	// NAND3
	widthNand3N = 3 * MIN_NMOS_SIZE * tech.featureSize;
	widthNand3P = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNand3N, &widthNand3P, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	initialized = true;
}

void Comparator::CalculateUnitArea(AreaModify _option) {
	if (!initialized) {
		cout << "[Comparator] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand2, wNand2, hNand3, wNand3;
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// NAND2
		CalculateGateArea(NAND, 2, widthNand2N, widthNand2P, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand2, &wNand2);
		// NAND3
		CalculateGateArea(NAND, 3, widthNand3N, widthNand3P, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand3, &wNand3);
		// Just sum up the area of all components...
		areaUnit = ((hInv * wInv) * 4 + (hNand2 * wNand2) * 4 + (hNand3 * wNand3) * 3) * numBit;
		
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
		// NAND2
		CalculateGateCapacitance(NAND, 2, widthNand2N, widthNand2P, hNand2, tech, &capNand2Input, &capNand2Output);
		// NAND3
		CalculateGateCapacitance(NAND, 3, widthNand3N, widthNand3P, hNand3, tech, &capNand3Input, &capNand3Output);
	}
}

void Comparator::CalculateArea(double widthArray){
	if (!initialized) {
		cout << "[Comparator] Error: Require initialization first!" << endl;
	} else {
		area = 0;
		height = 0;
		width = 0;
		
		width= widthArray;
		area = areaUnit * numComparator;
		height = area/width;
		
	}
}


void Comparator::CalculateLatency(double _rampInput, double _capLoad, double numRead){
	if (!initialized) {
		cout << "[Comparator] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		readLatency = 0;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double readLatencyIntermediate = 0;
		double rampNand2Output, rampNand3Output;
		
		// Just use the delay path from Gin to Gout for simplicity
		// 1st bit comparator
		// NAND2
		resPullDown = CalculateOnResistance(widthNand2N, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNand2Output + capNand3Input);
		gm = CalculateTransconductance(widthNand2N, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, rampInput, &rampNand2Output);
		// NAND3
		resPullUp = CalculateOnResistance(widthNand3P, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNand3Output + capNand2Input);
		gm = CalculateTransconductance(widthNand3P, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, rampNand2Output, &rampNand3Output);
		
		// 2nd bit to the second last bit comparator
		// NAND2
		resPullDown = CalculateOnResistance(widthNand2N, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNand2Output + capNand3Input);
		gm = CalculateTransconductance(widthNand2N, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatencyIntermediate += horowitz(tr, beta, rampNand3Output, &rampNand2Output);
		// NAND3
		resPullUp = CalculateOnResistance(widthNand3P, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNand3Output + capNand2Input);
		gm = CalculateTransconductance(widthNand3P, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatencyIntermediate += horowitz(tr, beta, rampNand2Output, &rampNand3Output);

		readLatency += readLatencyIntermediate * (numBit - 2);

		// Last bit comparator
		// NAND2
		resPullDown = CalculateOnResistance(widthNand2N, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNand2Output + capNand3Input);
		gm = CalculateTransconductance(widthNand2N, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, rampNand3Output, &rampNand2Output);
		// NAND3
		resPullUp = CalculateOnResistance(widthNand3P, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNand3Output + capLoad);
		gm = CalculateTransconductance(widthNand3P, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, rampNand2Output, &rampNand3Output);

		readLatency *= numRead;
	}
}

void Comparator::CalculatePower(double numRead, int numComparatorPerOperation) {
	if (!initialized) {
		cout << "[Comparator] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		/* Leakage power */
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 4 * numBit * numComparator;
		leakage += CalculateGateLeakage(NAND, 2, widthNand2N, widthNand2P, inputParameter.temperature, tech) * tech.vdd * 4 * numBit * numComparator;
		leakage += CalculateGateLeakage(NAND, 3, widthNand3N, widthNand3P, inputParameter.temperature, tech) * tech.vdd * 3 * numBit * numComparator;

		/* Read Dynamic energy */
		// INV
		readDynamicEnergy += ((capInvInput + capInvOutput) * 4) * tech.vdd * tech.vdd;
		// NAND2
		readDynamicEnergy += ((capNand2Input + capNand2Output) * 4) * tech.vdd * tech.vdd;
		// NAND3
		readDynamicEnergy += ((capNand3Input + capNand3Output) * 3) * tech.vdd * tech.vdd;

		readDynamicEnergy *= numBit * MIN(numComparatorPerOperation, numComparator) * numRead;
		
		if(param->validated){
			readDynamicEnergy *= param->epsilon; 	// switching activity of control circuits, epsilon = 0.05 by default
		}
	}
}

void Comparator::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void Comparator::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}


