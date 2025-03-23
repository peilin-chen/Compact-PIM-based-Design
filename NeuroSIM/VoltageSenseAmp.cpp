#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "VoltageSenseAmp.h"

using namespace std;

VoltageSenseAmp::VoltageSenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void VoltageSenseAmp::Initialize(int _numReadCol, double _clkFreq) {
	if (initialized)
		cout << "[VoltageSenseAmp] Warning: Already initialized!" << endl;
	
	voltageSenseDiff = 0.1;
	numReadCol = _numReadCol;
	clkFreq = _clkFreq;
    
	widthNmos = MIN_NMOS_SIZE * tech.featureSize;
	widthPmos = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void VoltageSenseAmp::CalculateUnitArea() {
	if (!initialized) {
		cout << "[VoltageSenseAmp] Error: Require initialization first!" << endl;
	} else {
		double hNmos, wNmos, hPmos, wPmos;
		
		CalculateGateArea(INV, 1, widthNmos, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNmos, &wNmos);
		CalculateGateArea(INV, 1, widthPmos, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hPmos, &wPmos);
		
		areaUnit = (hNmos * wNmos) * 6 + (hPmos * wPmos) * 5;
		
		// Resistance
		resPrecharge = CalculateOnResistance(widthPmos, PMOS, inputParameter.temperature, tech);
		// Capacitance
		CalculateGateCapacitance(INV, 1, widthNmos, 0, hNmos, tech, &capNmosGate, &capNmosDrain);
		CalculateGateCapacitance(INV, 1, widthPmos, 0, hPmos, tech, &capPmosGate, &capPmosDrain);

		capS1 = capNmosGate + capNmosDrain + capPmosDrain;
	}
}

void VoltageSenseAmp::CalculateArea(double _widthVoltageSenseAmp) {	// Just add up the area of all the components
	if (!initialized) {
		cout << "[VoltageSenseAmp] Error: Require initialization first!" << endl;
	} else {
		widthVoltageSenseAmp = _widthVoltageSenseAmp;

		double x = sqrt(areaUnit/HEIGHT_WIDTH_RATIO_LIMIT); // area = HEIGHT_WIDTH_RATIO_LIMIT * x^2
		if (widthVoltageSenseAmp > x)   // Limit W/H <= HEIGHT_WIDTH_RATIO_LIMIT
			widthVoltageSenseAmp = x;
		
		area = areaUnit * numReadCol;
		width = widthVoltageSenseAmp * numReadCol;
		height = areaUnit/widthVoltageSenseAmp;
	}
}

void VoltageSenseAmp::CalculateLatency(double capInputLoad, double numRead) {
	if (!initialized) {
		cout << "[VoltageSenseAmp] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		readLatency += 2.3 * resPrecharge * capS1 + voltageSenseDiff * (capS1 + capNmosDrain + capInputLoad) / (cell.readVoltage/cell.resMemCellOn - cell.readVoltage/cell.resMemCellOff);
		// Anni update: 
		// readLatency += 1/clkFreq * 2;	// Clock time for precharge and S/A enable
		
		readLatency *= numRead;
	}
}

void VoltageSenseAmp::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[VoltageSenseAmp] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;

		// Leakage (assume connection to the cell is floating, no leakage on the precharge side, but in S/A it's roughly like 2 NAND2)
		leakage += CalculateGateLeakage(NAND, 2, widthNmos, widthPmos, inputParameter.temperature, tech) * tech.vdd * 2;
		
		// Dynamic energy
		readDynamicEnergy = 9.845e-15 * (tech.vdd / 1.1) * (tech.vdd / 1.1);	// 65nm tech node
		readDynamicEnergy *= numReadCol;
		readDynamicEnergy *= numRead;
	}
}

void VoltageSenseAmp::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

