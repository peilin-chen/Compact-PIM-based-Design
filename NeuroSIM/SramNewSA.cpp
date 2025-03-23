#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "SramNewSA.h"

using namespace std;

SramNewSA::SramNewSA(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void SramNewSA::Initialize(int _numCol, int _levelOutput) {
	if (initialized)
		cout << "[SramNewSA] Warning: Already initialized!" << endl;

	numCol = _numCol;
	levelOutput = _levelOutput;   
	
	initialized = true;
}

void SramNewSA::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SramNewSA] Error: Require initialization first!" << endl;
	} else {
		height = width = area = 0;
		
		double heightUnit = 1.072e-5;
		double widthUnit = 1.011e-5;
		double areaUnit;
		
		areaUnit = (levelOutput-1) * heightUnit * widthUnit;
		
		width = _newWidth;
		area = areaUnit * numCol;
		height = area/width;
		
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
			default:	// NONE
				break;
		}

	}
}

void SramNewSA::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[SramNewSA] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
        
		readLatency = 0.1e-9;
		
		readLatency *= numRead;
	}
}

void SramNewSA::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[SramNewSA] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		readDynamicEnergy = 0.092e-12;
		readDynamicEnergy *= (levelOutput-1);
		
		readDynamicEnergy *= numCol;
		readDynamicEnergy *= numRead;
	}
}

void SramNewSA::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

