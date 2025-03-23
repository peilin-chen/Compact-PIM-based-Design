#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "BitShifter.h"

using namespace std;

BitShifter::BitShifter(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void BitShifter::Initialize(int _numUnit, int _numBit, double _clkFreq) {
	if (initialized)
		cout << "[BitShifter] Warning: Already initialized!" << endl;
	
	numUnit = _numUnit;
	numBit = _numBit;
	clkFreq = _clkFreq;
	numDff = numBit * numUnit;	
	
	dff.Initialize(numDff, clkFreq);
	
	initialized = true;
}

void BitShifter::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[BitShifter] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand;
		area = 0;
		height = 0;
		width = 0;
		dff.CalculateArea(NULL, NULL, NONE);
		area = dff.area;
		
		if (_newWidth && _option==NONE) {
			width = _newWidth;
			height = area/width;
		} else {
			height = _newHeight;
            width = area/height;
		}
		
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

	}
}

void BitShifter::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[BitShifter] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		dff.CalculateLatency(1e20, 1);
		readLatency += dff.readLatency; // read out parallely
		readLatency *= numRead;    
	}
}

void BitShifter::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[BitShifter] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		dff.CalculatePower(numRead, numDff, false);	
		readDynamicEnergy += dff.readDynamicEnergy;
		leakage += dff.leakage;
	}
}

void BitShifter::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

