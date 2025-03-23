#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "SenseAmp.h"

using namespace std;

SenseAmp::SenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void SenseAmp::Initialize(int _numCol, bool _currentSense, double _senseVoltage, double _pitchSenseAmp, double _clkFreq, int _numReadCellPerOperationNeuro) {
	if (initialized)
		cout << "[SenseAmp] Warning: Already initialized!" << endl;

	numCol = _numCol;
	currentSense = _currentSense;
	senseVoltage = _senseVoltage;
	pitchSenseAmp = _pitchSenseAmp;
	clkFreq = _clkFreq;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;

	if (pitchSenseAmp <= tech.featureSize * 6) {
		puts("[SenseAmp] Error: pitch too small, cannot do the layout");
	}

	initialized = true;
}

void SenseAmp::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		height = width = area = 0;
		if (currentSense) {
			// TODO
			puts("[SenseAmp] Do not support current sensing yet");
			exit(-1);
		} 
		double hSenseP, wSenseP, hSenseN, wSenseN, hSenseIso, wSenseIso, hSenseEn, wSenseEn;
		area = 0;
		height = 0;
		width = 0;
		// Exchange width and height as in the original code

		// 1.4 update: for < 14nm compatibility
		CalculateGateArea(INV, 1, 0, ((tech.featureSize <= 14*1e-9)? 2:1)* W_SENSE_P * tech.featureSize, pitchSenseAmp, tech, &wSenseP, &hSenseP);
		CalculateGateArea(INV, 1, 0, ((tech.featureSize <= 14*1e-9)? 2:1)* W_SENSE_ISO * tech.featureSize, pitchSenseAmp, tech, &wSenseIso, &hSenseIso);
		CalculateGateArea(INV, 1, ((tech.featureSize <= 14*1e-9)? 2:1)* W_SENSE_N * tech.featureSize, 0, pitchSenseAmp, tech, &wSenseN, &hSenseN);
		CalculateGateArea(INV, 1, ((tech.featureSize <= 14*1e-9)? 2:1)* W_SENSE_EN * tech.featureSize, 0, pitchSenseAmp, tech, &wSenseEn, &hSenseEn);
		
		// Just sum up the area of all components
		area += (wSenseP * hSenseP) * 2 + (wSenseN * hSenseN) * 2 + wSenseIso * hSenseIso + wSenseEn * hSenseEn;
		area *= numCol;
		
		if (_newWidth && _option==NONE) {
			width = _newWidth;
			height = area/width;
		} else if (_newHeight && _option==NONE) {
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
			default:	// NONE
				break;
		}

		// Capacitance
		// 1.4 update: for < 14nm compatibility
		capLoad = CalculateGateCap(((tech.featureSize <= 14*1e-9)? 2:1)*(W_SENSE_P + W_SENSE_N) * tech.featureSize, tech)
				+ CalculateDrainCap(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_N * tech.featureSize, NMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_P * tech.featureSize, PMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_ISO * tech.featureSize, PMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_MUX * tech.featureSize, NMOS, pitchSenseAmp, tech);
	}
}

void SenseAmp::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;

		/* Voltage sense amplifier */

		// 1.4 update: for < 14nm compatibility
		double gm = CalculateTransconductance(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_N * tech.featureSize, NMOS, tech)
				+ CalculateTransconductance(((tech.featureSize <= 14*1e-9)? 2:1)*W_SENSE_P * tech.featureSize, PMOS, tech);
		double tau = capLoad / gm;
		readLatency += tau * log(tech.vdd / senseVoltage);
		// readLatency += 1/clkFreq;   // Clock time for S/A enable

		readLatency *= numRead;
		//cout << "senseAmp readLatency: " << readLatency << endl;
	}
}

void SenseAmp::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		/* Voltage sense amplifier */
		// Leakage

		// 1.4 update: for < 14nm compatibility
		// 230920 update
		double idleCurrent =  CalculateGateLeakage(INV, 1, ((tech.featureSize <= 14*1e-9)? 2:1)* W_SENSE_EN * tech.featureSize, 0, inputParameter.temperature, tech);
		leakage += idleCurrent * tech.vdd * numCol;
		
		// Dynamic energy
		readDynamicEnergy += capLoad * tech.vdd * tech.vdd;
		readDynamicEnergy *= MIN(numReadCellPerOperationNeuro, numCol) * numRead;
	}
}

void SenseAmp::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

