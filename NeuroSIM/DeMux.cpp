#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "DeMux.h"

using namespace std;

DeMux::DeMux(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void DeMux::Initialize(int _numInput, int numRow){
	if (initialized)
		cout << "[DeMux] Warning: Already initialized!" << endl;

	numInput = _numInput;
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
    widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	// TG
	resTg = cell.resistanceOn / numRow * IR_DROP_TOLERANCE;
	widthTgN = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech)
							* tech.featureSize / (resTg*2);
	widthTgP = CalculateOnResistance(tech.featureSize, PMOS, inputParameter.temperature, tech)
							* tech.featureSize / (resTg*2);

	initialized = true;
}

void DeMux::CalculateArea(double _newHeight, double _newWidth, AreaModify _option){
	if (!initialized) {
		cout << "[DeMux] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hTg, wTg;
		area = 0;
		height = 0;
		width = 0;
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);

		// TG
		int numTgPair = numInput;
		int numRowTgPair;

		if (_newWidth && _option==NONE) {
			numRowTgPair = 1;
			double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize; // min standard cell width
			
			// 1.4 update: new cell dimension
			if (tech.featureSize == 14 * 1e-9)
			minCellWidth  *= ( (double)CPP_14nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 10 * 1e-9)
			minCellWidth  *= ( (double)CPP_10nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 7 * 1e-9)
			minCellWidth  *= ( (double)CPP_7nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 5 * 1e-9)
			minCellWidth  *= ( (double)CPP_5nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 3 * 1e-9)
			minCellWidth  *= ( (double)CPP_3nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 2 * 1e-9)
			minCellWidth  *= ( (double)CPP_2nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else if (tech.featureSize == 1 * 1e-9)
			minCellWidth  *= ( (double)CPP_1nm/(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			else
			minCellWidth  *= 1;			
			
			if (minCellWidth > _newWidth) {
				cout << "[DeMux] Error: pass gate width is even larger than the array width" << endl;
			}

			int numTgPairPerRow = (int)(_newWidth / (minCellWidth*2));    // Get max # Tg pair per row (this is not the final # Tg pair per row because the last row may have less # Tg pair)
			numRowTgPair = (int)ceil((double)numTgPair / numTgPairPerRow); // Get min # rows based on this max # Tg pair per row
			numTgPairPerRow = (int)ceil((double)numTgPair / numRowTgPair);     // Get # Tg pair per row based on this min # rows
			int TgWidth = _newWidth / numTgPairPerRow / 2;	// Division of 2 because there are 2 Tg per pair
			int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // Get the max number of folding

			// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
			CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

			// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
			CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);
			width = _newWidth;
			height = MAX(hTg * numRowTgPair, hInv);

		} else {
			// Default (just use pass gate without folding)
			CalculatePassGateArea(widthTgN, widthTgP, tech, 1, &hTg, &wTg);
			height = MAX(hTg, hInv);
			width = (wTg * 2 * numTgPair) + wInv;
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
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}

void DeMux::CalculateLatency(double _rampInput, double numRead) {	// rampInput actually is not used
	if (!initialized) {
		cout << "[DeMux] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		double tr;  /* time constant */
		readLatency = 0;
		
		// 1.4 update: needs check capTgDrain?
		
		tr = resTg * capTgDrain;
		readLatency += 2.3 * tr;

		readLatency *= numRead;
	}
}

void DeMux::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[DeMux] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		// Leakage power
		// INV
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;
		
		// Dynamic energy for both memory and neuro modes
		// INV
		readDynamicEnergy += capInvInput * tech.vdd * tech.vdd;
		// TG gates
		readDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numInput;
		// TG drain input
		readDynamicEnergy += (capTgDrain * 2) * cell.readVoltage * cell.readVoltage * numInput;
		// TG drain output
		readDynamicEnergy += capTgDrain * cell.readVoltage * cell.readVoltage * numInput;


		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[DeMux] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void DeMux::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void DeMux::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

