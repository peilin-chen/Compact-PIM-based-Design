#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Mux.h"

using namespace std;

Mux::Mux(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Mux::Initialize(int _numInput, int _numSelection, double _resTg, bool _FPGA){
	if (initialized)
		cout << "[Mux] Warning: Already initialized!" << endl;

	numInput = _numInput;
	numSelection = _numSelection;	/* Typically numColMuxed */
	FPGA = _FPGA;

	// Mux
	if (FPGA) {	// Assume digital Mux has standard NMOS and PMOS width
		widthTgN = MIN_NMOS_SIZE * tech.featureSize;
		widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
		EnlargeSize(&widthTgN, &widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
		
		// 1.4 update
		resTg = 1 / (1/CalculateOnResistance_normal(widthTgN, NMOS, inputParameter.temperature, tech) 
					+ 1/CalculateOnResistance_normal(widthTgP, PMOS, inputParameter.temperature, tech));
	} else {
		resTg = _resTg * IR_DROP_TOLERANCE;

		// 1.4 update: for compatibility with < 14 nm
		widthTgN = CalculateOnResistance_normal(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, NMOS, inputParameter.temperature, tech)
								* tech.featureSize / (resTg*2);
		widthTgP = CalculateOnResistance_normal(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, PMOS, inputParameter.temperature, tech)
								* tech.featureSize/ (resTg*2);
		
		// 1.4 update: no enlarge
		if (tech.featureSize <= 14*1e-9){
			widthTgN = 2* ceil(widthTgN/tech.featureSize) * tech.featureSize;
			widthTgP = 2* ceil(widthTgP/tech.featureSize) * tech.featureSize;
		} else {
			if (widthTgN < tech.featureSize)
				widthTgN = tech.featureSize;
			if (widthTgP < tech.featureSize)
				widthTgP = tech.featureSize;
		}
		// EnlargeSize(&widthTgN, &widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
		resTg = 1 / (1/CalculateOnResistance_normal(widthTgN, NMOS, inputParameter.temperature, tech)
					+ 1/CalculateOnResistance_normal(widthTgP, PMOS, inputParameter.temperature, tech));
	}
	initialized = true;
}

void Mux::CalculateArea(double _newHeight, double _newWidth, AreaModify _option){
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		int numTg = numInput * numSelection;
		area = 0;
		height = 0;
		width = 0;
		// TG
		if (FPGA) {	// Digital Mux
			CalculateGateArea(INV, 1, widthTgN, widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hTg, &wTg);
			if (_newWidth && _option==NONE) { // Tg in multiple rows given the total width
				// Calculate the number of Tg per row
				if (_newWidth < wTg) {
					cout << "[Mux] Error: Mux width is even larger than the assigned width !" << endl;
				} else {
					int numTgPerRow = (int)(_newWidth/wTg);
					if (numTgPerRow > numTg) {
						numTgPerRow = numTg;
					}
					numRowTg = (int)ceil((double)numTg / numTgPerRow);
					width = _newWidth;
					height = hTg * numRowTg;
				}
			} else {    // Assume one row of Tg by default
				width = wTg * numTg;
				height = hTg;
			}
		} else {	// Analog Mux
			if (_newWidth && _option==NONE) {
				numRowTg = 1;
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
					cout << "[Mux] Error: Mux width is even larger than the assigned width !" << endl;
				}

				int numTgPerRow = (int)(_newWidth / minCellWidth);	// Get max # Tg per row (this is not the final # Tg per row because the last row may have less # Tg)
				numRowTg = (int)ceil((double)numTg / numTgPerRow);	// Get min # rows based on this max # Tg per row
				numTgPerRow = (int)ceil((double)numTg / numRowTg);	// Get # Tg per row based on this min # rows
				TgWidth = _newWidth / numTgPerRow;
				int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // Get the max number of folding

				// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
				CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

				width = _newWidth;
				height = hTg * numRowTg;

			} else {
				// Default (just use pass gate without folding)
				CalculatePassGateArea(widthTgN, widthTgP, tech, 1, &hTg, &wTg);
				height = hTg;
				width = wTg * numTg;
			}
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

		widthTgShared = width/numInput;

		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}

void Mux::CalculateLatency(double _rampInput, double _capLoad, double numRead) {  // rampInput is from SL/BL, not fron EN signal
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		double resPullDown;
		double tr;  	/* time constant */
		double gm;  	/* transconductance */
		double beta;    /* for horowitz calculation */
		double rampNandOutput;
		readLatency = 0;


		// 1.4 update : needs check the 0.5 coefficients
		// TG
		tr = resTg * (capTgDrain + 0.5*capTgGateN + 0.5*capTgGateP + capLoad);	

		readLatency += 2.3 * tr;	// 2.3 means charging from 0% to 90%
		readLatency *= numRead;
	}
}

void Mux::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		// TG gates only
		readDynamicEnergy += capTgGateN * numInput * tech.vdd * tech.vdd;	// Selected pass gates (OFF to ON)
		readDynamicEnergy += (capTgDrain * 2) * numInput * cell.readVoltage * cell.readVoltage;	// Selected pass gates (OFF to ON)
		readDynamicEnergy *= numRead;
	}
}

void Mux::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void Mux::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

