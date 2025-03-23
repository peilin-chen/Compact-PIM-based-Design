///////////////////////////////////////////////////////////////////////////////////////////////////////
////// This NewMux is used for switch S/A and Vbl during reading and writing ...  Column connect //////
///////////////////////////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "NewMux.h"

using namespace std;

NewMux::NewMux(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	// TODO Auto-generated constructor stub
	initialized = false;
}

void NewMux::Initialize(int _numInput){
	if (initialized)
		cout << "[NewMux] Warning: Already initialized!" << endl;
	
	numInput = _numInput;
	
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthTgN, &widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	initialized = true;
}

void NewMux::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[NewMux] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		area = 0;
		height = 0;
		width = 0;
		if (_newWidth && _option==NONE) {
			numRowTgPair = 1;
			double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize; // min standard cell width for 1 Tg
			

			// 1.4 update: cell dimension update
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
				cout << "[NewMux] Error: NewMux width is even larger than the assigned width !" << endl;
			}

			int numTgPairPerRow = (int)(_newWidth / (minCellWidth*2));    // Get max # Tg pair per row (this is not the final # Tg pair per row because the last row may have less # Tg)
			///////////////////// numInput*3 because there are 3 Tg in each single mux //////////////////////////
			numRowTgPair = (int)ceil((double)numInput*3 / numTgPairPerRow); // Get min # rows based on this max # Tg pair per row
			numTgPairPerRow = (int)ceil((double)numInput*3 / numRowTgPair);     // Get # Tg pair per row based on this min # rows
			TgWidth = _newWidth / numTgPairPerRow / 2;	// division of 2 because there are 2 Tg in one pair
			int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // get the max number of folding

			// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
			CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

			width = _newWidth;
			height = hTg * numRowTgPair;

		} else {
			// Default (pass gate with folding=1)
			CalculatePassGateArea(widthTgN, widthTgP, tech, 1, &hTg, &wTg);
			width = wTg * 2 * numInput * 3;
			height = hTg;
		}
		
	    area = height * width;

	    // Modify layout
	    newHeight = _newHeight;
	    newWidth = _newWidth;
	    switch (_option) {
		    case MAGIC:
			    MagicLayout();       // if MAGIC, call Magiclayout() in FunctionUnit.cpp
			    break;
		    case OVERRIDE:
			    OverrideLayout();    // if OVERRIDE, call Overridelayout() in FunctionUnit.cpp
			    break;
		    default:    // NONE
			    break;
	    }

	    // Capacitance
	    // TG
	    capTgGateN = CalculateGateCap(widthTgN, tech);
	    capTgGateP = CalculateGateCap(widthTgP, tech);
	    CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}


void NewMux::CalculateLatency(double _rampInput, double _capLoad, double numRead, double numWrite) {	// For simplicity, assume shift register is ideal
	if (!initialized) {
		cout << "[NewMux] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		double tr;  /* time constant */
		readLatency = 0;

		// TG

		// 1.4 update: coefficients 0.5 neeed check
		tr = resTg * (capTgDrain + 0.5*capTgGateN + 0.5*capTgGateP + capLoad);	// Calibration: use resTg*2 (only one transistor is transmitting signal in the pass gate) may be more accurate, and include gate cap because the voltage at the source of NMOS and drain of PMOS is changing (assuming Cg = 0.5Cgs + 0.5Cgd)
		readLatency += 2.3 * tr;	// 2.3 means charging from 0% to 90%
		readLatency *= numRead;
		writeLatency = cell.writePulseWidth;     // write latency determined by write pulse width
		writeLatency *= numWrite;
	}
}

void NewMux::CalculatePower(double numRead, double numWrite, double numWritePulse, int mode_1T1R, double activityRowRead, double activityColWrite) {      
	if (!initialized) {
		cout << "[NewMux] Error: Require initialization first!" << endl;
	} else {
		
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		// Read dynamic energy
		readDynamicEnergy += capTgDrain*3 * cell.readVoltage * cell.readVoltage * numInput;    // 2 TG pass Vread to BL, total loading is 3 Tg Drain capacitance
		readDynamicEnergy += (capTgGateN + capTgGateP) * 2 * tech.vdd * tech.vdd * numInput;    // open 2 TG when selected
		readDynamicEnergy *= numRead;
		readDynamicEnergy *= activityRowRead;
		
		// Write dynamic energy (2-step write and average case half SET and half RESET)
		if (mode_1T1R) {
			// LTP
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWritePulse * numInput*activityColWrite/2;   // Selected columns, '/2' means half of the writing cells are LTP
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * (numInput - numInput*activityColWrite/2);   // Unselected columns 
			// LTD
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWritePulse * numInput*activityColWrite/2;   // Selected columns	
			writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numInput;
		}else {
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWritePulse * numInput*activityColWrite/2;   // Selected columns in LTP
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWritePulse * numInput*activityColWrite/2;   // Selected columns in LTD
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage/2 * cell.writeVoltage/2 * numInput * (1-activityColWrite);   // Total unselected columns in LTP and LTD within the 2-step write
			writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numInput;
		}
	}
}


void NewMux::PrintProperty(const char* str) {
	//cout << "NewMux Properties:" << endl;
	FunctionUnit::PrintProperty(str);
}

