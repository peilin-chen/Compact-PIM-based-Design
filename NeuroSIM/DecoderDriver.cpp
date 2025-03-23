#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "DecoderDriver.h"

using namespace std;

DecoderDriver::DecoderDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void DecoderDriver::Initialize(int _mode, int _numOutput /* # of array rows/columns */, int numLoad) {
	if (initialized)
		cout << "[Decoder Driver] Warning: Already initialized!" << endl;

	mode = _mode;
	numOutput = _numOutput;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	// TG

	// 1.4 update : Onresistance calculation update for <14 nm
	resTg = cell.resMemCellOn / numLoad * IR_DROP_TOLERANCE;
	widthTgN = CalculateOnResistance_normal( ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech)
				* tech.featureSize / (resTg*2);
	widthTgP = CalculateOnResistance_normal( ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, PMOS, inputParameter.temperature, tech)
				* tech.featureSize / (resTg*2);
	resTg = 1 / (1/CalculateOnResistance_normal(widthTgN, NMOS, inputParameter.temperature, tech)
			+ 1/CalculateOnResistance_normal(widthTgP, PMOS, inputParameter.temperature, tech));

	initialized = true;
}

void DecoderDriver::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hTg, wTg;
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;
		double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;

		// 1.4 update: new cell dimenstion
		
		if (tech.featureSize == 14 * 1e-9)
		minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 10 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 7 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 5 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 3 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 2 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
    	else if (tech.featureSize == 1 * 1e-9)
    	minCellHeight *= ( (double)MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
    	else
    	minCellHeight *= 1;

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
		
		
		
		area = 0;
		height = 0;
		width = 0;
		// TG
		if (_newHeight && _option==NONE) {
			if (_newHeight < minCellHeight) {
				cout << "[DecoderDriver] Error: pass gate height is even larger than the assigned height" << endl;
			}

			int numTgPerCol = (int)(_newHeight / minCellHeight);    // Get max # Tg per column (this is not the final # Tg per column because the last column may have less # Tg)
			numColTg = (int)ceil((double)numOutput / numTgPerCol); // Get min # columns based on this max # Tg per column
			numTgPerCol = (int)ceil((double)numOutput / numColTg);     // Get # Tg per column based on this min # columns
			TgHeight = _newHeight / numTgPerCol;
			CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);

		} else if (_newWidth && _option==NONE) {
			if (minCellWidth > _newWidth) {
				cout << "[DecoderDriver] Error: pass gate width is even larger than the assigned width" << endl;
			}

			int numTgPerRow = (int)(_newWidth / minCellWidth);    // Get max # Tg per row (this is not the final # Tg per row because the last row may have less # Tg)
			numRowTg = (int)ceil((double)numOutput / numTgPerRow); // Get min # rows based on this max # Tg per row
			numTgPerRow = (int)ceil((double)numOutput / numRowTg);     // Get # Tg per row based on this min # rows
			TgWidth = _newWidth / numTgPerRow;
			int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // get the max number of folding
			CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

		} else {
			CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg);
		}
		
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, minCellHeight, tech, &hInv, &wInv);

		double hUnit, wUnit;
		if (cell.accessType == CMOS_access) {	// 1T1R
			if (mode == ROW_MODE) {
				hUnit = MAX(hInv, hTg);
				wUnit = wInv + wTg * 3;
			} else {
				hUnit = hInv + hTg * 3;
				wUnit = MAX(wInv, wTg);
			}
		} else {	// Cross-point
			if (mode == ROW_MODE) {
				hUnit = MAX(hInv, hTg);
				wUnit = wInv + wTg * 2;
			} else {
				hUnit = hInv + hTg * 2;
				wUnit = MAX(wInv, wTg);
			}
		}

		if (mode == ROW_MODE) {	// Connect to rows
			if (_newHeight && _option==NONE) {
				int numColUnit, numUnitPerCol;
				numUnitPerCol = (int)(_newHeight/hUnit);
				numColUnit = (int)ceil((double)numOutput/numUnitPerCol);
				if (numColUnit > numOutput) {
					numColUnit = numOutput;
				}
				height = _newHeight;
				width = wUnit * numColUnit;
			} else {
				height = hUnit * numOutput;
				width = wUnit;
			}
		} else {	// Connect to columns
			if (_newWidth && _option==NONE) {
				int numRowUnit, numUnitPerRow;
				numUnitPerRow = (int)(_newWidth/wUnit);
				numRowUnit = (int)ceil((double)numOutput/numUnitPerRow);
				if (numRowUnit > numOutput) {
					numRowUnit = numOutput;
				}
				height = hUnit * numRowUnit;
				width = _newWidth;
			} else {
				height = hUnit;
				width = wUnit * numOutput;
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
		
		// Capacitance
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}

void DecoderDriver::CalculateLatency(double _rampInput, double _capLoad1, double _capLoad2, double _resLoad, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		capLoad1 = _capLoad1;   // worst-case load (1T1R SL, crosspoint WL/BL)
		capLoad2 = _capLoad2;   // 1T1R BL (which does not include transistor drain cap), doesn't matter for crosspoint
		resLoad = _resLoad;
		
		rampInput = _rampInput;
		double capOutput;
		double tr;	/* time constant */
		double gm;	/* transconductance */
		double beta;	/* for horowitz calculation */
		
		// TG
		// 1.4 update : capOutput change from capTgDrain * 2; neeeds check
		capOutput = capTgDrain*4 + capTgGateN*0.5 + capTgGateP*0.5;
		tr = resTg * (capOutput + capLoad1) + resLoad * capLoad1 / 2;
		readLatency += horowitz(tr, 0, rampInput, &rampOutput); // get from chargeLatency in the original SubArray.cpp
		readLatency *= numRead;

		//writeLatency = cell.writePulseWidth;
		writeLatency += horowitz(tr, 0, rampInput, &rampOutput);
		writeLatency *= numWrite;
	}
}

void DecoderDriver::CalculatePower(double numReadCellPerOp, double numWriteCellPerOp, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;

		// Leakage power
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numOutput;
		
		// Read dynamic energy
		if (cell.accessType == CMOS_access) {  // 1T1R
			// Selected SLs and BLs are floating
			// Unselected SLs and BLs are GND
			readDynamicEnergy += (capInvInput + capTgGateN * 2 + capTgGateP) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capInvOutput + capTgGateP * 2 + capTgGateN) * tech.vdd * tech.vdd * numReadCellPerOp;
		} else {	// Crosspoint
			// For WL decoder driver, the selected WLs are GND
			// For BL decoder driver, the selected BLs are floating
			// No matter which one, the unselected WLs/BLs are read voltage
			readDynamicEnergy += (capInvInput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capInvOutput + capTgGateP + capTgGateN) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capTgDrain * 2) * cell.readVoltage * cell.readVoltage * (numOutput-numReadCellPerOp);
		}
		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[Decoder Driver] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy
		if (cell.accessType == CMOS_access) {  // 1T1R
			// Worst case: RESET operation (because SL cap is larger than BL cap)
			writeDynamicEnergy += (capInvInput + capTgGateN * 2 + capTgGateP) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capInvOutput + capTgGateP * 2 + capTgGateN) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWriteCellPerOp;
		} else {    // Crosspoint
			writeDynamicEnergy += (capInvInput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capInvOutput + capTgGateP + capTgGateN) * tech.vdd * tech.vdd * numWriteCellPerOp;
			if (mode == ROW_MODE) {	// Connects to rows
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWriteCellPerOp;
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage/2 * cell.writeVoltage/2 * (numOutput-numWriteCellPerOp);
			} else {	// Connects to columns
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage/2 * cell.writeVoltage/2 * (numOutput-numWriteCellPerOp);
			}
		}
		writeDynamicEnergy *= numWrite;
		if (!writeLatency) {
			//cout << "[Decoder Driver] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency;
		}

	}
}

void DecoderDriver::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void DecoderDriver::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

