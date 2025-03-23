#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "ShiftAdd.h"

using namespace std;

extern Param *param;

ShiftAdd::ShiftAdd(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), adder(_inputParameter, _tech, _cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void ShiftAdd::Initialize(int _numUnit, int _numAdderBit, double _clkFreq, SpikingMode _spikingMode, int _numReadPulse) {
	if (initialized)
		cout << "[ShiftAdd] Warning: Already initialized!" << endl;
	
	numUnit = _numUnit;
	numAdderBit = _numAdderBit;		
	numAdder = numUnit;
	clkFreq = _clkFreq;
	spikingMode = _spikingMode;
	numReadPulse = _numReadPulse;	
	
	if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
		numDff = (numAdderBit+1 + numReadPulse-1) * numUnit;	// numAdderBit+1 because the adder output is 1 bit more than the input, and numReadPulse-1 is for shift-and-add extension (shift register)
		dff.Initialize(numDff, clkFreq);
		adder.Initialize(numAdderBit, numAdder, clkFreq);
	} else {	// SPIKING: count spikes
		numBitPerDff = pow(2, numAdderBit);
		numDff = numBitPerDff * numUnit;	// numUnit shift registers in total
		dff.Initialize(numDff, clkFreq);
	}

	/* Currently ignore INV and NAND in shift-add circuit */
	// PISO shift register (https://en.wikipedia.org/wiki/Shift_register)
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numInv = numUnit;
	// NAND2
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numNand = 3 * (numDff/numUnit-1) * numUnit;	// numDff/numUnit means the number of DFF for each shift register

	initialized = true;
}

void ShiftAdd::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand;

		// 1.4 update: GAA special layout
		if ((tech.featureSize <= 2e-9) && param->speciallayout) {		// INV
			CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
			// NAND2
			CalculateGateArea(NAND, 2, widthNandN/2.0, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		}
		else {		// INV
			CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
			// NAND2
			CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		}
		area = 0;
		height = 0;
		width = 0;
		// Adder
		if (_newWidth && _option==NONE) {
			if (spikingMode == NONSPIKING) {   // NONSPIKING: binary format
				adder.CalculateArea(NULL, _newWidth, NONE);
				dff.CalculateArea(NULL, _newWidth, NONE);
			} else {    // SPIKING: count spikes
				dff.CalculateArea(NULL, _newWidth, NONE);
			}
			// Assume the INV and NAND2 are on the same row and the total width of them is smaller than the adder or DFF
			if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
				double NEW_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;

				if (tech.featureSize == 14 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 10 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 7 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 5 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 3 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 2 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 1 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
				else
				NEW_CELL_HEIGHT *= 1;	

				height = adder.height + tech.featureSize*NEW_CELL_HEIGHT /* INV and NAND2 */ + dff.height;
				width = _newWidth;
			} else {	// SPIKING: count spikes
				double NEW_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;

				if (tech.featureSize == 14 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 10 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 7 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 5 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 3 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 2 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
				else if (tech.featureSize == 1 * 1e-9)
				NEW_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
				else
				NEW_CELL_HEIGHT *= 1;	
				
				height = tech.featureSize*NEW_CELL_HEIGHT /* INV and NAND2 */ + dff.height;
				width = _newWidth;
			}
			area = height * width;
		} else {
			if (spikingMode == NONSPIKING) {   // NONSPIKING: binary format
				adder.CalculateArea(_newHeight, NULL, NONE);
				dff.CalculateArea(_newHeight, NULL, NONE);
			} else {    // SPIKING: count spikes
				dff.CalculateArea(_newHeight, NULL, NONE);
			}
			// Assume the INV and NAND2 are on the same row and the total width of them is smaller than the adder or DFF
			if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
				height = _newHeight;
				width = adder.width + wInv + wNand /* INV and NAND2 */ + dff.width;
			} else {	// SPIKING: count spikes
				height = _newHeight;
				width = wInv + wNand /* INV and NAND2 */ + dff.width;
			}
			area = height * width;
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

void ShiftAdd::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		// Assume the delay of INV and NAND2 are negligible
		if (spikingMode == NONSPIKING) {   // NONSPIKING: binary format
			// We can shift and add the weighted sum data in the next vector pulse integration cycle
			// Thus the shift-and-add time can be partially hidden by the vector pulse integration time at the next cycle
			// But there is at least one time of shift-and-add, which is at the last vector pulse cycle		

			// Anni update	
			adder.CalculateLatency(1e20, dff.capTgDrain, 1);
			dff.CalculateLatency(1e20, 1);
			double shiftAddLatency = adder.readLatency + dff.readLatency;
			if (shiftAddLatency > cell.readPulseWidth)    // Completely hidden in the vector pulse cycle if smaller
				readLatency += (shiftAddLatency - cell.readPulseWidth) * (numRead - 1);
			readLatency += shiftAddLatency;    // At least need one time of shift-and-add
			if (param->synchronous) {
				readLatency = numRead; 	// #cycles
			}

		} else {	// SPIKING: count spikes
			// We can shift out the weighted sum data in the next vector pulse integration cycle
			// Thus the shiftout time can be partially hidden by the vector pulse integration time at the next cycle
			// But there is at least one time of shiftout, which is at the last vector pulse cycle

			// Anni update
			dff.CalculateLatency(1e20, numBitPerDff);	// Need numBitPerDff cycles to shift out the weighted sum data
			double shiftLatency = dff.readLatency;
			if (shiftLatency > cell.readPulseWidth)	// Completely hidden in the vector pulse cycle if smaller
				readLatency += (shiftLatency - cell.readPulseWidth) * (numRead - 1);
			readLatency += shiftLatency;	// At least need one time of shiftout
			if (param->synchronous) {
				readLatency = numBitPerDff * numRead;	// #cycles
			}
		}
	}
}

void ShiftAdd::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
			adder.CalculatePower(numRead, numAdder);
			dff.CalculatePower(numRead, numDff, param->validated);
			
			readDynamicEnergy += adder.readDynamicEnergy;
			readDynamicEnergy += dff.readDynamicEnergy;

			leakage += adder.leakage;
			leakage += dff.leakage;
		} else {	// SPIKING: count spikes
			dff.CalculatePower(numRead, numDff, param->validated);
			readDynamicEnergy += dff.readDynamicEnergy;
			leakage += dff.leakage;
		}
	}
}

void ShiftAdd::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void ShiftAdd::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

