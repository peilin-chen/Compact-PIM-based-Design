#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "RowDecoder.h"

using namespace std;

extern Param *param;

RowDecoder::RowDecoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit(){
	initialized = false;
}

void RowDecoder::Initialize(DecoderMode _mode, int _numAddrRow, bool _MUX, bool _parallel) {
	if (initialized)
		cout << "[Row Decoder] Warning: Already initialized!" << endl;
	
	mode = _mode;
	numAddrRow = _numAddrRow;
	MUX = _MUX;
	parallel = _parallel;

	if (parallel) {         // increase MUX Decoder by 8 times  
	    // Use 2-bit predecoding
	    // INV
	    widthInvN = 8 * MIN_NMOS_SIZE * tech.featureSize;
	    widthInvP = 8 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	    numInv = numAddrRow;	// The INV at outpur driver stage does not count here

	    // NAND2
	    widthNandN = 8 * 2 * MIN_NMOS_SIZE * tech.featureSize;
	    widthNandP = 8 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	    numNand = 4 * (int)(floor(numAddrRow/2));

	    // NOR (ceil(N/2) inputs)
	    widthNorN = 8 * MIN_NMOS_SIZE * tech.featureSize;
	    widthNorP = 8 * (int)ceil((double)numAddrRow/2) * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	    if (numAddrRow > 2)
		    numNor = pow(2, numAddrRow);
	    else
		    numNor = 0;

	    // Number of M3 for connection between NAND2 and NOR stages (if numAddrRow > 2)
	    if (numAddrRow > 2)
		    numMetalConnection = numNand + (numAddrRow%2) * 2;
	    else
		    numMetalConnection = 0;
	
	    // Output driver INV
	    widthDriverInvN = 8 * 3 * MIN_NMOS_SIZE * tech.featureSize;
	    widthDriverInvP = 8 * 3 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	}
	
	else {
	    // Use 2-bit predecoding
	    // INV
	    widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	    widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
		EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	    numInv = numAddrRow;	// The INV at outpur driver stage does not count here

	    // NAND2
	    widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	    widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
		EnlargeSize(&widthNandN, &widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	    numNand = 4 * (int)(floor(numAddrRow/2));

	    // NOR (ceil(N/2) inputs)
	    widthNorN = MIN_NMOS_SIZE * tech.featureSize;
	    widthNorP = (int)ceil((double)numAddrRow/2) * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
		EnlargeSize(&widthNorN, &widthNorP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	    if (numAddrRow > 2)
		    numNor = pow(2, numAddrRow);
	    else
		    numNor = 0;

	    // Number of M3 for connection between NAND2 and NOR stages (if numAddrRow > 2)
	    if (numAddrRow > 2)
		    numMetalConnection = numNand + (numAddrRow%2) * 2;
	    else
		    numMetalConnection = 0;
	
		// 1.4 update: final driver width 

	    // Output driver INV

		if (MUX){
	    widthDriverInvN = 2 * MIN_NMOS_SIZE * tech.featureSize * param->sizingfactor_MUX;
	    widthDriverInvP = 2 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * param->sizingfactor_MUX;
		}
		
		else {
	    widthDriverInvN = 2 * MIN_NMOS_SIZE * tech.featureSize * param->sizingfactor_WLdecoder;
	    widthDriverInvP = 2 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * param->sizingfactor_WLdecoder;			
		}
		// EnlargeSize(&widthDriverInvN, &widthDriverInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
    }

	initialized = true;
}

void RowDecoder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Row Decoder Area] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand, hNor, wNor, hDriverInv, wDriverInv;
		area = 0;
		height = 0;
		width = 0;

		// 1.4 update: Metal pitch update
		// Metal Pitch Update
		double Metal2_pitch = M2_PITCH;
		double Metal3_pitch = M3_PITCH;

		if (tech.featureSize == 14 * 1e-9){
		Metal2_pitch  *= ( (double)M2_PITCH_14nm/M2_PITCH);
		Metal3_pitch *=  ( (double)M3_PITCH_14nm/M3_PITCH);}
		else if (tech.featureSize == 10 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_10nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_10nm/M3_PITCH);}
		else if (tech.featureSize == 7 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_7nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_7nm/M3_PITCH);}
		else if (tech.featureSize == 5 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_5nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_5nm/M3_PITCH);}
		else if (tech.featureSize == 3 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_3nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_3nm/M3_PITCH);}
		else if (tech.featureSize == 2 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_2nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_2nm/M3_PITCH);}
		else if (tech.featureSize == 1 * 1e-9){
		Metal2_pitch *= ( (double)M2_PITCH_1nm /M2_PITCH);
		Metal3_pitch *= ( (double)M3_PITCH_1nm/M3_PITCH);}
		else{
		Metal2_pitch *= 1;
		Metal3_pitch *=1;}

		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// NOR (ceil(N/2) inputs)
		CalculateGateArea(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNor, &wNor);
		// Output Driver INV
		CalculateGateArea(INV, 1, widthDriverInvN, widthDriverInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hDriverInv, &wDriverInv);

		// 1.4 update: Metal pitch constant replaced with newly defined variables
		if (mode == REGULAR_ROW) {	// Connect to rows
			if (_newHeight && _option==NONE) {
				if ((hNor > _newHeight) || (hNand > _newHeight) || (hInv > _newHeight)) {
					cout << "[Row Decoder] Error: Row Decoder width is even larger than the assigned width !" << endl;
				}
				// NOR
				numColNor = 0;	// Number of columns of NOR
				numNorPerCol = (int)(_newHeight/hNor);
				if (numNorPerCol > numNor) {
					numNorPerCol = numNor;
				}
				if (numNorPerCol > 0) {	// Prevent division by 0
					numColNor = (int)ceil((double)numNor/numNorPerCol);
				}
				// NAND2
				numColNand = 0;  // Number of columns of NAND2
				numNandPerCol = (int)(_newHeight/hNand);
				if (numNandPerCol > numNand) {
					numNandPerCol = numNand;
				}
				if (numNandPerCol > 0) {	// Prevent division by 0
					numColNand = (int)ceil((double)numNand/numNandPerCol);
				}
				// INV
				numColInv = 0;  // Number of columns of INV
				numInvPerCol = (int)(_newHeight/hInv);
				if (numInvPerCol > numInv) {
					numInvPerCol = numInv;
				}
				if (numColInv > 0) {	// Prevent division by 0
					numColInv = (int)ceil((double)numInv/numInvPerCol);
				}
				
				height = _newHeight;
				width = wInv * numColInv + wNand * numColNand + Metal3_pitch * numMetalConnection * tech.featureSize + wNor * numColNor;
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					width += (wNand + wDriverInv * 2) * numColNor; // 1.4 update: fixed
				} else {    // REGULAR: 2 INV as output driver
					width += (wDriverInv * 2) * numColNor;
				}
			} else {
				height = MAX(hNor*numNor, hNand*numNand);
				width = wInv + wNand + Metal3_pitch * numMetalConnection * tech.featureSize + wNor;
				if (MUX) {	// Mux enable circuit (NAND + INV) + INV
					width += wNand + wDriverInv * 2;
				} else {	// REGULAR: 2 INV as output driver
					width += wDriverInv * 2;
				}
			}
		} else {	// mode==REGULAR_COL
			if (_newWidth && _option==NONE) {
				
				if ((wNor > _newWidth) || (wNand > _newWidth) || (wInv > _newWidth)) {
					cout << "[Row Decoder] Error: Row Decoder width is even larger than the assigned width !" << endl;
				}
				
				// NOR
				numRowNor = 0;  // Number of rows of NOR
				numNorPerRow = (int)(_newWidth/wNor);
				if (numNorPerRow > numNor) {
					numNorPerRow = numNor;
				}
				if (numNorPerRow > 0) {	// Prevent division by 0
					numRowNor = (int)ceil((double)numNor/numNorPerRow);
				}
				// NAND2
				numRowNand = 0;  // Number of rows of NAND2
				numNandPerRow = (int)(_newWidth/wNand);
				if (numNandPerRow > numNand) {
					numNandPerRow = numNand;
				}
				if (numNandPerRow > 0) {	// Prevent division by 0
					numRowNand = (int)ceil((double)numNand/numNandPerRow);
				}
				// INV
				numRowInv = 0;  // Number of rows of INV
				numInvPerRow = (int)(_newWidth/wInv);
				if (numInvPerRow > numInv) {
					numInvPerRow = numInv;
				}
				numRowInv = (int)ceil((double)numInv/numInvPerRow);

				width = _newWidth;
				height = hInv * numRowInv + hNand * numRowNand + Metal2_pitch * numMetalConnection * tech.featureSize + hNor * numRowNor;
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					height += (hNand + hDriverInv * 2) * numRowNor; // 1.4 update : fixed
				} else {    // REGULAR: 2 INV as output driver
					height += (hDriverInv * 2) * numRowNor;
				}
			} else {
				height = hInv + hNand + Metal2_pitch * numMetalConnection * tech.featureSize + hNor;
				width = MAX(wNor*numNor, wNand*numNand);
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					height += hNand + hDriverInv * 2;
				} else {    // REGULAR: 2 INV as output driver
					height += hDriverInv * 2;
				}
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
			default:	// NONE
				break;
		}
		
		// Capacitance
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// NAND2
		if (numNand) {
			CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		} else {
			capNandInput = capNandOutput = 0;
		}
		// NOR (ceil(N/2) inputs)
		if (numNor) {
			CalculateGateCapacitance(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, hNor, tech, &capNorInput, &capNorOutput);
		} else {
			capNorInput = capNorOutput = 0;
		}
		// Output Driver INV
		CalculateGateCapacitance(INV, 1, widthDriverInvN, widthDriverInvP, hDriverInv, tech, &capDriverInvInput, &capDriverInvOutput);

		
	}
}

// 1.4 update: update the arguments of the latency function

void RowDecoder::CalculateLatency(double _rampInput, double _capLoad1, double _capLoad2, double resLoad, double colnum, double numRead, double numWrite){
	if (!initialized) {
		cout << "[Row Decoder Latency] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad1 = _capLoad1;   // REGULAR: general capLoad, MUX: the NMOS Tg gates
		capLoad2 = _capLoad2;   // MUX: the PMOS Tg gates
		readLatency = 0;
		writeLatency = 0;

		double resPullDown, resPullUp;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double rampInvOutput = 1e20;
		double rampNandOutput = 1e20;
		double rampNorOutput = 1e20;

		// INV
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech);	// doesn't matter pullup/pulldown?
		if (numNand)
			tr = resPullDown * (capInvOutput + capNandInput * 2);	// one address line connects to 2 NAND inputs
		else
			tr = resPullDown * (capInvOutput + capLoad1);
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, rampInput, &rampInvOutput);
		writeLatency += horowitz(tr, beta, rampInput, &rampInvOutput);
		if (!numNand)
			rampOutput = rampInvOutput;

		// NAND2
		if (numNand) {
			resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
			if (numNor)
				tr = resPullDown * (capNandOutput + capNorInput * numNor/4);
			else
				tr = resPullDown * (capNandOutput + capLoad1);
			gm = CalculateTransconductance(widthNandN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampNandOutput);
			writeLatency += horowitz(tr, beta, rampInvOutput, &rampNandOutput);
			if (!numNor)
				rampOutput = rampNandOutput;
		}
	
		// NOR (ceil(N/2) inputs)
		if (numNor) {
			resPullUp = CalculateOnResistance(widthNorP, PMOS, inputParameter.temperature, tech) * 2;
			if (MUX)
				tr = resPullUp * (capNorOutput + capNandInput);
			else
				tr = resPullUp * (capNorOutput + capInvInput);
			gm = CalculateTransconductance(widthNorP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampNandOutput, &rampNorOutput);
			writeLatency += horowitz(tr, beta, rampNandOutput, &rampNorOutput);
			rampOutput = rampNorOutput;
		}

		// Output driver or Mux enable circuit
		if (MUX) {	// Mux enable circuit (NAND + INV) + INV

			// 1.4 update : new latency to allow driver sizing

			// 1st NAND
			resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capNandOutput + capInvInput);
			gm = CalculateTransconductance(widthNandN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampNorOutput, &rampNandOutput);
			writeLatency += horowitz(tr, beta, rampNorOutput, &rampNandOutput);
			
			// &capInvInput_final, &capInvOutput_final

			// 2nd INV
			resPullUp = CalculateOnResistance(widthDriverInvP, PMOS, inputParameter.temperature, tech);
			tr = resPullUp * (capDriverInvOutput + capDriverInvInput + capLoad1);
			gm = CalculateTransconductance(widthDriverInvP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampNandOutput, &rampInvOutput);
			writeLatency += horowitz(tr, beta, rampNandOutput, &rampInvOutput);

			// 3rd INV
			resPullDown = CalculateOnResistance(widthDriverInvN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capDriverInvOutput + capLoad2);
			gm = CalculateTransconductance(widthDriverInvN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
			writeLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
			rampOutput = rampInvOutput;



		} else {	// REGULAR: 2 INV as output driver
			// 1st INV
			resPullDown = CalculateOnResistance(widthDriverInvN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capDriverInvOutput + capDriverInvInput);
			gm = CalculateTransconductance(widthDriverInvN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampNorOutput, &rampInvOutput);
			writeLatency += horowitz(tr, beta, rampNorOutput, &rampInvOutput);
			// 2nd INV
			resPullUp = CalculateOnResistance(widthDriverInvP, PMOS, inputParameter.temperature, tech);

			// 1.4 update : resisitve load at the output is considered
			tr = resPullUp * (capDriverInvOutput + capLoad1) + (resLoad * (capLoad1 + capLoad1/colnum))/2;
			gm = CalculateTransconductance(widthDriverInvP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
			writeLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
			rampOutput = rampInvOutput;
		}
		
		readLatency *= numRead;

		writeLatency *= numWrite;
	}
}

void RowDecoder::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Row Decoder] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		// Leakage power
		// INV
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numInv;
		// NAND2
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numNand;
		// NOR (ceil(N/2) inputs)
		leakage += CalculateGateLeakage(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, inputParameter.temperature, tech) * tech.vdd * numNor;
		// Output driver or Mux enable circuit
		if (MUX) {
			leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numNor;
			leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 2 * numNor;
		} else {
			leakage += CalculateGateLeakage(INV, 1, widthDriverInvN, widthDriverInvP, inputParameter.temperature, tech) * tech.vdd * 2 * numNor;
		}

		// Read dynamic energy for both memory and neuro modes (rough calculation assuming all addr from 0 to 1)
		// INV
		readDynamicEnergy += (capInvInput + capNandInput * 2) * tech.vdd * tech.vdd * (int)floor(numAddrRow/2)*2;
		readDynamicEnergy += (capInvInput + capNorInput * numNor/2) * tech.vdd * tech.vdd * (numAddrRow - (int)floor(numAddrRow/2)*2);	// If numAddrRow is odd number
		// NAND2
		readDynamicEnergy += (capNandOutput + capNorInput * numNor/4) * tech.vdd * tech.vdd * numNand/4;	// every (NAND * 4) group has one NAND output activated
		
		// INV
		writeDynamicEnergy += (capInvInput + capNandInput * 2) * tech.vdd * tech.vdd * (int)floor(numAddrRow/2)*2;
		writeDynamicEnergy += (capInvInput + capNorInput * numNor/2) * tech.vdd * tech.vdd * (numAddrRow - (int)floor(numAddrRow/2)*2);	// If numAddrRow is odd number
		// NAND2
		writeDynamicEnergy += (capNandOutput + capNorInput * numNor/4) * tech.vdd * tech.vdd * numNand/4;	// every (NAND * 4) group has one NAND output activated
		
		// NOR (ceil(N/2) inputs)
		if (MUX)
			readDynamicEnergy += (capNorOutput + capNandInput) * tech.vdd * tech.vdd;	// one NOR output activated
		else
			readDynamicEnergy += (capNorOutput + capInvInput) * tech.vdd * tech.vdd;	// one NOR output activated
		
		// NOR (ceil(N/2) inputs)
		if (MUX)
			writeDynamicEnergy += (capNorOutput + capNandInput) * tech.vdd * tech.vdd;	// one NOR output activated
		else
			writeDynamicEnergy += (capNorOutput + capInvInput) * tech.vdd * tech.vdd;	// one NOR output activated
		
		// Output driver or Mux enable circuit
		if (MUX) {
			readDynamicEnergy += (capNandOutput + capDriverInvInput) * tech.vdd * tech.vdd;
			readDynamicEnergy += (capDriverInvOutput + capDriverInvInput) * tech.vdd * tech.vdd;
			readDynamicEnergy += capDriverInvOutput * tech.vdd * tech.vdd;
			
			writeDynamicEnergy += (capNandOutput + capDriverInvInput) * tech.vdd * tech.vdd;
			writeDynamicEnergy += (capDriverInvOutput + capDriverInvInput) * tech.vdd * tech.vdd;
			writeDynamicEnergy += capDriverInvOutput * tech.vdd * tech.vdd;
		} else {
			readDynamicEnergy += (capDriverInvInput + capDriverInvOutput) * tech.vdd * tech.vdd * 2;
			writeDynamicEnergy += (capDriverInvInput + capDriverInvOutput) * tech.vdd * tech.vdd * 2;
		}
		readDynamicEnergy *= numRead;
		writeDynamicEnergy *= numWrite;
		
		if(param->validated){
			readDynamicEnergy *= param->epsilon; 	// switching activity of control circuits, epsilon = 0.05 by default
		}
	}
}

void RowDecoder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void RowDecoder::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}
