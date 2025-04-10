//////////////////////////////////////////////////////////////////////////////////////////////////////////
////// This WLNewDecoderDriver is used for new BNN Row-by-Row RRAM mode, only one row activatied... //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "WLNewDecoderDriver.h"

using namespace std;

extern Param *param;

WLNewDecoderDriver::WLNewDecoderDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit(){
	initialized = false;
	invalid = false;
	//inputParameter = _inputParameter;
	//tech = _tech;
	//cell = _cell;
}

WLNewDecoderDriver::~WLNewDecoderDriver() {
	// TODO Auto-generated destructor stub
}

void WLNewDecoderDriver::Initialize(int _numWLRow) {
	if (initialized)
		cout << "[WL New Decoder Driver] Warning: Already initialized!" << endl;

	//inputCap = _inputCap;
	numWLRow = _numWLRow;

	// NAND2
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
    widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthNandN, &widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthInvN, &widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);
	
	// Transmission Gate
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthTgN, &widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

	initialized = true;
}

void WLNewDecoderDriver::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[WL New Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		double hNand, wNand, hInv, wInv, hTg, wTg;
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;
		
		// 1.4 update: new cell dimension
		if (tech.featureSize == 14 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 10 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 7 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 5 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 3 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 2 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
		else if (tech.featureSize == 1 * 1e-9)
		minCellHeight *= (MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
		else
		minCellHeight *= 1;		

		
		area = 0;
		height = 0;
		width = 0;
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, minCellHeight, tech, &hNand, &wNand);
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, minCellHeight, tech, &hInv, &wInv);
		// TG
		CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg);
		

		if (_newHeight && _option==NONE) {        // no modify area option
			int numColTg = 1;
			// Adjust TG height based on _newHeight if necessary
			double TgHeight = _newHeight / numWLRow;       // define transmission gate height to be (module height/number of WL row)
			if (TgHeight < minCellHeight) {               
				numColTg = (int)ceil(minCellHeight/TgHeight);       // Calculate the number of columns
				if (numColTg > numWLRow) {
					cout << "[WLNewDecoderDriver] Error: pass gate height is even larger than the array height" << endl;
				}
				TgHeight = _newHeight / (int)ceil((double)numWLRow/numColTg);         // (numWLRow/numColTg)is the actual row number of TG 
			}
			CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);           //calculate TG layout height and width
			
			double hUnit = MAX(MAX(hInv,hNand), hTg);        // find the largest height to be the unit height
			double wUnit = 3 * wNand + wInv + 2 * wTg;         // find the largest width to be the unit width
			int numUnitPerCol = (int)(_newHeight / hUnit);            // find how many DecoderOutput units in one column
			int numColUnit = (int)ceil((double)numWLRow / numUnitPerCol);             // how many column in total
			if (numUnitPerCol > numWLRow) {                // if all units in one single column ... 
				numUnitPerCol = numWLRow;
			}
			
			height = _newHeight;          // module height is the pre-defined height
			width = (3 * wNand + wInv + 2 * wTg) * numColUnit;            // unit width * number of unit column

		} else {              // MAGIC or OVERRIDE ...
			height = MAX(MAX(hInv,hNand), hTg) * numWLRow;
			width = 3 * wNand + wInv + 2 * wTg;
		}
		area = height * width;
	
	
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

		// Resistance
		// TG
		double resTgN, resTgP;

		// 1.4 update: change to CalculateOnResistance_normal
		resTgN = CalculateOnResistance_normal(widthTgN, NMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO;
		resTgP = CalculateOnResistance_normal(widthTgP, PMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO;
		resTg = 1/(1/resTgN + 1/resTgP);

		// Capacitance
		// NAND2
		CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
		
	}
}

void WLNewDecoderDriver::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[WL New Decoder Driver] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readLatency = writeLatency = 1e41;
	} else {
		readLatency = 0;
		writeLatency = 0;
		
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		double resPullDown;      // NAND2 pulldown resistance
		double resPullUp;     // INV pullup resistance
		double capOutput;
		double trnand;	/* NAND time constant */
		double gmnand;	/* NAND transconductance */
		double betanand;	/* for NAND horowitz calculation */
		double trinv;	/* INV time constant */
		double gminv;	/* INV transconductance */
		double betainv;	/* for INV horowitz calculation */
		double trtg;	/* TG time constant */
		//double rampNorOutput;
		
		// 1st stage: NAND2
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;      // pulldown 2 NMOS in series
		trnand = resPullDown * (capNandOutput + capInvInput);          // connect to INV
		gmnand = CalculateTransconductance(widthNandN, NMOS, tech);  
		betanand = 1 / (resPullDown * gmnand);
		readLatency += horowitz(trnand, betanand, rampInput, NULL);
		writeLatency += horowitz(trnand, betanand, rampInput, NULL);
		
		// 2ed stage: INV
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		trinv = resPullUp * (capInvOutput + 2 * capNandInput);       // connect to 2 NAND2 gate
		gminv = CalculateTransconductance(widthNandP, PMOS, tech);  
		betainv = 1 / (resPullUp * gminv);
		readLatency += horowitz(trinv, betainv, rampInput, NULL);
		writeLatency += horowitz(trinv, betainv, rampInput, NULL);
		
		// 3ed stage: NAND2
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;      
		trnand = resPullDown * (capNandOutput + capTgGateP + capTgGateN);      // connect to 2 transmission gates
		gmnand = CalculateTransconductance(widthNandN, NMOS, tech);  
		betanand = 1 / (resPullDown * gmnand);
		readLatency += horowitz(trnand, betanand, rampInput, NULL);
		writeLatency += horowitz(trnand, betanand, rampInput, NULL);
		
		// 4th stage: TG
		capOutput = 2 * capTgDrain;      
		trtg = resTg * (capOutput + capLoad) + resLoad * capLoad / 2;        // elmore delay model
		readLatency += horowitz(trtg, 0, 1e20, &rampOutput);	// get from chargeLatency in the original SubArray.cpp
		writeLatency += horowitz(trtg, 0, 1e20, &rampOutput);
		
		readLatency *= numRead;
		writeLatency *= numWrite;
	}
}


void WLNewDecoderDriver::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[WL New Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;

		// Leakage power
		// NAND2
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numWLRow * 2;
		// INV
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numWLRow * 2;
		// assuming no leakge in TG
		
		// Read dynamic energy (only one row activated)
		readDynamicEnergy += capNandInput * tech.vdd * tech.vdd;                           // NAND2 input charging ( 0 to 1 )
		readDynamicEnergy += (capInvOutput + capTgGateN) * tech.vdd * tech.vdd;            // INV output charging ( 0 to 1 )
		readDynamicEnergy += (capNandOutput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd;               // NAND2 output charging ( 0 to 1 )
		readDynamicEnergy += capTgDrain * cell.readVoltage * cell.readVoltage;         // TG gate energy
		readDynamicEnergy *= numRead;          // multiply reading operation times
		
		// Write dynamic energy (only one row activated)
		writeDynamicEnergy += capNandInput * tech.vdd * tech.vdd;
		writeDynamicEnergy += (capInvOutput + capTgGateN) * tech.vdd * tech.vdd;
		writeDynamicEnergy += (capNandOutput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd;     
		writeDynamicEnergy += capTgDrain * cell.writeVoltage * cell.writeVoltage;    
		writeDynamicEnergy *= numWrite;
	}
}

void WLNewDecoderDriver::PrintProperty(const char* str) {
	//cout << "WLNewDecoderDriver Properties:" << endl;
	FunctionUnit::PrintProperty(str);
	//cout << "Number of inverter stage: " << numStage << endl;
}

void WLNewDecoderDriver::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

