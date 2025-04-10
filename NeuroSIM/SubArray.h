#ifndef SUBARRAY_H_
#define SUBARRAY_H_

#include <vector>
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "formula.h"
#include "FunctionUnit.h"
#include "Adder.h"
#include "RowDecoder.h"
#include "Mux.h"
#include "WLDecoderOutput.h"
#include "DFF.h"
#include "DeMux.h"
#include "Precharger.h"
#include "SenseAmp.h"
#include "DecoderDriver.h"
#include "SRAMWriteDriver.h"
#include "ReadCircuit.h"
#include "SwitchMatrix.h"
#include "ShiftAdd.h"
#include "WLNewDecoderDriver.h"
#include "NewSwitchMatrix.h"
#include "CurrentSenseAmp.h"
#include "MultilevelSenseAmp.h"
#include "MultilevelSAEncoder.h"
#include "SarADC.h"
#include "LevelShifter.h"

using namespace std;

class SubArray: public FunctionUnit {
public:
	SubArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell);
	virtual ~SubArray() {}
	InputParameter& inputParameter;
	Technology& tech;
	MemCell& cell;

	/* Functions */
	void PrintProperty();
	void Initialize(int _numRow, int _numCol, double _unitWireRes);
	void CalculateArea();
	void CalculateLatency(double _rampInput, const vector<double> &columnResistance, bool CalculateclkFreq);
	void CalculatePower(const vector<double> &columnResistance);

	/* Properties */	
	bool initialized;	   // Initialization flag
	int numRow;			   // Number of rows
	int numCol;			   // Number of columns
	
	int numColMuxed;	   // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	int numWriteColMuxed;	// // How many columns share 1 write column decoder driver (for memory or neuro mode with digital RRAM)
	int totalNumWritePulse;
	int numWritePulseAVG;
	double heightArray;
	double widthArray;
	double areaArray;
	double readDynamicEnergyArray, writeDynamicEnergyArray;
	double writeLatencyArray;
	
	double lengthRow;	// Length of rows, Unit: m
	double lengthCol;	// Length of columns, Unit: m
	double unitWireRes;	// Unit wire resistance, Unit ohm/m
	double capRow1;		// Capacitance of row1 (BL for 1T1R, WL for Cross-point), Unit: F
	double capRow2;		// Capacitance of row2 (WL for 1T1R), Unit: F
	double capCol;		// Capacitance of column, Unit: F
	double resRow;		// Row resistance, Unit: ohm
	double resCol;		// Column resistance, Unit: ohm
	double resCellAccess;	// Resistance of access device, Unit: ohm
	double capCellAccess;	// Capacitance of access device, Unit: ohm
	double colDelay;	// Column delay, Unit: s

	double activityRowWrite;	// Activity for # of rows in the write
	double activityColWrite;	// Activity for # of columns in the write
	double activityRowRead, activityColRead;		// Activity for # of rows in the read
	int numReadPulse;		// # of read pulses for the input vector
	double numWritePulse;	// Average number of write pulse
	int maxNumWritePulse;	// Max # of write pulses for the device
	int maxNumIntBit;		// Max # bits for the integrate-and-fire neuron
	
	bool neuro;				// Neuro mode
	bool neuroSimReadSimulation;
	bool multifunctional;	// Multifunctional memory
	bool conventionalSequential;	
	bool conventionalParallel;
	bool BNNsequentialMode;     
	bool BNNparallelMode;      
	bool XNORsequentialMode;      
	bool XNORparallelMode;    
	bool SARADC;                // true: use sar adc; false: use MLSA	
	bool currentMode;
	bool validated;

	int levelOutput;
	
	ReadCircuitMode readCircuitMode;
	int numWriteCellPerOperationFPGA;   // Parameter for SRAM
	int numWriteCellPerOperationMemory;
	int numWriteCellPerOperationNeuro;  // For multifunctional mode or neuro only mode
	double clkFreq;			// Clock frequency (Hz)
	int avgWeightBit;		// Average weight for each synapse (value can range from 0 to numCellPerSynapse)
	int numCellPerSynapse;
	int numReadCellPerOperationFPGA;    // Parameter for SRAM
	int numReadCellPerOperationMemory;  // Parameter for SRAM
	int numReadCellPerOperationNeuro;   // Parameter for SRAM (use numBitPerCell cells to represent one D)
	bool parallelWrite; // Parameter for crossbar RRAM in neuro mode
	bool FPGA;
	bool LUT_dynamic;   // Parameter for FPGA
	bool backToBack;    // Parameter for FPGA Mux (2 Mux triangles can combine into a Mux rectangle)
	int numLut; // Parameter for FPGA
				// numLut means the number of LUTs per row for dynamic LUT
	int numReadLutPerOperationFPGA; // parameter for FPGA
	SpikingMode spikingMode;	// NONSPIKING: input data using pulses in binary representation
								// SPIKING: input data using # of pulses
	bool shiftAddEnable;    // 0 for partition because the shift-and-add circuit will be after the last A&R stage

	bool relaxArrayCellHeight;	// true: relax the memory cell height to match the height of periperal circuit unit that connects to the row (ex: standard cell height in the last stage of row decoder) if the latter is larger
	bool relaxArrayCellWidth;	// true: relax the memory cell width to match the width of periperal circuit unit that connects to the column (ex: pass gate width in the column mux) if the latter is larger

	double areaADC, areaAccum, areaOther, readLatencyADC, readLatencyAccum, readLatencyOther, readDynamicEnergyADC, readDynamicEnergyAccum, readDynamicEnergyOther;
	
	bool trainingEstimation, parallelTrans;
	int levelOutputTrans, numRowMuxedTrans, numReadPulseTrans;
	
	// 1.4 update : parameters for buffer insertion
	double widthInvN, widthInvP;
	double wInv, hInv, drivecapin, drivecapout, sectionres, targetdriveres;
	
	// Anni update
	int numRowParallel, numAdd; 
	double leakageSRAMInUse;

	/* Circuit Modules */
	LevelShifter			 wllevelshifter;
	LevelShifter			 sllevelshifter;
	LevelShifter			 bllevelshifter;
	RowDecoder               wlDecoder;
	WLDecoderOutput          wlDecoderOutput;
	WLNewDecoderDriver       wlNewDecoderDriver;
	NewSwitchMatrix          wlNewSwitchMatrix;
	CurrentSenseAmp          rowCurrentSenseAmp;
	Mux                      mux;
	RowDecoder               muxDecoder;
	SwitchMatrix             slSwitchMatrix;
	SwitchMatrix             blSwitchMatrix;
	SwitchMatrix             wlSwitchMatrix;
	DeMux                    deMux;
	ReadCircuit              readCircuit;
	Precharger               precharger;
	SenseAmp                 senseAmp;
	DecoderDriver            wlDecoderDriver;
	SRAMWriteDriver          sramWriteDriver;
	DFF                      dff;
	Adder                    adder;
	ShiftAdd                 shiftAddInput;
	ShiftAdd                 shiftAddWeight;
	MultilevelSenseAmp       multilevelSenseAmp;
	MultilevelSAEncoder      multilevelSAEncoder;
	SarADC                   sarADC;
};

#endif /* SUBARRAY_H_ */
