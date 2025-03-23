#ifndef PARAM_H_
#define PARAM_H_

class Param {
public:
	Param();

	int operationmode, operationmodeBack, memcelltype, accesstype, transistortype, deviceroadmap;      		
	
	double heightInFeatureSizeSRAM, widthInFeatureSizeSRAM, widthSRAMCellNMOS, widthSRAMCellPMOS, widthAccessCMOS, minSenseVoltage;
	 
	double heightInFeatureSize1T1R, widthInFeatureSize1T1R, heightInFeatureSizeCrossbar, widthInFeatureSizeCrossbar;
	
	int relaxArrayCellHeight, relaxArrayCellWidth;

	bool globalBusType, globalBufferType, tileBufferType, peBufferType, reLu, SARADC, currentMode, validated, synchronous;
	int globalBufferCoreSizeRow, globalBufferCoreSizeCol, tileBufferCoreSizeRow, tileBufferCoreSizeCol;																								
	
	double clkFreq, featuresize, readNoise, resistanceOn, resistanceOff, maxConductance, minConductance, gateCapFeFET, polarization;
	int temp, technode, wireWidth, multipleCells;
	double maxNumLevelLTP, maxNumLevelLTD, readVoltage, readPulseWidth, writeVoltage;
	double accessVoltage, resistanceAccess;
	double nonlinearIV, nonlinearity;
	double writePulseWidth, numWritePulse;
	double globalBusDelayTolerance, localBusDelayTolerance;
	double treeFoldedRatio, maxGlobalBusWidth;
	double algoWeightMax, algoWeightMin;
	
	int neuro, multifunctional, parallelWrite, parallelRead;
	int numlut, numColMuxed, numWriteColMuxed, levelOutput, avgWeightBit, numBitInput;
	int numRowSubArray, numColSubArray;
	int cellBit, synapseBit;
	int numTileRow, numTileCol;
	int numPERowPerTile, numPEColPerTile;
	int numSubArrayRowPerPE, numSubArrayColPerPE;
	bool Scenariotype;
	bool Baseline;
	double memory_bus_bandwidth;
	bool mode;
	bool pipeline_dynamic_dup;
	int batch_size;

	int XNORparallelMode, XNORsequentialMode, BNNparallelMode, BNNsequentialMode, conventionalParallel, conventionalSequential; 
	int numRowPerSynapse, numColPerSynapse;
	double AR, Rho, wireLengthRow, wireLengthCol, unitLengthWireResistance, wireResistanceRow, wireResistanceCol;
	
	double alpha, beta, gamma, delta, epsilon, zeta;
	
	double Metal0=0;
	double Metal1=0;
	double AR_Metal0=0;
	double AR_Metal1=0;
	double Rho_Metal0=0;
	double Rho_Metal1=0;
	double Metal0_unitwireresis=0;
	double Metal1_unitwireresis=0;

	bool Activationtype; // true: SRAM, False: RRAM

	// 1.4 update: Final driver sizing for row decoder conventional parallel mode (SRAM, RRAM)
	// multiplied by the driver width
	double sizingfactor_MUX= 1; 
	double sizingfactor_WLdecoder= 1; 

	// 1.4 update: switchmatrix parameter tuning
	double newswitchmatrixsizeratio=6;
	double switchmatrixsizeratio=1;
	
	// 1.4 update: Special layout
	double speciallayout;
	
	// 1.4 update: added parameters for buffer insertion
	double unitcap;
	double unitres;
	double drivecapin; 
	double buffernumber=0;
	double buffersizeratio=0;
	
	// 1.4 update: barrier thickness
	double barrierthickness= 0;
	
	// 1.4 update: new ADC modeling related parameters
	double dumcolshared;
	double columncap;
	double reference_energy_peri=0;
	
	// 1.4 update: array dimension/SRAM access resistance for multilevelsenseamp
	double arrayheight;
	double arraywidthunit;
	double resCellAccess;

	// 1.4 update 
	double inputtoggle;
	double outputtoggle;

	// 1.4 debug
	double ADClatency;
	double rowdelay;
	double muxdelay;
	
	// 1.4 update: technology node
	int technologynode;

	// Anni update: partial parallel mode
	int numRowParallel;

	// 230920 update
	double totaltile_num;
};

#endif