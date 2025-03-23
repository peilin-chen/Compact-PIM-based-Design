#ifndef CHIP_H_
#define CHIP_H_

/*** Functions ***/
void ChipDesignInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell);

vector<double> Mapping(bool findsubArrayDup, bool findpeDup, bool findtileDup, bool findNumTile, bool find_no_tileDup_NumTile, bool findUtilization, bool findSpeedUp, const vector<vector<double> > &netStructure, 
					double *desiredNumTileCM, double *desiredTileSizeCM, double *desiredPESizeCM, int *numTileRow, int *numTileCol, int raw_processing, int layer_by_layer, 
					int full_pipeline, int partial_pipeline);

void ChipInitialize(InputParameter& inputParameter, Technology& tech, MemCell& cell, const vector<vector<double> > &netStructure, const vector<double> &numTileEachLayer,
					double desiredNumTileCM, double desiredTileSizeCM, double desiredPESizeCM, int numTileRow, int numTileCol, int full_pipeline, const vector<double> &tileDup);
		
vector<double> ChipCalculateArea(InputParameter& inputParameter, Technology& tech, MemCell& cell, double desiredNumTileCM, double desiredTileSizeCM, double desiredPESizeCM, 
						int numTileRow, double *CMTileheight, double *CMTilewidth);

double ChipCalculatePerformance(InputParameter& inputParameter, Technology& tech, MemCell& cell, int layerNumber, const string &newweightfile, const string &oldweightfile, const string &inputfile, bool followedByMaxPool, const vector<vector<double> > &netStructure, 
							const vector<double> &numTileEachLayer, const vector<double> &speedUpEachLayer, 
							const vector<vector<double> > &tileLocaEachLayer, double desiredTileSizeCM, double desiredPESizeCM, 
							double CMTileheight, double CMTilewidth, double *readLatency, double *readDynamicEnergy, double *writeLatency, double *writeDynamicEnergy,
							double *leakage, double *leakageSRAMInUse, double *bufferLatency, double *bufferDynamicEnergy, double *icLatency, double *icDynamicEnergy,
							double *coreLatencyADC, double *coreLatencyAccum, double *coreLatencyOther, double *coreEnergyADC, double *coreEnergyAccum, double *coreEnergyOther, 
							bool CalculateclkFreq, double *clkPeriod, int *partition_col, int *partition_row, const vector<double> &tileDup, const vector<double> &peDup, 
							const vector<double> &subArrayDup, double *act_poolLatency, double *computationLatency);

vector<double> TileDup(double tileSize, const vector<double> &subArrayDup, const vector<double> &peDup, const vector<vector<double> > &netStructure, int numRowPerSynapse, int numColPerSynapse, 
						int raw_processing, int layer_by_layer, int full_pipeline, int partial_pipeline);

vector<double> PeDup(const vector<double> &subArrayDup, double peSize, double desiredTileSize,
								const vector<vector<double> > &netStructure, int numRowPerSynapse, int numColPerSynapse);

vector<double> SubArrayDup(double desiredPESizeCM, const vector<vector<double> > &netStructure, int numRowPerSynapse, int numColPerSynapse);

vector<double> OverallEachLayer(bool utilization, bool speedUp, bool numTile, const vector<double> &tileDup, const vector<double> &peDup, 
										const vector<double> &subArrayDup, double desiredTileSizeCM, const vector<vector<double> > &netStructure, 
										int numRowPerSynapse, int numColPerSynapse);

vector<vector<double> > LoadInWeightData(const string &weightfile, int numRowPerSynapse, int numColPerSynapse, double maxConductance, double minConductance);
vector<vector<double> > CopyArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol);
vector<vector<double> > ReshapeArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol, int numPE, int weightMatrixRow);
vector<vector<double> > LoadInInputData(const string &inputfile);
vector<vector<double> > CopyInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow);
vector<vector<double> > ReshapeInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow, int numPE, int weightMatrixRow);

#endif /* CHIP_H_ */