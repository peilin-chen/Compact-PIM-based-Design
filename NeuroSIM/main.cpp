#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "Tile.h"
#include "Chip.h"
#include "ProcessingUnit.h"
#include "SubArray.h"
#include "Definition.h"

using namespace std;

vector<vector<double> > getNetStructure(const string &inputfile);

string generateNextAddress(int &address_counter) {
    stringstream ss;
    
    ss << "0x" << setw(8) << setfill('0') << hex << address_counter;
    address_counter++;
    return ss.str();
}

int main(int argc, char * argv[]) {   

    auto start = chrono::high_resolution_clock::now();
    
    gen.seed(0);
    
    vector<vector<double> > netStructure;
    netStructure = getNetStructure(argv[1]); // get neural network structure
    
    // define weight/input/memory precision from wrapper
    param->synapseBit     = atoi(argv[2]);          // precision of weight
    param->numBitInput    = atoi(argv[3]);          // precision of input activation
    param->numRowSubArray = atoi(argv[4]);          // number row of subarray
    param->numRowParallel = atoi(argv[5]);          // number of enabled rows of subarray

    if (param->cellBit > param->synapseBit) {
        cout << "ERROR!: Memory precision is even higher than synapse precision, please modify 'cellBit' in Param.cpp!" << endl;
        param->cellBit = param->synapseBit;
    }

    // warning for the incompatible modes for a given device technology
    if ( param->memcelltype == 2 && param->technode <=14) {
        cout << "ERROR!: RRAM-CIM is not supported beyond 22 nm!" << endl;
        exit(-1);
    }
    if ( param->deviceroadmap == 1 && param->technode <=14) {
        cout << "ERROR!: HP technology is not supported for 14 nm and beyond!" << endl;
        exit(-1);
    }
    if ( param->temp != 300 && param->technode <=14) {
        cout << "ERROR!: only 300K operation is supported for 14 nm and beyond!" << endl;
        exit(-1);
    }

    /*** initialize operationMode as default ***/
    param->conventionalParallel   = 0;
    param->conventionalSequential = 0;
    param->BNNparallelMode        = 0; // parallel BNN
    param->BNNsequentialMode      = 0; // sequential BNN
    param->XNORsequentialMode     = 0; // Use several multi-bit RRAM as one synapse
    param->XNORparallelMode       = 0; // Use several multi-bit RRAM as one synapse

    switch(param->operationmode) {
        case 6:     param->XNORparallelMode = 1;               break;     
        case 5:     param->XNORsequentialMode = 1;             break;     
        case 4:     param->BNNparallelMode = 1;                break;     
        case 3:     param->BNNsequentialMode = 1;              break;    
        case 2:     param->conventionalParallel = 1;           break;     
        case 1:     param->conventionalSequential = 1;         break;    
        case -1:    break;
        default:    exit(-1);
    }
    
    if (param->XNORparallelMode || param->XNORsequentialMode) {
        param->numRowPerSynapse = 2;
    } 
    else {
        param->numRowPerSynapse = 1;
    }

    if (param->BNNparallelMode) {
        param->numColPerSynapse = 2;
    } 
    else if (param->XNORparallelMode || param->XNORsequentialMode || param->BNNsequentialMode) {
        param->numColPerSynapse = 1;
    } 
    else {
        param->numColPerSynapse = ceil((double)param->synapseBit/(double)param->cellBit); 
    }
    
    if ((param->conventionalSequential == 1 || param->conventionalSequential == 3 || param->conventionalSequential == 5 ) && (param->memcelltype==1))
    {
        param->numColMuxed = param->numColPerSynapse;
    }

    ChipDesignInitialize(inputParameter, tech, cell);
   
    double desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM;
    int numTileRow, numTileCol;

    vector<vector<double> > tileLocaEachLayer;
    
    vector<double> numTileEachLayer           ;
    vector<double> no_tileDup_numTileEachLayer;
    vector<double> speedUpEachLayer           ;
    vector<double> utilizationEachLayer       ;
    vector<double> subArrayDup                ;
    vector<double> peDup                      ;
    vector<double> tileDup                    ;

    double totalNumTile_original = 0;

    for (int i=0; i<netStructure.size(); i++) {
        double numtileEachLayerRow_original = ceil((double) netStructure[i][2]*(double) netStructure[i][3]*(double) netStructure[i][4]*(double) param->numRowPerSynapse/(param->numPERowPerTile*param->numSubArrayRowPerPE*param->numRowSubArray));
	    double numtileEachLayerCol_original = ceil((double) netStructure[i][5]*(double) param->numColPerSynapse/(double) (param->numPERowPerTile*param->numSubArrayRowPerPE*param->numRowSubArray));

        totalNumTile_original = totalNumTile_original + numtileEachLayerRow_original*numtileEachLayerCol_original;
    }

    int raw_processing   = 0;
    int layer_by_layer   = 0;
    int full_pipeline    = 0;
    int partial_pipeline = 0;

    // define using which processing method
    if (param->Baseline) {
        raw_processing = 1;
    }
    else {
        if (param->Scenariotype) {
            layer_by_layer = 1;
        }
        else {
            if (totalNumTile_original <= (param->numTileRow*param->numTileCol)) {
                full_pipeline = 1;
            }
            else {
                partial_pipeline = 1;
            }
        }
    }
    //cout << 1 << endl;
    subArrayDup                 = Mapping(true , false, false, false, false, false, false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    peDup                       = Mapping(false, true , false, false, false, false, false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    tileDup                     = Mapping(false, false, true , false, false, false, false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    numTileEachLayer            = Mapping(false, false, false, true , false, false, false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);    
    no_tileDup_numTileEachLayer = Mapping(false, false, false, false, true , false, false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    utilizationEachLayer        = Mapping(false, false, false, false, false, true , false, netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    speedUpEachLayer            = Mapping(false, false, false, false, false, false, true , netStructure, &desiredNumTileCM, &desiredTileSizeCM, &desiredPESizeCM, &numTileRow, &numTileCol, raw_processing, layer_by_layer, full_pipeline, partial_pipeline);
    
    //cout << 2 << endl;

    cout << "------------------------------ Hardware Resources --------------------------------" <<  endl;
    cout << endl;
    cout << "Number of Tile in the Accelerator: " << param->numTileRow << "x" << param->numTileCol << endl;
    cout << "Tile Storage Size: " << desiredTileSizeCM << "x" << desiredTileSizeCM << endl;
    cout << "PE Storage Size: " << desiredPESizeCM << "x" << desiredPESizeCM << endl; 
    cout << "SubArray Size: " << param->numRowSubArray << "x" << param->numColSubArray << endl;
    cout << endl;

    vector<vector<int> > record;
    vector<int> recordRow;

    if (raw_processing || partial_pipeline) {
        for (int i=0; i<netStructure.size(); i++) {
            if ((no_tileDup_numTileEachLayer[i]) <= (param->numTileRow*param->numTileCol)) {

                double p = no_tileDup_numTileEachLayer[i];
                recordRow.push_back(i);
                
                for(int j=i+1; j<netStructure.size(); j++) {
                    p = p + no_tileDup_numTileEachLayer[j];
                    if (p <= (param->numTileRow*param->numTileCol)) {
                        recordRow.push_back(j);
                    }
                    else {
                        break;
                    }
                }
                record.push_back(recordRow);
                i = recordRow.size()-1;
            }
            else {
                recordRow.push_back(i);
                record.push_back(recordRow);
                i = recordRow.size()-1;
            }
        }
    }

    vector<double> tileLocaEachLayerRow;
	vector<double> tileLocaEachLayerCol;

	int temp = 0;

	if ((!full_pipeline) && (!raw_processing) && (!partial_pipeline)) {
		for (int i=0; i<netStructure.size(); i++) {
			tileLocaEachLayerRow.push_back(0);
			tileLocaEachLayerCol.push_back(0);
		}
	}
	else if (full_pipeline) {
		double thisTileTotal;
		for (int i=0; i<netStructure.size(); i++) {
			if (i==0) {
				tileLocaEachLayerRow.push_back(0);
				tileLocaEachLayerCol.push_back(0);
			} 
			else {
				thisTileTotal += numTileEachLayer[i-1];
				tileLocaEachLayerRow.push_back((int)thisTileTotal/(numTileRow));
				tileLocaEachLayerCol.push_back((int)thisTileTotal%(numTileRow));
			}
		}
	}
	else if (raw_processing) {
		if (record.size() > 1) {
			for (int i=0; i<record.size(); i++) {
				double thisTileTotal;
				for (int j=temp; j<record[i].size(); j++) {
					if (j==0 || (j==record[i-1].size() && i!=0)) {
						tileLocaEachLayerRow.push_back(0);
						tileLocaEachLayerCol.push_back(0);
					}
					else {
						thisTileTotal += numTileEachLayer[j-1];
						tileLocaEachLayerRow.push_back((int)thisTileTotal/(numTileRow));
						tileLocaEachLayerCol.push_back((int)thisTileTotal%(numTileRow));
					}
					temp = j + 1;
				}
				thisTileTotal = 0;
			}
		}
		else {
			double thisTileTotal;
			for (int i=0; i<netStructure.size(); i++) {
				if (i==0) {
					tileLocaEachLayerRow.push_back(0);
					tileLocaEachLayerCol.push_back(0);
				} 
				else {
					thisTileTotal += numTileEachLayer[i-1];
					tileLocaEachLayerRow.push_back((int)thisTileTotal/(numTileRow));
					tileLocaEachLayerCol.push_back((int)thisTileTotal%(numTileRow));
				}
			}
		}
	}
	else {
		if (record.size() > 1) {
			for (int i=0; i<record.size(); i++) {
				double thisTileTotal;
				for (int j=temp; j<record[i].size(); j++) {
					if (j==0 || (j==record[i-1].size() && i!=0)) {
						tileLocaEachLayerRow.push_back(0);
						tileLocaEachLayerCol.push_back(0);
					}
					else {
						thisTileTotal += numTileEachLayer[j-1];
						tileLocaEachLayerRow.push_back((int)thisTileTotal/(numTileRow));
						tileLocaEachLayerCol.push_back((int)thisTileTotal%(numTileRow));
					}
					temp = j + 1;
				}
				thisTileTotal = 0;
			}
		}
	}

	tileLocaEachLayer.push_back(tileLocaEachLayerRow);
	tileLocaEachLayer.push_back(tileLocaEachLayerCol);

    cout << "----------------- # of tile used for each layer -----------------" <<  endl;
    cout << endl;
    
    if (raw_processing || partial_pipeline) {
        if (record.size() > 1) {
            cout << "Note: accelerator cannot store all the weights of this NN!" << endl;
            cout << "This NN is divided into " << record.size() << " parts to process!" << endl;
            for (int i=0; i<record.size(); i++){
                if (i == 0) {
                    cout << "Part " << (i+1) << " includes " << record[i].size() <<" layers!" << endl;
                }
                else {
                    cout << "Part " << (i+1) << " includes " << record[i].size() - record[i-1].size() <<" layers!" << endl;
                }
            }
            cout << endl;
        }
    }

    vector<int> partition_col_array;
    vector<int> partition_row_array;

    for (int i=0; i<netStructure.size(); i++) {

        int weightMatrixRow = netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*param->numRowPerSynapse;
	    int weightMatrixCol = netStructure[i][5]*param->numColPerSynapse;

        int col = 0;
        int row = 0;
    
        if (numTileEachLayer[i] > param->numTileRow*param->numTileCol) {
            int flag = ceil(numTileEachLayer[i]/(param->numTileRow*param->numTileCol));
			
            flag = ceil(flag/2);

            for (int i=0; i<flag; i++) {
                if (weightMatrixCol > weightMatrixRow) {
                    weightMatrixCol = weightMatrixCol/2;
                    col = col+1;
                }
                else {
                    weightMatrixRow = weightMatrixRow/2;
                    row = row+1;
                }
            }
        }

        int partition_Col = pow(2, col);
        int partition_Row = pow(2, row);

        partition_col_array.push_back(partition_Col);
        partition_row_array.push_back(partition_Row);

        if (raw_processing) {
            if (partition_Col>1 || partition_Row>1) {
                cout << "layer" << i+1 << ": " << (numTileEachLayer[i]/(partition_Row*partition_Col)) << " (Row partition of this layer: " << partition_Row << 
                " Col partition of this layer: " << partition_Col << ")" << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
            else {
                cout << "layer" << i+1 << ": " << numTileEachLayer[i] << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
        }
        else if (layer_by_layer) {
            if (partition_Col>1 || partition_Row>1) {
                cout << "layer" << i+1 << ": " << (numTileEachLayer[i]/(partition_Row*partition_Col)) << " (Row partition of this layer: " << partition_Row << 
                " Col partition of this layer: " << partition_Col << ")" << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
            else {
                cout << "layer" << i+1 << ": " << numTileEachLayer[i] << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
        }
        else if (full_pipeline) {
            cout << "layer" << i+1 << ": " << numTileEachLayer[i] << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
        }
        else if (partial_pipeline) {
            if (partition_Col>1 || partition_Row>1) {
                cout << "layer" << i+1 << ": " << (numTileEachLayer[i]/(partition_Row*partition_Col)) << " (Row partition of this layer: " << partition_Row << 
                " Col partition of this layer: " << partition_Col << ")" << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
            else {
                cout << "layer" << i+1 << ": " << numTileEachLayer[i] << " (Original Tile requirements: " << no_tileDup_numTileEachLayer[i] << ")"<< endl;
            }
        }
    }
    cout << endl;

    cout << "----------------- Speed-up of each layer ------------------" <<  endl;
    cout << endl;

    for (int i=0; i<netStructure.size(); i++) {
       
        int maxDup = netStructure[i][0]*netStructure[i][1];

        if (raw_processing) {
            cout << "layer" << i+1 << ": " << speedUpEachLayer[i] << "(subarray: "<< subArrayDup[i] << ")" << 
            " (pe: "<< peDup[i] << ")" << " (tile: "<< tileDup[i] << ")" << " (Ideal speed-up: " << maxDup << ")" << endl;
        }
        else if (layer_by_layer) {
            cout << "layer" << i+1 << ": " << speedUpEachLayer[i] << "(subarray: "<< subArrayDup[i] << ")" << 
            " (pe: "<< peDup[i] << ")" << " (tile: "<< tileDup[i] << ")" << " (Ideal speed-up: " << maxDup << ")" << endl;
        }
        else if (full_pipeline) {
            cout << "layer" << i+1 << ": " << speedUpEachLayer[i] << "(subarray: "<< subArrayDup[i] << ")" << 
            " (pe: "<< peDup[i] << ")" << " (tile: "<< tileDup[i] << ")" << " (Ideal speed-up: " << maxDup << ")" << endl;
        }
        else if (partial_pipeline) {
            cout << "layer" << i+1 << ": " << speedUpEachLayer[i] << "(subarray: "<< subArrayDup[i] << ")" << 
            " (pe: "<< peDup[i] << ")" << " (tile: "<< tileDup[i] << ")" << " (Ideal speed-up: " << maxDup << ")" << endl;
        }
    }
    cout << endl;

    cout << "----------------- Utilization of each layer ------------------" <<  endl;
    cout << endl;

    double realMappedMemory = 0;

    for (int i=0; i<netStructure.size(); i++) {

        if (raw_processing) {
            double utilization = (tileDup[i]*peDup[i]*subArrayDup[i]*netStructure[i][2]*netStructure[i][3]*netStructure[i][4]
									*param->numRowPerSynapse*netStructure[i][5]*param->numColPerSynapse)/(numTileEachLayer[i]*desiredTileSizeCM*desiredTileSizeCM);
            cout << "layer" << i+1 << ": " << utilization << endl;
        }
        else if (layer_by_layer) {
            cout << "layer" << i+1 << ": " << utilizationEachLayer[i]/(partition_col_array[i]*partition_row_array[i]) << endl;
        }
        else if (full_pipeline) {
            double utilization = (tileDup[i]*peDup[i]*subArrayDup[i]*netStructure[i][2]*netStructure[i][3]*netStructure[i][4]
									*param->numRowPerSynapse*netStructure[i][5]*param->numColPerSynapse)/(numTileEachLayer[i]*desiredTileSizeCM*desiredTileSizeCM);
            cout << "layer" << i+1 << ": " << utilization << endl;
            realMappedMemory = realMappedMemory + utilization*numTileEachLayer[i];
        }
        else if (partial_pipeline) {
            double utilization = (tileDup[i]*peDup[i]*subArrayDup[i]*netStructure[i][2]*netStructure[i][3]*netStructure[i][4]
									*param->numRowPerSynapse*netStructure[i][5]*param->numColPerSynapse)/(numTileEachLayer[i]*desiredTileSizeCM*desiredTileSizeCM);
            cout << "layer" << i+1 << ": " << utilization << endl;
        }
    }
    
    if (full_pipeline) {
        cout << "Memory Utilization of Whole Chip: " << realMappedMemory/(numTileRow*numTileCol)*100 << " % " << endl;
    }
    cout << endl;

    if (raw_processing) {
        ChipInitialize(inputParameter, tech, cell, netStructure, numTileEachLayer, desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM, numTileRow, numTileCol, full_pipeline, tileDup);
    }
    else if (layer_by_layer) {
        ChipInitialize(inputParameter, tech, cell, netStructure, numTileEachLayer, desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM, numTileRow, numTileCol, full_pipeline, tileDup);
    }
    else if (full_pipeline) {
        ChipInitialize(inputParameter, tech, cell, netStructure, numTileEachLayer, desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM, numTileRow, numTileCol, full_pipeline, tileDup);
    }
    else if (partial_pipeline) {
        ChipInitialize(inputParameter, tech, cell, netStructure, numTileEachLayer, desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM, numTileRow, numTileCol, full_pipeline, tileDup);
    }

    double chipArea, chipAreaIC, chipAreaADC, chipAreaAccum, chipAreaOther, chipAreaArray;
    double CMTileheight = 0;
    double CMTilewidth  = 0;

    vector<double> chipAreaResults;
    
    chipAreaResults = ChipCalculateArea(inputParameter, tech, cell, desiredNumTileCM, desiredTileSizeCM, desiredPESizeCM, numTileRow, &CMTileheight, &CMTilewidth);
    
    chipArea      = chipAreaResults[0];
    chipAreaIC    = chipAreaResults[1];
    chipAreaADC   = chipAreaResults[2];
    chipAreaAccum = chipAreaResults[3];
    chipAreaOther = chipAreaResults[4];
    chipAreaArray = chipAreaResults[5];

    cout << "------------------------------ Area Overhead --------------------------------" <<  endl;
    cout << endl;
    cout << "ChipArea : " << chipArea*1e12 << "um^2" << endl;
    cout << "Chip total CIM array : " << chipAreaArray*1e12 << "um^2" << endl;
    cout << "Total IC Area on chip (Global and Tile/PE local): " << chipAreaIC*1e12 << "um^2" << endl;
    cout << "Total ADC (or S/As and precharger for SRAM) Area on chip : " << chipAreaADC*1e12 << "um^2" << endl;
    cout << "Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << chipAreaAccum*1e12 << "um^2" << endl;
    cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << chipAreaOther*1e12 << "um^2" << endl;
    cout << endl;
    
    // using by ChipCalculatePerformance function unit
    double clkPeriod                = 0;
    double layerclkPeriod           = 0;
    double layerReadLatency         = 0;
    double layerReadDynamicEnergy   = 0;
    double layerWriteLatency        = 0;
    double layerWriteDynamicEnergy  = 0;
    double tileLeakage              = 0;
    double tileLeakageSRAMInUse     = 0;
    double layerbufferLatency       = 0;
    double layerbufferDynamicEnergy = 0;
    double layericLatency           = 0;
    double layericDynamicEnergy     = 0;
    double coreLatencyADC           = 0;
    double coreLatencyAccum         = 0;
    double coreLatencyOther         = 0;
    double coreEnergyADC            = 0;
    double coreEnergyAccum          = 0;
    double coreEnergyOther          = 0;

    double layer_act_poolLatency   = 0;
    double layerComputationLatency = 0;
    
    int partition_col = 0;
    int partition_row = 0;
    // use to record overall performance 
    double chipReadLatency             = 0;
    double chipReadDynamicEnergy       = 0;
    double chipWriteLatency            = 0;
    double chipWriteDynamicEnergy      = 0;
    double chipLeakageEnergy           = 0;
    double chipLeakage                 = 0;
    double chipbufferLatency           = 0;
    double chipbufferReadDynamicEnergy = 0;
    double chipicLatency               = 0;
    double chipicReadDynamicEnergy     = 0;
    double chipLatencyADC              = 0;
    double chipLatencyAccum            = 0;
    double chipLatencyOther            = 0;
    double chipEnergyADC               = 0;
    double chipEnergyAccum             = 0;
    double chipEnergyOther             = 0;

    double off_chip_dramLatency = 0;
    double off_chip_dramEnergy  = 0;
    
    // netStructure[i][6]: represent followed by pooling or not 
    // netStructure[i][7]: represent the strides of convolutional kernel

    double numComputationPerLayer = 0;
    double numComputation         = 0;

    double data_transfer_num = 0;

    cout << "-------------------------------------- Hardware Performance --------------------------------------" <<  endl;  
    cout << endl;
    if (raw_processing) {

        int numRead_transactions = 0;
        int numWrite_transactions = 0;
        int dram_bus_width = 128; // this param needs to be changed based on the chosen memory

        // open a file to write data movement trace
        ofstream outfile; 
        outfile.open("data_movement-RP.trace", ios::out | ios::trunc);

        if (!outfile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return 1;
        }

        int timestamp=-1;
        int address_counter = 0;
        int flag = 0;

        for (int i=0; i<record.size(); i++) {

            timestamp = timestamp+1;

            // record loading IFM from off-chip memory
            numRead_transactions = param->batch_size*netStructure[flag][0]*netStructure[flag][1]*netStructure[flag][2]*netStructure[flag][3]*netStructure[flag][4]*(flag==0 ? param->numBitInput:(param->numBitInput+param->synapseBit));
            numRead_transactions = ceil(numRead_transactions/dram_bus_width);

            address_counter = address_counter - numWrite_transactions;

            if(i==0) {
                cout << "Access number of IFM: " << numRead_transactions << endl;
            }

            for (int j=0; j<numRead_transactions; j++) {

                string transaction_type = "READ";

                string memory_address = generateNextAddress(address_counter);

                outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
            }

            data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;

            for (int j=flag; j<record[i].size(); j++) {
                flag = j+1;

                if (partition_col_array[j]>1 || partition_row_array[j]>1) {
                    // if partition_row_array[j]>1, there are additional (read(IFM))*(partition_row_array[j]-1) operations.
                    for (int k=0; k<(partition_row_array[j]-1); k++) {

                        address_counter = address_counter - numRead_transactions;

                        for (int z=0; z<numRead_transactions; z++) {

                            string transaction_type = "READ";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }
                        data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
                    }

                    // if partition_col_array[j]>1, there are additional (write+read(this layer computation result))*(partition_col_array[j]-1) operations.
                    for (int k=0; k<(partition_col_array[j]-1); k++) {
                        numWrite_transactions = param->batch_size*netStructure[j][0]*netStructure[j][1]*netStructure[j][2]*(param->numBitInput+param->synapseBit);
                        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                        for (int z=0; z<numWrite_transactions; z++) {
                            string transaction_type = "WRITE";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }
                        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

                        address_counter = address_counter - numWrite_transactions;

                        for (int z=0; z<numWrite_transactions; z++) {
                            string transaction_type = "READ";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }
                        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
                    }
                }

                // record loading weights from off-chip memory
                numRead_transactions =  netStructure[j][2]*netStructure[j][3]*netStructure[j][4]*netStructure[j][5]*param->synapseBit;
                numRead_transactions = ceil(numRead_transactions/dram_bus_width);

                for (int k=0; k<speedUpEachLayer[j]; k++) {
                    
                    if (k != 0) {
                        address_counter = address_counter-numRead_transactions;
                    }

                    for (int z=0; z<numRead_transactions; z++) {

                        string transaction_type = "READ";

                        string memory_address = generateNextAddress(address_counter);

                        outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                    }
                    data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
                }
            }

            // record writing intermediate data to DRAM
            if (i != (record.size()-1)) {
                numWrite_transactions = param->batch_size*netStructure[flag-1][0]*netStructure[flag-1][1]*netStructure[flag-1][2]*(param->numBitInput+param->synapseBit);
                numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                timestamp = timestamp+1;

                for (int j=0; j<numWrite_transactions; j++) {
                    string transaction_type = "WRITE";

                    string memory_address = generateNextAddress(address_counter);

                    outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                }
                data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
            }
        }

        // record writing computation result to DRAM
        int last_layer = netStructure.size()-1;

        numWrite_transactions = param->batch_size*netStructure[last_layer][5]*(param->numBitInput+param->synapseBit);
        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

        cout << "Access number of Final result: " << numWrite_transactions << endl;

        timestamp = timestamp+1; // Finally, write computation result into DRAM 

        for (int i=0; i<numWrite_transactions; i++) {
            string transaction_type = "WRITE";

            string memory_address = generateNextAddress(address_counter);

            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
        }

        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

        outfile.close();

        // show the detailed hardware performance for each layer
        for (int i=0; i<netStructure.size(); i++) {
            numComputationPerLayer = ((record.size()>=2)?param->batch_size:1)*2*(netStructure[i][0] * netStructure[i][1] * netStructure[i][2] * netStructure[i][3] * netStructure[i][4] * netStructure[i][5]);
            numComputation = numComputation + numComputationPerLayer;
            cout << "-----------------------------------------------------------------------------" << endl;
            cout << "-------------------- Estimation of Layer " << i+1 << " ----------------------" << endl;
            //cout << "-----------------------------------------------------------------------------" << endl;
            //cout << endl;
            ChipCalculatePerformance(inputParameter, tech, cell, i, argv[2*i+6], argv[2*i+6], argv[2*i+7], netStructure[i][6], 
                        netStructure, numTileEachLayer, speedUpEachLayer, tileLocaEachLayer, desiredTileSizeCM, desiredPESizeCM, 
                        CMTileheight, CMTilewidth, &layerReadLatency, &layerReadDynamicEnergy, &layerWriteLatency, &layerWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &layerbufferLatency, 
                        &layerbufferDynamicEnergy, &layericLatency, &layericDynamicEnergy, &coreLatencyADC, &coreLatencyAccum, &coreLatencyOther, &coreEnergyADC, &coreEnergyAccum, 
                        &coreEnergyOther, false, &layerclkPeriod, &partition_col, &partition_row, tileDup, peDup, subArrayDup, &layer_act_poolLatency, &layerComputationLatency);

            if (partition_col>1 || partition_row>1) {
                cout << "--------Note: accelerator cannot store all the weights of this layer!--------" << endl;
                cout << "------ Row partition of this layer: " << partition_row << " Col partition of this layer: " << partition_col << "------" << endl;
                cout << "-----------------------------------------------------------------------------" << endl;
                cout << endl;
            }
            else {
                cout << "-----------------------------------------------------------------------------" << endl;
                cout << endl;
            }

            double numTileOtherLayer  = 0;
            double layerLeakageEnergy = 0;      
            
            if (record.size() > 1) {
                for (int j=0; j<record.size(); j++) {
                    if (i < record[j].size()) {
                        for (int k=i; k<record[j].size(); k++) {
                            if (k != i) {
                                numTileOtherLayer += numTileEachLayer[k];
                            }
                        }
                    }
                }
            }
            else {
                for (int j=0; j<netStructure.size(); j++) {
                    if (j != i) {
                        numTileOtherLayer += numTileEachLayer[j];
                    }
                }
            }
            
            layerLeakageEnergy = (numTileOtherLayer * tileLeakage + numTileEachLayer[i] * tileLeakageSRAMInUse) * layerReadLatency;

            cout << "layer" << i+1 << "'s inference time is: " << (layerReadLatency+layerWriteLatency)*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readLatency is: " << layerReadLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readDynamicEnergy is: " << layerReadDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s writeLatency is: " << layerWriteLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s writeDynamicEnergy is: " << layerWriteDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s leakagePower is: " << numTileEachLayer[i] * tileLeakage*1e6 << "uW" << endl;
            cout << "layer" << i+1 << "'s leakageEnergy is: " << layerLeakageEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s buffer latency is: " << layerbufferLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s buffer readDynamicEnergy is: " << layerbufferDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s ic latency is: " << layericLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s ic readDynamicEnergy is: " << layericDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s computation latency is: " << layerComputationLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s activation and pool latency is: " << layer_act_poolLatency*1e9 << "ns" << endl;
            
            cout << endl;
            cout << "************ Breakdown of Latency and Dynamic Energy *************" << endl;
            cout << endl;

            cout << "ADC (or S/As and precharger for SRAM) readLatency is : " << coreLatencyADC*1e9 << "ns" << endl;
            cout << "Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << coreLatencyAccum*1e9 << "ns" << endl;
            cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << coreLatencyOther*1e9 << "ns" << endl;
            cout << "ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << coreEnergyADC*1e12 << "pJ" << endl;
            cout << "Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << coreEnergyAccum*1e12 << "pJ" << endl;
            cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << coreEnergyOther*1e12 << "pJ" << endl;
            cout << endl;

            chipReadLatency += layerReadLatency;
            chipReadDynamicEnergy += layerReadDynamicEnergy;
            chipWriteLatency += layerWriteLatency;
            chipWriteDynamicEnergy += layerWriteDynamicEnergy;
            chipLeakageEnergy += layerLeakageEnergy;
            chipLeakage += tileLeakage*numTileEachLayer[i];
            chipbufferLatency += layerbufferLatency;
            chipbufferReadDynamicEnergy += layerbufferDynamicEnergy;
            chipicLatency += layericLatency;
            chipicReadDynamicEnergy += layericDynamicEnergy;
            
            chipLatencyADC += coreLatencyADC;
            chipLatencyAccum += coreLatencyAccum;
            chipLatencyOther += coreLatencyOther;
            chipEnergyADC += coreEnergyADC;
            chipEnergyAccum += coreEnergyAccum;
            chipEnergyOther += coreEnergyOther;
        }
    }
    else if (layer_by_layer) {

        int numRead_transactions = 0;
        int numWrite_transactions = 0;
        int dram_bus_width = 128; // this param needs to be changed based on the chosen memory

        // open a file to write data movement trace
        ofstream outfile; 
        outfile.open("data_movement-LP.trace", ios::out | ios::trunc);

        if (!outfile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return 1;
        }

        int timestamp=-1;
        int address_counter = 0;

        // show the detailed hardware performance for each layer
        for (int i=0; i<netStructure.size(); i++) {
            numComputationPerLayer = 2*(netStructure[i][0] * netStructure[i][1] * netStructure[i][2] * netStructure[i][3] * netStructure[i][4] * netStructure[i][5]);
            numComputation = numComputation + numComputationPerLayer;

            timestamp = timestamp+1;

            // record loading IFM from off-chip memory
            numRead_transactions = netStructure[i][0]*netStructure[i][1]*netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*(i==0 ? param->numBitInput:(param->numBitInput+param->synapseBit));
            numRead_transactions = ceil(numRead_transactions/dram_bus_width);

            address_counter = address_counter - numWrite_transactions;

            for (int j=0; j<numRead_transactions; j++) {

                string transaction_type = "READ";

                string memory_address = generateNextAddress(address_counter);

                outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
            }
            
            data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;

            if (partition_col_array[i]>1 || partition_row_array[i]>1) {
                // if partition_row_array[j]>1, there are additional (read(IFM))*(partition_row_array[j]-1) operations.
                for (int j=0; j<(partition_row_array[i]-1); j++) {

                    address_counter = address_counter - numRead_transactions;

                    for (int k=0; k<numRead_transactions; k++) {

                        string transaction_type = "READ";

                        string memory_address = generateNextAddress(address_counter);

                        outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                    }
                    data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
                }

                // if partition_col_array[j]>1, there are additional (write+read(this layer computation result))*(partition_col_array[j]-1) operations.
                for (int j=0; j<(partition_col_array[i]-1); j++) {
                    numWrite_transactions = netStructure[i][0]*netStructure[i][1]*netStructure[i][2]*(param->numBitInput+param->synapseBit);
                    numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                    for (int k=0; k<numWrite_transactions; k++) {
                        string transaction_type = "WRITE";

                        string memory_address = generateNextAddress(address_counter);

                        outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                    }
                    data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

                    address_counter = address_counter - numWrite_transactions;

                    for (int k=0; k<numWrite_transactions; k++) {
                        string transaction_type = "READ";

                        string memory_address = generateNextAddress(address_counter);

                        outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                    }
                    data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
                }
            }

            // record loading weights from off-chip memory
            numRead_transactions =  netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*netStructure[i][5]*param->synapseBit;
            numRead_transactions = ceil(numRead_transactions/dram_bus_width);

            for (int k=0; k<speedUpEachLayer[i]; k++) {
                
                if (k != 0) {
                    address_counter = address_counter-numRead_transactions;
                }

                for (int j=0; j<numRead_transactions; j++) {

                    string transaction_type = "READ";

                    string memory_address = generateNextAddress(address_counter);

                    outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                }
                data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
            }

            cout << "-----------------------------------------------------------------------------" << endl;
            cout << "-------------------- Estimation of Layer " << i+1 << " ----------------------" << endl;
            //cout << "-----------------------------------------------------------------------------" << endl;
            //cout << endl;
            ChipCalculatePerformance(inputParameter, tech, cell, i, argv[2*i+6], argv[2*i+6], argv[2*i+7], netStructure[i][6],
                        netStructure, numTileEachLayer, speedUpEachLayer, tileLocaEachLayer, desiredTileSizeCM, desiredPESizeCM, 
                        CMTileheight, CMTilewidth, &layerReadLatency, &layerReadDynamicEnergy, &layerWriteLatency, &layerWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &layerbufferLatency, 
                        &layerbufferDynamicEnergy, &layericLatency, &layericDynamicEnergy, &coreLatencyADC, &coreLatencyAccum, &coreLatencyOther, &coreEnergyADC, &coreEnergyAccum, 
                        &coreEnergyOther, false, &layerclkPeriod, &partition_col, &partition_row, tileDup, peDup, subArrayDup, &layer_act_poolLatency, &layerComputationLatency);

            if (partition_col>1 || partition_row>1) {
                cout << "--------Note: accelerator cannot store all the weights of this layer!--------" << endl;
                cout << "------ Row partition of this layer: " << partition_row << " Col partition of this layer: " << partition_col << "------" << endl;
                cout << "-----------------------------------------------------------------------------" << endl;
                cout << endl;
            }
            else {
                cout << "-----------------------------------------------------------------------------" << endl;
                cout << endl;
            }

            // record writing intermediate data to DRAM
            if (i != (netStructure.size()-1)) {
                numWrite_transactions = netStructure[i][0]*netStructure[i][1]*netStructure[i][2]*(param->numBitInput+param->synapseBit);
                numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                timestamp = timestamp+1;

                for (int j=0; j<numWrite_transactions; j++) {
                    string transaction_type = "WRITE";

                    string memory_address = generateNextAddress(address_counter);

                    outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                }
                data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
            }
            
            double layerLeakageEnergy = 0;      
            
            layerLeakageEnergy = (numTileEachLayer[i] * tileLeakageSRAMInUse) * layerReadLatency;

            cout << "layer" << i+1 << "'s inference time is: " << (layerReadLatency+layerWriteLatency)*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readLatency is: " << layerReadLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readDynamicEnergy is: " << layerReadDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s writeLatency is: " << layerWriteLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s writeDynamicEnergy is: " << layerWriteDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s leakagePower is: " << numTileEachLayer[i] * tileLeakage*1e6 << "uW" << endl;
            cout << "layer" << i+1 << "'s leakageEnergy is: " << layerLeakageEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s buffer latency is: " << layerbufferLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s buffer readDynamicEnergy is: " << layerbufferDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s ic latency is: " << layericLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s ic readDynamicEnergy is: " << layericDynamicEnergy*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s computation latency is: " << layerComputationLatency*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s activation and pool latency is: " << layer_act_poolLatency*1e9 << "ns" << endl;
            
            cout << endl;
            cout << "************ Breakdown of Latency and Dynamic Energy *************" << endl;
            cout << endl;

            cout << "ADC (or S/As and precharger for SRAM) readLatency is : " << coreLatencyADC*1e9 << "ns" << endl;
            cout << "Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << coreLatencyAccum*1e9 << "ns" << endl;
            cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << coreLatencyOther*1e9 << "ns" << endl;
            cout << "ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << coreEnergyADC*1e12 << "pJ" << endl;
            cout << "Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << coreEnergyAccum*1e12 << "pJ" << endl;
            cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << coreEnergyOther*1e12 << "pJ" << endl;
            cout << endl;

            chipReadLatency += layerReadLatency;
            chipReadDynamicEnergy += layerReadDynamicEnergy;
            chipWriteLatency += layerWriteLatency;
            chipWriteDynamicEnergy += layerWriteDynamicEnergy;
            chipLeakageEnergy += layerLeakageEnergy;
            chipLeakage += tileLeakage*numTileEachLayer[i];
            chipbufferLatency += layerbufferLatency;
            chipbufferReadDynamicEnergy += layerbufferDynamicEnergy;
            chipicLatency += layericLatency;
            chipicReadDynamicEnergy += layericDynamicEnergy;
            
            chipLatencyADC += coreLatencyADC;
            chipLatencyAccum += coreLatencyAccum;
            chipLatencyOther += coreLatencyOther;
            chipEnergyADC += coreEnergyADC;
            chipEnergyAccum += coreEnergyAccum;
            chipEnergyOther += coreEnergyOther;
        }

        // record writing computation result to DRAM
        int last_layer = netStructure.size()-1;

        numWrite_transactions = netStructure[last_layer][5]*(param->numBitInput+param->synapseBit);
        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

        timestamp = timestamp+1; // Finally, write computation result into DRAM 

        for (int j=0; j<numWrite_transactions; j++) {
            string transaction_type = "WRITE";

            string memory_address = generateNextAddress(address_counter);

            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
        }

        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

        outfile.close();

    }
    else if (full_pipeline) {
        // pipeline system
        // firstly define system clock
        double systemClock = 0;
        
        vector<double> readLatencyPerLayer;
        vector<double> readDynamicEnergyPerLayer;
        vector<double> writeLatencyPerLayer;
        vector<double> writeDynamicEnergyPerLayer;
        vector<double> leakagePowerPerLayer;
        vector<double> bufferLatencyPerLayer;
        vector<double> bufferEnergyPerLayer;
        vector<double> icLatencyPerLayer;
        vector<double> icEnergyPerLayer;
        
        vector<double> coreLatencyADCPerLayer;
        vector<double> coreEnergyADCPerLayer;
        vector<double> coreLatencyAccumPerLayer;
        vector<double> coreEnergyAccumPerLayer;
        vector<double> coreLatencyOtherPerLayer;
        vector<double> coreEnergyOtherPerLayer;

        vector<double> act_poolLatencyPerLayer;
        vector<double> computationLatencyPerLayer;

        int numRead_transactions = 0;
        int numWrite_transactions = 0;
        int dram_bus_width = 128; // this param needs to be changed based on the chosen memory

        // open a file to write data movement trace
        ofstream outfile; 
        outfile.open("data_movement-FPP.trace", ios::out | ios::trunc);

        if (!outfile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return 1;
        }
        
        int timestamp;
        int address_counter = 0;

        for (int i=0; i<netStructure.size(); i++) {

            numComputationPerLayer = 2*(netStructure[i][0] * netStructure[i][1] * netStructure[i][2] * netStructure[i][3] * netStructure[i][4] * netStructure[i][5]);
            numComputation = numComputation + numComputationPerLayer;

            // record loading weights from off-chip memory
            numRead_transactions = netStructure[i][2]*netStructure[i][3]*netStructure[i][4]*netStructure[i][5]*param->synapseBit;
            numRead_transactions = ceil(numRead_transactions/dram_bus_width);

            timestamp = 0; // At first, read all weigths from DRAM

            for (int k=0; k<speedUpEachLayer[i]; k++) {
                
                if (k != 0) {
                    address_counter = address_counter-numRead_transactions;
                }

                for (int j=0; j<numRead_transactions; j++) {

                    string transaction_type = "READ";

                    string memory_address = generateNextAddress(address_counter);

                    outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                }
                data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
            }

            ChipCalculatePerformance(inputParameter, tech, cell, i, argv[2*i+6], argv[2*i+6], argv[2*i+7], netStructure[i][6], netStructure, numTileEachLayer, speedUpEachLayer, tileLocaEachLayer, 
                                    desiredTileSizeCM, desiredPESizeCM, CMTileheight, CMTilewidth, &layerReadLatency, &layerReadDynamicEnergy, &layerWriteLatency, &layerWriteDynamicEnergy, &tileLeakage, 
                                    &tileLeakageSRAMInUse, &layerbufferLatency, &layerbufferDynamicEnergy, &layericLatency, &layericDynamicEnergy, &coreLatencyADC, &coreLatencyAccum, &coreLatencyOther, 
                                    &coreEnergyADC, &coreEnergyAccum, &coreEnergyOther, false, &layerclkPeriod, &partition_col, &partition_row, tileDup, peDup, subArrayDup, &layer_act_poolLatency, 
                                    &layerComputationLatency);
            
            systemClock = MAX(systemClock, layerReadLatency);
            
            readLatencyPerLayer.push_back(layerReadLatency);
            readDynamicEnergyPerLayer.push_back(layerReadDynamicEnergy);
            writeLatencyPerLayer.push_back(layerWriteLatency);
            writeDynamicEnergyPerLayer.push_back(layerWriteDynamicEnergy);
            leakagePowerPerLayer.push_back(numTileEachLayer[i] * (tileLeakage * (systemClock-readLatencyPerLayer[i]) + tileLeakageSRAMInUse * readLatencyPerLayer[i]) / systemClock);
            bufferLatencyPerLayer.push_back(layerbufferLatency);
            bufferEnergyPerLayer.push_back(layerbufferDynamicEnergy);
            icLatencyPerLayer.push_back(layericLatency);
            icEnergyPerLayer.push_back(layericDynamicEnergy);
            
            coreLatencyADCPerLayer.push_back(coreLatencyADC);
            coreEnergyADCPerLayer.push_back(coreEnergyADC);
            coreLatencyAccumPerLayer.push_back(coreLatencyAccum);
            coreEnergyAccumPerLayer.push_back(coreEnergyAccum);
            coreLatencyOtherPerLayer.push_back(coreLatencyOther);
            coreEnergyOtherPerLayer.push_back(coreEnergyOther);

            act_poolLatencyPerLayer.push_back(layer_act_poolLatency);
            computationLatencyPerLayer.push_back(layerComputationLatency);
        }

        // record loading IFM from off-chip memory
        numRead_transactions = param->batch_size*netStructure[0][0]*netStructure[0][1]*netStructure[0][2]*netStructure[0][3]*netStructure[0][4]*param->numBitInput;
        numRead_transactions = ceil(numRead_transactions/dram_bus_width);

        for (int j=0; j<numRead_transactions; j++) {

            string transaction_type = "READ";

            string memory_address = generateNextAddress(address_counter);

            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
        }

        data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;

        // record writing computation result to DRAM
        int last_layer = netStructure.size()-1;

        numWrite_transactions = param->batch_size*netStructure[last_layer][5]*(param->numBitInput+param->synapseBit);
        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

        timestamp = 1; // Finally, write computation result into DRAM 

        for (int j=0; j<numWrite_transactions; j++) {
            string transaction_type = "WRITE";

            string memory_address = generateNextAddress(address_counter);

            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
        }

        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

        outfile.close();

        for (int i=0; i<netStructure.size(); i++) {
            
            cout << "-----------------------------------------------------------------------------" << endl;
            cout << "-------------------- Estimation of Layer " << i+1 << " ----------------------" << endl;
            cout << "-----------------------------------------------------------------------------" << endl;
            cout << endl;

            cout << "layer" << i+1 << "'s inference time is: " << (readLatencyPerLayer[i]+writeLatencyPerLayer[i])*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readLatency is: " << readLatencyPerLayer[i]*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s readDynamicEnergy is: " << readDynamicEnergyPerLayer[i]*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s writeLatency is: " << writeLatencyPerLayer[i]*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s writeDynamicEnergy is: " << writeDynamicEnergyPerLayer[i]*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s leakagePower is: " << leakagePowerPerLayer[i]*1e6 << "uW" << endl;
            cout << "layer" << i+1 << "'s leakageEnergy is: " << leakagePowerPerLayer[i] * systemClock *1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s buffer latency is: " << bufferLatencyPerLayer[i]*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s buffer readDynamicEnergy is: " << bufferEnergyPerLayer[i]*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s ic latency is: " << icLatencyPerLayer[i]*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s ic readDynamicEnergy is: " << icEnergyPerLayer[i]*1e12 << "pJ" << endl;
            cout << "layer" << i+1 << "'s computation latency is: " << computationLatencyPerLayer[i]*1e9 << "ns" << endl;
            cout << "layer" << i+1 << "'s activation and pool latency is: " << act_poolLatencyPerLayer[i]*1e9 << "ns" << endl;

            cout << endl;
            cout << "************ Breakdown of Latency and Dynamic Energy *************" << endl;
            cout << endl;
            cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << coreLatencyADCPerLayer[i]*1e9 << "ns" << endl;
            cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << coreLatencyAccumPerLayer[i]*1e9 << "ns" << endl;
            cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << coreLatencyOtherPerLayer[i]*1e9 << "ns" << endl;
            cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << coreEnergyADCPerLayer[i]*1e12 << "pJ" << endl;
            cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << coreEnergyAccumPerLayer[i]*1e12 << "pJ" << endl;
            cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << coreEnergyOtherPerLayer[i]*1e12 << "pJ" << endl;
            cout << endl;
            
            chipReadLatency = systemClock;
            chipReadDynamicEnergy += readDynamicEnergyPerLayer[i];
            chipWriteLatency += writeLatencyPerLayer[i];
            chipWriteDynamicEnergy += writeDynamicEnergyPerLayer[i];
            chipLeakageEnergy += leakagePowerPerLayer[i] * systemClock;
            chipLeakage += leakagePowerPerLayer[i];
            chipbufferLatency = MAX(chipbufferLatency, bufferLatencyPerLayer[i]);
            chipbufferReadDynamicEnergy += bufferEnergyPerLayer[i];
            chipicLatency = MAX(chipicLatency, icLatencyPerLayer[i]);
            chipicReadDynamicEnergy += icEnergyPerLayer[i];
            
            chipLatencyADC = MAX(chipLatencyADC, coreLatencyADCPerLayer[i]);
            chipLatencyAccum = MAX(chipLatencyAccum, coreLatencyAccumPerLayer[i]);
            chipLatencyOther = MAX(chipLatencyOther, coreLatencyOtherPerLayer[i]);
            chipEnergyADC += coreEnergyADCPerLayer[i];
            chipEnergyAccum += coreEnergyAccumPerLayer[i];
            chipEnergyOther += coreEnergyOtherPerLayer[i];
        }
    }
    else if (partial_pipeline) {
        
        int numRead_transactions = 0;
        int numWrite_transactions = 0;
        int dram_bus_width = 128; // this param needs to be changed based on the chosen memory

        // open a file to write data movement trace
        ofstream outfile; 
        outfile.open("data_movement-PPP.trace", ios::out | ios::trunc);

        if (!outfile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return 1;
        }

        int timestamp=-1;
        int address_counter = 0;
        int flag = 0;

        for (int i=0; i<record.size(); i++) {

            timestamp = timestamp+1;

            // record loading IFM from off-chip memory
            numRead_transactions = param->batch_size*netStructure[flag][0]*netStructure[flag][1]*netStructure[flag][2]*netStructure[flag][3]*netStructure[flag][4]*(flag==0 ? param->numBitInput:(param->numBitInput+param->synapseBit));
            numRead_transactions = ceil(numRead_transactions/dram_bus_width);

            address_counter = address_counter - numWrite_transactions;

            for (int j=0; j<numRead_transactions; j++) {

                string transaction_type = "READ";

                string memory_address = generateNextAddress(address_counter);

                outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
            }

            data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;

            for (int j=flag; j<record[i].size(); j++) {
                flag = j+1;

                if (partition_col_array[j]>1 || partition_row_array[j]>1) {
                    // if partition_row_array[j]>1, there are additional (read(IFM))*(partition_row_array[j]-1) operations.
                    for (int k=0; k<(partition_row_array[j]-1); k++) {

                        address_counter = address_counter - numRead_transactions;

                        for (int z=0; z<numRead_transactions; z++) {

                            string transaction_type = "READ";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }
                        data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
                    }

                    // if partition_col_array[j]>1, there are additional (write+read(this layer computation result))*(partition_col_array[j]-1) operations.
                    for (int k=0; k<(partition_col_array[j]-1); k++) {
                        numWrite_transactions = param->batch_size*netStructure[j][0]*netStructure[j][1]*netStructure[j][2]*(param->numBitInput+param->synapseBit);
                        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                        for (int z=0; z<numWrite_transactions; z++) {
                            string transaction_type = "WRITE";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }

                        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

                        address_counter = address_counter - numWrite_transactions;

                        for (int z=0; z<numWrite_transactions; z++) {
                            string transaction_type = "READ";

                            string memory_address = generateNextAddress(address_counter);

                            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                        }
                        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
                    }
                }

                // record loading weights from off-chip memory
                numRead_transactions =  netStructure[j][2]*netStructure[j][3]*netStructure[j][4]*netStructure[j][5]*param->synapseBit;
                numRead_transactions = ceil(numRead_transactions/dram_bus_width);

                for (int k=0; k<speedUpEachLayer[j]; k++) {
                    if (k != 0) {
                        address_counter = address_counter-numRead_transactions;
                    }

                    for (int z=0; z<numRead_transactions; z++) {

                        string transaction_type = "READ";

                        string memory_address = generateNextAddress(address_counter);

                        outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                    }
                    data_transfer_num = data_transfer_num + numRead_transactions*dram_bus_width;
                }
            }

            // record writing intermediate data to DRAM
            if (i != (record.size()-1)) {
                numWrite_transactions = param->batch_size*netStructure[flag-1][0]*netStructure[flag-1][1]*netStructure[flag-1][2]*(param->numBitInput+param->synapseBit);
                numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

                timestamp = timestamp+1;

                for (int j=0; j<numWrite_transactions; j++) {
                    string transaction_type = "WRITE";

                    string memory_address = generateNextAddress(address_counter);

                    outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
                }
                data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;
            }
        }

        // record writing computation result to DRAM
        int last_layer = netStructure.size()-1;

        numWrite_transactions = param->batch_size*netStructure[last_layer][5]*(param->numBitInput+param->synapseBit);
        numWrite_transactions = ceil(numWrite_transactions/dram_bus_width);

        timestamp = timestamp+1; // Finally, write computation result into DRAM 

        for (int i=0; i<numWrite_transactions; i++) {
            string transaction_type = "WRITE";

            string memory_address = generateNextAddress(address_counter);

            outfile << timestamp << "," << transaction_type << "," << memory_address << endl;
        }

        data_transfer_num = data_transfer_num + numWrite_transactions*dram_bus_width;

        outfile.close();

        // partial pipeline system
        // firstly define system clock
        double systemClock = 0;
        
        vector<double> readLatencyPerLayer;
        vector<double> readDynamicEnergyPerLayer;
        vector<double> writeLatencyPerLayer;
        vector<double> writeDynamicEnergyPerLayer;
        vector<double> leakagePowerPerLayer;
        vector<double> bufferLatencyPerLayer;
        vector<double> bufferEnergyPerLayer;
        vector<double> icLatencyPerLayer;
        vector<double> icEnergyPerLayer;
        
        vector<double> coreLatencyADCPerLayer;
        vector<double> coreEnergyADCPerLayer;
        vector<double> coreLatencyAccumPerLayer;
        vector<double> coreEnergyAccumPerLayer;
        vector<double> coreLatencyOtherPerLayer;
        vector<double> coreEnergyOtherPerLayer;

        vector<double> act_poolLatencyPerLayer;
        vector<double> computationLatencyPerLayer;

        double bufferLatency    = 0;
        double icLatency        = 0;
        double LatencyADC       = 0;
        double LatencyAccum     = 0;
        double LatencyOther     = 0;

        int temp1 = 0;
        int temp2 = 0;

        if (record.size() > 1) {
            for (int i=0; i<record.size(); i++) {
                for (int k=temp1; k<record[i].size(); k++) {
                    numComputationPerLayer = param->batch_size*2*(netStructure[k][0] * netStructure[k][1] * netStructure[k][2] * netStructure[k][3] * netStructure[k][4] * netStructure[k][5]);
                    numComputation = numComputation + numComputationPerLayer;

                    ChipCalculatePerformance(inputParameter, tech, cell, k, argv[2*k+6], argv[2*k+6], argv[2*k+7], netStructure[k][6],
                                netStructure, numTileEachLayer, speedUpEachLayer, tileLocaEachLayer, desiredTileSizeCM, desiredPESizeCM, CMTileheight, CMTilewidth, &layerReadLatency, 
                                &layerReadDynamicEnergy, &layerWriteLatency, &layerWriteDynamicEnergy, &tileLeakage, &tileLeakageSRAMInUse, &layerbufferLatency, &layerbufferDynamicEnergy, 
                                &layericLatency, &layericDynamicEnergy, &coreLatencyADC, &coreLatencyAccum, &coreLatencyOther, &coreEnergyADC, &coreEnergyAccum, &coreEnergyOther, false, 
                                &layerclkPeriod, &partition_col, &partition_row, tileDup, peDup, subArrayDup, &layer_act_poolLatency, &layerComputationLatency);
                    if (param->synchronous) {
                        layerReadLatency *= clkPeriod;
                        layerbufferLatency *= clkPeriod;
                        layericLatency *= clkPeriod;
                        coreLatencyADC *= clkPeriod;
                        coreLatencyAccum *= clkPeriod;
                        coreLatencyOther *= clkPeriod;
                    }           
                    
                    systemClock = MAX(systemClock, layerReadLatency);
                    
                    readLatencyPerLayer.push_back(layerReadLatency);
                    readDynamicEnergyPerLayer.push_back(layerReadDynamicEnergy);
                    writeLatencyPerLayer.push_back(layerWriteLatency);
                    writeDynamicEnergyPerLayer.push_back(layerWriteDynamicEnergy);
                    leakagePowerPerLayer.push_back(numTileEachLayer[k] * (tileLeakage * (systemClock-readLatencyPerLayer[k]) + tileLeakageSRAMInUse * readLatencyPerLayer[k]) / systemClock);
                    bufferLatencyPerLayer.push_back(layerbufferLatency);
                    bufferEnergyPerLayer.push_back(layerbufferDynamicEnergy);
                    icLatencyPerLayer.push_back(layericLatency);
                    icEnergyPerLayer.push_back(layericDynamicEnergy);
                    
                    coreLatencyADCPerLayer.push_back(coreLatencyADC);
                    coreEnergyADCPerLayer.push_back(coreEnergyADC);
                    coreLatencyAccumPerLayer.push_back(coreLatencyAccum);
                    coreEnergyAccumPerLayer.push_back(coreEnergyAccum);
                    coreLatencyOtherPerLayer.push_back(coreLatencyOther);
                    coreEnergyOtherPerLayer.push_back(coreEnergyOther);

                    act_poolLatencyPerLayer.push_back(layer_act_poolLatency);
                    computationLatencyPerLayer.push_back(layerComputationLatency);

                    temp1 = k+1;
                }

                for (int j=temp2; j<record[i].size(); j++) {
                    cout << "-----------------------------------------------------------------------------" << endl;
                    cout << "-------------------- Estimation of Layer " << j+1 << " ----------------------" << endl;
                    //cout << "-----------------------------------------------------------------------------" << endl;
                    //cout << endl;   
                    if (partition_col_array[j]>1 || partition_row_array[j]>1) {
                        cout << "--------Note: accelerator cannot store all the weights of this layer!--------" << endl;
                        cout << "------ Row partition of this layer: " << partition_row_array[j] << " Col partition of this layer: " << partition_col_array[j] << "------" << endl;
                        cout << "-----------------------------------------------------------------------------" << endl;
                        cout << endl;
                    }
                    else {
                        cout << "-----------------------------------------------------------------------------" << endl;
                        cout << endl;
                    }
                    
                    cout << "layer" << j+1 << "'s inference time is: " << (readLatencyPerLayer[j]+writeLatencyPerLayer[j])*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s readLatency is: " << readLatencyPerLayer[j]*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s readDynamicEnergy is: " << readDynamicEnergyPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "layer" << j+1 << "'s writeLatency is: " << writeLatencyPerLayer[j]*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s writeDynamicEnergy is: " << writeDynamicEnergyPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "layer" << j+1 << "'s leakagePower is: " << leakagePowerPerLayer[j]*1e6 << "uW" << endl;
                    cout << "layer" << j+1 << "'s leakageEnergy is: " << leakagePowerPerLayer[j] * systemClock *1e12 << "pJ" << endl;
                    cout << "layer" << j+1 << "'s buffer latency is: " << bufferLatencyPerLayer[j]*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s buffer readDynamicEnergy is: " << bufferEnergyPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "layer" << j+1 << "'s ic latency is: " << icLatencyPerLayer[j]*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s ic readDynamicEnergy is: " << icEnergyPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "layer" << j+1 << "'s computation latency is: " << computationLatencyPerLayer[j]*1e9 << "ns" << endl;
                    cout << "layer" << j+1 << "'s activation and pool latency is: " << act_poolLatencyPerLayer[j]*1e9 << "ns" << endl;

                    cout << endl;
                    cout << "************ Breakdown of Latency and Dynamic Energy *************" << endl;
                    cout << endl;
                    cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << coreLatencyADCPerLayer[j]*1e9 << "ns" << endl;
                    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << coreLatencyAccumPerLayer[j]*1e9 << "ns" << endl;
                    cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << coreLatencyOtherPerLayer[j]*1e9 << "ns" << endl;
                    cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << coreEnergyADCPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << coreEnergyAccumPerLayer[j]*1e12 << "pJ" << endl;
                    cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << coreEnergyOtherPerLayer[j]*1e12 << "pJ" << endl;
                    cout << endl;
                    
                    //chipReadLatency = systemClock;
                    chipReadDynamicEnergy += readDynamicEnergyPerLayer[j];
                    chipWriteLatency += writeLatencyPerLayer[j];
                    chipWriteDynamicEnergy += writeDynamicEnergyPerLayer[j];
                    chipLeakageEnergy += leakagePowerPerLayer[j] * systemClock;
                    chipLeakage += leakagePowerPerLayer[j];
                    bufferLatency = MAX(bufferLatency, bufferLatencyPerLayer[j]);
                    chipbufferReadDynamicEnergy += bufferEnergyPerLayer[j];
                    icLatency = MAX(icLatency, icLatencyPerLayer[j]);
                    chipicReadDynamicEnergy += icEnergyPerLayer[j];
                    
                    LatencyADC = MAX(LatencyADC, coreLatencyADCPerLayer[j]);
                    LatencyAccum = MAX(LatencyAccum, coreLatencyAccumPerLayer[j]);
                    LatencyOther = MAX(LatencyOther, coreLatencyOtherPerLayer[j]);
                    chipEnergyADC += coreEnergyADCPerLayer[j];
                    chipEnergyAccum += coreEnergyAccumPerLayer[j];
                    chipEnergyOther += coreEnergyOtherPerLayer[j];

                    temp2 = j+1;
                }

                //cout << "systemClock: " << systemClock << endl;
                chipReadLatency   += systemClock  ;
                chipbufferLatency += bufferLatency;
                chipicLatency     += icLatency    ;
                chipLatencyADC    += LatencyADC   ;
                chipLatencyAccum  += LatencyAccum ;
                chipLatencyOther  += LatencyOther ;

                systemClock   = 0; 
                bufferLatency = 0;
                icLatency     = 0;
                LatencyADC    = 0;
                LatencyAccum  = 0;
                LatencyOther  = 0;

            }
        }
    }

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "--------------------------------- Summary -----------------------------------" <<  endl;
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << endl;

    //cout << "Computation number: " << numComputation << endl;

    if (param->mode) {
        off_chip_dramLatency = data_transfer_num/((param->memory_bus_bandwidth)*1e9*8);
        off_chip_dramEnergy = 2950035190.37; // unit: pJ
    }
    else {
        off_chip_dramLatency = 0;
        off_chip_dramEnergy  = 0;
    }

    if (raw_processing) {
        if (record.size()>=2) {
            cout << "Chip raw-processing-system-total-time (batch) is: " << (param->batch_size*chipReadLatency+param->batch_size*chipWriteLatency+off_chip_dramLatency)*1e9 << "ns" << endl;
        }
        else {
            cout << "Chip raw-processing-system-total-time (per image) is: " << (chipReadLatency+chipWriteLatency+off_chip_dramLatency/param->batch_size)*1e9 << "ns" << endl;                                                                                                     
        }
        cout << "Chip raw-processing-system-total-energy (per image) is: " << (((record.size()>=2)?param->batch_size:1)*chipReadDynamicEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipLeakageEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipWriteDynamicEnergy*1e12+((record.size()>=2)?off_chip_dramEnergy:(off_chip_dramEnergy/param->batch_size)))/(((record.size()>=2)?param->batch_size:1)) << "pJ" << endl;
        cout << "Chip raw-processing-system readLatency (per image) is: " << chipReadLatency*1e9 << "ns" << endl;
        cout << "Chip raw-processing-system readDynamicEnergy (per image) is: " << chipReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip raw-processing-system writeLatency (per image) is: " << chipWriteLatency*1e9 << "ns" << endl;
        cout << "Chip raw-processing-system writeDynamicEnergy (per image) is: " << chipWriteDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip raw-processing-system leakage Energy (per image) is: " << chipLeakageEnergy*1e12 << "pJ" << endl;
        cout << "Chip raw-processing-system leakage Power (per image) is: " << chipLeakage*1e6 << "uW" << endl;
        cout << "Chip raw-processing-system buffer readLatency (per image) is: " << chipbufferLatency*1e9 << "ns" << endl;
        cout << "Chip raw-processing-system buffer readDynamicEnergy (per image) is: " << chipbufferReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip raw-processing-system ic readLatency (per image) is: " << chipicLatency*1e9 << "ns" << endl;
        cout << "Chip raw-processing-system ic readDynamicEnergy (per image) is: " << chipicReadDynamicEnergy*1e12 << "pJ" << endl;
        if (param->mode) {
            cout << "Chip raw-processing-system off-chip DRAM latency (per image) is: " << off_chip_dramLatency*1e9/(param->batch_size) << "ns" << endl;
            cout << "Chip raw-processing-system off-chip DRAM energy (per image) is: " << off_chip_dramEnergy/(param->batch_size) << "pJ" << endl;
        }

        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << chipLatencyADC*1e9 << "ns" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << chipLatencyAccum*1e9 << "ns" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << chipLatencyOther*1e9 << "ns" << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << chipEnergyADC*1e12 << "pJ" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << chipEnergyAccum*1e12 << "pJ" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << chipEnergyOther*1e12 << "pJ" << endl;
        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;

        if (param->validated) {
            cout << "Energy Efficiency TOPS/W (raw_processing): " << numComputation/(((record.size()>=2)?param->batch_size:1)*chipReadDynamicEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipLeakageEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipWriteDynamicEnergy*1e12+((record.size()>=2)?off_chip_dramEnergy:(off_chip_dramEnergy/param->batch_size)))/param->zeta << endl; // post-layout energy increase, zeta = 1.23 by default
        }
        else {
            cout << "Energy Efficiency TOPS/W (raw_processing): " << numComputation/(((record.size()>=2)?param->batch_size:1)*chipReadDynamicEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipLeakageEnergy*1e12+((record.size()>=2)?param->batch_size:1)*chipWriteDynamicEnergy*1e12+((record.size()>=2)?off_chip_dramEnergy:(off_chip_dramEnergy/param->batch_size))) << endl;
        }

        cout << "Throughput TOPS (raw_processing): " << numComputation/(((record.size()>=2)?param->batch_size:1)*chipReadLatency*1e12+((record.size()>=2)?param->batch_size:1)*chipWriteLatency*1e12+((record.size()>=2)?off_chip_dramLatency:(off_chip_dramLatency/param->batch_size))*1e12) << endl; 
        cout << "Throughput FPS (raw_processing): " << ((record.size()>=2)?param->batch_size:1)*(1/(((record.size()>=2)?param->batch_size:1)*chipReadLatency+((record.size()>=2)?param->batch_size:1)*chipWriteLatency+((record.size()>=2)?off_chip_dramLatency:(off_chip_dramLatency/param->batch_size)))) << endl;
        cout << "Compute efficiency TOPS/mm^2 (raw_processing): " << numComputation/(((record.size()>=2)?param->batch_size:1)*chipReadLatency*1e12+((record.size()>=2)?param->batch_size:1)*chipWriteLatency*1e12+((record.size()>=2)?off_chip_dramLatency:(off_chip_dramLatency/param->batch_size))*1e12)/(chipArea*1e6) << endl;
        cout << endl;
    }
    else if (layer_by_layer) {
        cout << "Chip layer-based-system-total-time (per image) is: " << (chipReadLatency+chipWriteLatency+off_chip_dramLatency)*1e9 << "ns" << endl;
        cout << "Chip layer-based-system-total-energy (per image) is: " << (chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy) << "pJ" << endl;
        cout << "Chip layer-based-system readLatency (per image) is: " << chipReadLatency*1e9 << "ns" << endl;
        cout << "Chip layer-based-system readDynamicEnergy (per image) is: " << chipReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip layer-based-system writeLatency (per image) is: " << chipWriteLatency*1e9 << "ns" << endl;
        cout << "Chip layer-based-system writeDynamicEnergy (per image) is: " << chipWriteDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip layer-based-system leakage Energy (per image) is: " << chipLeakageEnergy*1e12 << "pJ" << endl;
        cout << "Chip layer-based-system leakage Power (per image) is: " << chipLeakage*1e6 << "uW" << endl;
        cout << "Chip layer-based-system buffer readLatency (per image) is: " << chipbufferLatency*1e9 << "ns" << endl;
        cout << "Chip layer-based-system buffer readDynamicEnergy (per image) is: " << chipbufferReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip layer-based-system ic readLatency (per image) is: " << chipicLatency*1e9 << "ns" << endl;
        cout << "Chip layer-based-system ic readDynamicEnergy (per image) is: " << chipicReadDynamicEnergy*1e12 << "pJ" << endl;
        if (param->mode) {
            cout << "Chip layer-based-system off-chip DRAM latency (per image) is: " << off_chip_dramLatency*1e9 << "ns" << endl;
            cout << "Chip layer-based-system off-chip DRAM energy (per image) is: " << off_chip_dramEnergy << "pJ" << endl;
        }

        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << chipLatencyADC*1e9 << "ns" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << chipLatencyAccum*1e9 << "ns" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << chipLatencyOther*1e9 << "ns" << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << chipEnergyADC*1e12 << "pJ" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << chipEnergyAccum*1e12 << "pJ" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << chipEnergyOther*1e12 << "pJ" << endl;
        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;

        if (param->validated) {
            cout << "Energy Efficiency TOPS/W (layer-based mapper): " << numComputation/(chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy)/param->zeta << endl; // post-layout energy increase, zeta = 1.23 by default
        }
        else {
            cout << "Energy Efficiency TOPS/W (layer-based mapper): " << numComputation/(chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy) << endl;
        }

        cout << "Throughput TOPS (layer-based mapper): " << numComputation/(chipReadLatency*1e12+chipWriteLatency*1e12+off_chip_dramLatency*1e12) << endl;
        cout << "Throughput FPS (layer-based mapper): " << (1/(chipReadLatency+chipWriteLatency+off_chip_dramLatency)) << endl;
        cout << "Compute efficiency TOPS/mm^2 (layer-based mapper): " << numComputation/(chipReadLatency*1e12+chipWriteLatency*1e12+off_chip_dramLatency*1e12)/(chipArea*1e6) << endl;
        cout << endl;
    }
    else if (full_pipeline) {
        cout << "Chip full-pipeline-system-total-time (per image) is: " << (chipReadLatency+chipWriteLatency+off_chip_dramLatency/param->batch_size)*1e9 << "ns" << endl;
        cout << "Chip full-pipeline-system-total-energy (per image) is: " << (chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy/param->batch_size) << "pJ" << endl;
        cout << "Chip full-pipeline-system readLatency (per image) is: " << chipReadLatency*1e9 << "ns" << endl;
        cout << "Chip full-pipeline-system readDynamicEnergy (per image) is: " << chipReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip full-pipeline-system writeLatency (per image) is: " << chipWriteLatency*1e9 << "ns" << endl;
        cout << "Chip full-pipeline-system writeDynamicEnergy (per image) is: " << chipWriteDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip full-pipeline-system leakage Energy (per image) is: " << chipLeakageEnergy*1e12 << "pJ" << endl;
        cout << "Chip full-pipeline-system leakage Power (per image) is: " << chipLeakage*1e6 << "uW" << endl;
        cout << "Chip full-pipeline-system buffer readLatency (per image) is: " << chipbufferLatency*1e9 << "ns" << endl;
        cout << "Chip full-pipeline-system buffer readDynamicEnergy (per image) is: " << chipbufferReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip full-pipeline-system ic readLatency (per image) is: " << chipicLatency*1e9 << "ns" << endl;
        cout << "Chip full-pipeline-system ic readDynamicEnergy (per image) is: " << chipicReadDynamicEnergy*1e12 << "pJ" << endl;
        if (param->mode) {
            cout << "Chip full-pipeline-system off-chip DRAM latency (per image) is: " << off_chip_dramLatency*1e9/param->batch_size << "ns" << endl;
            cout << "Chip full-pipeline-system off-chip DRAM energy (per image) is: " << off_chip_dramEnergy/param->batch_size << "pJ" << endl;
        }

        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << chipLatencyADC*1e9 << "ns" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << chipLatencyAccum*1e9 << "ns" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << chipLatencyOther*1e9 << "ns" << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << chipEnergyADC*1e12 << "pJ" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << chipEnergyAccum*1e12 << "pJ" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << chipEnergyOther*1e12 << "pJ" << endl;
        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;

        if (param->validated) {
            cout << "Energy Efficiency TOPS/W (full pipelined mapper): " << numComputation/(chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy/param->batch_size)/param->zeta << endl; // post-layout energy increase, zeta = 1.23 by default
        }
        else {
            cout << "Energy Efficiency TOPS/W (full pipelined mapper): " << numComputation/(chipReadDynamicEnergy*1e12+chipLeakageEnergy*1e12+chipWriteDynamicEnergy*1e12+off_chip_dramEnergy/param->batch_size) << endl;
        }

        cout << "Throughput TOPS (full pipelined mapper): " << numComputation/(chipReadLatency*1e12+chipWriteLatency*1e12+off_chip_dramLatency*1e12/param->batch_size) << endl;
        cout << "Throughput FPS (full pipelined mapper): " << (1/(chipReadLatency+chipWriteLatency+off_chip_dramLatency/param->batch_size)) << endl;
        cout << "Compute efficiency TOPS/mm^2 (full pipelined mapper): " << numComputation/(chipReadLatency*1e12+chipWriteLatency*1e12+off_chip_dramLatency*1e12/param->batch_size)/(chipArea*1e6) << endl;
        cout << endl;
    }
    else if (partial_pipeline) {
        cout << "Chip partial-pipeline-system-total-time (batch) is: " << (param->batch_size*chipReadLatency+param->batch_size*chipWriteLatency+off_chip_dramLatency)*1e9 << "ns" << endl;
        cout << "Chip partial-pipeline-system-total-energy (batch) is: " << (param->batch_size*chipReadDynamicEnergy*1e12+param->batch_size*chipLeakageEnergy*1e12+param->batch_size*chipWriteDynamicEnergy*1e12+off_chip_dramEnergy) << "pJ" << endl;
        cout << "Chip partial-pipeline-system readLatency (per image) is: " << chipReadLatency*1e9 << "ns" << endl;
        cout << "Chip partial-pipeline-system readDynamicEnergy (per image) is: " << chipReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip partial-pipeline-system writeLatency (per image) is: " << chipWriteLatency*1e9 << "ns" << endl;
        cout << "Chip partial-pipeline-system writeDynamicEnergy (per image) is: " << chipWriteDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip partial-pipeline-system leakage Energy (per image) is: " << chipLeakageEnergy*1e12 << "pJ" << endl;
        cout << "Chip partial-pipeline-system leakage Power (per image) is: " << chipLeakage*1e6 << "uW" << endl;
        cout << "Chip partial-pipeline-system buffer readLatency (per image) is: " << chipbufferLatency*1e9 << "ns" << endl;
        cout << "Chip partial-pipeline-system buffer readDynamicEnergy (per image) is: " << chipbufferReadDynamicEnergy*1e12 << "pJ" << endl;
        cout << "Chip partial-pipeline-system ic readLatency (per image) is: " << chipicLatency*1e9 << "ns" << endl;
        cout << "Chip partial-pipeline-system ic readDynamicEnergy (per image) is: " << chipicReadDynamicEnergy*1e12 << "pJ" << endl;
        if (param->mode) {
            cout << "Chip partial-pipeline-system off-chip DRAM latency (batch) is: " << off_chip_dramLatency*1e9 << "ns" << endl;
            cout << "Chip partial-pipeline-system off-chip DRAM energy (batch) is: " << off_chip_dramEnergy << "pJ" << endl;
        }

        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << chipLatencyADC*1e9 << "ns" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << chipLatencyAccum*1e9 << "ns" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readLatency is : " << chipLatencyOther*1e9 << "ns" << endl;
        cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << chipEnergyADC*1e12 << "pJ" << endl;
        cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << chipEnergyAccum*1e12 << "pJ" << endl;
        cout << "----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, IC, pooling and activation units) readDynamicEnergy is : " << chipEnergyOther*1e12 << "pJ" << endl;
        cout << endl;
        cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
        cout << endl;

         if (param->validated) {
            cout << "Energy Efficiency TOPS/W (partial pipelined mapper-batch): " << numComputation/(param->batch_size*chipReadDynamicEnergy*1e12+param->batch_size*chipLeakageEnergy*1e12+param->batch_size*chipWriteDynamicEnergy*1e12+off_chip_dramEnergy)/param->zeta << endl; // post-layout energy increase, zeta = 1.23 by default
        }
        else {
            cout << "Energy Efficiency TOPS/W (partial pipelined mapper-batch): " << numComputation/(param->batch_size*chipReadDynamicEnergy*1e12+param->batch_size*chipLeakageEnergy*1e12+param->batch_size*chipWriteDynamicEnergy*1e12+off_chip_dramEnergy) << endl;
        }

        cout << "Throughput TOPS (partial pipelined mapper-batch): " << numComputation/(param->batch_size*chipReadLatency*1e12+param->batch_size*chipWriteLatency*1e12+off_chip_dramLatency*1e12) << endl;
        cout << "Throughput FPS (partial pipelined mapper-batch): " << param->batch_size*(1/(param->batch_size*chipReadLatency+param->batch_size*chipWriteLatency+off_chip_dramLatency)) << endl;
        cout << "Compute efficiency TOPS/mm^2 (partial pipelined mapper-batch): " << numComputation/(param->batch_size*chipReadLatency*1e12+param->batch_size*chipWriteLatency*1e12+off_chip_dramLatency*1e12)/(chipArea*1e6) << endl;
        cout << endl;
    }

    cout << "-------------------------------------- Hardware Performance Done --------------------------------------" <<  endl;
    cout << endl;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop-start);

    cout << "------------------------------ Simulation Performance --------------------------------" <<  endl;
    cout << "Total Run-time of NeuroSim: " << duration.count() << " seconds" << endl;
    cout << endl;

    return 0;
}

vector<vector<double> > getNetStructure(const string &inputfile) {
    ifstream infile(inputfile.c_str());      
    string inputline;
    string inputval;
    
    int ROWin=0, COLin=0;      
    if (!infile.good()) {        
        cerr << "Error: the input file cannot be opened!" << endl;
        exit(1);
    }
    else {
        while (getline(infile, inputline, '\n')) {       
            ROWin++;                                
        }
        infile.clear();
        infile.seekg(0, ios::beg);      
        if (getline(infile, inputline, '\n')) {        
            istringstream iss (inputline);      
            while (getline(iss, inputval, ',')) {       
                COLin++;
            }
        }   
    }
    infile.clear();
    infile.seekg(0, ios::beg);          

    vector<vector<double> > netStructure;               
    for (int row=0; row<ROWin; row++) { 
        vector<double> netStructurerow;
        getline(infile, inputline, '\n');             
        istringstream iss;
        iss.str(inputline);
        for (int col=0; col<COLin; col++) {       
            while(getline(iss, inputval, ',')){ 
                istringstream fs;
                fs.str(inputval);
                double f=0;
                fs >> f;                
                netStructurerow.push_back(f);           
            }           
        }       
        netStructure.push_back(netStructurerow);
    }
    infile.close();
    
    return netStructure;
}   




