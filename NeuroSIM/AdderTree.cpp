#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "AdderTree.h"

using namespace std;
// Anni update
extern Param *param;

AdderTree::AdderTree(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), adder(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void AdderTree::Initialize(int _numSubcoreRow, int _numAdderBit, int _numAdderTree, double _clkFreq) {
	if (initialized)
		cout << "[AdderTree] Warning: Already initialized!" << endl;
	
	numSubcoreRow = _numSubcoreRow;                  // # of row of subcore in the synaptic core
	numStage = ceil(log2(numSubcoreRow));            // # of stage of the adder tree, used for CalculateLatency ...
	numAdderBit = _numAdderBit;                      // # of input bits of the Adder
	numAdderTree = _numAdderTree;                    // # of Adder Tree
	clkFreq = _clkFreq;
	
	initialized = true;
}

void AdderTree::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[AdderTree] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand;
		area = 0;
		height = 0;
		width = 0;
		// Adder
		int numAdderEachStage = 0;                          // define # of adder in each stage
		int numBitEachStage = numAdderBit;                  // define # of bits of the adder in each stage
		int numAdderEachTree = 0;                           // define # of Adder in each Adder Tree
		int i = ceil(log2(numSubcoreRow));
		int j = numSubcoreRow;
		
		while (i != 0) {  // calculate the total # of full adder in each Adder Tree
			numAdderEachStage = ceil(j/2);
			numAdderEachTree += numBitEachStage*numAdderEachStage;
			numBitEachStage += 1;
			j = ceil(j/2);
			i -= 1;
		}
		adder.Initialize(numAdderEachTree, numAdderTree, clkFreq);   
		// dff.Initialize((numAdderBit+numStage)*numAdderTree, clkFreq);
		
		if (_newWidth && _option==NONE) {
			adder.CalculateArea(NULL, _newWidth, NONE);
			width = _newWidth;
			height = adder.area/width;
		} else if (_newHeight && _option==NONE) {
			adder.CalculateArea(_newHeight, NULL, NONE);
			height = _newHeight;
			width = adder.area/height;
		} else {
			cout << "[AdderTree] Error: No width assigned for the adder tree circuit" << endl;
			exit(-1);
		}
		area = height*width;
		adder.initialized = false;
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

void AdderTree::CalculateLatency(double numRead, int numUnitAdd, double _capLoad) {
	if (!initialized) {
		cout << "[AdderTree] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		int numAdderEachStage = 0;                          // define # of adder in each stage
		int numBitEachStage = numAdderBit;                  // define # of bits of the adder in each stage
		int numAdderEachTree = 0;                           // define # of Adder in each Adder Tree
		int i = 0;
		int j = 0;
		
		if (!numUnitAdd) {
			i = ceil(log2(numSubcoreRow));
			j = numSubcoreRow;
		} else {
			i = ceil(log2(numUnitAdd));
			j = numUnitAdd;
		}

		// 1.4 update: adder tree delay updated
		numAdderEachStage = ceil(j/2);
		adder.Initialize(numBitEachStage, numAdderEachStage, clkFreq); 
		adder.CalculateLatency(1e20, _capLoad, 1);
		readLatency += adder.readLatency;
		numBitEachStage += 1;
		j = ceil(j/2);
		i -= 1;
		adder.initialized = false;

		if (i>0) {
			while (i != 0) {   // calculate the total # of full adder in each Adder Tree
				numAdderEachStage = ceil(j/2);
				adder.Initialize(2, numAdderEachStage, clkFreq);   
				adder.CalculateLatency(1e20, _capLoad, 1);
				readLatency += adder.readLatency;
				numBitEachStage += 1;
				j = ceil(j/2);
				i -= 1;
				
				adder.initialized = false;
			}
		}

		//Anni update
		if (param->synchronous) {
			readLatency  = ceil(readLatency*clkFreq);	//#cycles
		}

        readLatency *= numRead;		
	}
}

void AdderTree::CalculatePower(double numRead, int numUnitAdd) {
	if (!initialized) {
		cout << "[AdderTree] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		int numAdderEachStage = 0;                          // define # of adder in each stage
		int numBitEachStage = numAdderBit;                  // define # of bits of the adder in each stage
		int numAdderEachTree = 0;                           // define # of Adder in each Adder Tree
		int i = 0;
		int j = 0;
		
		if (!numUnitAdd) {
			i = ceil(log2(numSubcoreRow));
			j = numSubcoreRow;
		} else {
			i = ceil(log2(numUnitAdd));
			j = numUnitAdd;
		}
		
		while (i != 0) {  // calculate the total # of full adder in each Adder Tree
			numAdderEachStage = ceil(j/2);
			adder.Initialize(numBitEachStage, numAdderEachStage, clkFreq);     
			adder.CalculatePower(1, numAdderEachStage);	
			readDynamicEnergy += adder.readDynamicEnergy;	
			leakage += adder.leakage;
			numBitEachStage += 1;
			j = ceil(j/2);
			i -= 1;
			
			adder.initialized = false;
		}
		readDynamicEnergy *= numAdderTree;	
		readDynamicEnergy *= numRead;
		leakage *= numAdderTree;
	}
}

void AdderTree::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

