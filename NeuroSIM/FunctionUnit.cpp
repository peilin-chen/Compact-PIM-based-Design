#include <iostream>
#include <fstream>
#include "FunctionUnit.h"

using namespace std;

FunctionUnit::FunctionUnit() {
	height = width = 0;
	area = 0;
	emptyArea = 0;
	usedArea = 0;
	totalArea = 0;	// may use it if the circuit units are scattered
	readLatency = writeLatency = 0;
	readDynamicEnergy = writeDynamicEnergy = 0;
	leakage = 0;

	newWidth = newHeight = 0;
	readPower = writePower = 0;
}

void FunctionUnit::PrintProperty(const char* str) {
	cout << "---------------------------------------------------------" << endl;
	cout << str << endl;
	cout << "Area = " << height*1e6 << "um x " << width*1e6 << "um = " << area*1e12 << "um^2" << endl;
	if (totalArea)
		cout << "Total Area = " << totalArea*1e12 << "um^2" << endl;
	cout << "Timing:" << endl;
	cout << " - Read Latency = " << readLatency*1e9 << "ns" << endl;
	cout << " - Write Latency = " << writeLatency*1e9 << "ns" << endl;
	cout << "Power:" << endl;
	cout << " - Read Dynamic Energy = " << readDynamicEnergy*1e12 << "pJ" << endl;
	cout << " - Write Dynamic Energy = " << writeDynamicEnergy*1e12 << "pJ" << endl;
	cout << " - Leakage Power = " << leakage*1e6 << "uW" << endl;
	// cout << " - Read Power = " << readPower*1e6 << "uW" << endl;
	// cout << " - Write Power = " << writePower*1e6 << "uW" << endl;
}

void FunctionUnit::SaveOutput(const char* str) {
	ofstream outfile;                                           
	outfile.open("SynapticCOREoutput.txt", ios::app);     
	outfile << "---------------------------------------------------------" << endl;
	outfile << str << endl;
	outfile << "Area = " << height*1e6 << "um x " << width*1e6 << "um = " << area*1e12 << "um^2" << endl;
	if (totalArea)
		outfile << "Total Area = " << totalArea*1e12 << "um^2" << endl;
	outfile << "Timing:" << endl;
	outfile << " - Read Latency = " << readLatency*1e9 << "ns" << endl;
	outfile << " - Write Latency = " << writeLatency*1e9 << "ns" << endl;
	outfile << "Power:" << endl;
	outfile << " - Read Dynamic Energy = " << readDynamicEnergy*1e12 << "pJ" << endl;
	outfile << " - Write Dynamic Energy = " << writeDynamicEnergy*1e12 << "pJ" << endl;
	outfile << " - Leakage Power = " << leakage*1e6 << "uW" << endl;
	// outfile << " - Read Power = " << readPower*1e6 << "uW" << endl;
	// outfile << " - Write Power = " << writePower*1e6 << "uW" << endl;
	outfile << '\n' << endl;
	outfile.close();
}

void FunctionUnit::MagicLayout() {
	if (newHeight) {
		width = area / newHeight;
		height = newHeight;
	} else if (newWidth) {
		height = area / newWidth;
		width = newWidth;
	}
}

void FunctionUnit::OverrideLayout() {
	if (newHeight && newWidth) {
		height = newHeight;
		width = newWidth;
	} else {
		puts("Need to provide both newHeight and newWidth for OverrideLayout()");
		exit(-1);
	}
	area = height * width;
}

