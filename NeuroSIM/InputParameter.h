#ifndef INPUTPARAMETER_H_
#define INPUTPARAMETER_H_

#include "typedef.h"

class InputParameter {
public:
	/* Properties */
	int processNode;				/* Process node (nm) */
	DeviceRoadmap deviceRoadmap;	/* ITRS roadmap: HP or LSTP */
	int temperature;				/* The ambient temperature, Unit: K */
	TransistorType transistorType;	/* Conventional CMOS, 2D FET, or TFET */
};

#endif /* INPUTPARAMETER_H_ */
