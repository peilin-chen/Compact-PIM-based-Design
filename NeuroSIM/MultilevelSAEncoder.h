/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
* 
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#ifndef MultilevelSAEncoder_H_
#define MultilevelSAEncoder_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class MultilevelSAEncoder: public FunctionUnit {
public:
	MultilevelSAEncoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~MultilevelSAEncoder() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numLevel, int _numEncoder);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double numRead);
	void CalculatePower(double numRead);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capNandInput, capNandOutput, capNandLgInput, capNandLgOutput, capInvInput, capInvOutput;
	double widthInvN, widthInvP, widthNandN, widthNandP;
	double rampInput, rampOutput;
	int numEncoder;     // number of encoder needed
	int numLevel;       // number of levels from MultilevelSA
	int numInput;       // number of NAND gate in encoder
	int numGate;        // number of NAND gate in encoder 

	// 1.4 update
	double capInvInput_Encoder, capInvOutput_Encoder;
	
};

#endif /* MultilevelSAEncoder_H_ */