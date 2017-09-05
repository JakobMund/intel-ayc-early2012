/*
 * AYC-machine.h
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

#ifndef AYC_MACHINE_H_
#define AYC_MACHINE_H_

typedef char __v16qi __attribute__ ((__vector_size__ (16)));

union v16qi	{	// 128 bit (char)
  __v16qi v;
  char c[16];
};

v16qi compare16_ss4_pcmpestrm (v16qi *refc, v16qi *testc);

#endif /* AYC_MACHINE_H_ */
