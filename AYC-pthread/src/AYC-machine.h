/*
 * AYC-machine.h
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

#include <smmintrin.h>

#ifndef AYC_MACHINE_H_
#define AYC_MACHINE_H_

//typedef char __v16qi __attribute__ ((__vector_size__ (16)));
//typedef __m128i __v16qi

union v16qi	{	// 128 bit (char)
  // __v16qi v;
  __m128i v;
  char c[16];
};

void compare16_ss4_pcmpestrm (char refc, char *testc, char **R, unsigned long ref_index, unsigned long test_index);

#endif /* AYC_MACHINE_H_ */
