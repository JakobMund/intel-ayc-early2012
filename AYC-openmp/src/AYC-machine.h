/*
 * AYC-machine.h
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

#ifndef AYC_MACHINE_H_
#define AYC_MACHINE_H_

#include <nmmintrin.h>	// SSE 4.2

#if __GNUC__

//  __v16qi v;
typedef char __v16qi __attribute__ ((__vector_size__ (16)));
union v16qi	{	// 128 bit (char)
  __m128i v;
  char c[16];
};

v16qi compare16_ss4_pcmpestrm (v16qi *refc, v16qi *testc);

#endif

#if __INTEL_COMPILER

union v16qi {
	__m128i v;
	char c[16];
};

v16qi compare16_ss4_pcmpestrm (v16qi *refc, v16qi *testc);

#endif

#endif /* AYC_MACHINE_H_ */
