/*
 * AYC-machine.h
 *
 * Introduces machine specific data structures;
 * in this case, a vectorized (aligned) memory block
 * of 16 bytes, accessable as a whole (__m128i) or
 * individually on a byte basis (char[16]).
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

#ifndef AYC_MACHINE_H_
#define AYC_MACHINE_H_

#include <nmmintrin.h>	// SSE 4.2

#if __GNUC__
typedef char __v16qi __attribute__ ((__vector_size__ (16)));	// definition for vectorized, memory-aligned types in GCC
union v16qi	{	// 128 bit (char)
	__m128i v;
  char c[16];
};
#endif

#if __INTEL_COMPILER

union v16qi {
	__m128i v;
	char c[16];
};

#endif

#endif /* AYC_MACHINE_H_ */
