/*
 * AYC-machine.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

//#define __DEBUG__

#include <stdio.h>
#include <memory.h>
#include <smmintrin.h>	// SSE 4.2

#include "AYC-machine.h"

void compare16_ss4_pcmpestrm (char refc, char *testc, char **R, unsigned long ref_index, unsigned long test_index) {
	union v16qi o1, o2, rmask;
	const unsigned short size = 16;
	const unsigned short chunksize = 4;

	memset(&o1, refc, size); // fill o1 (match against refc)
	memcpy(&o2, testc, size);

#ifdef __DEBUG__
	printf("DEBUG: Comparing %c to ", refc);
	for (short int i = 0; i < size; i++) {
		printf("%c", testc[i]);
	}
	printf("\n");
#endif __DEBUG__

	const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10, signed byte=10)

	//rmask.v = __builtin_ia32_pcmpestrm128(o1.v, 16, o2.v, 16, mode);
	rmask.v = _mm_cmpestrm(o1.v, 16, o2.v, 16, mode);

	memcpy( (R[ref_index] + test_index) , &rmask.v, 16 * sizeof(char));

	//printf("DEBUG: Result is ");
	/*for (unsigned short i = 0; i < size; i++ ) { // TODO: unnecessary, if would we use a bit representation for R
			R[ref_index][test_index+i] = (rmask.c[i] == 0xFFFFFFFF ? '+' : '-');
	}*/

}

