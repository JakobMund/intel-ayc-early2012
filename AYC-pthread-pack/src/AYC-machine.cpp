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

v16qi compare16_ss4_pcmpestrm (v16qi *o1, v16qi *o2) {
	union v16qi rmask;
	const unsigned short size = 16;

#ifdef __DEBUG__
	printf("DEBUG: Comparing %c to ", refc);
	for (short int i = 0; i < size; i++) {
		printf("%c", testc[i]);
	}
	printf("\n");
#endif __DEBUG__

	//const int mode = 0x0A; // 0000 1010b (scatter=0, mode=10 (cmp-each), signed byte=10)
	const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10 (cmp-each), signed byte=10)
	rmask.v = __builtin_ia32_pcmpestrm128((*o1).v, 16, (*o2).v, 16, mode);

	//memcpy( (R[ref_index] + test_index) , &rmask.v, 16 * sizeof(char));
	return rmask;

	//printf("DEBUG: Result is ");
	/*for (unsigned short i = 0; i < size; i++ ) { // TODO: unnecessary, if would we use a bit representation for R
			R[ref_index][test_index+i] = (rmask.c[i] == 0xFFFFFFFF ? '+' : '-');
	}*/

}

