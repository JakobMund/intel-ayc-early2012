/*
 * AYC-machine.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: mund
 */

//#define __DEBUG__

#include <nmmintrin.h>	// SSE 4.1

#include "AYC-machine.h"

#ifdef __GNUC__
/*v16qi compare16_ss4_pcmpestrm (v16qi *o1, v16qi *o2) {
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

	return rmask;

}*/
#endif

#if __INTEL_COMPILER
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
	rmask.v =  _mm_cmpestrm ((*o1).v, 16, (*o2).v, 16, mode);


	return rmask;

}
#endif
