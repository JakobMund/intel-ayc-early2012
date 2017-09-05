/*
 * AYC-comparison.cpp
 *
 *  Created on: Mai 04, 2012
 *      Author: mund
 */

#define CMP_HIT 0xFFFFFFFF
#define CMP_MISS 0x00

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "AYC-types.h"
#include "AYC-machine.h"

#include <nmmintrin.h>	// SSE 4.2

unsigned long lin_scan_left(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y) {
	for (unsigned long i = begin_pos; i >= 0; i--) {
		if (refseq[x+i] != testseq[y+i]) {
			return i+1;
		}
	}
	return 0;
}

unsigned long lin_scan_right(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long max_steps) {
	for (unsigned long i = begin_pos; i <= max_steps; i++) {
		if (refseq[x+i] != testseq[y+i]) {
			return i-1;
		}
	}
	return max_steps;
}

void compute_diagonal(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned long lb = 0, ub = 0;
	int gapsize = min_size >= 8 ? min_size - 8 : -3; // TODO: test if this is correct
	unsigned int gapsize0 = min_size >= 8 ? min_size-8 : 0;
	unsigned int gapsize4 = gapsize+4;
	 v16qi R;
	//v16qi v_rss, v_tss;
	//__v4si test;
	for (unsigned long i = 0; i <= max_steps; i += (gapsize4)*4 ) {
		//long t_offset = i+y+gapsize;
		//memcpy(&v_tss.c[0],  &testseq[t_offset], 4);
		//memcpy(&v_tss.c[4],  &testseq[t_offset + gapsize4], 4);
		//memcpy(&v_tss.c[8],  &testseq[t_offset + 2*gapsize4], 4);
		//memcpy(&v_tss.c[12], &testseq[t_offset + 3*gapsize4], 4);

		//long r_offset = i+x+gapsize;
		//memcpy(&v_rss.c[0],  &refseq[r_offset], 4);
		//memcpy(&v_rss.c[4],  &refseq[r_offset + gapsize4], 4);
		//memcpy(&v_rss.c[8],  &refseq[r_offset + 2*gapsize4], 4);
		//memcpy(&v_rss.c[12], &refseq[r_offset + 3*gapsize4], 4);

		//R = compare16_ss4_pcmpestrm(&v_rss, &v_tss);
		//v_tss.v = _mm_set_epi32(testseq[t_offset[0]], testseq[t_offset[1]], testseq[t_offset[2]], testseq[t_offset[3]];
		//v_rss.v = _mm_set_epi32(refseq[t_offset[0]], refseq[t_offset[1]], refseq[t_offset[2]], refseq[t_offset[3]]);

		/////////////////////////////////////////////////////////

		//const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10 (cmp-each), signed byte=10)
		//R.v = _mm_cmpestrm(v_rss.v, 16, v_tss.v, 16, mode);
		/*R.v = _mm_cmpestrm(
					_mm_setr_epi32(testseq[t_offset[0]], testseq[t_offset[1]], testseq[t_offset[2]], testseq[t_offset[3]]), 16,
					_mm_setr_epi32(refseq[r_offset[0]],  refseq[r_offset[1]],  refseq[r_offset[2]],  refseq[r_offset[3]]),  16, mode);*/

		/*printf("DEBUG: CMP %c%c%c%c to %c%c%c%c\n", testseq[t_offset[0]], testseq[t_offset[0]+1], testseq[t_offset[0]+2], testseq[t_offset[0]+3],
				refseq[r_offset[0]], refseq[r_offset[0]+1], refseq[r_offset[0]+2], refseq[r_offset[0]+3]); */

		//printf("DEBUG: checking REF=%u TEST=%u...", i+x+gapsize0, i+y+gapsize0);

		__asm__ __volatile__ (
				//"cli;"
		        "pinsrd $0x00,%1,%%xmm0;"      		// xmm0[31:0]   <- t_0
				"pinsrd $0x01,%2,%%xmm0;"      		// xmm0[63:32]  <- t_1
				"pinsrd $0x02,%3,%%xmm0;"      		// xmm0[95:64]  <- t_2
				"pinsrd $0x03,%4,%%xmm0;"      		// xmm0[127:96] <- t_3
				"pinsrd $0x00,%5,%%xmm1;"      		// xmm0[31:0]   <- t_0
				"pinsrd $0x01,%6,%%xmm1;"      		// xmm0[63:32]  <- t_1
				"pinsrd $0x02,%7,%%xmm1;"      		// xmm0[95:64]  <- t_2
				"pinsrd $0x03,%8,%%xmm1;"      		// xmm0[127:96] <- t_3
				"mov $16,%%rax;"					// length of input in EAX (RAX 64B)
				"mov $16,%%rdx;"					// length of input in EDX (RDX 64B)
		        "pcmpestrm $0x4A,%%xmm1,%%xmm0;"    	// xmm0 <- xmm0 STRCMP xmm1
		        "movdqa %%xmm0, %0;"      				// R <- xmm0
				//"sti;"

		        :"=X"(R)          // output operand, %0
		        :"m"(testseq[i+y + gapsize0]), // %1
		         "m"(testseq[i+y + gapsize0 + gapsize4]),
		         "m"(testseq[i+y + gapsize0 + 2*gapsize4]),
		         "m"(testseq[i+y + gapsize0 + 3*gapsize4]),
		         "m"(refseq[i+x + gapsize0]), // %2
		         "m"(refseq[i+x + gapsize0 + gapsize4]), // %2
		         "m"(refseq[i+x + gapsize0 + 2*gapsize4]), // %2
		         "m"(refseq[i+x + gapsize0 + 3*gapsize4]) // %2
		        : "%xmm0","%xmm1", "%rax", "%rdx"  // clobbered registers
		    );

		//printf("OK\n");


		/////////////////////////////////////////////////////////

		if (R.c[0]  == CMP_HIT && R.c[1]  == CMP_HIT && R.c[2]  == CMP_HIT && R.c[3]  == CMP_HIT) { // 1. quadruple
			unsigned long index = i + gapsize+4; // TODO: not sure about these limits
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[7]  == CMP_HIT && R.c[4]  == CMP_HIT && R.c[5] == CMP_HIT && R.c[6]  == CMP_HIT) { // 2. quadruple
			unsigned long index = i + (gapsize+4)*2;
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[8]  == CMP_HIT && R.c[9]  == CMP_HIT && R.c[10]  == CMP_HIT && R.c[11]  == CMP_HIT) { // 3. quadruple
			unsigned long index = i + (gapsize+4)*3;
			lb = lin_scan_left (index -4, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[12]  == CMP_HIT && R.c[13]  == CMP_HIT && R.c[14]  == CMP_HIT && R.c[15]  == CMP_HIT) { // 4. quadruple
			// get upper and lower bounds
			unsigned long index = i + (gapsize+4)*4;
			lb = lin_scan_left (index-5, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}
	}
	return;
}

void compute_diagonal_4(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned long lb = 0, ub = 0;
	int gapsize = min_size >= 8 ? min_size - 8 : -3; // TODO: test if this is correct
	unsigned int gapsize0 = min_size >= 8 ? min_size-8 : 0;
	unsigned int gapsize4 = gapsize+4;
	//volatile v16qi R;
	v16qi v_rss, v_tss, R;
	//__v4si test;
	for (unsigned long i = 0; i <= max_steps; i += (gapsize4)*4 ) {
		long t_offset = i+y+gapsize0;
		memcpy(&v_tss.c[0],  &testseq[t_offset], 4);
		memcpy(&v_tss.c[4],  &testseq[t_offset + gapsize4], 4);
		memcpy(&v_tss.c[8],  &testseq[t_offset + 2*gapsize4], 4);
		memcpy(&v_tss.c[12], &testseq[t_offset + 3*gapsize4], 4);

		long r_offset = i+x+gapsize0;
		memcpy(&v_rss.c[0],  &refseq[r_offset], 4);
		memcpy(&v_rss.c[4],  &refseq[r_offset + gapsize4], 4);
		memcpy(&v_rss.c[8],  &refseq[r_offset + 2*gapsize4], 4);
		memcpy(&v_rss.c[12], &refseq[r_offset + 3*gapsize4], 4);

		//R = compare16_ss4_pcmpestrm(&v_rss, &v_tss);
		//v_tss.v = _mm_set_epi32(testseq[t_offset[0]], testseq[t_offset[1]], testseq[t_offset[2]], testseq[t_offset[3]];
		//v_rss.v = _mm_set_epi32(refseq[t_offset[0]], refseq[t_offset[1]], refseq[t_offset[2]], refseq[t_offset[3]]);

		/////////////////////////////////////////////////////////

		const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10 (cmp-each), signed byte=10)
		R.v = _mm_cmpestrm(v_rss.v, 16, v_tss.v, 16, mode);
		/*R.v = _mm_cmpestrm(
					_mm_setr_epi32(testseq[t_offset[0]], testseq[t_offset[1]], testseq[t_offset[2]], testseq[t_offset[3]]), 16,
					_mm_setr_epi32(refseq[r_offset[0]],  refseq[r_offset[1]],  refseq[r_offset[2]],  refseq[r_offset[3]]),  16, mode);*/

		/*printf("DEBUG: CMP %c%c%c%c to %c%c%c%c\n", testseq[t_offset[0]], testseq[t_offset[0]+1], testseq[t_offset[0]+2], testseq[t_offset[0]+3],
				refseq[r_offset[0]], refseq[r_offset[0]+1], refseq[r_offset[0]+2], refseq[r_offset[0]+3]); */

		//printf("DEBUG: checking REF=%u TEST=%u...", i+x+gapsize0, i+y+gapsize0);

		/*__asm__ __volatile__ (
				//"cli;"
		        "pinsrd $0x00,%1,%%xmm0;"      		// xmm0[31:0]   <- t_0
				"pinsrd $0x01,%2,%%xmm0;"      		// xmm0[63:32]  <- t_1
				"pinsrd $0x02,%3,%%xmm0;"      		// xmm0[95:64]  <- t_2
				"pinsrd $0x03,%4,%%xmm0;"      		// xmm0[127:96] <- t_3
				"pinsrd $0x00,%5,%%xmm1;"      		// xmm0[31:0]   <- t_0
				"pinsrd $0x01,%6,%%xmm1;"      		// xmm0[63:32]  <- t_1
				"pinsrd $0x02,%7,%%xmm1;"      		// xmm0[95:64]  <- t_2
				"pinsrd $0x03,%8,%%xmm1;"      		// xmm0[127:96] <- t_3
				"mov $16,%%rax;"					// length of input in EAX (RAX 64B)
				"mov $16,%%rdx;"					// length of input in EDX (RDX 64B)
		        "pcmpestrm $0x4A,%%xmm1,%%xmm0;"    	// xmm0 <- xmm0 STRCMP xmm1
		        "movdqa %%xmm0, %0;"      				// R <- xmm0
				//"sti;"

		        :"=X"(R)          // output operand, %0
		        :"m"(testseq[i+y + gapsize0]), // %1
		         "m"(testseq[i+y + gapsize0 + gapsize4]),
		         "m"(testseq[i+y + gapsize0 + 2*gapsize4]),
		         "m"(testseq[i+y + gapsize0 + 3*gapsize4]),
		         "m"(refseq[i+x + gapsize0]), // %2
		         "m"(refseq[i+x + gapsize0 + gapsize4]), // %2
		         "m"(refseq[i+x + gapsize0 + 2*gapsize4]), // %2
		         "m"(refseq[i+x + gapsize0 + 3*gapsize4]) // %2
		        : "%xmm0","%xmm1", "%rax", "%rdx"  // clobbered registers
		    );*/

		//printf("OK\n");


		/////////////////////////////////////////////////////////

		if (R.c[0]  == CMP_HIT && R.c[1]  == CMP_HIT && R.c[2]  == CMP_HIT && R.c[3]  == CMP_HIT) { // 1. quadruple
			unsigned long index = i + gapsize+4; // TODO: not sure about these limits
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[7]  == CMP_HIT && R.c[4]  == CMP_HIT && R.c[5] == CMP_HIT && R.c[6]  == CMP_HIT) { // 2. quadruple
			unsigned long index = i + (gapsize+4)*2;
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[8]  == CMP_HIT && R.c[9]  == CMP_HIT && R.c[10]  == CMP_HIT && R.c[11]  == CMP_HIT) { // 3. quadruple
			unsigned long index = i + (gapsize+4)*3;
			lb = lin_scan_left (index -4, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[12]  == CMP_HIT && R.c[13]  == CMP_HIT && R.c[14]  == CMP_HIT && R.c[15]  == CMP_HIT) { // 4. quadruple
			// get upper and lower bounds
			unsigned long index = i + (gapsize+4)*4;
			lb = lin_scan_left (index-5, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}
	}
	return;
}

void compute_diagonal_2(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned long lb = 0, ub = 0;
	unsigned long t_offset[8], r_offset[8];
	const int gapsize = min_size > 4 ? min_size - 4 : -1; // TODO: test if this is correct
	v16qi v_rss, v_tss, R;
	//__v4si test;
	for (unsigned long i = 0; i <= max_steps; i += (gapsize+2)*8 ) {

		t_offset[0] = i+y + gapsize;					memcpy(&v_tss.c[0],  &testseq[t_offset[0]], 2);
		t_offset[1] = i+y + ((gapsize +2) * 2) -2;		memcpy(&v_tss.c[2],  &testseq[t_offset[1]], 2);
		t_offset[2] = i+y + ((gapsize +2) * 3) -2;		memcpy(&v_tss.c[4],  &testseq[t_offset[2]], 2);
		t_offset[3] = i+y + ((gapsize +2) * 4) -2;		memcpy(&v_tss.c[6],  &testseq[t_offset[3]], 2);
		t_offset[4] = i+y + ((gapsize +2) * 5) -2;		memcpy(&v_tss.c[8],  &testseq[t_offset[4]], 2);
		t_offset[5] = i+y + ((gapsize +2) * 6) -2;		memcpy(&v_tss.c[10], &testseq[t_offset[5]], 2);
		t_offset[6] = i+y + ((gapsize +2) * 7) -2;		memcpy(&v_tss.c[12], &testseq[t_offset[6]], 2);
		t_offset[7] = i+y + ((gapsize +2) * 8) -2;		memcpy(&v_tss.c[14], &testseq[t_offset[7]], 2);

		r_offset[0] = i+x + gapsize;					memcpy(&v_rss.c[0],  &refseq[r_offset[0]], 2);
		r_offset[1] = i+x + ((gapsize +2) * 2) -2;		memcpy(&v_rss.c[2],  &refseq[r_offset[1]], 2);
		r_offset[2] = i+x + ((gapsize +2) * 3) -2;		memcpy(&v_rss.c[4],  &refseq[r_offset[2]], 2);
		r_offset[3] = i+x + ((gapsize +2) * 4) -2;		memcpy(&v_rss.c[6],  &refseq[r_offset[3]], 2);
		r_offset[4] = i+x + ((gapsize +2) * 5) -2;		memcpy(&v_rss.c[8],  &refseq[r_offset[4]], 2);
		r_offset[5] = i+x + ((gapsize +2) * 6) -2;		memcpy(&v_rss.c[10], &refseq[r_offset[5]], 2);
		r_offset[6] = i+x + ((gapsize +2) * 7) -2;		memcpy(&v_rss.c[12], &refseq[r_offset[6]], 2);
		r_offset[7] = i+x + ((gapsize +2) * 8) -2;		memcpy(&v_rss.c[14], &refseq[r_offset[7]], 2);

		const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10 (cmp-each), signed byte=10)
		R.v = _mm_cmpestrm(v_rss.v, 16, v_tss.v, 16, mode);
		/*R.v = _mm_cmpestrm(
				_mm_set_epi32(testseq[t_offset[0]], testseq[t_offset[1]], testseq[t_offset[2]], testseq[t_offset[3]]), 16,
				_mm_set_epi32(refseq[r_offset[0]],  refseq[r_offset[1]],  refseq[r_offset[2]],  refseq[r_offset[3]]),  16, mode); */

		if (R.c[0]  == CMP_HIT && R.c[1]  == CMP_HIT ) { // 1. tuple
			unsigned long index = i + gapsize+2; // TODO: not sure about these limits
			// get upper and lower bounds
			lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if (R.c[2] == CMP_HIT && R.c[3]  == CMP_HIT) { // 2. tuple
					unsigned long index = i + (gapsize+2)*2; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[4] == CMP_HIT && R.c[5]  == CMP_HIT) { // 3. tuple
					unsigned long index = i + (gapsize+2)*3; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[6] == CMP_HIT && R.c[7]  == CMP_HIT) { // 4. tuple
					unsigned long index = i + (gapsize+2)*4; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[8] == CMP_HIT && R.c[9]  == CMP_HIT) { // 5. tuple
					unsigned long index = i + (gapsize+2)*5; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[10] == CMP_HIT && R.c[11]  == CMP_HIT) { // 6. tuple
					unsigned long index = i + (gapsize+2)*6; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[12] == CMP_HIT && R.c[13]  == CMP_HIT) { // 7. tuple
					unsigned long index = i + (gapsize+2)*7; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

		if (R.c[14] == CMP_HIT && R.c[15]  == CMP_HIT) { // 8. tuple
					unsigned long index = i + (gapsize+2)*8; // TODO: not sure about these limits
					// get upper and lower bounds
					lb = lin_scan_left (index -2 -1, refseq, testseq, x, y); // TODO: not sure about these limits
					ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
					if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
						i = ub; continue;
					}
				}

	}
	return;
}
