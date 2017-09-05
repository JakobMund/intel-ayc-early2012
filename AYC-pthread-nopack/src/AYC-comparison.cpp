/*
 * AYC-comparison.cpp
 *
 *	Comparison-based operations. The core function in this file is the "compute_diagonal"
 *	function, which generates results for individual alignments of (offsets between)
 *	the reference and test sequences.
 *
 *  Created on: Mai 04, 2012
 *      Author: mund
 */

#define CMP_HIT 0xFFFFFFFF
#define CMP_HIT8 0xFFFFFFFF
#define CMP_HITL 0x0F
#define CMP_HITLL 0x0F00
#define CMP_HITM 0xF0
#define CMP_HITMM 0xF000
#define CMP_MISS 0x00

// definition of the linear scanning functions to use..
#define LIN_SCAN_LEFT lin_scan_left
#define LIN_SCAN_RIGHT lin_scan_right

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "AYC-types.h"
#include "AYC-machine.h"

#include <nmmintrin.h>	// SSE 4.2

unsigned long lin_scan_left16(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y);
unsigned long lin_scan_right16(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long max_steps);

/*
 * Linear scanning in left (negative) direction
 *
 * Linear scan and compare a reference sequence vs a test-sequence in reverse order from begin_pos,
 * comparing one-by-one, and returning the index of the first non-matching position
 */
long lin_scan_left(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y) {
	for (unsigned long i = begin_pos; i >= 0; i--) {
		if (refseq[x+i] != testseq[y+i]) {
			return i+1;
		}
	}
	return 0;
}

/*
 * Linear scanning in right (positive) direction
 *
 * Linear scan and compare a reference sequence vs a test-sequence in forward order from begin_pos,
 * comparing one-by-one, and returning the index of the first non-matching position
 */
long lin_scan_right(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long max_steps) {
	for (unsigned long i = begin_pos; i <= max_steps; i++) {
		if (refseq[x+i] != testseq[y+i]) {
			return i-1;
		}
	}
	return max_steps;
}

/**
 * Fast comparison with minimal operations
 *
 * This functions computes what was called "diagonals" in the original program. A diagonal is a
 * fixed alignment (offsets) for the reference and test sequences. The algorithm is minimal in the
 * number of comparisons done in the way that it only compares chunks of 4 chars (bytes), and skips
 * min_size - 2*chunk_size (here: 8) chars in-between. Only when a match of a complete chunk is
 * encountered, linear scans from this position are started.
 *
 * The comparison itself is done using STTNI (SSE 4.2) instructions to compare 4 chunks a time, i.e.,
 * 4*4 Byte. This SIMD techniques enables to scan large parts of non-matching data very fast and
 * efficiently.
 */
void compute_diagonal(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	const unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned long lb, ub;
	const unsigned int gapsize = min_size - 8;
	const unsigned int gapsize4 = gapsize+4;
	if (max_steps < min_size) { // if length of diagonal is smaller than min_size, no matching is possible
			return;
		}
	int z;
	for (unsigned long i = 0; i <= max_steps; i += (gapsize4)*4 ) {
		// Load 4 * 4 byte chunks into xmm0,xmm1 registers and do a STTNI (SSE 4.2) comparison
		__asm__ __volatile__ (
				        "pinsrd $0x00,%1,%%xmm0;"      		// xmm0[31:0]   <- t_0
						"pinsrd $0x01,%2,%%xmm0;"      		// xmm0[63:32]  <- t_1
						"pinsrd $0x02,%3,%%xmm0;"      		// xmm0[95:64]  <- t_2
						"pinsrd $0x03,%4,%%xmm0;"      		// xmm0[127:96] <- t_3
						"pinsrd $0x00,%5,%%xmm1;"      		// xmm1[31:0]   <- r_0
						"pinsrd $0x01,%6,%%xmm1;"      		// xmm1[63:32]  <- r_1
						"pinsrd $0x02,%7,%%xmm1;"      		// xmm1[95:64]  <- r_2
						"pinsrd $0x03,%8,%%xmm1;"      		// xmm1[127:96] <- r_3
				        "pcmpistrm $0x08,%%xmm1,%%xmm0;"    // xmm0 <- xmm0 STRCMP xmm1 (no scatter/bitwise mask)
						"pextrw $0x00,%%xmm0,%0;"			// xmm0[31:0] -> z
				        :"=r"(z)          					// output operand holding result %0
				        :"m"(testseq[i+y + gapsize]), 		// %1..
				         "m"(testseq[i+y + gapsize + gapsize4]),
				         "m"(testseq[i+y + gapsize + 2*gapsize4]),
				         "m"(testseq[i+y + gapsize + 3*gapsize4]),
				         "m"(refseq[i+x + gapsize]), 		// ..%5..
				         "m"(refseq[i+x + gapsize + gapsize4]),
				         "m"(refseq[i+x + gapsize + 2*gapsize4]),
				         "m"(refseq[i+x + gapsize + 3*gapsize4]) // ..%8
				        : "%xmm0","%xmm1" 	 // clobbered registers
				    );

		if ((z & CMP_HITL) == CMP_HITL) { // (complete) matching found in 1. quadruple, needs further investigation
			unsigned long index = i + gapsize4; // TODO: not sure about these limits
			// get upper and lower bounds
			lb = LIN_SCAN_LEFT (index -4 -1, refseq, testseq, x, y);
			ub = LIN_SCAN_RIGHT (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: verify why -2 ???? TODO: whats wrong with dt?
				// add to potential result list
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if ((z & CMP_HITM) == CMP_HITM) { // (complete) matching found in 2. quadruple, needs further investigation
			unsigned long index = i + (gapsize4)*2;
			// get upper and lower bounds
			lb = LIN_SCAN_LEFT (index -4 -1, refseq, testseq, x, y); // TODO: not sure about these limits
			ub = LIN_SCAN_RIGHT (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				// add to potential result list
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if ((z & CMP_HITLL) == CMP_HITLL) { // (complete) matching found in 3. quadruple, needs further investigation
			unsigned long index = i + (gapsize4)*3;
			lb = LIN_SCAN_LEFT (index -4, refseq, testseq, x, y);
			ub = LIN_SCAN_RIGHT (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				// add to potential result list
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
				hd->mfr = false;
				hd->next = (*rlist); (*rlist) = hd;
				i = ub; continue;
			}
		}

		if ((z & CMP_HITMM) == CMP_HITMM) { // (complete) matching found in 4. quadruple, needs further investigation
			// get upper and lower bounds
			unsigned long index = i + (gapsize4)*4;
			lb = LIN_SCAN_LEFT (index-5, refseq, testseq, x, y);
			ub = LIN_SCAN_RIGHT (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // equals UB - LB +1 >= min_size
				// add to potential result list
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

/*
 * UNUSED alternative for compute_diagonal
 *
 * Not based on skipping parts of the input, but rather make all comparisons, but as fast as possible using
 * STTNI instructions. This function can be used for small min_sizes (6 <= N <= 8), but scales poorly beyond
 * this compared to compute_diagonal.
 */
void compute_diagonal_16(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
	// lower bound for iteration
	const unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x) : (dim_y-y);
	long lb = 0, ub = 0;
	if (max_steps < min_size) {
			return;
		}
	long l_lb = 0, l_ub = 0;
	bool matching = false;

	for (unsigned long i = 0; i <= max_steps; i += 16 ) {

		__asm__ __volatile__ (
				        "movdqu %2, %%xmm0;"      			// xmm0[127:0]  <- testseq
						//"movdqu %3, %%xmm1;"				// xmm1[127:0]	<- refseq
				        "pcmpistri $0x18,%3,%%xmm0;"    // xmm0 STRCMP xmm1 (0x08 = uns. bytes, equal-each, LSB, inverted)
						"movl %%ecx, %0;"
						"pcmpistri $0x58,%3,%%xmm0;"    // xmm0 STRCMP xmm1 (0x48 = uns. bytes, equal each, MSB, inverted)
						"movl %%ecx, %1;"
				        :"=m"(l_lb),          	// output operands, local lower/upper bounds
				         "=m"(l_ub)          	// output operands, local lower/upper bounds
				        :"m"(refseq[i+x]), 		// %2
				         "m"(testseq[i+y]) 		// %3
				        : "%xmm0","%xmm1", "%ecx"  // clobbered registers
				    );


		if (matching && (l_lb != 16 || i+16 >= max_steps)) { // end of match
			ub += l_lb -1; // OK
			if (ub >= max_steps ) ub = max_steps-1;

			if (ub-lb > min_size-2) {
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
						hd->val.t_lb = y+lb+1-dt; hd->val.t_ub = y+ub+1-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
					} // OK

			if (l_ub != 15) {
				lb = i+l_ub+1; // ...start new matching
				ub = i+16;
			} else {
				matching = false;
			}
		} // OK
		else if (matching && l_lb == 16) { // all matching!
			ub += 16;
		} // OK
		else if (l_lb == 16) { // not in matching mode, but all match!
			lb = i;
			ub = i+16;
			matching = true;
		} // OK
		else {
			if (l_ub != 15) {
				lb = i+l_ub+1;
				ub = i+16;
				matching = true;
			} else {
				matching = false;
			}
		}

		if (l_ub - l_lb >= min_size-1) { // check for intra-matchings

			if (l_lb == 0 || l_ub == 0) {
					printf("Matching: %c%c%c%c.....%c%c%c%c with\n", refseq[i+x+0],refseq[i+x+1],refseq[i+x+2],refseq[i+x+3],  refseq[i+x+12],refseq[i+x+13],refseq[i+x+14],refseq[i+x+15]);
					printf("          %c%c%c%c.....%c%c%c%c\n", testseq[i+y+0],testseq[i+y+1],testseq[i+y+2],testseq[i+y+3],  testseq[i+y+12],testseq[i+y+13],testseq[i+y+14],testseq[i+y+15]);
					printf("Result:    LB=%u,       UB=%u\n\n",l_lb,l_ub);
					}
			int z = 0;
			for (int j= i+l_lb; j < i+l_ub+2; j++) {
				if (refseq[x+j] == testseq[y+j]) {
					z++;
				} else {
					if (z >= min_size) {
						RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
						hd->val.s_lb = x+j+1-z; hd->val.s_ub = x+j;
						hd->val.t_lb = y+j+1-z-dt; hd->val.t_ub = y+j-dt;
						hd->mfr = false;
						hd->next = (*rlist); (*rlist) = hd;
					}
					z = 0;
				}
			}
		}
	}
	return;
}

unsigned long lin_scan_left16(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y) {
	int ub = 0;
	for (unsigned long i = begin_pos; i >= 16; i -= 16) {
		__asm__ __volatile__ (
			"movdqu %1, %%xmm0;"      		// xmm0[127:0]  <- refseq
			"pcmpistri $0x58,%2,%%xmm0;"    // xmm0 STRCMP xmm1 (0x48 = uns. bytes, equal each, MSB, inverted)
			"movl %%ecx, %0;"
			:"=X"(ub)          			// output operands, local upper bound
			:"m"(refseq[i+x]), 				// %1
			 "m"(testseq[i+y]) 				// %2
			: "%xmm0","%xmm1", "%ecx"  // clobbered registers
		);
		if (ub != 16) {
			return i-ub+2;
		}
	}
	return lin_scan_left(16, refseq, testseq, x, y);
}

unsigned long lin_scan_right16(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long max_steps) {
	int lb = 0;
	for (unsigned long i = begin_pos; i <= max_steps; i += 16) {
		__asm__ (
			"movdqu %1, %%xmm0;"      			// xmm0[127:0]  <- refseq
			"pcmpistri $0x18,%2,%%xmm0;"    // xmm0 STRCMP xmm1 (0x08 = uns. bytes, equal-each, LSB, inverted)
			"movl %%ecx, %0;"
			:"=X"(lb)          	// output operands, local lower bound
			:"m"(refseq[i+x]), 		// %1
			 "m"(testseq[i+y]) 		// %2
			: "%xmm0","%xmm1", "%ecx"  // clobbered registers
		);
		if (lb != 16) {
			return i+lb-1;
		}
	}
	return max_steps;
}

/*
 * version for 6 <= N <= 7,
 * low performance though
 */
void compute_diagonal_4(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	long lb = 0, ub = 0;
	int gapsize = min_size >= 8 ? min_size - 8 : -3; // TODO: test if this is correct
	int gapsize0 = min_size >= 8 ? min_size-8 : 0;
	int gapsize4 = gapsize+4;
	//volatile v16qi R;
	v16qi v_rss, v_tss, R;
	//__v4si test;
	for (long i = 0; i <= max_steps; i += (gapsize4)*4 ) {
		long t_offset = i+y+gapsize;
		memcpy(&v_tss.c[0],  &testseq[t_offset], 4);
		memcpy(&v_tss.c[4],  &testseq[t_offset + gapsize4], 4);
		memcpy(&v_tss.c[8],  &testseq[t_offset + 2*gapsize4], 4);
		memcpy(&v_tss.c[12], &testseq[t_offset + 3*gapsize4], 4);

		long r_offset = i+x+gapsize;
		memcpy(&v_rss.c[0],  &refseq[r_offset], 4);
		memcpy(&v_rss.c[4],  &refseq[r_offset + gapsize4], 4);
		memcpy(&v_rss.c[8],  &refseq[r_offset + 2*gapsize4], 4);
		memcpy(&v_rss.c[12], &refseq[r_offset + 3*gapsize4], 4);

		const int mode = 0x4A; // 0100 1010b (scatter=1, mode=10 (cmp-each), signed byte=10)
		R.v = _mm_cmpestrm(v_rss.v, 16, v_tss.v, 16, mode);

		if (R.c[0]  == CMP_HIT && R.c[1]  == CMP_HIT && R.c[2]  == CMP_HIT && R.c[3]  == CMP_HIT) { // 1. quadruple
			unsigned long index = i + gapsize+4;
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size-2 && index+y < dim_y && index+x < dim_x && y+ub+1-dt >= min_size) { // TODO: whats wrong with dt?
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				if (lb < 0) lb = 0;
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
				if (lb < 0) lb = 0;
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
				if (lb < 0) lb = 0;
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
				if (lb < 0) lb = 0;
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
