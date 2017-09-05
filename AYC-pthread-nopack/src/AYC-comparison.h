/*
 * AYC-comparison.h
 *
 *  Created on: May 4, 2012
 *      Author: mund
 */

#ifndef AYC_COMPARISON_H_
#define AYC_COMPARISON_H_

#include "AYC-types.h"

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
void compute_diagonal(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);

/*
 * UNUSED alternative for compute_diagonal
 *
 * Not based on skipping parts of the input, but rather make all comparisons, but as fast as possible using
 * STTNI instructions. This function can be used for small min_sizes (6 <= N <= 8), but scales poorly beyond
 * this compared to compute_diagonal.
 */
void compute_diagonal_16(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);

/*
 * USED for 6 <= N <= 8 only currently
 */
void compute_diagonal_4(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);

#endif /* AYC_COMPARISON_H_ */
