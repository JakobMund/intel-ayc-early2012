/*
 * AYC-comparison.h
 *
 *  Created on: May 4, 2012
 *      Author: mund
 */

#ifndef AYC_COMPARISON_H_
#define AYC_COMPARISON_H_

#include "AYC-types.h"

void compute_diagonal(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);
void compute_diagonal_2(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);
void compute_diagonal_4(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist, unsigned long dt);

#endif /* AYC_COMPARISON_H_ */
