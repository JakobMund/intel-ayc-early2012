/*
 * AYC-pthread.h
 *
 *  Created on: May 6, 2012
 *      Author: mund
 */

#ifndef AYC_PTHREAD_H_
#define AYC_PTHREAD_H_

#include <pthread.h>
#include "AYC-types.h"

/*
 * Individual job computing a _subset_ of all diagonals (alignments)
 * between reference and test sequences.
 */
void par_compute_diagonal(unsigned int num_threads, unsigned long refseq_max, unsigned long testseq_min,
		unsigned long testseq_max, unsigned char *refseq, unsigned char *testseq, unsigned int min_size,
		RESLIST **rlist, unsigned long dt);

/*
 * par_process_files is the main parallel control
 * flow of the algorithm. First, it reads in all test files and the
 * reference file in parallel.
 * Afterwards, it distributes the computation of diagonals (matchings)
 * across a thread pool.
 * When the computation of the individual diagonals is finished, the
 * results are concatened and substrings (sub-matches) are removed. This
 * step happens in parallel for all test sequences.
 */
void par_process_files(unsigned short num_threads, char *reffname, char **argv, SEQINFO **si,
		unsigned short num_test_files, unsigned int min_match_size, SEQINFO_LIST **sil);

#endif /* AYC_PTHREAD_H_ */
