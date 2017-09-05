/*
 * AYC-pthread.cpp
 *
 *  Created on: May 6, 2012
 *      Author: mund
 */

//#define __DEBUG__
#define MAX_TEST_SEQUENCES 64
#define COMPUTE_DIAGONAL compute_diagonal
#define COMPUTE_DIAGONAL_SHORT compute_diagonal_4

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <omp.h>

#include "AYC-types.h"
#include "AYC-sequences.h"
#include "AYC-comparison.h"
#include "AYC-result-list.h"
#include "AYC-openmp.h"

void par_process_files(unsigned short num_threads, char **argv, SEQINFO **si, unsigned short num_test_files, unsigned int min_match_size, unsigned char *refseq, unsigned long refseq_l, SEQINFO_LIST **sil)
{
	unsigned char *testseq[num_test_files];

	omp_set_num_threads(num_threads);
	omp_set_nested(1);


	for (unsigned short f = 0; f < num_test_files; f++) { // for each reference file

		unsigned dt = 0;
		si[f] = (SEQINFO *) malloc (MAX_TEST_SEQUENCES * sizeof(SEQINFO));
		unsigned long testseq_l = 0;
		unsigned char *testseq;
		unsigned int seq_no = get_sequence(argv[f+4], &testseq, &testseq_l, min_match_size, si[f]);

		for (unsigned int s = 0; s < seq_no; s++) {
			SEQINFO_LIST *c_sil;
			RESLIST *rlist = NULL;


			unsigned long chunk_r=refseq_l / (num_threads * 16);
			unsigned long chunk_t=(si[f][s].max - si[f][s].min +1) / (num_threads * 16);
			// compute diagonals in parallel
			//#pragma omp parallel default(shared)

			#pragma omp parallel default(shared)
			#pragma omp for nowait schedule(dynamic,chunk_r)
			for (unsigned long i = 0; i < refseq_l; i++) {
					if (min_match_size >= 8) {
						COMPUTE_DIAGONAL(refseq,testseq,i,si[f][s].min,refseq_l,si[f][s].max,min_match_size, &rlist, dt);
					} else {
						COMPUTE_DIAGONAL_SHORT(refseq,testseq,i,si[f][s].min,refseq_l,si[f][s].max,min_match_size, &rlist, dt);
					}
			}

			// ... and all "lower" diagonals (start vertically with second element)
			//#pragma omp for
			#pragma omp for schedule(dynamic,chunk_t)
			for (unsigned long j = si[f][s].min+1; j < si[f][s].max; j++) {
				if (min_match_size >= 8) {
					COMPUTE_DIAGONAL(refseq,testseq,0,j,refseq_l,si[f][s].max,min_match_size, &rlist, dt);
				} else {
					COMPUTE_DIAGONAL_SHORT(refseq,testseq,0,j,refseq_l,si[f][s].max,min_match_size, &rlist, dt);
				}
			}

			// remove duplicates, and sort output
			rlist = remove_substrings (rlist);
			rlist = remove_marked (rlist);

			si[f][s].rl = rlist;
			dt = si[f][s].max;
			// buffer results until evaluation completes
			// serialize access to the SIList
			c_sil = (SEQINFO_LIST *) malloc (sizeof(SEQINFO_LIST));
			c_sil->si = &(si[f][s]);

			c_sil->next = *sil;
			*sil = c_sil;
		}
	}
}
