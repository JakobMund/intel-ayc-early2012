/*
 * AYC-pthread.h
 *
 *  Created on: May 6, 2012
 *      Author: mund
 */

#ifndef AYC_OPENMP_H_
#define AYC_OPENMP_H_

#include <omp.h>
#include "AYC-types.h"

void par_process_files(unsigned short num_threads, char **argv, SEQINFO **si, unsigned short num_test_files,
			unsigned int min_match_size, unsigned char *refseq, unsigned long refseq_l, SEQINFO_LIST **sil);

#endif /* AYC_OPENMP_H_ */
