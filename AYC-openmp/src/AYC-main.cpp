//============================================================================
// Name        : AYC-pthread.cpp
// Author      : Jakob Mund, Sebastian Eder
// Version     : 1.0
// Copyright   : 
// Description : AYC contest, new implementation based on Basti's Idea
// Contest-ID  : 4a09b31365e51a609800f4660b1bccaf
//============================================================================

//#define __DEBUG__	// Enable debug outputs
//#define __PROFILE__ // Enable basic profiling of core routines

#define MAX_TEST_SEQUENCES 32

#ifdef __DEBUG__
#include <mcheck.h>
#endif

#ifdef __PROFILE__
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#define __PROFILE_START() 		gettimeofday(&start, NULL);
#define __PROFILE_END( MSG ) 	gettimeofday(&end, NULL);\
			seconds  = end.tv_sec  - start.tv_sec; useconds = end.tv_usec - start.tv_usec;\
			mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;\
			printf("PROFILE: %s %ld ms\n", #MSG, mtime );
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <omp.h>

#include "AYC-types.h"
#include "AYC-sequences.h"
#include "AYC-result-list.h"
#include "AYC-machine.h"
#include "AYC-comparison.h"
#include "AYC-openmp.h"

int main(int argc, char **argv) {
#ifdef __DEBUG__
	mcheck(NULL);
#endif

#ifdef __PROFILE__
	struct timeval start, end;
	long mtime, seconds, useconds;
#ifdef __DEBUG__
	printf("WARNING: Profiling results may be flawed, as __DEBUG__ is active\n");
#endif
#endif

	// parse arguments
	if (argc < 5) {
		printf("USAGE: ayc-pthread <max_threads> <min_match_size> <ref_seq> <test_seq1> [<test_seq2> <...>]\n");
		return 0;
	}
	unsigned int max_threads = atoi(argv[1]);		// maximum number of worker threads, as given as an argument
	//pthread_t threads[max_threads];
	unsigned int min_match_size = atoi(argv[2]);	// minimum length for matching
	unsigned short num_test_files = (argc-4);		// number of files containing test sequences
	char *reffname = argv[3];						// file names of the input files
	//char *testfname[num_test_files];
#ifdef __DEBUG__
	printf("DEBUG: %u test sequence files submitted\n", num_test_files);
#endif
	//unsigned int seq_no[num_test_files];			// number of test sequences (in each file)
	SEQINFO *si[num_test_files], ri[1];
	SEQINFO_LIST *sil = NULL;
	//SEQINFO_LIST *c_sil = NULL;

	unsigned char *refseq; // *testseq[num_test_files];
	unsigned long refseq_l = 0; // testseq_l[num_test_files];

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	get_sequence(reffname, &refseq, &refseq_l, min_match_size, ri);
#ifdef __PROFILE__
	__PROFILE_END(Loading Reference Sequence)
#endif

#ifdef __DEBUG__
	printf("DEBUG: size of reference sequence is %u\n", refseq_l);
#endif

	// process the test sequences in parallel
	par_process_files(max_threads, argv, si, num_test_files, min_match_size, refseq, refseq_l, &sil);

	// Finally, print the results in alphabetical order
#ifdef __PROFILE__
	__PROFILE_START()
#endif
	results(sil);
#ifdef __PROFILE__
	__PROFILE_END(Printing_Results)
#endif

#ifdef __DEBUG__
	printf("DEBUG: Program successfully terminated.\n");
#endif

	return 0;
}
