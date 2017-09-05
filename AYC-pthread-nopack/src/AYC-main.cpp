//============================================================================
//
// Name        : AYC-pthread.cpp
// Author      : Jakob Mund (mund@in.tum.de), Sebastian Eder (eders@in.tum.de)
// Affiliation : Technical University, Munich (TUM)
// Version     : 1.0 (FINAL)
// Team name   : Team Zubrowka
// Description : AYC early 2012; Genome detection based on an approach scaling
//				 efficiently for large amounts of data and cores, as well as
//				 for larger MINIMUM_MATCHING_SIZES. See enclosed readme.txt
//				 for further informations about the approach and the
//				 implementation itself.
// Contest-ID  : 4a09b31365e51a609800f4660b1bccaf
//
//============================================================================

//#define __DEBUG__	// Enable debug outputs
//#define __PROFILE__ // Enable basic profiling of core routines

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

#include "AYC-types.h"
#include "AYC-sequences.h"
#include "AYC-result-list.h"
#include "AYC-machine.h"
#include "AYC-comparison.h"
#include "AYC-pthread.h"

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
	unsigned int min_match_size = atoi(argv[2]);	// minimum length for matching
	unsigned short num_test_files = (argc-4);		// number of files containing test sequences
	char *reffname = argv[3];						// file names of the input files
#ifdef __DEBUG__
	printf("DEBUG: %u test sequence files submitted\n", num_test_files);
#endif
	SEQINFO *si[num_test_files];
	SEQINFO_LIST *sil = NULL;

	// begin parallel processing (entry point for the pthreads-based implementation)
	par_process_files(max_threads, reffname, argv, si, num_test_files, min_match_size, &sil);

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
