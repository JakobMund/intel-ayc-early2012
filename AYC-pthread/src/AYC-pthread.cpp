//============================================================================
// Name        : AYC-pthread.cpp
// Author      : Jakob Mund, Sebastian Eder
// Version     : 1.0
// Copyright   : 
// Description : AYC contest, implementation based on POSIX threads
//============================================================================

//#define __DEBUG__	// Enable debug outputs
//#define __PROFILE__ // Enable basic profiling of core routines

#define CMP_HIT 0xFFFFFFFF
#define CMP_MISS 0x00

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

#include "AYC-sequences.h"
#include "AYC-result-list.h"
#include "AYC-machine.h"

inline char __cmp_result (char c1, char c2) {
	if (c1 == c2)
		return CMP_HIT;
	else
		return CMP_MISS;
}

void compare(unsigned char prc, unsigned char *ptc, char **R, int ref_index, int test_index) {
	char refc[4];
	char *testc = (char *) malloc (16 * sizeof(char));
	__unpack(prc, refc);
	/*__unpack(ptc[0], testc);
	__unpack(ptc[1], testc+4);
	__unpack(ptc[2], testc+8);
	__unpack(ptc[3], testc+12);	// TODO: __unpack16 */
	__unpack16(ptc, testc);

//#ifdef __SSE4.2__
	compare16_ss4_pcmpestrm(refc[0], testc, R, ref_index, test_index);
	compare16_ss4_pcmpestrm(refc[1], testc, R, ref_index+1, test_index);
	compare16_ss4_pcmpestrm(refc[2], testc, R, ref_index+2, test_index);
	compare16_ss4_pcmpestrm(refc[3], testc, R, ref_index+3, test_index);
//#endif

#ifndef __SSE4.2__
	/*R[ref_index][test_index] = __cmp_result (refc[0], testc[0]);
	R[ref_index+1][test_index] = __cmp_result (refc[1], testc[0]);
	R[ref_index+2][test_index] = __cmp_result (refc[2], testc[0]);
	R[ref_index+3][test_index] = __cmp_result (refc[3], testc[0]);

	R[ref_index][test_index+1]   = __cmp_result (refc[0], testc[1]);
	R[ref_index+1][test_index+1] = __cmp_result (refc[1], testc[1]);
	R[ref_index+2][test_index+1] = __cmp_result (refc[2], testc[1]);
	R[ref_index+3][test_index+1] = __cmp_result (refc[3], testc[1]);

	R[ref_index][test_index+2]   = __cmp_result (refc[0], testc[2]);
	R[ref_index+1][test_index+2] = __cmp_result (refc[1], testc[2]);
	R[ref_index+2][test_index+2] = __cmp_result (refc[2], testc[2]);
	R[ref_index+3][test_index+2] = __cmp_result (refc[3], testc[2]);

	R[ref_index][test_index+3]   = __cmp_result (refc[0], testc[3]);
	R[ref_index+1][test_index+3] = __cmp_result (refc[1], testc[3]);
	R[ref_index+2][test_index+3] = __cmp_result (refc[2], testc[3]);
	R[ref_index+3][test_index+3] = __cmp_result (refc[3], testc[3]); */
#endif
}

unsigned long get_diagonal (char **R, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned long min_size, RESLIST **rlist) {
	// lower bound for iteration
	unsigned long c = 0;
	unsigned int max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned int res = 0;
	for (unsigned long i = 0; i <= max_steps; i++) {
		if (R[i+x][i+y] == CMP_HIT) {
			res++;
			if (res >= min_size && (i+x+1 >= dim_x || i+y+1 >= dim_y || R[i+x+1][i+y+1] != CMP_HIT)) {	// check if exceeding array, and if not, is the next char a match?
				// matching found, add to result list
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST)); // TODO: alloc chunkwise, then free afterwards!
				hd->val.s_lb = (i+1+x-res+1); hd->val.t_lb = (i+1+y-res+1); hd->val.s_ub = (i+1+x); hd->val.t_ub = (i+1+y);
				hd->next = (*rlist);
				(*rlist) = hd;
				c++;
			}
		} else if (res != 0){
			res = 0;
		}
	}
	return c;
}

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
	// use char for the data type for the result matrix
	// as char is equal to CHAR_BIT (which is 8 Bit on almost every platform).
	// hence, each char will hold the result of 8 comparisons.
	// usage is R[number of test sequence][index of ref. seq][index of test seq]
	char **R = NULL;

	// parse arguments
	if (argc != 5) {
		printf("USAGE: ayc-pthread <max_threads> <min_match_size> <ref_seq> <test_seqs>\n");
		return 0;
	}
	int max_threads = atoi(argv[1]);	// maximum number of worker threads, as given as an argument
	int min_match_size = atoi(argv[2]);	// minimum length for matching
	char *reffname = argv[3];			// file names of the input files
	char *testfname = argv[4];
	FILE *reffile = NULL, *testfile = NULL;

	unsigned char *refseq, *testseq;
	unsigned long refseq_l = 0, testseq_l = 0;

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	get_sequence(reffname, reffile, &refseq, &refseq_l);
	get_sequence(testfname, testfile, &testseq, &testseq_l);
#ifdef __PROFILE__
	__PROFILE_END(Loading)
#endif

#ifdef __DEBUG__
	mprobe(testseq);
	mprobe(refseq);
	printf("DEBUG: size of input is %lu / %lu characters\n", refseq_l, testseq_l);
#endif

	// allocate memory for R (TODO: pack the results, saving memory (and potentially matching more at once))
	size_t alloc_x = (refseq_l + (refseq_l%4) +1 );
	size_t alloc_y = (testseq_l + (testseq_l%16) +1);
	R = (char**) malloc(alloc_x * sizeof (char*));
	for (unsigned long i = 0; i < alloc_x; i++) {
	    R[i] = (char*) malloc(alloc_y * sizeof (char));
	  }

#ifdef __DEBUG__
	mprobe(R); // TODO: what's the problem here?
#endif

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	// compute the comparison matrix
	for (unsigned long i = 0; i < refseq_l; i += 4) {
		unsigned char prc = refseq[i/4];
		for (unsigned long j = 0; j < testseq_l; j += 16) {
			unsigned char ptc[4];
			ptc[0] = testseq[j/4];		// 4 *
			ptc[1] = testseq[j/4 + 1]; // 4 *
			ptc[2] = testseq[j/4 + 2]; // 4 *
			ptc[3] = testseq[j/4 + 3]; // 4 = 16 chars
			compare(prc, ptc, R, i, j);
		}
	}
#ifdef __PROFILE__
	__PROFILE_END(Comparison_matrix)
#endif

#ifdef __DEBUG__
	mprobe(R);
#endif

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	unsigned long ccc = 0;
	RESLIST *res_list = NULL;
	// compute and output all "upper" diagonals (start horizontally with first element) ...
	for (unsigned long i = 0; i < refseq_l; i++) {
		ccc += get_diagonal(R,i,0,refseq_l,testseq_l,min_match_size, &res_list);
	}
	// ... and all "lower" diagonals (start vertically with second element)
	for (unsigned long j = 1; j < testseq_l; j++) {
		ccc += get_diagonal(R,0,j,refseq_l,testseq_l,min_match_size, &res_list);
	}
#ifdef __PROFILE__
	__PROFILE_END(Computing_diagonals)
#endif

#ifdef __DEBUG__
	//mprobe(res_list);
	printf("DEBUG: Found %lu matches in diagonal-phase\n", ccc);
#endif

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	// remove duplicates
	res_list = remove_substrings (res_list);
#ifdef __PROFILE__
	__PROFILE_END(Inclusion_Removal)
#endif

	// finally, print the results
	format_reslist(res_list);
	return 0;
}
