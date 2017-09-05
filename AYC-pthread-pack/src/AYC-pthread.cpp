//============================================================================
// Name        : AYC-pthread.cpp
// Author      : Jakob Mund, Sebastian Eder
// Version     : 1.0
// Copyright   : 
// Description : AYC contest, new implementation based on Basti's Idea
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

unsigned long lin_scan_left(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y) {
	for (unsigned long i = begin_pos; i > 0; i--) {
		char rss[4];
		char tss[4];
		__unpack(refseq[(int) floor((x+i) / 4)], rss);
		__unpack(testseq[(int) floor((y+i) / 4)], tss);
		if (rss[(x+i)%4] != tss[(y+i)%4]) {
			return i+1;
		}
	}
	return 0;
}

unsigned long lin_scan_right(unsigned long begin_pos, unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long max_steps) {
	for (unsigned long i = begin_pos; i <= max_steps; i++) {
		char rss[4];
		char tss[4];
		__unpack(refseq[(int) floor((x+i) / 4)], rss);
		__unpack(testseq[(int) floor((y+i) / 4)], tss);
		if (rss[(x+i)%4] != tss[(y+i)%4]) {
			return i-1;
		}
	}
	return max_steps;
}

unsigned long compute_diagonal(unsigned char *refseq, unsigned char *testseq, unsigned long x, unsigned long y, unsigned long dim_x, unsigned long dim_y, unsigned int min_size, RESLIST **rlist) {
#ifdef __DEBUG__
	unsigned long c = 0;
#endif
	// lower bound for iteration
	unsigned long max_steps = (dim_x-x) <= (dim_y-y) ? (dim_x-x-1) : (dim_y-y-1);
	unsigned char _rss[4], _rss2[4], _tss[4], _tss2[4];	// buffers for reference/test sub-sequences
	long lb = 0, ub = 0;
	unsigned short t_offset[4], r_offset[4];
	const unsigned int gapsize = min_size > 8 ? min_size - 8 : 9;
	v16qi v_rss, v_tss, R;
	for (unsigned long i = 0; i <= max_steps; i += (gapsize+4) * 4 ) {

		t_offset[0] = i+y+gapsize;
		t_offset[1] = i+y + ((gapsize +4) * 2) -4;
		t_offset[2] = i+y + ((gapsize +4) * 3) -4;
		t_offset[3] = i+y + ((gapsize +4) * 4) -4;

		r_offset[0] = i+x+gapsize;
		r_offset[1] = i+x + ((gapsize +4) * 2) -4;
		r_offset[2] = i+x + ((gapsize +4) * 3) -4;
		r_offset[3] = i+x + ((gapsize +4) * 4) -4;

		/*int d1 = (i+x+gapsize)%4 - (i+y+gapsize)%4,
			d2 = (i+x + ((gapsize +4) * 2) -4)%4 - (i+y + ((gapsize +4) * 2) -4)%4,
			d3 = (i+x + ((gapsize +4) * 4) -4)%4 - (i+y + ((gapsize +4) * 4) -4)%4,
			d4 = (i+x + ((gapsize +4) * 4) -4)%4 - (i+y + ((gapsize +4) * 4) -4)%4;

		if (d1 != d2 || d2 != d3 || d3 != d4 || d4 != d1) {
			printf("Delta(1) = %i. Delta (2) = %i, Delta (3) = %i, Delta (4) = %i\n", d1, d2, d3, d4);
		}*/

		_rss[0] = refseq[(unsigned long) ceil (r_offset[0] / 4) +1];
		_rss[1] = refseq[(unsigned long) ceil (r_offset[1] / 4) +1];
		_rss[2] = refseq[(unsigned long) ceil (r_offset[2] / 4) +1];
		_rss[3] = refseq[(unsigned long) ceil (r_offset[3] / 4) +1];

		_tss[0] = testseq[(unsigned long) ceil (t_offset[0] / 4) +1];
		_tss[1] = testseq[(unsigned long) ceil (t_offset[1] / 4) +1];
		_tss[2] = testseq[(unsigned long) ceil (t_offset[2] / 4) +1];
		_tss[3] = testseq[(unsigned long) ceil (t_offset[3] / 4) +1];

		_rss2[0] = refseq[(unsigned long) ceil (r_offset[0] / 4) ];
		_rss2[1] = refseq[(unsigned long) ceil (r_offset[1] / 4) ];
		_rss2[2] = refseq[(unsigned long) ceil (r_offset[2] / 4) ];
		_rss2[3] = refseq[(unsigned long) ceil (r_offset[3] / 4) ];
		__unpack16vs(_rss, _rss2, &v_rss, r_offset);

		_tss2[0] = testseq[(unsigned long) ceil (t_offset[0] / 4) ];
		_tss2[1] = testseq[(unsigned long) ceil (t_offset[1] / 4) ];
		_tss2[2] = testseq[(unsigned long) ceil (t_offset[2] / 4) ];
		_tss2[3] = testseq[(unsigned long) ceil (t_offset[3] / 4) ];
		__unpack16vs(_tss, _tss2, &v_tss, t_offset);


		R = compare16_ss4_pcmpestrm(&v_rss, &v_tss);

		if(y==770 && x == 0 && i == 0) {
			printf("t_offset[0]=%i, r_offset[0]=%i \n", t_offset[0], r_offset[0]);
			printf("matching %c%c%c%c (ref) vs. %c%c%c%c (test) (after shift)\n",
					v_rss.c[0], v_rss.c[1], v_rss.c[2], v_rss.c[3],
					v_tss.c[0], v_tss.c[1], v_tss.c[2], v_tss.c[3]
				);
		}

		if (R.c[0]  == CMP_HIT && R.c[1]  == CMP_HIT && R.c[2]  == CMP_HIT && R.c[3]  == CMP_HIT) { // 1. quadruple
			unsigned long index = i + gapsize+4;
			// get upper and lower bounds
			if(y==770 && x == 0 && i == 0) {
				printf("match! \n", t_offset[0], r_offset[0]);
			}
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size && index+y < dim_y && index+x < dim_x) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1; hd->val.t_ub = y+ub+1;
				hd->next = (*rlist); (*rlist) = hd;
#ifdef __DEBUG__
				c++;
#endif
				i = ub; continue;
			}
		}

		if (R.c[7]  == CMP_HIT && R.c[4]  == CMP_HIT && R.c[5] == CMP_HIT && R.c[6]  == CMP_HIT) { // 2. quadruple
			unsigned long index = i + (gapsize+4)*2;
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size && index+y < dim_y && index+x < dim_x) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1; hd->val.t_ub = y+ub+1;
				hd->next = (*rlist); (*rlist) = hd;
#ifdef __DEBUG__
				c++;
#endif
				i = ub; continue;
			}
		}

		if (R.c[8]  == CMP_HIT && R.c[9]  == CMP_HIT && R.c[10]  == CMP_HIT && R.c[11]  == CMP_HIT) { // 3. quadruple
			unsigned long index = i + (gapsize+4)*3;
			// get upper and lower bounds
			lb = lin_scan_left (index -4 -1, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size && index+y < dim_y && index+x < dim_x) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1; hd->val.t_ub = y+ub+1;
				hd->next = (*rlist); (*rlist) = hd;
#ifdef __DEBUG__
				c++;
#endif
				i = ub; continue;
			}
		}

		if (R.c[12]  == CMP_HIT && R.c[13]  == CMP_HIT && R.c[14]  == CMP_HIT && R.c[15]  == CMP_HIT) { // 4. quadruple
			// get upper and lower bounds
			unsigned long index = i + (gapsize+4)*4;
			lb = lin_scan_left (index-5, refseq, testseq, x, y);
			ub = lin_scan_right (index, refseq, testseq, x, y, max_steps);
			if (ub-lb > min_size && index+y < dim_y && index+x < dim_x) { // equals UB - LB +1 >= min_size
				RESLIST *hd = (RESLIST*) malloc (sizeof(RESLIST));
				hd->val.s_lb = x+lb+1; hd->val.s_ub = x+ub+1;
				hd->val.t_lb = y+lb+1; hd->val.t_ub = y+ub+1;
				hd->next = (*rlist); (*rlist) = hd;
#ifdef __DEBUG__
				c++;
#endif
				i = ub; continue;
			}
		}
	}
#ifdef __DEBUG__
	return c;
#endif
#ifndef __DEBUG__
	return 0;
#endif
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

	// parse arguments
	if (argc != 5) {
		printf("USAGE: ayc-pthread <max_threads> <min_match_size> <ref_seq> <test_seqs>\n");
		return 0;
	}
	unsigned int max_threads = atoi(argv[1]);		// maximum number of worker threads, as given as an argument
	unsigned int min_match_size = atoi(argv[2]);	// minimum length for matching
	char *reffname = argv[3];						// file names of the input files
	char *testfname = argv[4];
	FILE *reffile = NULL, *testfile = NULL;

	unsigned char *refseq, *testseq;
	unsigned long refseq_l = 0, testseq_l = 0;

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	get_sequence(reffname, reffile, &refseq, &refseq_l, min_match_size);
	get_sequence(testfname, testfile, &testseq, &testseq_l, min_match_size);
#ifdef __PROFILE__
	__PROFILE_END(Loading)
#endif

#ifdef __DEBUG__
	mprobe(testseq);
	mprobe(refseq);
	printf("DEBUG: size of input is %lu / %lu characters\n", refseq_l, testseq_l);
#endif

#ifdef __PROFILE__
	__PROFILE_START()
#endif
	unsigned long ccc = 0;
	RESLIST *res_list = NULL;
	// compute and output all "upper" diagonals (start horizontally with first element) ...
	for (unsigned long i = 0; i < refseq_l; i++) {
		ccc += compute_diagonal(refseq,testseq,i,0,refseq_l,testseq_l,min_match_size, &res_list);
	}
	// ... and all "lower" diagonals (start vertically with second element)
	for (unsigned long j = 1; j < testseq_l; j++) {
		ccc += compute_diagonal(refseq,testseq,0,j,refseq_l,testseq_l,min_match_size, &res_list);
	}
#ifdef __PROFILE__
	__PROFILE_END(Computing_diagonals)
#endif

#ifdef __DEBUG__
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
	printf("Done!");
	return 0;
}
