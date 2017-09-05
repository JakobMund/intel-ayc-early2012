/*
 * AYC-pthread.cpp
 *
 *  Created on: May 6, 2012
 *      Author: mund
 */

//#define __DEBUG__
//#define __PROFILE__

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

#define MAX_TEST_SEQUENCES 128							// maximal number of test_sequences per file
#define COMPUTE_DIAGONAL compute_diagonal				// diagonal function for min_match_size above a given threshold
#define COMPUTE_DIAGONAL_SHORT compute_diagonal_4		// diagonal function to use for inputs with min_match_size below a given threshold
#define TASK_MULTFAC 64									// number of jobs := number of threads * TASK_MULTFAC

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#include "AYC-types.h"
#include "AYC-sequences.h"
#include "AYC-comparison.h"
#include "AYC-result-list.h"
#include "AYC-pthread.h"
#include "AYC-pool.h"

#ifndef MIN
	#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/*
 * Individual job computing a _subset_ of all diagonals (alignments)
 * between reference and test sequences.
 */
void job_compute_diagonal(void *targ) {
	struct TD_compute_diagonal *td;
	td = (struct TD_compute_diagonal *) targ;
	RESLIST *last = NULL, *curr;

	// compute this thread's partition offsets
	unsigned long offset_x_lb = td->i * td->length_x;
	unsigned long offset_x_ub = MIN(offset_x_lb + td->length_x, td->refseq_max);

	unsigned long offset_y_lb = td->testseq_min + (td->i * td->length_y) +1;
	unsigned long offset_y_ub = MIN(offset_y_lb + td->length_y, td->testseq_max);

#ifdef __DEBUG__
	printf("DEBUG (Thread %i): Range(X)=%u..%u (%u), Range(Y)=%u..%u (%u) \n", td->i,
			offset_x_lb, offset_x_ub-1, td->length_x,
			offset_y_lb, offset_y_ub-1, td->length_y);
#endif

	// base on the matching size, use the appropriate functions
	if (td->min_size >= 8) {
		// compute the thread' part of "upper" diagonals ...
		for (unsigned long i = offset_x_lb; i < offset_x_ub; i++) {
			COMPUTE_DIAGONAL(td->refseq,td->testseq,i, td->testseq_min,td->refseq_max,td->testseq_max,td->min_size, td->rlist, td->dt);
		}
		// ... and all "lower" diagonals (start vertically with second element)
		for (unsigned long j = offset_y_lb; j < offset_y_ub; j++) {
			COMPUTE_DIAGONAL(td->refseq,td->testseq,0,j,td->refseq_max,td->testseq_max,td->min_size, td->rlist, td->dt);
		}
	} else {
		// compute the thread' part of "upper" diagonals ...
		for (unsigned long i = offset_x_lb; i < offset_x_ub; i++) {
			COMPUTE_DIAGONAL_SHORT(td->refseq,td->testseq,i, td->testseq_min,td->refseq_max,td->testseq_max,td->min_size, td->rlist, td->dt);
		}
		// ... and all "lower" diagonals (start vertically with second element)
		for (unsigned long j = offset_y_lb; j < offset_y_ub; j++) {
			COMPUTE_DIAGONAL_SHORT(td->refseq,td->testseq,0,j,td->refseq_max,td->testseq_max,td->min_size, td->rlist, td->dt);
		}
	}

	curr = *td->rlist;
	while (curr) {
		last = curr;
		curr = curr->next;
	}
	*td->last = last;
	return;
}

/*
 * jobs for POSTprocessing a test sequence, i.e.,
 * for concatenating results from the individual jobs
 * computing subsets of the diagonals, and afterwards,
 * removing substrings (contained strings, like
 * implemented in the reference implementation).
 *
 * The results are storted in a SEQINFO_LIST data structure.
 */
void job_process_sequence(void *arg)
{
	struct TD_process_sequence *td;
	td = (struct TD_process_sequence *) arg;
	SEQINFO **si = td->si;
	SEQINFO_LIST *c_sil;

	RESLIST *hd = NULL, *last = NULL, *ilast = NULL;
	for (unsigned int i = 0; i < td->res_count; i++) {
		ilast = td->res_last[td->dn+i];
			if( td->res_hd[td->dn+i] != NULL) {
				if(hd == NULL)
					hd = td->res_hd[td->dn+i];
				if (last != NULL)
					last->next = td->res_hd[td->dn+i];
				last = ilast;
			}
		}

	si[td->f][td->s].rl = hd;

	// remove duplicates, and sort output
	si[td->f][td->s].rl = remove_substrings (si[td->f][td->s].rl);
	si[td->f][td->s].rl = remove_marked (si[td->f][td->s].rl);

	// buffer results until evaluation completes
	// serialize access to the SIList
	c_sil = (SEQINFO_LIST *) malloc (sizeof(SEQINFO_LIST));
	c_sil->si = &si[td->f][td->s];
	pthread_spin_lock(td->sil_lock);
	c_sil->next = *(td->sil);
	*(td->sil) = c_sil;
	pthread_spin_unlock(td->sil_lock);

	return;
}

/*
 * job for reading and parsing a file containing
 * test sequences.
 */
void job_process_file(void *arg)
{
	struct TD_process_file *td;
	td = (struct TD_process_file *) arg;
	SEQINFO **si = td->si;
	si[td->f] = (SEQINFO *) malloc (MAX_TEST_SEQUENCES * sizeof(SEQINFO));
	(*td->seq_no) = get_testsequence(td->testfname, td->testseq, si[td->f], td->seq_total);
#ifdef __DEBUG__
	printf("DEBUG: Size of test-sequences from file %s: %lu\n", td->testfname, si[td->f]->max - si[td->f]->min);
#endif
	return;
}

/*
 * job read reading and parsing a file containing
 * the reference sequence.
 */
void job_process_ref(void *arg)
{
	struct TD_process_ref *td;
	td = (struct TD_process_ref *) arg;
	get_refsequence(td->filename, td->seq, td->length);
#ifdef __DEBUG__
	printf("DEBUG: size of reference sequence is %u\n", td->length);
#endif
	return;
}

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
void par_process_files(unsigned short num_threads, char *reffname, char **argv, SEQINFO **si, unsigned short num_test_files, unsigned int min_match_size, SEQINFO_LIST **sil)
{
	// read as many test-sequences in parallel as possible
	tjoblist_t joblist_rf;
	tjob_t job_rf[num_test_files+1];
	unsigned int seq_no[num_test_files];
	struct TD_process_file td[num_test_files];
	struct TD_process_ref tdr;
	pthread_spinlock_t sil_lock;
	pthread_spin_init(&sil_lock, PTHREAD_PROCESS_PRIVATE);
	unsigned char **testseq = (unsigned char **) malloc (sizeof(unsigned char*) * num_test_files);
	unsigned char *refseq;
	unsigned long refseq_l;
	unsigned int num_seq_total = 0;

	// for each file containing test sequences => create a job
	for (unsigned short f = 0; f < num_test_files; f++) { // for each test file
		td[f].f = f;
		td[f].testfname = argv[f+4];
		td[f].testseq = &testseq[f];
		td[f].seq_no = &seq_no[f];
		td[f].seq_total = &num_seq_total;
		td[f].si = si;
		job_rf[f+1].job 	= &job_process_file;
		job_rf[f+1].jarg 	= (void*) &td[f];
	}
	tdr.filename = reffname; // add one job for the reference file (to the beginning of the job-list)
	tdr.length = &refseq_l;
	tdr.seq = &refseq;
	job_rf[0].job 	= &job_process_ref;
	job_rf[0].jarg = (void*) &tdr;

	// create and execute the job list, then wait for completion
	joblist_rf.jobs = job_rf;
	joblist_rf.num_jobs = num_test_files+1;
	run_parallel(&joblist_rf, num_threads);

#ifdef __DEBUG__
	printf("DEBUG: Total #sequences is %u\n", num_seq_total);
#endif

	tjoblist_t joblist_cd, joblist_rs;
	unsigned long num_cd_jobs = num_threads * TASK_MULTFAC;
	unsigned int cd_per_seq = ceil (num_cd_jobs / num_seq_total) > 1 ? ceil (num_cd_jobs / num_seq_total) : 1;
#ifdef __DEBUG__
	printf("DEBUG: num_cd_jobs=%u, cd_per_seq*num_seq_total=%u\n", num_cd_jobs, cd_per_seq*num_seq_total);
#endif
	unsigned long length_x = ceil ((double) refseq_l / cd_per_seq);
	tjob_t job_cd[num_cd_jobs];				// jobs for computing diagonals
	tjob_t job_rs[num_seq_total];			// jobs for	removing subsequences
											// TODO: parallelize this further (by partitioning the list)
	unsigned int sn = 0, dn = 0;
	unsigned long dt;
	struct TD_compute_diagonal cd_arg[num_cd_jobs];
	RESLIST *res_hd[num_cd_jobs], *res_last[num_cd_jobs];

	// partitioning into num_threads parts, and start execution
	for (unsigned short f = 0; f < num_test_files; f++) {
		dt = 0;
		for (unsigned int s = 0; s < seq_no[f]; s++) {
			unsigned long length_y = ceil ( (si[f][s].max - si[f][s].min) / cd_per_seq);
			for (unsigned int i = 0; i < cd_per_seq; i++) {
				cd_arg[dn+i].i 				= i;
				cd_arg[dn+i].length_x		= length_x;
				cd_arg[dn+i].length_y		= length_y;
				cd_arg[dn+i].refseq_max		= refseq_l;
				cd_arg[dn+i].testseq_min	= si[f][s].min;
				cd_arg[dn+i].testseq_max	= si[f][s].max;
				cd_arg[dn+i].refseq 		= refseq;
				cd_arg[dn+i].testseq 		= testseq[f];
				cd_arg[dn+i].min_size		= min_match_size;
				cd_arg[dn+i].dt				= dt;
				res_hd[dn+i]				= NULL;
				cd_arg[dn+i].rlist			= &res_hd[dn+i];
				cd_arg[dn+i].last			= &res_last[dn+i];

				job_cd[dn+i].job 	= &job_compute_diagonal;
				job_cd[dn+i].jarg 	= (void*) &cd_arg[dn+i];
			}
			dn += cd_per_seq;
			dt = si[f][s].max;
		}
	}

	// run the jobs (computing diagonals, removing substrings)
	joblist_cd.jobs = job_cd; joblist_cd.num_jobs = dn; run_parallel(&joblist_cd, num_threads);

	struct TD_process_sequence tds[num_seq_total];
	dn = 0; sn = 0;

	for (unsigned short f = 0; f < num_test_files; f++) {
		dt = 0;
		for (unsigned int s = 0; s < seq_no[f]; s++) {
			// when we're at it, add the jobs for the (parallel) removal of substrings
			tds[sn].res_hd		= res_hd;
			tds[sn].res_last	= res_last;
			tds[sn].res_count 	= cd_per_seq;
			tds[sn].dn 			= dn;
			tds[sn].f = f;
			tds[sn].s = s;
			tds[sn].si = si;
			tds[sn].sil_lock = &sil_lock;
			tds[sn].sil = sil;

			job_rs[sn].job = &job_process_sequence;
			job_rs[sn].jarg = (void*) &tds[sn];

			dn += cd_per_seq;
			dt = si[f][s].max;
			sn++;
		}
	}
	// run the jobs (computing diagonals, removing substrings)
	joblist_rs.jobs = job_rs; joblist_rs.num_jobs = sn; run_parallel(&joblist_rs, num_threads);
}
