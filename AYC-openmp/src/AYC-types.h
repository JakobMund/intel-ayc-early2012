/*
 * AYC-types.h
 *
 *  Created on: May 4, 2012
 *      Author: mund
 */

#ifndef AYC_TYPES_H_
#define AYC_TYPES_H_

#include <semaphore.h>

// double interval list, used to store results (including duplicates)
struct dinv_l {
	struct _dinv {
		unsigned long s_lb;
		unsigned long s_ub;
		unsigned long t_lb;
		unsigned long t_ub;
	} val;
	bool mfr; // marked for removal
	struct dinv_l *next;
};
typedef dinv_l RESLIST;

struct __seq_info {
	char name[128];
	unsigned short nl;
	unsigned long min;
	unsigned long max;
	RESLIST *rl;
};
typedef __seq_info SEQINFO;

struct _si_dllist {
	SEQINFO *si;
	struct _si_dllist *next;
};
typedef _si_dllist SEQINFO_LIST;

struct TD_compute_diagonal {
	unsigned int i;
	unsigned long length_x;
	unsigned long length_y;
	unsigned long refseq_max;
	unsigned long testseq_min;
	unsigned long testseq_max;
	RESLIST *rlist;
	RESLIST *last;

	unsigned char *refseq;
	unsigned char *testseq;
	unsigned int min_size;
	unsigned long dt;

	pthread_mutex_t *jobs_lock;
	pthread_cond_t *jobs_done;
	unsigned long *num_jobs;

	sem_t *jobs_sem;
};

struct TD_process_file {
	unsigned int f;
	unsigned int available_threads;
	char *testfname;
	SEQINFO **si;
	unsigned int min_match_size;
	unsigned char *refseq;
	unsigned long refseq_l;
	SEQINFO_LIST **sil;
	pthread_mutex_t *sil_mutex;
};

struct TD_process_sequence {
	unsigned int s;
	unsigned int f;
	unsigned long dt;
	unsigned int available_threads;
	unsigned int min_match_size;
	unsigned char *refseq;
	unsigned char *testseq;
	unsigned long refseq_l;
	SEQINFO **si;
	SEQINFO_LIST **sil;
	pthread_mutex_t *sil_mutex;
};

#endif /* AYC_TYPES_H_ */
