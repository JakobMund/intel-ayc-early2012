/*
 * AYC-types.h
 *
 * Definition of data types (mainly structs),
 * including representation of results,
 * and data structures to be passed to threads
 * and created for jobs passed to threads
 *
 *  Created on: May 4, 2012
 *      Author: mund
 */

#ifndef AYC_TYPES_H_
#define AYC_TYPES_H_

#include <pthread.h>

/*
 * double interval list, used to store results (including duplicates)
 */
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

/*
 * information about a test sequence
 */
struct __seq_info {
	char name[110];		// 110 chars, because __seq_info is 128 byte, i.e., 2^7
	unsigned short nl;
	unsigned long min;
	unsigned long max;
	RESLIST *rl;
};
typedef __seq_info SEQINFO;

/*
 * linked list of test sequence information
 */
struct _si_dllist {
	SEQINFO *si;
	struct _si_dllist *next;
};
typedef _si_dllist SEQINFO_LIST;

/*
 * thread data - computing subsets of diagonals
 */
struct TD_compute_diagonal {
	unsigned int i;				// in
	unsigned long length_x;
	unsigned long length_y;
	unsigned long refseq_max;
	unsigned long testseq_min;
	unsigned long testseq_max;
	unsigned char *refseq;
	unsigned char *testseq;
	unsigned int min_size;
	unsigned long dt;

	RESLIST **rlist;			// out
	RESLIST **last;
};

/*
 * Thread Data - processing test sequence files
 */
struct TD_process_file {
	unsigned int f;
	char *testfname;
	unsigned char **testseq;
	unsigned int *seq_no;
	unsigned int *seq_total;
	SEQINFO **si;
};

/*
 * Thread Data - processing the reference sequence file
 */
struct TD_process_ref {
	char *filename;
	unsigned char **seq;
	unsigned long *length;
	unsigned long *__DUMMY__;
};

/*
 * Thread Data - post-processing individual sequences
 */
struct TD_process_sequence {
	unsigned int s;
	unsigned int f;
	unsigned int res_count;
	RESLIST **res_hd;
	RESLIST **res_last;
	unsigned int dn;

	SEQINFO **si;
	SEQINFO_LIST **sil;
	pthread_spinlock_t *sil_lock;
};

#endif /* AYC_TYPES_H_ */
