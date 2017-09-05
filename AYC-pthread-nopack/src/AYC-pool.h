/*
 * AYC-pool.h
 *
 *  Created on: May 11, 2012
 *      Author: mund
 */

#ifndef AYC_POOL_H_
#define AYC_POOL_H_

#include <pthread.h>

typedef struct _tpool_t {
	unsigned int num_threads;
	pthread_t *threads;
	bool keep_alive;
} tpool_t;

typedef struct _tjob_t {
	void (*job)(void *);
	void *jarg;
} tjob_t;

typedef struct _tjoblist_t {
	unsigned int num_jobs;
	tjob_t *jobs;
} tjoblist_t;

typedef struct _tworker_t {
	tjoblist_t *jobs;
	unsigned int *remaining;
	bool keep_alive;
	pthread_mutex_t	*start_mutex;
	pthread_spinlock_t *jobs_lock;
	pthread_cond_t *jobs_complete;
} tworker_t;

/*
 * Create a worker pool containing num_threads
 * threads, which shall be kept alive after they
 * finished all jobs iff keep_alive is set.
 * (not implemented yet).
 */
tpool_t create_pool(unsigned int num_threads, bool keep_alive);

/*
 * Destroy a given thread pool by freeing memory.
 */
void destroy_pool(tpool_t *tp);

/*
 * Run a set of jobs in a thread pool, and wait for
 * completion of all jobs when done (implicit barrier).
 * Destroy the pool afterwards only if keep_alive is
 * true (not implemented yet).
 */
void run_jobs(tjoblist_t *jobs, tpool_t *tp);

/*
 * Run a set of jobs in parallel using num_threads
 * threads.
 */
void run_parallel(tjoblist_t *jobs, unsigned int num_threads);

#endif /* AYC_POOL_H_ */
