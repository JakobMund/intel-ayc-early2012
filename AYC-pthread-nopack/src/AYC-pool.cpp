/*
 * AYC-pool.cpp
 *
 * Implementation of a POSIX treads-based thread
 * (worker) pool, with bundles of jobs submitted at
 * a time, and an implicit barrier after all jobs have
 * been completed.
 *
 * The implementation does NOT use persistent (reuse of)
 * threads, but creates them when necessary. This provides
 * accurate %CPU-usage numbers, and does not influence
 * the main execution thread at all. The number of threads
 * created is linear in the number of cores available,
 * and the threads are directly mapped to physical cores
 * (using CPU affinity hints).
 *
 *  Created on: May 11, 2012
 *      Author: mund
 */

//#define __DEBUG__

#include <stdlib.h>
#include <pthread.h>
#include <errno.h>

#include <stdio.h>

#include "AYC-pool.h"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __DEBUG__
#include <stdio.h>
static unsigned long __DEBUG_THREAD_ID = 0;
#endif

#ifndef MIN
	#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/*
 * Create a worker pool containing num_threads
 * threads, which shall be kept alive after they
 * finished all jobs iff keep_alive is set.
 * (not implemented yet).
 */
tpool_t create_pool(unsigned int num_threads, bool keep_alive) {
	tpool_t tp;
	tp.threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
	tp.num_threads = num_threads;
	tp.keep_alive = keep_alive;
	return tp;
}

/*
 * Destroy a given thread pool by freeing memory.
 */
void destroy_pool(tpool_t *tp) {
	free(tp->threads);
}

/*
 * Each thread is created as a worker, which waits for jobs
 * until no more jobs are queued, then exits.
 * A spin lock is used for synchronization, as the time spend
 * in the critical section is not worth doing a context switch
 * (as in an ordinary mutex).
 */
void *worker(void *arg) {
	tworker_t *winfo = (tworker_t *) arg;
	int c = -1;
	tjob_t *curr_job = NULL;
	pthread_mutex_lock(winfo->start_mutex);
	pthread_mutex_unlock(winfo->start_mutex);
	while(1) {
		// check if all jobs done, and if not, dequeue a job
		pthread_spin_lock(winfo->jobs_lock);
		c = *(winfo->remaining);
		if (c <= 0 && !winfo->keep_alive) {
#ifdef __DEBUG__
			printf("DEBUG: Worker Thread (TPool=%u) shutting down\n", __DEBUG_THREAD_ID);
#endif
			pthread_spin_unlock(winfo->jobs_lock);
			pthread_exit (NULL);
		} else {
			curr_job = &winfo->jobs->jobs[c-1];
			*winfo->remaining = c - 1; // dequeue
			pthread_spin_unlock(winfo->jobs_lock);
		}

		// run the job
		(*curr_job->job)(curr_job->jarg);

		// check if we were the last thread
		if (c <= 1) {
#ifdef __DEBUG__
			printf("DEBUG: Worker Thread (TPool=%u) signals completion\n", __DEBUG_THREAD_ID);
#endif
			pthread_cond_signal(winfo->jobs_complete);
			break;
		}
	}
	return NULL;
}

/*
 * Run a set of jobs in a thread pool, and wait for
 * completion of all jobs when done (implicit barrier).
 * Destroy the pool afterwards only if keep_alive is
 * true (not implemented yet).
 */
void run_jobs(tjoblist_t *jobs, tpool_t *tp) {
	pthread_mutex_t start_mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t jobs_complete = PTHREAD_COND_INITIALIZER;
	pthread_spinlock_t jobs_lock;
	pthread_spin_init(&jobs_lock, PTHREAD_PROCESS_PRIVATE);
	pthread_attr_t worker_attr;	// attributes of the worker threads
	pthread_attr_init(&worker_attr);
	pthread_attr_setdetachstate(&worker_attr, PTHREAD_CREATE_JOINABLE);
	tworker_t warg[tp->num_threads];
	unsigned int jobs_remaining = jobs->num_jobs;
	pthread_mutex_lock(&start_mutex);

#ifdef __DEBUG__
	__DEBUG_THREAD_ID++;
	printf("DEBUG: Creating TPool(ID=%u) spawning %u threads for %u jobs\n", __DEBUG_THREAD_ID, tp->num_threads, jobs->num_jobs);
#endif
    cpu_set_t cpuset;
	for (unsigned int i = 0; i < MIN(tp->num_threads, jobs->num_jobs); i++) {
		warg[i].jobs = jobs;
		warg[i].jobs_complete = &jobs_complete;
		warg[i].start_mutex = &start_mutex;
		warg[i].jobs_lock = &jobs_lock;
		warg[i].keep_alive = tp->keep_alive;
		warg[i].remaining = &jobs_remaining;
		pthread_create(&tp->threads[i], &worker_attr, &worker, (void*) &warg[i]);
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        pthread_setaffinity_np(tp->threads[i], sizeof(cpu_set_t), &cpuset);
	}
	pthread_attr_destroy(&worker_attr);
#ifdef __DEBUG__
	printf("DEBUG: TPool(ID=%u) waiting for jobs to complete...\n", __DEBUG_THREAD_ID);
#endif
	pthread_cond_wait(&jobs_complete, &start_mutex);
	pthread_mutex_unlock(&start_mutex);
#ifdef __DEBUG__
	printf("DEBUG: TPool(ID=%u) waiting for threads to finish..\n", __DEBUG_THREAD_ID);
#endif
	for (unsigned int i = 0; i < MIN(tp->num_threads, jobs->num_jobs); i++) {
		pthread_join(tp->threads[i], NULL);
	}
#ifdef __DEBUG__
	printf("DEBUG: TPool(ID=%u) finished all jobs\n", __DEBUG_THREAD_ID);
#endif
	if (!tp->keep_alive) {
		pthread_cond_destroy(&jobs_complete);
		pthread_mutex_destroy(&start_mutex);
		pthread_spin_destroy(&jobs_lock);
		destroy_pool (tp);
	}
}

/*
 * Run a set of jobs in parallel using num_threads
 * threads.
 */
void run_parallel(tjoblist_t *jobs, unsigned int num_threads) {
	tpool_t tp = create_pool(num_threads, false);
	run_jobs(jobs, &tp);
}
