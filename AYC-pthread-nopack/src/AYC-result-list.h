/*
 * AYC-result-list.h
 *
 *  Created on: Apr 22, 2012
 *      Author: mund
 */

#include "AYC-sequences.h"

#ifndef AYC_RESULT_LIST_H_
#define AYC_RESULT_LIST_H_

/*
 * pop_min removes the smallest element from the list
 * and returns it in top. Furthermore, the return value
 * determines whether a the list is empty (false) or
 * not (true). The sorting is based on the reference
 * implementation, i.e., order by the reference
 * upper bound, than by the test sequence upper bound.
 */
bool pop_min (RESLIST **curr, RESLIST **top);

/*
 * Checks whether some matchings are contained in other
 * alignments (offsets) of the reference and test sequence,
 * and shall thus be removed (to generate the exact same
 * output as the reference implementation).
 *
 * Note that this can be easily be parallelized by a variable
 * number of threads, as the result list can be partitioned
 * among them. Only the update of the "mark-for-removal" flag
 * must be protected, e.g. by a mutex.
 *
 * However, running this function in parallel was not beneficial
 * in the given benchmark suite, probably because the result set was
 * too small.
 */
RESLIST *remove_substrings (RESLIST *rl);

/*
 * Remove items from the results list which were previously marked
 * as "to be removed" (mfr == true). This task can be run on each
 * test sequence individually, but can not be directy parallelized
 * further.
 */
RESLIST *remove_marked (RESLIST *rl);

/*
 * Formats the output according to the original
 * reference implementation.
 */
void format_results(SEQINFO *si);

/*
 * Sort the output by the name of the test sequences,
 * and then print the corresponding results of each sequence
 * using the same ordering and format as the reference
 * implementation.
 */
void results(SEQINFO_LIST *sil);

#endif /* AYC_RESULT_LIST_H_ */
