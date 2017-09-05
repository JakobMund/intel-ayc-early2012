/*
 * AYC-result-list.h
 *
 *  Created on: Apr 22, 2012
 *      Author: mund
 */

#include "AYC-sequences.h"

#ifndef AYC_RESULT_LIST_H_
#define AYC_RESULT_LIST_H_

bool pop_min (RESLIST **curr, RESLIST **top);
RESLIST *remove_substrings (RESLIST *rl);
RESLIST *remove_marked (RESLIST *rl);
void format_results(SEQINFO *si);
void results(SEQINFO_LIST *sil);

#endif /* AYC_RESULT_LIST_H_ */
