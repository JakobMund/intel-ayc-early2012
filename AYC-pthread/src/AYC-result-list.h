/*
 * AYC-result-list.h
 *
 *  Created on: Apr 22, 2012
 *      Author: mund
 */

#ifndef AYC_RESULT_LIST_H_
#define AYC_RESULT_LIST_H_

// double interval list, used to store results (including duplicates)
struct dinv_l {
	struct _dinv {
		unsigned long s_lb;
		unsigned long s_ub;
		unsigned long t_lb;
		unsigned long t_ub;
	} val;
	struct dinv_l *next;
};

typedef dinv_l RESLIST;

RESLIST *remove_substrings (RESLIST *outer_curr);
void format_reslist(RESLIST *rl);

#endif /* AYC_RESULT_LIST_H_ */
