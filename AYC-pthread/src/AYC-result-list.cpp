#define __DEBUG__

#include <stdlib.h>
#include <stdio.h>

#include "AYC-result-list.h"

RESLIST *remove_substrings (RESLIST *outer_curr) {
	RESLIST *inner_curr = outer_curr, *pred = NULL, *first = outer_curr;

	unsigned long c = 0;
	unsigned long cc = 0;
	// O(n²), where "n" is the number of "maximal" matches (>= min_size) w.r.t. _ONE_ diagonal
	// TODO: Optimization needed (parallelization possible? Idea: sort beforehand, than split amongst processors)
	while(outer_curr) {
			while (inner_curr) {
				cc++;
				// check whether outer_curr includes inner_curr
				if (outer_curr->val.s_lb <= inner_curr->val.s_lb && outer_curr->val.s_ub >= inner_curr->val.s_ub &&
						outer_curr->val.t_lb <= inner_curr->val.t_lb && outer_curr->val.t_ub >= inner_curr->val.t_ub && inner_curr != outer_curr) {
					c++;

					// outer_curr includes inner_curr => remove inner_curr from list
					// TODO: fix memory leak
					if (pred == NULL) {		// remove the head of the list
						first = first->next;
						inner_curr = inner_curr->next;
						continue;
					} else {				// remove 'inner_cur_ from the middle of the list
						pred->next = inner_curr->next;
						inner_curr = inner_curr->next;
						continue;
					}
				}
				pred = inner_curr;
				inner_curr = inner_curr->next;
			}
			inner_curr = first; pred = NULL;
			outer_curr = outer_curr->next;
	   }
#ifdef __DEBUG__
	printf("Removed %lu contained words (using %lu iterations)\n", c, cc);
#endif

	return first;
}

void format_reslist(RESLIST *rl) {
	RESLIST *curr = rl;
	while(curr) {
		printf("%lu %lu %lu %lu\n", curr->val.s_lb, curr->val.s_ub, curr->val.t_lb, curr->val.t_ub);
		curr = curr->next;
	}
}
