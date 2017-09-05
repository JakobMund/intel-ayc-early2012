//#define __DEBUG__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "AYC-sequences.h"
#include "AYC-result-list.h"

RESLIST *remove_marked (RESLIST *rl) {
	RESLIST *pred = NULL, *hd = rl, *curr = rl;
	int c = 0;
	while (curr) {
		if (curr->mfr) {
			c++;
			if (pred != NULL) {
				pred->next = curr->next;
			} else {
				hd = curr->next;
			}
		} else {
			pred = curr;
		}
		curr = curr->next;
	}
	return hd;
}

RESLIST *remove_substrings (RESLIST *outer_curr) {
	RESLIST *inner_curr = outer_curr, *first = outer_curr;

#ifdef _DEBUG__
	unsigned long c = 0;
	unsigned long cc = 0;
#endif
	// O(n²), where "n" is the number of "maximal" matches (>= min_size) w.r.t. _ONE_ diagonal
	// TODO: Optimization needed (parallelization possible? Idea: sort beforehand, than split amongst processors)
	while(outer_curr) {
			while (inner_curr) {
#ifdef __DEBUG__
				cc++;
#endif
				// check whether outer_curr includes inner_curr


				long i_no 		= inner_curr->val.s_ub - inner_curr->val.s_lb +1;	// horizontal
				long o_no 		= outer_curr->val.s_ub - outer_curr->val.s_lb +1;
				long dt			= outer_curr->val.t_ub - inner_curr->val.t_ub;
				long o_no_norm_t 	= o_no - dt;
				long ds			= outer_curr->val.s_ub - inner_curr->val.s_ub;
				long o_no_norm_s 	= o_no - ds;

					/*if (inner_curr != outer_curr &&									// remove duplicates (imprecise bounds in PAR evaluation)
							inner_curr->val.s_ub == outer_curr->val.s_ub &&
							inner_curr->val.s_lb == outer_curr->val.s_lb &&
							inner_curr->val.t_ub == outer_curr->val.t_ub &&
							inner_curr->val.t_lb == outer_curr->val.t_lb
							)
					{
						inner_curr->mfr = true;
#ifdef __DEBUG__
						c++;
#endif
					} */
					 if ( inner_curr->val.s_ub +1 == outer_curr->val.s_ub - dt
							&& inner_curr->val.s_ub >= outer_curr->val.s_lb
							&& inner_curr->val.t_ub >= outer_curr->val.t_lb
							&& o_no_norm_t >= i_no
							&& dt >= 0
					) {
						inner_curr->mfr = true;
#ifdef __DEBUG__
						c++;
#endif
					}
					 if ( inner_curr->val.t_ub +1 == outer_curr->val.t_ub - ds
							&& inner_curr->val.s_ub >= outer_curr->val.s_lb
							&& inner_curr->val.t_ub >= outer_curr->val.t_lb
							&& o_no_norm_s >= i_no
							&& ds >= 0
					) {
						inner_curr->mfr = true;
#ifdef __DEBUG__
						c++;
#endif
					}
				inner_curr = inner_curr->next;
			}
			inner_curr = first;
			outer_curr = outer_curr->next;
	   }
#ifdef __DEBUG__
	printf("DEBUG: Removed %lu contained words (using %lu iterations)\n", c, cc);
#endif

	return first;
}

bool pop_min (RESLIST **rl, RESLIST **top) {
	RESLIST *curr = (*rl);
	RESLIST *local_min = curr, *local_pred = NULL, *pred = NULL, *hd = curr;
	if(curr == NULL) {
		return false;
	}
	while (curr) {
				if(curr->val.s_ub < local_min->val.s_ub	 // find the minimum in the remaining list
						|| (curr->val.s_ub == local_min->val.s_ub
								&& curr->val.t_ub < local_min->val.t_ub)) {
					local_min = curr;
					local_pred = pred;
				}
				pred = curr;
				curr = curr->next;
			}
	if(local_pred == NULL) { // head of list was minimal element
		(*top) = hd;
		(*rl) = hd->next;
		return true;
	} else {
		(*top) = local_min;
		local_pred->next = local_min->next;
		return true;
	}

}

void format_results(SEQINFO *si) {
	RESLIST *top;
	printf("%s\n", si->name);
	while (pop_min(&si->rl, &top)) {
		printf("%lu %lu %lu %lu\n", top->val.s_lb, top->val.s_ub, top->val.t_lb, top->val.t_ub);
	}
}

void results(SEQINFO_LIST *sil) {
	// get minimal (alph. ordering) test-sequence, and remove it
	SEQINFO_LIST *curr = sil, *inner_curr = sil, *hd = sil;
	SEQINFO *candidate_si = NULL;
	char candidate_name[128], cutoff_name[128];
	candidate_name[0] = '\0'; cutoff_name[0] = '\0';
	while(curr) {
		candidate_name[0] = '\0';
		while (inner_curr) {
			if (   (cutoff_name[0] == '\0' ||
					strcmp(cutoff_name, inner_curr->si->name) < 0)
				&& (candidate_name[0] == '\0'
						|| strcmp(inner_curr->si->name, candidate_name) < 0) ) {

					candidate_si = inner_curr->si;
					sprintf(candidate_name, "%s", inner_curr->si->name );
				}
			inner_curr = inner_curr->next;
		}
		sprintf(cutoff_name, "%s", candidate_name);
		if(candidate_si) {
			format_results(candidate_si);
		}
		inner_curr = hd;
		curr = curr->next;
	}
}
