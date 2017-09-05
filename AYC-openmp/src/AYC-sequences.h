/*
 * AYC-sequences.h
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

//#include <string>

#ifndef AYC_SEQUENCES_H_
#define AYC_SEQUENCES_H_

#include "AYC-types.h"

unsigned int get_sequence(char *filename, unsigned char **seq, unsigned long *length, unsigned int max_size, SEQINFO *si);

#endif /* AYC_SEQUENCES_H_ */
