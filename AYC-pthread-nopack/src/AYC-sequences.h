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
#include <nmmintrin.h>	// SSE 4.2

/*
 * get_refsequences reads filename and parses it as a reference sequence,
 * returning the data in **seq and its length in *length
 *
 * the parsing itself is accomplished using STTNI instructions, where 16 bytes
 * are read at once, and if no Line-Feed (\n, ASCII 0x0A) character is present,
 * the chunk is copied at once (MOVDQU).
 */
void get_refsequence(char *filename, unsigned char **seq, unsigned long *length);

/*
 * get_refsequences reads filename and parses it as a reference sequence,
 * returning the data in **seq. Further sequence information (name, starting
 * index, ending index) is stored in *si, and *seq_total is increased by the total
 * number of sequences contained.
 *
 * the parsing itself is accomplished using STTNI instructions, where 16 bytes
 * are read at once, and if no Line-Feed (\n, ASCII 0x0A) character is present,
 * the chunk is copied at once (MOVDQU).
 */
unsigned int get_testsequence(char *filename, unsigned char **seq, SEQINFO *si, unsigned int* seq_total);

#endif /* AYC_SEQUENCES_H_ */
