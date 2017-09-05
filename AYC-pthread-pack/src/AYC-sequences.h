/*
 * AYC-sequences.h
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

#ifndef AYC_SEQUENCES_H_
#define AYC_SEQUENCES_H_

#include "AYC-machine.h"

unsigned char __pack(char buf[4]);
void __unpack(unsigned char pc, char *buf);
void __unpack16(unsigned char *pc, char *buf);
void __unpack16v(unsigned char *pc, v16qi *buf);
void __unpack16vs(unsigned char *pc1, unsigned char *pc2, v16qi *buf, unsigned short *shift);

unsigned short __to_bin (char c);
char __to_char (unsigned short c);

void get_sequence(char *filename, FILE *fp, unsigned char **seq, unsigned long *length, unsigned int max_size);

#endif /* AYC_SEQUENCES_H_ */
