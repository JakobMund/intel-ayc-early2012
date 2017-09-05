/*
 * AYC-sequences.h
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

#ifndef AYC_SEQUENCES_H_
#define AYC_SEQUENCES_H_

unsigned char __pack(char buf[4]);
void __unpack(unsigned char pc, char *buf);
void __unpack16(unsigned char *pc, char *buf);
unsigned short __to_bin (char c);
char __to_char (unsigned short c);
void get_sequence(char *filename, FILE *fp, unsigned char **seq, unsigned long *length);

#endif /* AYC_SEQUENCES_H_ */
