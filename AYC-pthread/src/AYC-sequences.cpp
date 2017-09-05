/*
 * AYC-sequences.cpp
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <math.h>
#include <mcheck.h>

#define BIN_A 0
#define BIN_T 1
#define BIN_G 2
#define BIN_C 3

inline unsigned short __to_bin (char c) {
	switch(c) {
	case 'A': return BIN_A; break;
	case 'T': return BIN_T; break;
	case 'G': return BIN_G; break;
	case 'C': return BIN_C; break;
	}
	return 0;
}

inline char __to_char (unsigned short c) {
	switch (c) {
	case BIN_A: return 'A'; break;
	case BIN_T: return 'T'; break;
	case BIN_G: return 'G'; break;
	case BIN_C: return 'C'; break;
	}
	return 'F';
}

inline unsigned char __pack(char buf[4]) {
	return ((char) (__to_bin(buf[0]) << 6 )
				 + (__to_bin(buf[1]) << 4 )
				 + (__to_bin(buf[2]) << 2 )
				 +  __to_bin(buf[3]));
}

void __unpack(unsigned char pc, char *buf) {
	// TODO: asm might be able to improve this part
	char lbuf[4];
	lbuf[0] = __to_char(pc >> 6);
	lbuf[1] = __to_char((pc >> 4) % 4);
	lbuf[2] = __to_char((pc >> 2) % 4);
	lbuf[3] = __to_char(pc % 4);
	memcpy(buf, lbuf, 4);
	return;
}

void __unpack16(unsigned char *pc, char *buf) {
	// TODO: same as for __unpack: asm might be able to improve this part
	char lbuf[16];
	lbuf[0] = __to_char(pc[0] >> 6);
	lbuf[1] = __to_char((pc[0] >> 4) % 4);
	lbuf[2] = __to_char((pc[0] >> 2) % 4);
	lbuf[3] = __to_char(pc[0] % 4);

	lbuf[4] = __to_char(pc[1] >> 6);
	lbuf[5] = __to_char((pc[1] >> 4) % 4);
	lbuf[6] = __to_char((pc[1] >> 2) % 4);
	lbuf[7] = __to_char(pc[1] % 4);

	lbuf[8] = __to_char(pc[2] >> 6);
	lbuf[9] = __to_char((pc[2] >> 4) % 4);
	lbuf[10] = __to_char((pc[2] >> 2) % 4);
	lbuf[11] = __to_char(pc[2] % 4);

	lbuf[12] = __to_char(pc[3] >> 6);
	lbuf[13] = __to_char((pc[3] >> 4) % 4);
	lbuf[14] = __to_char((pc[3] >> 2) % 4);
	lbuf[15] = __to_char(pc[3] % 4);

	memcpy(buf, lbuf, 16);
	return;
}

void get_sequence(char *filename, FILE *fp, unsigned char **seq, unsigned long *length) {
	int c;
	int i = 0;
	unsigned int iarray_size;	// initial and actual size of the dynamic array
	unsigned int array_size;
	char buf[4];

	fp = fopen(filename, "r");
	if (fp == NULL)
		return;

	// try to determine the size from the file, and use it as the initial array size
	int fd = fileno(fp);
    struct stat fileinfo;
    fstat(fd, &fileinfo);
    iarray_size = ceil(fileinfo.st_size / 4);

    // initial size of the array
    (*seq) = (unsigned char*) malloc(iarray_size);
    array_size = iarray_size;

    *length = 0;
    bool comment = false;
	while( (c = fgetc(fp)) != EOF ) {
		// ignore comments and white-spaces
		if (c == '>') {
			comment = true; continue;
		} else if (c == '\n' && comment) {
			comment = false; continue;
		}
		if (comment || (c != 'T' && c != 'A' && c != 'G' && c != 'C')) {
			continue;
		}
		buf[*length % 4] = c;
		(*length)++;
		// when buffer is filled with 4 chars, pack them
		if (*length % 4 == 0) {
			// but first check whether array is sufficiently large
			// TODO: check if this really works
			if ( array_size < (*length/4)+1 ) {
				array_size *= 2;
#ifdef __DEBUG__
				printf("DEBUG: updating array size to %lu\n", array_size);
#endif
				(*seq) = (unsigned char*) realloc(*seq, array_size * sizeof(unsigned char));
			}

			(*seq)[i] = __pack(buf);
			i++;
		}
	}
	// pack the last 0-3 chars (may include some dirty bits), if necessary
	if (*length % 4 != 0) {
		(*seq)[i] = __pack(buf);
	}

	fclose(fp);
	return;
}
