/*
 * AYC-sequences.cpp
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

//#define __DEBUG__

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <math.h>

#include "AYC-sequences.h"



unsigned int get_sequence(char *filename, unsigned char **seq, unsigned long *length, unsigned int max_size, SEQINFO *si) {
	int c;
	unsigned int nl = 0, seq_no = -1;
	unsigned int iarray_size;	// initial and actual size of the dynamic array
	unsigned int array_size;
	FILE *fp = NULL;

	fp = fopen(filename, "r");
	if (fp == NULL)
		return 0;

	// try to determine the size from the file, and use it as the initial array size
	int fd = fileno(fp);
    struct stat fileinfo;
    fstat(fd, &fileinfo);
    iarray_size = ceil(fileinfo.st_size) + 256;	// allocate 256 more bytes, as we may exceed the borders when comparing 16 chars at once

    // initial size of the array
    (*seq) = (unsigned char*) malloc(iarray_size);
    //posix_memalign((void **) seq, sizeof(void*), iarray_size);
    array_size = iarray_size;

    *length = 0;
    bool comment = false;
	while( (c = fgetc(fp)) != EOF ) {
		// ignore comments and white-spaces
		if (c == '>') {
			if(seq_no != -1) {
							si[seq_no].max = *length;
						}
			++seq_no;
			si[seq_no].min = (*length);
			si[seq_no].rl = NULL;	// TODO: Initialization required?
			si[seq_no].nl = 0;
			comment = true;
			continue;	// TODO: special case: two comments in a row
		} else if (c == '\n' && comment) {
			si[seq_no].name[si[seq_no].nl] = '\0';
			comment = false;
			continue;
		} else if (comment) {
			sprintf(&si[seq_no].name[si[seq_no].nl], "%c", c);
			si[seq_no].nl = si[seq_no].nl + 1;
		} else if (c != 'T' && c != 'A' && c != 'G' && c != 'C') {
			continue;
		} else {
		(*seq)[*length] = c;
		(*length)++;
		}
	}
	si[seq_no].max = *length;
	fclose(fp);
	return seq_no+1;
}
