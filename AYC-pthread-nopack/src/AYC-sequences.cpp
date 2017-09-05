/*
 * AYC-sequences.cpp
 *
 * Reading and parsing input files, i.e.,
 * reference sequences and test sequences,
 * and storing them in data structures that
 * are henceforth used.
 *
 *  Created on: Apr 21, 2012
 *      Author: mund
 */

//#define __DEBUG__
#define FREAD_CHUNK 8192			// 8 * 1024 (8KByte)
#define FREAD_BUFF 8388608

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "AYC-machine.h"
#include "AYC-sequences.h"

//#define __PROFILE__

#ifdef __PROFILE__
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#define __PROFILE_START() 		gettimeofday(&start, NULL);
#define __PROFILE_END( MSG ) 	gettimeofday(&end, NULL);\
			seconds  = end.tv_sec  - start.tv_sec; useconds = end.tv_usec - start.tv_usec;\
			mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;\
			printf("PROFILE: %s %ld ms\n", #MSG, mtime );
	struct timeval start, end;
	long mtime, seconds, useconds;
#ifdef __DEBUG__
	printf("WARNING: Profiling results may be flawed, as __DEBUG__ is active\n");
#endif
#endif

/*
 * get_refsequences reads filename and parses it as a reference sequence,
 * returning the data in **seq and its length in *length
 *
 * the parsing itself is accomplished using STTNI instructions, where 16 bytes
 * are read at once, and if no Line-Feed (\n, ASCII 0x0A) character is present,
 * the chunk is copied at once (MOVDQU).
 */
void get_refsequence(char *filename, unsigned char **seq, unsigned long *length) {
	unsigned int array_size;
	FILE *fp = NULL;
	fp = fopen(filename, "r");
	//setvbuf (fp, NULL, _IOFBF, FREAD_BUFF);

	// try to determine the size from the file, and use it as the initial array size
	int fd = fileno(fp);
    struct stat fileinfo;
    fstat(fd, &fileinfo);
    array_size = fileinfo.st_size + 256;	// allocate 256 more bytes, as we may exceed the borders when comparing 16 chars at once

    (*seq) = (unsigned char*) malloc(array_size);
    char *buffer = (char*) malloc(array_size);
    unsigned int read_sum = 0;
    unsigned int read = 1;
    while(read > 0) {
        read = fread(buffer+read_sum, 1, array_size, fp);
        read_sum += read;
    }
    fclose(fp);

    unsigned int inval_chars = 0;
    unsigned int i = 0;
    unsigned int pos_nl = 0;

#ifdef __PROFILE__
    __PROFILE_START()
#endif

    if (buffer[0] == '>') { // ignore the name for the reference sequence
    	while(i < read_sum) {
    		if (buffer[i] == '\n') {
    			i++; inval_chars++; break;
    		}
    		i++; inval_chars++;
    	}
	}

	while (i < read_sum) {
    	__asm__ (
    			"movdqu %2,%%xmm0;" 				// move unaligned 16 byte from the read-buffer into xmm0
    			"mov $0x0a, %%ecx;"					// move '\n' (0x0a) into xmm1
    			"pinsrb $0x00,%%ecx,%%xmm1;"
    			"xor %%ecx,%%ecx;"					// clear ECX
    			"pcmpistri $0x00, %%xmm0, %%xmm1;"	// find first occurence of '\n' in the buffer at position i
    			"cmp $16, %%ecx;"
    			"movdqu %%xmm0, %3;"
    			"jne nocpy;"
    			"add $16, %1;"
    "nocpy: 	 mov %%ecx,%0;"
    			:"=r"(pos_nl),					// result %0: index of
    			 "=m"(i)
    			:"m"(buffer[i]), 				// operand %1: address of the buffer at position i
    			 "m"((*seq)[i-inval_chars])		// operand %2: destination buffer (i minus number of invalid chars so far)
    			: "%xmm0", "%xmm1", "%ecx"
    	);

    	if (pos_nl != 16) {
    		i += pos_nl + 1; inval_chars++;
    	}
   }

	/*__asm__ __volatile__ (
					"movl %0, %%eax;"					// i
					"movl %2, %%edx;"					// BASE_ADDR(buffer)
					"movl %3, %%ebx;"					// BASE_ADDR(*seq)
					"subl %1, %%ebx;"					// substract inval_chars from BASE_ADDR(*seq)

		"loop: 		 movdqu (%%edx,%%eax,),%%xmm0;" 	// move unaligned 16 byte from the read-buffer into xmm0
					"movl $0x0a, %%ecx;"				// move '\n' (0x0a)  into xmm1
					"pinsrb $0x00,%%ecx,%%xmm1;"
					"xor %%ecx,%%ecx;"					// clear ECX
					"pcmpistri $0x00, %%xmm0, %%xmm1;"	// find first occurence of '\n' in the buffer at position i
					"movdqu %%xmm0, (%%ebx,%%eax,);"

					"cmp $16, %%ecx;"
					"je no_nl;"
					"add $1, %%ecx;"
					"sub $1, %%ebx;"
		"no_nl:		 add %%ecx, %%eax;"

					"cmp %%eax, %4;"
					"jb loop;"
					:
					:"m"(i),
					 "m"(inval_chars),
					 "m"(buffer), 						// operand %1: address of the buffer at position i
					 "m"((*seq)),						// operand %2: destination buffer (i minus number of invalid chars so far)
					 "m"(read_sum)
					: "%xmm0", "%xmm1", "%ecx", "%eax", "%ebx", "%edx"
			); */

#ifdef __PROFILE__
    __PROFILE_END(REF-File loading..)
#endif

    *length = read_sum-inval_chars;

   /* if (*length > 1000000000) {
    	printf("ALL YOUR BASE BELONG TO BISON!\n");
    } */
	return;
}

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
unsigned int get_testsequence(char *filename, unsigned char **seq, SEQINFO *si, unsigned int *seq_total) {
	//int c;
	int seq_no = 0;
	unsigned int array_size;
	FILE *fp = NULL;
	fp = fopen(filename, "r");

	// try to determine the size from the file, and use it as the initial array size
	int fd = fileno(fp);
    struct stat fileinfo;
    fstat(fd, &fileinfo);
    array_size = fileinfo.st_size + 256;	// allocate 256 more bytes, as we may exceed the borders when comparing 16 chars at once

    // initial size of the array
    (*seq) = (unsigned char*) malloc(array_size);

    int inval_chars = 0;
    char* buffer = (char*) malloc(array_size);
    unsigned int read_sum = 0;
    unsigned int read = 1;
    while(read > 0) {
        read = fread(buffer+read_sum, 1, array_size, fp);
        read_sum += read;
    }
    fclose(fp);

    unsigned int i = 0;
    unsigned int pos_nl;

    if(buffer[0] == '>') { // read first sequence name, if beginning with >
		while(i < read_sum) {
			if (buffer[i] == '\n') {
				i++; inval_chars++;
				break;
			} else if (buffer[i] != '>') {
				si[seq_no].name[si[seq_no].nl] = buffer[i];
				si[seq_no].nl += 1;
			}
			i++; inval_chars++;
		}
    }

    while (i < read_sum) {
		__asm__ (
				"movdqu %2,%%xmm0;" 				// move unaligned 16 byte from the read-buffer into xmm0
				"mov $0x0a, %%ecx;"					// move '\n' (0x0a) into xmm1
				"pinsrb $0x00,%%ecx,%%xmm1;"
				"xor %%ecx,%%ecx;"					// clear ECX
				"pcmpistri $0x00, %%xmm0, %%xmm1;"	// find first occurence of '\n' in the buffer at position i
				"cmp $16, %%ecx;"
				"movdqu %%xmm0, %3;"
				"jne nocpyt;"
				"add $16, %1;"
	"nocpyt: 	 mov %%ecx,%0;"
				:"=r"(pos_nl),					// result %0: index of
				 "=m"(i)
				:"m"(buffer[i]), 				// operand %1: address of the buffer at position i
				 "m"((*seq)[i-inval_chars])		// operand %2: destination buffer (i minus number of invalid chars so far)
				: "%xmm0", "%xmm1", "%ecx"
		);

		if (pos_nl != 16) {
			if (buffer[i+pos_nl+1] != '>') {
				i += pos_nl + 1; inval_chars++;
			} else { // next sequence starting next line
				i += pos_nl + 1; inval_chars++;
				si[seq_no].max = i-inval_chars;
				++seq_no;
				si[seq_no].min = (i-inval_chars);
				si[seq_no].nl = 0;

				while(i < read_sum) {
							if (buffer[i] == '\n') {
								i++; inval_chars++;
								break;
							} else if (buffer[i] != '>') {
								si[seq_no].name[si[seq_no].nl] = buffer[i];
								si[seq_no].nl += 1;
							}
							i++; inval_chars++;
						}
			}
		}
    }

	si[seq_no].max = read_sum-inval_chars;
	if (seq_total != NULL) {
		(*seq_total) += (seq_no+1);
	}
	return seq_no+1;
}
