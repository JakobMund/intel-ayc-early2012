//============================================================================
// Name        : AYC-sample-gen.cpp
// Author      : Jakob Mund
// Version     : 1.0
// Description : Generate sample data for the AYC using a RNG
//============================================================================

#define MAX_FILENAME_LENGTH 50
#define DEFAULT_REF_FILENAME "refseq.txt"
#define DEFAULT_TAR_FILENAME "input.txt"
#define MAX_LINE_SIZE 70

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <fstream>

using namespace std;

void inline arg_to_size (char* str, int *size) {
	*size = atoi(str);
	if(strchr(str,'K') != NULL || strchr(str,'k') != NULL) {
		*size *= 1000;	// 10^3
	}
	if(strchr(str,'M') != NULL || strchr(str,'m') != NULL) {
		*size *= 1000000; // 10^6
	}
	if(strchr(str,'G') != NULL || strchr(str,'g') != NULL) {
		*size *= 1000000000; // 10^9
	}
}

int main(int argc, char** argv) {

	int refsize, numtargets, tarsize;
	char reffilename[MAX_FILENAME_LENGTH], tarfilename[MAX_FILENAME_LENGTH];
	char *refarray, *tararray;
	ofstream reffile, tarfile;

	unsigned int iseed;
	short rng_char;

	// correct number of arguments is three (resp. four including bin name)
	if (argc < 4 || argc > 6) {
		cout << "USAGE: ayc-sample-gen <length of reference seq> <number of target seq> <length of target seq> [<filename ref.seq.>] [<filename tar.seq.>]\n";
	} else {
		// extract CLI arguments
		arg_to_size(argv[1], &refsize);
		arg_to_size(argv[2], &numtargets);
		arg_to_size(argv[3], &tarsize);

		// set filenames (from arguments, or use defaults)
		if (argc > 4)	strcpy(reffilename, argv[3]);
		else			strcpy(reffilename, DEFAULT_REF_FILENAME);
		if (argc > 5)	strcpy(tarfilename, argv[4]);
		else			strcpy(tarfilename, DEFAULT_TAR_FILENAME);

		// alloc mem for output
		refarray = (char*) malloc(refsize*sizeof(char));
		tararray = (char*) malloc(numtargets*tarsize*sizeof(char));

		// generate random sequences...
		iseed = (unsigned int) time (NULL);
		srand (iseed);
		// ... for the reference sequence
		for(int i = 0; i < refsize; i++) {
			rng_char = rand () % 4;
			switch (rng_char) {
			case 0:  *(refarray+i) = 'C'; break;
			case 1:  *(refarray+i) = 'A'; break;
			case 2:  *(refarray+i) = 'G'; break;
			default: *(refarray+i) = 'T'; break;
			}
		}
		// ... and for the test sequences
		for(int j = 0; j < numtargets; j++) {
			for(int i = 0; i < tarsize; i++) {
				rng_char = rand () % 4;
				switch (rng_char) {
				case 0:  *(tararray+i+(j*tarsize)) = 'C'; break;
				case 1:  *(tararray+i+(j*tarsize)) = 'A'; break;
				case 2:  *(tararray+i+(j*tarsize)) = 'G'; break;
				default: *(tararray+i+(j*tarsize)) = 'T'; break;
				}
			}
		}

		// write the reference sequence ...
		reffile.open(reffilename);
		reffile << "> Reference sequence (auto-generated, " << refsize << " characters)\n";
		for(int i = 0; i < refsize; i++) {
			reffile.put(*(refarray+i));
			if ((i+1) % MAX_LINE_SIZE == 0) {		// new line each MAX_LINE_SIZE characters
				reffile.put('\n');
			}
		}
		reffile.close();
		// ... and the test sequences to files
		tarfile.open(tarfilename);
		tarfile << "> Test sequences (auto-generated, " << numtargets << " sequences, " << tarsize << " characters each)\n";
		for(int j = 0; j < numtargets; j++) {
			tarfile << "> Test_sequence_" << j << "\n";
			for(int i = 0; i < tarsize; i++) {
				tarfile.put(*(tararray+i+(j*tarsize)));
				if ((i+1) % MAX_LINE_SIZE == 0) {		// new line each MAX_LINE_SIZE characters
					tarfile.put('\n');
				}
			}
			tarfile.put('\n');
		}
		tarfile.close();

		cout << "Sample written to \n   >" << reffilename << "< [reference sequence, size=" <<  argv[1] << "] and\n   >" << tarfilename << "< [test sequences, size=" << numtargets << "x" << argv[3] << "]\n";
	}

	return 0;
}
