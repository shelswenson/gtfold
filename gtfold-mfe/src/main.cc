/*
 GTfold: compute minimum free energy of RNA secondary structure
 Copyright (C) 2008  David A. Bader
 http://www.cc.gatech.edu/~bader

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Authored by Amrita Mathuriya August 2007 - January 2009.*/
/* Modified by Sonny Hernandez May 2007 - Aug 2007. All comments added marked by "SH: "*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "loader.h"
#include "algorithms.h"
#include "traceback.h"
#include "main.h"
#include "main-c.h"

using namespace std;

/* GLOBAL VARIABLES */
enum BOOL ILSA; /* A boolean variable to know if we are executing with Internal loop speedup algorithm (ILA) or not. ILSA finds the optimal internal loop by exploring all possibilities. */
enum BOOL NOISOLATE;
enum BOOL USERDATA;
#ifdef DYNALLOC
int LENGTH;
unsigned char *RNA; /* Contains RNA string in terms of 0, 1, 2, 3 for A, C, G and U respectively*/
int *structure; /* An array to contain the optimal structure */
int *V; /* int V[LENGTH][LENGTH]; */
int *W;
int **VBI; /* VBI(i,j) will contain the energy of optimal internal loop closed with (i,j) base pair */
int **VM; /* VM(i, j) will contain the energy of optimla multiloop closed with (i,j) base pair */
int **WM; /* This array is introduced to help multiloop calculations. WM(i,j) contains the optimal energy of string segment from si to sj if this forms part of a multiloop */
int *indx; /* This array is used to index V array. Here V array is mapped from 2D to 1D and indx array is used to get the mapping back.*/
int *constraints;
#else
/* This are previously used variables, now they are not used. */
unsigned char RNA[LENGTH];
int structure[LENGTH];
int VBI[LENGTH][LENGTH];
int VM[LENGTH][LENGTH];
int V[(LENGTH-1)*(LENGTH)/2 + 1]; /* int V[LENGTH][LENGTH]; */
int WM[LENGTH][LENGTH];
int W[LENGTH];
int indx [LENGTH];
#endif

/* Function for displaying help */
void help() {
	fprintf(
			stderr,
			"Usage: gtfold [-ilsa] [-noisolate] [-constraints filename] [-datadir datadirloc] filename(sequence)\n\n-ilsa = Calculation of all possible internal loops using the speedup algorithm\n-noisolate=prevents isolated base pairs from forming\nSequence file has to be in one of the two formats: Single line or FASTA\nConstraint filename is a optional parameter.\n\nSyntax for giving constraints is:\n\t\tfor forcing (i,j) base pair, F i j and\n\t\tto prohibit (i,j) base pair, P i j.\n\n");
	exit(-1);
}

/* Function for calculating time */
double get_seconds() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

/*Function for printing the sequence*/
void printSequence(int len) {
	int i = 1;
	for (i = 1; i <= len; i++) {
		if (RNA[i] == 0)
			printf("A");
		else if (RNA[i] == 1)
			printf("C");
		else if (RNA[i] == 2)
			printf("G");
		else
			printf("T");
	}
	printf("\n");
}

/*Function for printing the input constraints*/
void printConstraints(int len) {
	int i = 1;
	for (i = 1; i <= len; i++) {
		if (constraints[i] > 0 && constraints[i] > i)
			printf("(");
		else if (constraints[i] > 0 && constraints[i] < i)
			printf(")");
		else if (constraints[i] < 0)
			printf("x");
		else
			printf(".");
	}
	printf("\n");
}

/*Function for printing the predicted structure*/
void printStructure(int len) {
	int i = 1;
	for (i = 1; i <= len; i++) {
		if (structure[i] > 0 && structure[i] > i)
			printf("(");
		else if (structure[i] > 0 && structure[i] < i)
			printf(")");
		else
			printf(".");
	}
	printf("\n");
}

/* Initialize global variables. */
void init_variables(int len) {

	int i;

#ifdef DYNALLOC
	LENGTH = len + 1;

	RNA = (unsigned char *) malloc(LENGTH * sizeof(unsigned char));
	if (RNA == NULL) {
		perror("Cannot allocate variable 'RNA'");
		exit(-1);
	}

	structure = (int *) malloc(LENGTH * sizeof(int));
	if (structure == NULL) {
		perror("Cannot allocate variable 'structure'");
		exit(-1);
	}

	V = (int *) malloc(((LENGTH - 1) * (LENGTH) / 2 + 1) * sizeof(int));
	if (V == NULL) {
		perror("Cannot allocate variable 'V'");
		exit(-1);
	}

	W = (int *) malloc(LENGTH * sizeof(int));
	if (W == NULL) {
		perror("Cannot allocate variable 'W'");
		exit(-1);
	}

	VBI = (int **) malloc(LENGTH * sizeof(int *));
	if (VBI == NULL) {
		perror("Cannot allocate variable 'VBI'");
		exit(-1);
	}
	for (i = 0; i < LENGTH; i++) {
		VBI[i] = (int *) malloc(LENGTH * sizeof(int));
		if (VBI[i] == NULL) {
			perror("Cannot allocate variable 'VBI[i]'");
			exit(-1);
		}
	}

	VM = (int **) malloc(LENGTH * sizeof(int *));
	if (VM == NULL) {
		perror("Cannot allocate variable 'VM'");
		exit(-1);
	}
	for (i = 0; i < LENGTH; i++) {
		VM[i] = (int *) malloc(LENGTH * sizeof(int));
		if (VM[i] == NULL) {
			perror("Cannot allocate variable 'VM[i]'");
			exit(-1);
		}
	}

	WM = (int **) malloc(LENGTH * sizeof(int *));
	if (WM == NULL) {
		perror("Cannot allocate variable 'WM'");
		exit(-1);
	}
	for (i = 0; i < LENGTH; i++) {
		WM[i] = (int *) malloc(LENGTH * sizeof(int));
		if (WM[i] == NULL) {
			perror("Cannot allocate variable 'WM[i]'");
			exit(-1);
		}
	}

	indx = (int *) malloc(LENGTH * sizeof(int));
	if (indx == NULL) {
		perror("Cannot allocate variable 'indx'");
		exit(-1);
	}

	constraints = (int*) malloc((len + 1) * sizeof(int));
	if (constraints == NULL) {
		perror("Cannot allocate variable 'constraints'");
		exit(-1);
	}

#endif
	return;
}

/* deallocate global variables */
void free_variables() {
	int i;

#ifdef DYNALLOC
	free(indx);
	for (i = 0; i < LENGTH; i++)
		free(WM[i]);
	free(WM);
	for (i = 0; i < LENGTH; i++)
		free(VM[i]);
	free(VM);
	for (i = 0; i < LENGTH; i++)
		free(VBI[i]);
	free(VBI);
	free(W);
	free(V);
	free(constraints);
	free(structure);
	free(RNA);

#endif

	return;

}

/* main function - This calls
 *  1) Read command line arguments.
 *  2) populate() from loader.cc to read the thermodynamic parameters defined in the files given in data directory.
 *  3) Initialize variables
 *  4) Calls calculate function defined in algorithms.c for filling up energy tables.
 *  5) Calls trace function defined in trace.c file for tracing the optimal secondary structure
 *  6) Then it generates .ct file from the 1D array structure.
 *  */
int main(int argc, char** argv) {
	int i;
	ifstream cf;
	int bases;
	string s, seq;
	int energy;
	double t1;
	ILSA = FALSE;
	NOISOLATE = FALSE;
	fprintf(stdout,
			"GTfold: A Scalable Multicore Code for RNA Secondary Structure Prediction\n");
	fprintf(
			stdout,
			"(c) 2007-2010  D.A. Bader, S. Mallidi, A. Mathuriya, C.E. Heitsch, S.C. Harvey\n");
	fprintf(stdout, "Georgia Institute of Technology\n\n");

	/* Reading command line arguments */
	if (argc < 2)
		help();

	int fileIndex = 0, consIndex = 0, dataIndex = 0;
	i = 1;
	while (i < argc) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i], "-ilsa") == 0) {
				ILSA = TRUE;
			} else if (strcmp(argv[i], "-noisolate") == 0) {
				NOISOLATE = TRUE;
			} else if (strcmp(argv[i], "-help") == 0) {
				help();
			} else if (strcmp(argv[i], "-constraints") == 0) {
				if (i < argc)
					consIndex = ++i;
				else
					help();
			} else if (strcmp(argv[i], "-datadir") == 0) {
				USERDATA = TRUE;
				if (i < argc)
					dataIndex = ++i;
				else
					help();
			}
		} else {
			fileIndex = i;
		}
		i++;
	}

	/*
	 for ( int i = 1; i < argc; i++ ) {
	 if ( argv[i][0] == '-' ) {
	 if ( strcmp(argv[i], "-ilsa")==0 ) {
	 ILSA = TRUE;
	 }
	 else if ( strcmp(argv[i], "-noisolate")==0 ){
	 NOISOLATE = TRUE;
	 }
	 else if ( strcmp(argv[i], "-help")==0 )  help();
	 else help();
	 } else {
	 if(fileIndex ==0) fileIndex = i;
	 else conIndex=i;
	 }
	 }
	 */

	if (fileIndex == 0)
		help();
	if (ILSA == TRUE)
		fprintf(stdout, "Running with Internal Loop Speedup Algorithm\n");
	if (NOISOLATE == TRUE)
		fprintf(stdout, "Not allowing isolated base pairs\n");
	else
		fprintf(stdout, "Allowing isolated base pairs\n");
	if (consIndex != 0)
		fprintf(stdout, "Constraint file index: %d\n", consIndex);

	fprintf(stdout, "Opening file: %s\n", argv[fileIndex]);
	cf.open(argv[fileIndex], ios::in);
	if (cf != NULL)
		fprintf(stdout, "File opened.\n\n");
	else {
		fprintf(stdout, "File open failed.\n\n");
		exit(-1);
	}

	seq = "";
	s = "";

	//Handle FASTA input
	char ss[10000];
	cf.getline(ss, 10000);
	if (ss[0] != '>') {
		char *fline;
		fline = strtok(ss, " ");
		while (fline != NULL) {
			seq.append(fline);
			fline = strtok(NULL, " ");
		}
	}

	while (!cf.eof()) {
		cf >> s;
		seq.append(s);
		s = "";
	}
	s = seq;

	bases = s.length();

	init_variables(bases);

	cout << "Sequence: " << s << endl;
	fprintf(stdout, "Sequence length: %5d\n\n", bases);

	cf.close();

	int fit = 0, pit = 0, it = 0;
	int **fbp = NULL, **pbp = NULL;
	int numfConstraints = 0, numpConstraints = 0;

	if (consIndex == 0)
		goto nocons;

	fprintf(stdout, "Running with constraints\n");
	fprintf(stdout, "Opening constraint file: %s\n", argv[consIndex]);
	cf.open(argv[consIndex], ios::in);
	if (cf != NULL)
		fprintf(stdout, "Constraint file opened.\n");
	else {
		fprintf(stderr, "Error opening constraint file\n\n");
		exit(-1);
	}

	char cons[10];

	while (!cf.eof()) {
		cf.getline(cons, 10);
		if (cons[0] == 'F' || cons[0] == 'f')
			numfConstraints++;
		if (cons[0] == 'P' || cons[0] == 'p')
			numpConstraints++;
	}

	cf.close();

	fprintf(stdout, "Number of Constraints given: %d\n\n", numfConstraints
			+ numpConstraints);
	if (numfConstraints + numpConstraints != 0)
		fprintf(stdout, "Reading Constraints.\n");
	else {
		fprintf(stderr, "No Constraints found.\n\n");
		goto nocons;
	}

	//int fbp[numfConstraints][2], pbp[numpConstraints][2];
	fbp = (int**) malloc(numfConstraints * sizeof(int*));
	pbp = (int**) malloc(numpConstraints * sizeof(int*));

	for (it = 0; it < numfConstraints; it++) {
		fbp[it] = (int*) malloc(2* sizeof (int));
	}
	for(it=0; it<numpConstraints; it++) {
		pbp[it] = (int*)malloc(2*sizeof(int));
	}
	cf.open(argv[consIndex], ios::in);

	while(!cf.eof()) {
		cf.getline(cons,10);
		char *p=strtok(cons, " ");
		p = strtok(NULL, " ");
		if(cons[0]=='F' || cons[0]=='f') {
			int fit1=0;
			while(p!=NULL) {
				fbp[fit][fit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			fit++;
		}
		if( cons[0]=='P' || cons[0]=='p') {
			int pit1=0;
			while(p!=NULL) {
				pbp[pit][pit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			pit++;
		}
	}

	fprintf(stdout, "Forced base pairs: ");
	for(it=0; it<numfConstraints; it++) {
		fprintf(stdout, "(%d,%d) ", fbp[it][0], fbp[it][1]);
	}
	fprintf(stdout, "\nProhibited base pairs: ");
	for(it=0; it<numpConstraints; it++) {
		fprintf(stdout, "(%d,%d) ", pbp[it][0], pbp[it][1]);
	}
	fprintf(stdout, "\n\n");

	nocons:
	/* SH: Conversion of the sequence to numerical values. */
	for(i = 1; i <= bases; i++) {
		RNA[i] = getBase(s.substr(i-1,1));
		if (RNA[i]=='X') {
			fprintf(stderr,"ERROR: Base unrecognized\n");
			exit(0);
		}
	}

	if(USERDATA==TRUE)
		populate(argv[dataIndex]);
	else
		populate(NULL); /* Defined in loader.cc file to read in the thermodynamic parameter values from the tables in the ../data directory. */

	initTables(bases); /* Initialize global variables */

	fprintf(stdout,"Computing minimum free energy structure. . . \n");
	fflush(stdout);

	t1 = get_seconds();
	energy = calculate(bases, fbp, pbp, numfConstraints, numpConstraints); /* Runs the Dynamic programming algorithm to calculate the optimal energy. Defined in algorithms.c file.*/
	t1 = get_seconds() - t1;

	fprintf(stdout," Done.\n");

	fprintf(stdout,"Minimum Free Energy = %12.2f\n\n", energy/100.00);
	fprintf(stdout,"MFE running time (in seconds): %9.6f\n\n", t1);

	t1 = get_seconds();
	trace(bases); /* Traces the optimal structure*/
	t1 = get_seconds() - t1;

	std::stringstream ss1, ss2;
	char suboptfile[1024];
	ss1 << bases;
	ss2 << energy/100.0;

	i = 0;

	strcpy (suboptfile, argv[fileIndex]);
	strcat ( suboptfile, ".ct");
	ofstream outfile;
	outfile.open ( suboptfile );

	fprintf(stdout, "Writing secondary structure to the file: %s\n", suboptfile);

#if 0
	outfile << bases << " " << energy/100.0;
	outfile << endl << s;

	for ( i = 1; i <= bases; i++ )
		outfile << "\n" <<  i << " " << structure[i] ;
#endif
	/* Generate the output file containing the optimal secondary structure in .ct format */
#if 1
	outfile << bases << "\t  dG = " << energy/100.0;
	i = 1;
	while ( i <= bases ) {
		outfile << endl << i << "\t" << s[i-1] << "\t" << i-1 << "\t" << (i+1)%(bases+1) << "\t" << structure[i] << "\t" << i;
		i++;
	}
	outfile << endl;
#endif


	outfile.close();

	fprintf(stdout,"\n\n");
	fprintf(stdout,"Traceback running time (in seconds): %9.6f\n", t1);

	fprintf(stdout, "\n\nFolding complete\n\n");
	printSequence(bases);
	printConstraints(bases);
	printStructure(bases);

	free_variables();

	return 0;

}