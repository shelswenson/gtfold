/*
GTfold: compute minimum free energy of RNA secondary structure
Copyright (C) 2008 David A. Bader
http://www.cc.gatech.edu/~bader
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/ 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "utils.h"
#include "energy.h"
#include "global.h"
#include "algorithms.h"
#include "constraints.h"
#include "shapereader.h"
#ifdef _OPENMP 
#include "omp.h"
#endif

//#define DEBUG 1

void initializeMatrix(int len) {
	int i, j;

	for (i = 1; i <= len; ++i) 
		for (j = len; j >= i; --j) 
			if (canPair(RNA[i],RNA[j]) && j-i > TURN) 
				PP[i][j]  = 1;
}

void prefilter(int len, int prefilter1, int prefilter2) {
	char** in;
	int i, j, k, count;

	in = (char**)malloc(len*sizeof(char*));
	for (i = 1; i <= len; ++i) in[i - 1] = (char*)malloc(len*sizeof(char));

	for (i = 1; i <= len - prefilter2 + 1; ++i)
		for (j = len; j >= prefilter2 && j >= i; --j) {
			count = 0;
			for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
				if (PP[i + k][j - k] == 1) ++count;
			if (count >= prefilter1)
				for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
					++in[i + k - 1][j - k - 1];
		}

	for (i = 1; i <= len; ++i) {
		for (j = len; j >= i; --j)
			if (!in[i - 1][j - 1]) PP[i][j] = 0;
		free(in[i - 1]);
	}

	free(in);
}

int calcVBI(int i, int j) {
	int p=0, q=0;
	int VBIij = INFINITY_;

	for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
		int minq = j-i+p-MAXLOOP-2;
		if (minq < p+1+TURN) minq = p+1+TURN;
		int maxq = (p==(i+1))?(j-2):(j-1);

		for (q = minq; q <= maxq; q++) {
			if (PP[p][q]==0) continue;
			if (!canILoop(i,j,p,q)) continue;
			VBIij = MIN(eL(i, j, p, q) + V(p,q), VBIij);
		}
	}

	return VBIij;
}

int calcVBI1(int i, int j) {
	int p=0, q=0;
	int VBIij = INFINITY_;

	for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
		int minq = j-i+p-MAXLOOP-2;
		if (minq < p+1+TURN) minq = p+1+TURN;
		int maxq = (p==(i+1))?(j-2):(j-1);

		for (q = minq; q <= maxq; q++) {
			if (PP[p][q]==0) continue;
			if (!canILoop(i,j,p,q)) continue;
			VBIij = MIN(eL1(i, j, p, q) + V(p,q), VBIij);
		}
	}

	return VBIij;
}

int calcVBI2(int i, int j, int  len) {
	int d, ii, jj; 
	int energy = INFINITY_;

	for (d = j-i-3; d >= TURN+1 && d >= j-i-2-MAXLOOP; --d) 
		for (ii = i + 1; ii < j - d && ii <= len; ++ii)
		{    
			jj = d + ii;
			if (PP[ii][jj]==1)
				energy = MIN(energy, eL1(i, j, ii, jj) + V(ii, jj));
		}    

	return energy;
}

int calculate(int len) { //, int nThreads, int unamode, int mismatch) { 
	int b, i, j;
#ifdef _OPENMP
	if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
#endif
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
	fprintf(stdout,"Thread count: %3d \n",omp_get_num_threads());
#endif

	initializeMatrix(len);
	if (g_unamode || g_prefilter_mode) {
		prefilter(len,g_prefilter1,g_prefilter2);
	}

	for (b = TURN+1; b <= len-1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
		for (i = 1; i <= len - b; i++) {
			j = i + b;
			int newWM = INFINITY_; 
			
			if (PP[i][j]==1) {
				int eh = canHairpin(i,j)?eH(i,j):INFINITY_; //hair pin
				int es = canStack(i,j)?eS(i,j)+V(i+1,j-1):INFINITY_; // stack

				// Internal Loop BEGIN
				if (g_unamode) 
					VBI(i,j) = calcVBI1(i,j);
				else
					VBI(i,j) = calcVBI(i,j);
				// Internal Loop END

				// Multi Loop BEGIN
				int VMij =  WMPrime[i+1][j-1]; 
				int VMidj = WMPrime[i+2][j-1]; 
				int VMijd = WMPrime[i+1][j-2]; 
				int VMidjd = WMPrime[i+2][j-2]; 

				int d3 = canSS(j-1)?Ed3(i,j,j-1):INFINITY_;
				int d5 = canSS(i+1)?Ed5(i,j,i+1):INFINITY_;
				VMij = MIN(VMij, (VMidj + d5 +Ec)) ;
				VMij = MIN(VMij, (VMijd + d3 +Ec));

				if (g_unamode || g_mismatch) {
					VMij = MIN(VMij, (VMidjd + Estackm(i,j) + 2*Ec));
				} else {
					VMij = MIN(VMij, (VMidjd + d5 + d3+ 2*Ec));
				}

				VMij = VMij + Ea + Eb + auPenalty(i,j);
				VM(i,j) = canStack(i,j)?VMij:INFINITY_;
				// Multi Loop END

				V(i,j) = MIN4(eh,es,VBI(i,j),VM(i,j));
			}
			else V(i,j) = INFINITY_;

			int h; 
			for (h = i+TURN+1 ; h <= j-TURN-2; h++) {
				// Added auxillary storage WMPrime to speedup multiloop calculations
				WMPrime[i][j] = MIN(WMPrime[i][j], WMU(i,h-1) + WML(h,j)); 
				//newWM = (!forcePair(i,j))?MIN(newWM, WMU(i,h-1) + WML(h,j)):newWM;
			}
			
			//ZS: This sum corresponds to when i,j are NOT paired with each other.
			//So we need to make sure only terms where i,j aren't pairing are considered. 
			newWM = (!forcePair(i,j))?MIN(newWM, WMPrime[i][j]):newWM;
			
			newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM); 
			newWM = canSS(i)?MIN(V(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM):newWM; //i dangle
			newWM = canSS(j)?MIN(V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM):newWM;  //j dangle

			if (g_unamode || g_mismatch) {
				if (i<j-TURN-2)
					newWM = (canSS(i)&&canSS(j))?MIN(V(i+1,j-1) + Estackm(j-1,i+1) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM):newWM; 
			}
			else {
				newWM = (canSS(i)&&canSS(j))?MIN(V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM):newWM; //i,j dangle
			}

			newWM = canSS(i)?MIN(WMU(i+1,j) + Ec, newWM):newWM; //i dangle
			newWM = canSS(j)?MIN(WML(i,j-1) + Ec, newWM):newWM; //j dangle
			WMU(i,j) = WML(i,j) = newWM;
		}
	}

	for (j = TURN+2; j <= len; j++) {
		int i, Wj, Widjd, Wijd, Widj, Wij, Wim1;
		Wj = INFINITY_;
		for (i = 1; i < j-TURN; i++) {
			Wij = Widjd = Wijd = Widj = INFINITY_;
			Wim1 = MIN(0, W[i-1]); 
			Wij = V(i, j) + auPenalty(i, j) + Wim1;
		
			if (g_unamode || g_mismatch) {
				Widjd = (canSS(i)&&canSS(j))?V(i+1,j-1) + auPenalty(i+1,j-1) + Estacke(j-1,i+1) + Wim1:Widjd;
			} else {
				Widjd = (canSS(i)&&canSS(j))?V(i+1,j-1) + auPenalty(i+1,j-1) + Ed3(j-1,i + 1,i) + Ed5(j-1,i+1,j) + Wim1:Widjd;
			}
			
			Wijd = canSS(j)?V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + Wim1:Wijd;
			Widj = canSS(i)?V(i+1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1:Widj;
			Wj = MIN(MIN4(Wij, Widjd, Wijd, Widj), Wj); 
		}
		W[j] = canSS(j)?MIN(Wj, W[j-1]):Wj;
	}

#ifdef DEBUG
	FILE* file = fopen("VM.txt", "w");
	int ii, jj;
	for (ii = 1; ii <= len; ++ii) {    
		for (jj = 1; jj <= len; ++jj) {
			fprintf(file, "%d %d %d\n",ii,jj,VM(ii,jj));
		}
	}    
	fclose(file);
	
 	file = fopen("WM.txt", "w");
	for (ii = 1; ii <= len; ++ii) {    
		for (jj = 1; jj <= len; ++jj) {
			fprintf(file, "%d %d %d\n",ii,jj,WM(ii,jj));
		}
	}    
	fclose(file);
#endif
	
	return W[len];
}
