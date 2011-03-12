/*
 A class to score an RNA structure reading in from a file.
 */

#ifndef RNASCROING_H
#define RNASCORING_H
#include "TreeScoring.h"
#include "data.h"

#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))
#define auPen(i, j) ((( (i)==BASE_U || (j)==BASE_U ) && ( (i)==BASE_A || (i)==BASE_G || (j)==BASE_A || (j)==BASE_G )) ? 1 : 0)
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

int eH(int i, int j, int* RNA, nndb_constants* param);
int eS(int i, int j, int* RNA, nndb_constants* param);
int eL(int i, int j, int ip, int jp, int* RNA, nndb_constants* param);
int eM(int nr_branches, int nr_unpaired, nndb_constants* param);
int eE(int* RNA, nndb_constants* param);

#endif
