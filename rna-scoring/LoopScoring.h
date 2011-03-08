/*
 A class to score an RNA structure reading in from a file.
 */

#ifndef RNASCROING_H
#define RNASCORING_H
#include "TreeScoring.h"
#include "data.h"

#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

int eH(BasePair* pair, UnpairedRegion* unpaired);
int eS(BasePair* exteriorPair, BasePair* interiorPair, nndb_constants* param);
int eL(BasePair* exteriorPair, BasePair* interiorPair, 
    UnpairedRegion* firstUnpaired, UnpairedRegion* secondUnpaired);
int eM(BasePair** pairs, UnpairedRegion** unpairedRegions);
int eE(BasePair** pairs, UnpairedRegion** unpairedRegions);

#endif
