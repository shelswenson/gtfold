/*
 A class to score an RNA structure reading in from a file.
 */

#ifndef RNASCROING_H
#define RNASCORING_H
#include "TreeScoring.h"

int eH(BasePair* pair, UnpairedRegion* unpaired);
int eS(BasePair* exteriorPair, BasePair* interiorPair);
int eL(BasePair* exteriorPair, BasePair* interiorPair, UnpairedRegion* firstUnpaired, UnpairedRegion* secondUnpaired);
int eM(BasePair** pairs, UnpairedRegion** unpairedRegions);
int eE(BasePair** pairs, UnpairedRegion** unpairedRegions);

#endif