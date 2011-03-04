#include <stdio.h>
#include "LoopScoring.h"

int eH(BasePair* pair, UnpairedRegion* unpaired) 
{
	return 1;
}

int eS(BasePair* exteriorPair, BasePair* interiorPair)
{
	return 1;
}

int eL(BasePair* exteriorPair, BasePair* interiorPair, UnpairedRegion* firstUnpaired, UnpairedRegion* secondUnpaired)
{
	return 1;
}

int eM(BasePair** pairs, UnpairedRegion** unpairedRegions)
{
	return 1;
}

int eE(BasePair** pairs, UnpairedRegion** unpairedRegions)
{
	return 0;
}
