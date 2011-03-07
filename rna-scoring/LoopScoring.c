#include <stdio.h>
#include "LoopScoring.h"
#include "TreeScoring.h"

int eH(BasePair* pair, UnpairedRegion* unpaired) 
{
	return 1;
}

int eS(BasePair *exteriorPair, BasePair *interiorPair)
{
    char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    printf("Bases were %c-%c stacked on %c-%c\n", 
                  bases[exteriorPair->highBase.base], bases[exteriorPair->lowBase.base],
                  bases[interiorPair->highBase.base], bases[exteriorPair->lowBase.base]); 
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
