/*
	A class to score a RNA sructure stored in a tree (node) data structure.
*/

#ifndef TREESCORING_H
#define TREESCORING_H
#include "StructureReader.h" 

/* Nucleotides in a base pair */
typedef	struct _BasePair
{
    BaseData lowBase; /* The low base. */
    BaseData highBase; /* The high base, if a pair. */
} BasePair;

/* Data for unparied region */
typedef	struct _UnpairedRegion
{
    BaseData firstBase; /* The first base in the unpaired region, if region is non-empty. */
    BaseData lastBase; /* The high base in the unpaired rebion, if region has length at least two. */
	int	numBases; /* The number of bases in the unpaired region (possibly zero).*/
	//struct _Base**; /* The actual bases in the unpaired region listed in 5' to 3' order.*/
} UnpairedRegion;

int ScoreTree(TreeNode* root);

#endif