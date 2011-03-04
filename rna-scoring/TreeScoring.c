#include "TreeScoring.h"
#include "StructureReader.h"
#include "LoopScoring.h"
//#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//Shel: Function for scoring a node (recursive)
int ScoreNode(TreeNode* node)
{
	char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};//remove
	int result;
	result = 0;
	int *pairedChildren;
	pairedChildren = NULL;
	int numPairedChildren;
	numPairedChildren = 0;
	int i;
	
	for (i = 0 ; i < node->numChildren ; i++) // find location and number of paired children and add scores of associated loops
	{
		if ((node->children[i])->isPair) 
		{
			result += ScoreNode(node->children[i]);
			numPairedChildren += 1;
			pairedChildren = realloc(pairedChildren, sizeof(int) * numPairedChildren);
			pairedChildren[numPairedChildren - 1] = i;
			printf("Child %d of (%d %c - %c %d) is paired\n", pairedChildren[numPairedChildren-1],//remove 
				   node->lowBase.index, bases[node->lowBase.base],
				   bases[node->highBase.base], node->highBase.index);
		}
	}
	
	if (numPairedChildren == 0)  // must be a hairpin
	{
		result += eH(NULL,NULL);
		printf("Found a Hairpin Loop\n");
	}
	else if (numPairedChildren == 1)  // must be stack, bulge, or internal
	{
		if (node->numChildren == 1)  // must be stack 
		{
			result += eS(NULL,NULL);
			printf("Found a Stacked Pair\n");
		}
		else 
		{  // must be bulge or internal 
			result += eL(NULL,NULL,NULL,NULL);
			 
			printf("Found a Bulge or Inernal Loop\n");
		}
	}
	else  // must be a multi-loop
	{	
		result += eM(NULL,NULL);
		printf("Found a Mulit-Loop\n");
	}
	
	return result;
}



//Shel: Function for scoring a tree from the root
int ScoreTree(TreeNode* root)
{
	int result;
	result = 0;
	int *pairedChildren;
	pairedChildren = NULL;
	int numPairedChildren = 0;
	int i;
	for (i = 0 ; i < root->numChildren ; i++)
	{
		if ((root->children[i])->isPair) 
		{
			result += ScoreNode(root->children[i]);
			numPairedChildren += 1;
			pairedChildren = realloc(pairedChildren, sizeof(int) * numPairedChildren);
			pairedChildren[numPairedChildren - 1] = i;
		}
	}
	printf("Root has %d paired children", numPairedChildren);
	result += eE(NULL,NULL);
	
	return result;
}

