#include "TreeScoring.h"
#include "StructureReader.h"
#include "LoopScoring.h"
//#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//Shel: Function for scoring a node (recursive)
int ScoreNode(TreeNode* node, nndb_constants* param){
	int result;
	result = 0;
	int *pairedChildren;
	pairedChildren = NULL;
	int numPairedChildren;
	numPairedChildren = 0;
	int i;
	
	for (i = 0 ; i < node->numChildren ; i++){ 
       // find location and number of paired children 
       //and add scores of associated loops
		if ((node->children[i])->isPair) {
			result += ScoreNode(node->children[i], param);
			numPairedChildren += 1;
			pairedChildren = realloc(pairedChildren, sizeof(int) * numPairedChildren);
			pairedChildren[numPairedChildren - 1] = i;
		}
	}
	
	if (node->lowBase.index != 0) 
	{
		
		if (numPairedChildren == 0)  // must be a hairpin
		{
			
            result += eH(NULL,NULL);
			//printf("Found a Hairpin Loop\n");
		}
		else if (numPairedChildren == 1)  // must be stack, bulge, or internal
		{
			if (node->numChildren == 1)  // must be stack 
			{
            
            //printf("Found a Stacked Pair\n");

				BasePair *first;
				BasePair *second;
				first = node;
				second = (node->children[0]);
				
             //result += eS(NULL,NULL);
            result += eS(first,second, param);
				
			}
			else 
			{  // must be bulge or internal 
				result += eL(NULL,NULL,NULL,NULL);
			 
				//printf("Found a Bulge or Inernal Loop\n");
			}
		}
		else  // must be a multi-loop
		{	
			result += eM(NULL,NULL);
			//printf("Found a Mulit-Loop\n");
		}
	}
	else { // must be external
		result += eE(NULL,NULL);
		//printf("Found a External Loop\n");
	}

	return result;
}


