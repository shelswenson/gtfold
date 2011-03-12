#include "TreeScoring.h"
#include "StructureReader.h"
#include "LoopScoring.h"
//#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int result;
//Shel: Function for scoring a node (recursive)
int ScoreNode(TreeNode* node, int* RNA, nndb_constants* param){
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
			result += ScoreNode(node->children[i], RNA, param);
			numPairedChildren += 1;
			pairedChildren = realloc(pairedChildren, sizeof(int) * numPairedChildren);
			pairedChildren[numPairedChildren - 1] = i;
		}
	}
	
	if (node->lowBase.index != 0) 
	{
		
		if (numPairedChildren == 0)  // must be a hairpin
		{
			
          //result += eH(0,0, RNA, param);
          //result += 0;

			printf("Found a Hairpin Loop\n");
		}
		else if (numPairedChildren == 1)  // must be stack, bulge, or internal
		{
			if (node->numChildren == 1)  // must be stack 
			{
			 //	printf("Found a stacked pair:");
          	int energy = eS(node->lowBase.index, node->highBase.index, RNA, param);
	       	result += energy; 
				printf("Found a Stacked Pair with energy %i\n", energy);
			}
			else 
			{  // must be bulge or internal 
			//	result += eL(NULL,NULL,NULL,NULL);
			 
				//printf("Found a Bulge or Inernal Loop\n");
			}
		}
		else  // must be a multi-loop
		{	
		//	result += eM(NULL,NULL);
			//printf("Found a Mulit-Loop\n");
		}
	}
	else { // must be external
	//	result += eE(NULL,NULL);
		//printf("Found a External Loop\n");
	}

	return result;
}


