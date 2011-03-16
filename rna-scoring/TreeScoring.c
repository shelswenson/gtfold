#include "TreeScoring.h"
#include "StructureReader.h"
#include "LoopScoring.h"
//#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Shel: Function for scoring a node (recursive)
int ScoreNode(TreeNode* node, int* RNA, nndb_constants* param){
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
			
			 int energy = eH(node->lowBase.index,node->highBase.index, RNA, param);
          result += energy;
			printf("%d \t %d: Hairpin Loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
		}
		else if (numPairedChildren == 1)  // must be stack, bulge, or internal
		{
			if (node->numChildren == 1)  // must be stack 
			{
          	int energy = eS(node->lowBase.index, node->highBase.index, RNA, param);
	       	result += energy; 
	     		printf("%d \t %d: Stacked pair with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
			}
			else 
			{  // must be bulge or internal 
			   
				int energy = eL(node->lowBase.index, node->highBase.index, 
				              node->children[pairedChildren[0]]->lowBase.index,
				              node->children[pairedChildren[0]]->highBase.index,
								  RNA, param);
            result += energy;
				printf("%d \t %d: Bulge or Inernal Loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
			}
		}
		else  // must be a multi-loop
		{	
			int nr_branches = node->numChildren; //ZS: I'm not sure about this????
			int nr_unpaired = node->numChildren - numPairedChildren;
			int energy = eM(nr_branches, nr_unpaired, param);
			result += energy;
			printf("%d \t %d: Multi-loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
		}
	}
	else { // must be external
      int energy = eE(RNA, param);
		result += energy; 
		printf("%d \t %d: External loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
	}

	return result;
}


