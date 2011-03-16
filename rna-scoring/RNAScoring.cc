extern "C" {
    #include "StructureReader.h"
    #include "RNAScoring.h"
    #include "TreeScoring.h"    
}
//#include <ctype.h>
//#include <limits.h>
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include "loader.h"
//#include "data.h"
//#include "constants.h"
#include <cstdlib>
#include <iostream>
#include "LoopScoring.h"
#include <math.h>
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
		fprintf(stderr, "USAGE: RNAScoring <filename>\n");
		return 1;
    }
	
	printf("hi");
	char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
                      'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    ResultBundle* resultBundle = CreateFromFile(argv[1]);
    TreeNode* tree = resultBundle->treenode;
    int* RNA = resultBundle->RNA_seq;
    
   // PrintTree(tree, 0);
	 nndb_constants* param = populate("data/Turner99", 1);
	
	int one; int two; int three; int four; int score; 
	
	int tree_score = ScoreNode(tree, resultBundle->RNA_seq, param);
	printf("Tree score is ");
	printf("%.2f\n", (double)tree_score/100);

    return 0;
}
