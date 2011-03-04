#include "RNAScoring.h"
#include "StructureReader.h"
#include "TreeScoring.h"
//#include <ctype.h>
//#include <limits.h>
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
		fprintf(stderr, "USAGE: RNAScoring <filename>\n");
		return 1;
    }
	
    TreeNode* tree = CreateFromFile(argv[1]);
    //PrintTree(tree, 0);
	int tree_score = ScoreTree(tree);
	printf("Tree score is ");
	printf("%d\n", tree_score);
    return 0;
}
