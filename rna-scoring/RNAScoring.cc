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

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
		fprintf(stderr, "USAGE: RNAScoring <filename>\n");
		return 1;
    }
	
    TreeNode* tree = CreateFromFile(argv[1]);
    PrintTree(tree, 0);
	nndb_constants* param = populate("./data", 0);
	
	int tree_score = ScoreNode(tree);
	printf("Tree score is ");
	printf("%d\n", tree_score);
	
    int user_id;
    printf ("Enter a number: %d\n", user_id);

    scanf ("%d", &user_id);

    return 0;
}
