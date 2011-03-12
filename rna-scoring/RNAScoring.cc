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
	
	 char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
                      'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    ResultBundle* resultBundle = CreateFromFile(argv[1]);
    TreeNode* tree = resultBundle->treenode;
    int* RNA = resultBundle->RNA_seq;
    
    //print out the whole RNA array. 
	int cnt;
	char basecode[4] = {'A','C','G','U'};
	int value = resultBundle->length;
	printf("Size %d\n", value);
	for(cnt = 1; cnt <= value; cnt++){
       printf("%d: %d\n", cnt, RNA[cnt]);
       //printf("%d: %d\n", cnt, *(resultBundle->RNA_seq+cnt));
   }
	printf("FINISHED READING\n");
   // PrintTree(tree, 0);
	 nndb_constants* param = populate("data/Turner99", 1);
	
	int one; int two; int three; int four; int score; 


//Output all stack values below: 
//	for(one = 0; one<4; one++){
//     for(two = 0; two<4; two++){
//        for(three=0; three<4; three++){
//           for(four=0; four<4; four++){
//              score = param->stack[fourBaseIndex(one, two, three, four)];
//               printf("%c-%c %c-%c: %d \n", bases[(int)pow(2,one)], bases[(int)pow(2,two)], 
//               bases[(int)pow(2,three)], bases[(int)pow(2,four)], score);       
//           }
//        }
//     }
//   }
	
//	int tree_score = ScoreNode(tree, resultBundle->RNA_seq, param);
//	printf("Tree score is ");
//	printf("%d\n", tree_score);

    return 0;
}
