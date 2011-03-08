#include <stdio.h>
#include "LoopScoring.h"
#include "TreeScoring.h"

char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
                      'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

int *decomposeToNucleotideBits(int n, int* array){
   //ZS: decomposes an ambiguous base into nucleotides. 
   //"array" is an array with 4 positions
   //0th position is 1 if A is included, 0 if not
   //1st position is 1 if C is included, 0 if not
   //eg. [1100] represents M = "A or C"
   //    [1111] represents N = "A or C or G or U"
   int A = 0x1;
   int C = 0x2;
   int G = 0x4;
   int U = 0x8;
   if (n & A){array[0]=1;}
   else      {array[0]=0;}
   if (n & C){array[1]=1;}
   else      {array[1]=0;}
   if (n & G){array[2]=1;}
   else      {array[2]=0;}
   if (n & U){array[3]=1;}
   else      {array[3]=0;}
   return array;
} 

int maskToGTfoldBases(int value){
    //ZS - this function uses bitmasking to efficiently convert "scoring" bases 
    //to "GTfold" bases that can be used in "fourBaseIndex"
    //GTfold bases: A=0, C=1, G=2, U=3
    //Scoring bases: A=1, C=2, G=4, U=8
    if(value==1){
       return 0;}
    else if(value==2){
       return 1;}
    else if(value==4){
       return 2;}
    else if(value==8){
       return 3;}
    else{
         printf("ERROR! Base cannot be masked to GTfold type: %i", value);
         return -1;
    }
}

int sum(int* array){
    //finds the sum of elements in an array of size 4.
    int result = 0;
    int i;
    for(i = 0; i<4; i++){
          if(array[i]==1){
             result++;
          }
    }
    return result; 
}

int scoreStackPossibilities(int* extBase, int* intBase, nndb_constants* param){
    //TODO!!
    //ZS: returns the score for the stack. 
    //extBase is a 2D array containing the arrays for the external basepair. 
    //intBase is a 2D array containing the arrays for the internal basepair. 
    int score = 0;
    return score;     
}


int eH(BasePair* pair, UnpairedRegion* unpaired) {
	return 1;
}

int _eS(BasePair* exteriorPair, BasePair* interiorPair, nndb_constants* param){
    //ZS: Score all combinations of possible pairs represented by the 
    //ambiguous base, and find the average of those scores. 
    
    //Create arrays for the "exterior" and "interior" pairs
    int *extBase[2][4];
    int *extHighBaseArray = extBase[0];
    extHighBaseArray = decomposeToNucleotideBits(exteriorPair->highBase.base,extHighBaseArray);
    int *extLowBaseArray = extBase[1];
    extLowBaseArray = decomposeToNucleotideBits(exteriorPair->lowBase.base,extLowBaseArray);
    int *intBase[2][4];
    int *intHighBaseArray = intBase[0];
    intHighBaseArray = decomposeToNucleotideBits(interiorPair->highBase.base,intHighBaseArray);
    int *intLowBaseArray = intBase[1];
    intLowBaseArray = decomposeToNucleotideBits(interiorPair->lowBase.base,intLowBaseArray);

    //printf("%d - A: %d C: %d G: %d U: %d sum: %d \n",exteriorPair->highBase.base,
    //       *extBase[0], *(extBase[0]+1), *(extBase[0]+2), 
    //       *(extBase[0]+3), sum(extBase[0]));        
    
    //Find a cumulative score for all the possibilities
    int score = scoreStackPossibilities(*extBase, *intBase, param);
    
    //Divide by the number of combinations we have gone through to get the average.
    return score;
}

int isAmbiguous(int number){
    return !(number==1||number==2||number==4||number==8);
}

int eS(BasePair* exteriorPair, BasePair* interiorPair, nndb_constants* param){
    //ZS: Only scores non-ambiguous nucleotides 
    //Ignore if ambiguous base 
    if(isAmbiguous(exteriorPair->highBase.base) || isAmbiguous(exteriorPair->lowBase.base)
             || isAmbiguous(interiorPair->highBase.base) || isAmbiguous(interiorPair->lowBase.base)){
        printf("Warning: ignoring ambiguous base in stack score. (Stacking pairs: %c-%c, %c-%c)\n",
        bases[exteriorPair->highBase.base], bases[exteriorPair->lowBase.base],
        bases[interiorPair->highBase.base], bases[interiorPair->lowBase.base]); 
        return 0;
    }
    
    int score = param->stack[fourBaseIndex(maskToGTfoldBases(exteriorPair->lowBase.base), 
              maskToGTfoldBases(exteriorPair->highBase.base),
              maskToGTfoldBases(interiorPair->lowBase.base), 
              maskToGTfoldBases(interiorPair->highBase.base))];
    printf("Score for stacking pairs %c-%c, %c-%c: %d\n",
        bases[exteriorPair->lowBase.base], bases[exteriorPair->highBase.base],
        bases[interiorPair->lowBase.base], bases[interiorPair->highBase.base], score);

    return score;
}


int eL(BasePair* exteriorPair, BasePair* interiorPair,
                 UnpairedRegion* firstUnpaired, UnpairedRegion* secondUnpaired){
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
