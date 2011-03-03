#include "RNAScoring.h"
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Read data on a single base, filling baseData. Returns false if nothing could be read. Puts any pair
** index in pairIndex, or -1 if no pair.
*/
unsigned char ReadBase(FILE* filePtr, int index, BaseData* baseData, int* pairIndex, unsigned char isBPSEQ)
{
    while (1) /* Read as much as we have to to get the next base */
    {
	int junk = index - 1;
	int nextChar = 0;

	/* Check for a number. If not, ignore the line, or maybe the file is done. */
	int numRead = fscanf(filePtr, "%d", &junk);
	if (numRead != 1)
	{
	    /* Ignore this line, if it exists */
	    nextChar = fgetc(filePtr);
	    while (1)
	    {
		if (nextChar == EOF || nextChar == '\n')
		    break;
	    	nextChar = fgetc(filePtr);
	    }

	    /* If we are out of input, just return NULL */
	    if (nextChar == EOF)
	    {
		return 0;
	    }

	    /* Try the next line. */
	    continue;
	}

	/* Check for the next ID. If not, ignore the line and try the next one. */
	if (junk != index)
	{
	    while (1)
	    {
		nextChar = fgetc(filePtr);
		if (nextChar == EOF || nextChar == '\n')
		    break;
	    }
	    continue;
	}

	/* We're here if we have a valid index. */
	baseData->index = index;

	nextChar = fgetc(filePtr);
	while (isspace(nextChar))
	    nextChar = fgetc(filePtr);

	switch (nextChar)
	{
	    case 'a':
	    case 'A':
		baseData->base = A;
		break;

	    case 'c':
	    case 'C':
		baseData->base = C;
		break;

	    case 'g':
	    case 'G':
		baseData->base = G;
		break;

	    case 't':
	    case 'T':
	    case 'u':
	    case 'U':
		baseData->base = U;
		break;


	    case 'm':
	    case 'M':
		baseData->base = M;
		break;
		
        case 'r':
	    case 'R':
		baseData->base = R;
		break;
		
        case 's':
	    case 'S':
		baseData->base = S;
		break;
		
		case 'v':
	    case 'V':
		baseData->base = V;
		break;
		
        case 'w':
	    case 'W':
		baseData->base = W;
		break;
		
        case 'y':
	    case 'Y':
		baseData->base = Y;
		break;
		
	    case 'h':
	    case 'H':
		baseData->base = H;
		break;
		
        case 'k':
	    case 'K':
		baseData->base = K;
		break;
		
        case 'd':
	    case 'D':
		baseData->base = D;
		break;
		
        case 'b':
	    case 'B':
		baseData->base = B;
		break;
		
        case 'n':
	    case 'N':
		baseData->base = N;
		break;


	    default:
		fprintf(stderr, "Bad base: index %d\n", index);
		return 0;
	}

	nextChar = fgetc(filePtr);
	while (isspace(nextChar))
	    nextChar = fgetc(filePtr);
	ungetc(nextChar, filePtr);

	if (!isBPSEQ)
	{
	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad prev id: id %d\n", index);
		return 0;
	    }

	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);

	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad next id: id %d\n", index);
		return 0;
	    }

	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);
	}

	if (!fscanf(filePtr, "%d", pairIndex))
	{
	    fprintf(stderr, "Bad pair %d\n", index);
	    return 0;
	}

	if (!isBPSEQ)
	{
	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);

	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad trailing id: id %d\n", index);
		return 0;
	    }
	}

	return 1;
    }
}

TreeNode* CreateNode(
    int lowIndex,
    int* nextIndex,
    FILE* filePtr,
    unsigned char isBPSEQ,
    BaseData* returnData,
    int* returnPair)
{
    if (!ReadBase(filePtr, lowIndex, returnData, returnPair, isBPSEQ))
    {
	return NULL;
    }

    TreeNode* result = (TreeNode*)malloc(sizeof(TreeNode));
    result->numChildren = 0;
    result->children = NULL;
    result->lowBase.index = returnData->index;
    result->lowBase.base = returnData->base;

    if (!(*returnPair))
    {
	result->isPair = 0;
	*nextIndex = lowIndex + 1;
	return result;
    }

    result->isPair = 1;
    result->highBase.index = *returnPair;
    *nextIndex = lowIndex + 1;
    TreeNode* child;
    do
    {
	child = CreateNode(*nextIndex, nextIndex, filePtr, isBPSEQ, returnData, returnPair);
	result->numChildren += 1;
	result->children = (TreeNode**)realloc(result->children, sizeof(TreeNode*) * result->numChildren);
	result->children[result->numChildren - 1] = child;
    } while (child && *nextIndex < result->highBase.index);

    if (*nextIndex == result->highBase.index)
    {
	ReadBase(filePtr, *nextIndex, returnData, returnPair, isBPSEQ);
	result->highBase.base = returnData->base;
	*nextIndex += 1;
    }

    return result;
}

TreeNode* CreateFromFile(char* filename)
{
    // Figure out what kind of file we have and try to load it.
    char* extension = strrchr(filename, '.');
    unsigned char isBPSEQ = 0;
    if (extension && !strncmp(extension, ".bpseq", 6))
    {
	isBPSEQ = 1;
    }
    else if (extension && !strncmp(extension, ".ct", 3))
    {
	isBPSEQ = 0;
    }
    else
    {
	fprintf(stderr, "Unknown file type: %s\n", filename);
	return NULL;
    }

    FILE* filePtr = fopen(filename, "r");
    if (!filePtr)
    {
	fprintf(stderr, "Unable to open file: %s\n", filename);
	return NULL;
    }

    TreeNode* result = (TreeNode*)malloc(sizeof(TreeNode));
    result->numChildren = 0;
    result->children = NULL;
    result->lowBase.index = 0;
    result->lowBase.base = 0;
    result->highBase.index = 0;
    result->highBase.base = 0;
    result->isPair = 0;
    int nextIndex = 1;
    int pairIndex;
    BaseData bData;
    TreeNode* child = CreateNode(nextIndex, &nextIndex, filePtr, isBPSEQ, &bData, &pairIndex);
    while (child)
    {
	result->numChildren += 1;
	result->children = (TreeNode**)realloc(result->children, sizeof(TreeNode*) * result->numChildren);
	result->children[result->numChildren - 1] = child;
	child = CreateNode(nextIndex, &nextIndex, filePtr, isBPSEQ, &bData, &pairIndex);
    }

    fclose(filePtr);

    if (!result)
    {
	fprintf(stderr, "Empty or malformed file: %s\n", filename);
    }

    return result;
}

#define PrintTabs(indent) for (i = 0; i < indent; ++i) printf("   ");

void PrintTree(TreeNode* tree, int indent)
{
    char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    int i;

    PrintTabs(indent);
    if (tree->lowBase.index == 0)
    {
	printf("Root\n");
    }
    else if (tree->isPair)
    {
	printf("(%d %c - %c %d)\n",
	       tree->lowBase.index, bases[tree->lowBase.base],
	       bases[tree->highBase.base], tree->highBase.index);
    }
    else
    {
	printf("%d %c\n", tree->lowBase.index, bases[tree->lowBase.base]);
    }

    for (i = 0 ; i < tree->numChildren; ++i)
    {
	PrintTree(tree->children[i], indent + 1);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
	fprintf(stderr, "USAGE: RNAScoring <filename>\n");
	return 1;
    }

    TreeNode* tree = CreateFromFile(argv[1]);
    PrintTree(tree, 0);

    return 0;
}

