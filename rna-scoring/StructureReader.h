/*
    A class to hold an RNA structure.
*/

#ifndef RNASTRUCTURE_H
#define RNASTRUCTURE_H

typedef enum _Base
{
    A = 0x1,
    C = 0x2,
    G = 0x4,
    U = 0x8,
    
    //ZS: Support for noncanonical bases using the natural binary code
    M = 0x3,  //3: A or C (1+2)
    R = 0x5,  //5: G or A (1+4) 
    S = 0x6,  //6: C or G (2+4) 
    V = 0x7,  //7: G or A or C (1+2+3) 
    W = 0x9,  //9: or A (1+8) 
    Y = 0xA,  //10: U or C (2+8) 
    H = 0xB,  //11: U or A or C (1+2+8) 
    K = 0xC,  //12: U or G (4+8) 
    D = 0xD,  //13: U or G or A (1+4+8)
    B = 0xE,  //14: U or G or C (2+4+8)
    N = 0xF   //15: any nucleotide (1+2+4+8)
    
} Base;

/* Data about a single base */
typedef struct _BaseData
{
    int index;
    Base base;
} BaseData;

/* The graph node data structure */
typedef struct _TreeNode
{
    BaseData lowBase; /* The low base. */
    BaseData highBase; /* The high base, if a pair. */
    unsigned char isPair; /* Non-zero if this is a pair. */
    int numChildren; /* The number of children nodes. */
    struct _TreeNode** children; /* The child nodes. */
} TreeNode;

#endif
