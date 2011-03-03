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
    U = 0x8
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
