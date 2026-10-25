#ifndef STRUCTS_H
#define STRUCTS_H

#include "graph.h"
#include "sets.h"
#include "time.h"

typedef struct _community {
     int id, n;
     int * nodeSet;
     GRAPH *G; // the graph we came from
     double score;
     int edgesIn, edgesOut;
} COMMUNITY;

/******************** Sets of non-overlapping Communities (partition) ***********/
typedef struct _communitySet {
    unsigned n; // current number of non-empty communities
    GRAPH *G; // the graph we came from
    COMMUNITY **C; // array of pointers to COMMUNITY
    int *whichCommunity; // Tells where each node belongs to which community
    int *whichMember; // Within the community, tells at which index the node is located at 
    //SET *common; // In case merge of 2 communities have overlap, record them (Only for useful for overlapping communities)
    double total; // Cumulative score of partition
    SET *visited; // A set for community update to avoid double counting (Only used in MoveOneNode)  
    SET *marked; // Marking which ones will be moved   
    int *toMove; // Marking, but only holds which ones to move instead of marking which nodes in the graph to move 
    int numMoved; // Number of nodes that will be moved
    
    // Everything below here is just for timing
    int iters;
    clock_t clk;
} PARTITION;


#endif // STRUCTS_H
