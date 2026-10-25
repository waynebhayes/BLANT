#ifndef COMMUNITIES_H
#define COMMUNITIES_H

#include <stdio.h>
#include "structs.h"

COMMUNITY *CommunityAlloc(GRAPH *G, int id);
void CommunityFree(COMMUNITY *C);
void PrintCommunity(COMMUNITY *C);
COMMUNITY *CommunityAddNode(COMMUNITY *C, PARTITION *P, int node);
COMMUNITY *CommunityDelNode(COMMUNITY *C, PARTITION *P, int node);
int CommunityEdgeCount(COMMUNITY *C);
int CommunityEdgeOutwards(PARTITION *P, COMMUNITY *C);
PARTITION *PartitionAlloc(GRAPH *G);
PARTITION *PartitionAddCommunity(PARTITION *P, COMMUNITY *C);
PARTITION *PartitionDelCommunity(PARTITION *P, int c);
void PartitionFree(PARTITION *P);
void PartitionRead(FILE *fp, PARTITION *P);
int NodeInDegree(PARTITION *P, COMMUNITY *C, int node);
double ScorePartition(Boolean global, foint f);
double PerturbPartition(foint f);
Boolean MaybeAcceptPerturb(Boolean accept, foint f);
void HillClimbing(PARTITION *P, int tries);
void SAR(int iters, foint f);

#endif 
