#include <stdio.h>
#include "misc.h"
#include "rand48.h"
#include "graph.h"
#include "sets.h"
#include "sim_anneal.h"

#define TARGET_EDGE_DENSITY 0.5

#define VERBOSE 0 // 0 = no noisy outpt, 3 = lots, 1..2 is intermediate

/************************** Community routines *******************/
typedef struct _community {
    int id, n;
    SET *nodeSet;
    GRAPH *G; // the graph we came from
    double score;
    int edgesIn, edgesOut;
} COMMUNITY;

COMMUNITY *CommunityAlloc(GRAPH *G){
    COMMUNITY *C = Calloc(sizeof(COMMUNITY), 1);
    C->id = 0; // not sure if this will be handy or not...
    C->score = 0;
    C->n = 0;
    C->G = G;
    C->edgesIn = 0;
    C->edgesOut = 0;
    C->nodeSet = SetAlloc(G->n);
    return C;
}

void CommunityFree(COMMUNITY *C) {
    SetFree(C->nodeSet);
    Free(C);
}

// It's not an error to add a member multiple times
COMMUNITY *CommunityAddNode(COMMUNITY *C, int i) {
    if(!SetIn(C->nodeSet, i)){
	SetAdd(C->nodeSet, i);
	C->n++;
	assert(C->n <= C->G->n);
	assert(SetCardinality(C->nodeSet) == C->n);
    }
    else if(VERBOSE > 2)
	printf("%d Already in\n", i);
    return C;
}

// It's not an error to delete a node that's already not there
COMMUNITY *CommunityDelNode(COMMUNITY *C, int i) {
    if(SetIn(C->nodeSet, i)) {
	assert(C->n > 0);
	SetDelete(C->nodeSet, i);
	C->n--;
    }
    else if(VERBOSE>2)
	printf("Cannot find %d\n", i);
    return C;
}

int NodeInDegree(COMMUNITY *C, int node){
   int i, inDeg = 0;
   for(i=0;i<C->G->degree[node];i++){
        if(SetIn(C->nodeSet, C->G->neighbor[node][i])) inDeg++;
   }
   return inDeg;
}



void CommunityUpdate(COMMUNITY * oldCom, COMMUNITY * newCom, SET * nodes){
    // Update oldCom and newCom inEdges and outEdges based on node
    // Make sure to call ComUpdate BEFORE you move the nodes 
    GRAPH * G = oldCom->G;
    int * memberList = Calloc(SetCardinality(nodes), sizeof(int));
    int i, j, n = SetToArray(memberList, nodes);
#if VERBOSE > 2
    printf("BEFORE: oc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut);
#endif
    int * visited = Calloc(G->n, sizeof(int));     
    for(i = 0; i < G->n; ++i){
	visited[i] = 0; 
    }
    int within = 0, old = 0, new = 0, neither = 0; 
    for(i = 0; i < n; ++i){
	int node = memberList[i];
	visited[node] = 1;
	for(j = 0; j < G->degree[node]; ++j){
	    int u = G->neighbor[node][j];
	    if(!visited[u]){
		if(SetIn(nodes, u)){
		    // Edges that are connected within the movedNodes. 
		    --oldCom->edgesIn;
		    ++newCom->edgesIn;
		    ++within;
		} 
		else if(SetIn(oldCom->nodeSet, u)){
		    // Edges that are connected to oldCom
		    --oldCom->edgesIn;
		    ++old;
		} 	
		else if(SetIn(newCom->nodeSet, u)){
		    // Edges that are connected to newCom
		    ++newCom->edgesIn;
		    ++new;
		}
		else{
		    // Edges that are not connected to either community
		    ++newCom->edgesOut;
		    --oldCom->edgesOut;
		    ++neither;
		}
	    }
	}
    }
#if VERBOSE > 2
    printf("\nwithin %d, old %d, new %d, neither %d\n", within, old, new, neither);
#endif
    int diff = old-new;
    oldCom->edgesOut += diff;
    newCom->edgesOut += diff;
#if VERBOSE > 2
    printf("AFTER: oc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut);
#endif
    Free(memberList);
    Free(visited);    
}

int CommunityEdgeCount(COMMUNITY *C) {
    int *memberList = Calloc(C->n, sizeof(int));
    int i, j, n = SetToArray(memberList, C->nodeSet);
    assert(n == C->n);
    int numEdges=0;
    for(i=0; i<n; i++) {
	for(j=i+1; j<n; j++) {
	    int u = memberList[i], v = memberList[j];
	    assert(u < C->G->n && v < C->G->n);
	    if(GraphAreConnected(C->G,u,v)) ++numEdges;
	}
    }
    Free(memberList);
    return numEdges;
}

int CommunityEdgeOutwards(COMMUNITY *C){
    int out = 0;
    int *memberList = Calloc(C->n, sizeof(int));
    int i, j, k = SetToArray(memberList, C->nodeSet);
    for(i = 0; i < C->n; i++){
        int node = memberList[i];
	out += C->G->degree[node] - NodeInDegree(C, node);
    }
    Free(memberList);
    return out;
}

void PrintCommunity(COMMUNITY * C){
    int * memberList = Calloc(sizeof(int), C->n);
    int i, j, n = SetToArray(memberList, C->nodeSet);
    char space = '\0';
    for(i = 0; i < n; ++i) {
	if(space) putchar(space);
	printf("%s", C->G->name[memberList[i]]);
	space=' ';
    }
    printf("\n");
    Free(memberList);
}

/******************** Sets of non-overlapping Communities (partition) ***********/


typedef struct _communitySet {
    unsigned n; // current number of non-empty communities
    GRAPH *G; // the graph we came from
    COMMUNITY **C; // array of pointers to COMMUNITY
    int *where; // Tells where each node belongs to which community
    SET *moved; // Records moved nodes in case of rejection 
    SET *common; // In case merge of 2 communities have overlap, record them (Only for useful for overlapping communities)
    double total; // Cumulative score of partition
} PARTITION;

// Returns an empty partition containing no communities (they're not even allocated)
PARTITION *PartitionAlloc(GRAPH *G) {
    PARTITION *P = Calloc(sizeof(PARTITION), 1);
    P->G=G;
    P->C = Calloc(sizeof(COMMUNITY**), G->n);
    P->where = Calloc(sizeof(int), G->n);
    P->moved = SetAlloc(G->n);
    P->common = SetAlloc(G->n);
    return P;
}

PARTITION *PartitionAddCommunity(PARTITION *P, COMMUNITY *C) {
    assert(P->n < P->G->n);
    int i;
    for(i=0; i<P->G->n; i++)
        if(!P->C[i]) {
            P->C[i] = C;
            P->n++;
	    int * memberList = Calloc(sizeof(int), C->n);
	    int j, k, n = SetToArray(memberList, C->nodeSet);
	    for(j = 0; j < n; ++j)
		P->where[memberList[j]] = i;
	    Free(memberList);
            return P;// find empty community slot
        }
        return NULL;
}

PARTITION *PartitionDelCommunity(PARTITION *P, int c){
    if(P->C[c]){
	P->total -= P->C[c]->score;
        CommunityFree(P->C[c]);
        P->n--;
        if(c != P->n){
            int *memberList = Calloc(P->C[P->n]->n, sizeof(int));
            int i, j, k = SetToArray(memberList, P->C[P->n]->nodeSet);
            for(i = 0; i < k; i++)
                P->where[memberList[i]] = c;
            P->C[c] = P->C[P->n];
            Free(memberList);
        }
        P->C[P->n] = NULL;
    }
    return P;
}

void PartitionFree(PARTITION *P) {
    int i;
    for(i=0; i<P->n; i++) if(P->C[i]) {CommunityFree(P->C[i]);}
    Free(P->C);
    Free(P->where);
    SetFree(P->moved);
    SetFree(P->common);
    Free(P);
}

// oldCom is where to move the nodes back in case of reject
// newCom is where to get the nodes from
// moveDel is if MoveRandomNode deletes a community if it moves a node where it is the only node in the community
static int _moveOption = -1, _oldCom = -1, _newCom = -1, _moveDel = 0;

void MoveRandomNode(PARTITION *P){ 
    int u = P->G->n * drand48();
    int oldCom = P->where[u];
    COMMUNITY * oc = P->C[oldCom]; 
    _oldCom = oldCom;
    int newCom;
    do{newCom = (int)(P->n * drand48());}
    while(newCom == P->where[u]); //|| SetIn(P->C[newCom]->nodeSet, u)); // Ensures that it will move to a community that u is not already in
    COMMUNITY * nc = P->C[newCom];					 // Only useful for overlapping communities
    SetAdd(P->moved, u); 
    CommunityUpdate(oc, nc, P->moved);
    CommunityDelNode(oc, u);    
    CommunityAddNode(nc, u);
    P->where[u] = newCom;
#if VERBOSE > 1
    printf("Mv(%d,%d->%d) ", u, oldCom, newCom);
#endif
    _newCom = newCom;

    if(oc->n == 0){
	PartitionDelCommunity(P, oldCom);
	if(_newCom != P->n){
	    _oldCom = P->n; 
	}
	// In case oc only had 1 node and we moved it away, delete oc
	_moveDel = 1;
    }
#if VERBOSE > 2
    printf("o %d, n %d, P->n %d\n", _oldCom, _newCom, P->n);
#endif
}

// Merges C1 and C2 into just C1
void MergeCommunities(PARTITION *P, int c1, int c2){   
    COMMUNITY *C1 = P->C[c1];
    COMMUNITY *C2 = P->C[c2];
#if VERBOSE > 1    
    printf("Mg(%d,%d) |%d,%d| ", c1, c2, C1->n, C2->n);
#endif    
    int *memberList = Calloc(C2->n, sizeof(int));
    int i, j, n = SetToArray(memberList, C2->nodeSet);
    for(i = 0; i < C2->n; i++){
	int u = memberList[i]; 
	P->where[u] = c1;
	SetAdd(P->moved, u);
    }
    Free(memberList);
    CommunityUpdate(C2, C1, P->moved);

    SetIntersect(P->common, C1->nodeSet, C2->nodeSet);
    SetUnion(C1->nodeSet, C2->nodeSet, C1->nodeSet);
        
    C1->n = SetCardinality(C1->nodeSet);
    SetEmpty(C2->nodeSet);
    // Automatically takes care of P->n--;
    PartitionDelCommunity(P, c2);
    if(c1 == P->n){
	// In the case C1 is the last community, to preseve the list of communities
	// C1 is placed where C2 was previosuly.
	_newCom = c2;
    }
    else{
	_newCom = c1;    
    }
    _oldCom = P->n; 
    
}

// Splits P->C[c_id] into a new community that has numNodes random nodes from the old community
void SplitCommunity(PARTITION *P, int c_id, int numNodes){
    COMMUNITY *oldCom = P->C[c_id];
    assert(numNodes > 0);
    assert(numNodes < oldCom->n);
#if VERBOSE > 1
    printf("Sp(%d,|%d|->%d) ", c_id, oldCom->n, numNodes);
#endif
    COMMUNITY *newCom = CommunityAlloc(P->G);

    int *memberList = Calloc(sizeof(int), oldCom->n);
    int i, j, n = SetToArray(memberList, oldCom->nodeSet);

    for(i = oldCom->n-1; i > 0; i--){
        j = drand48() * (i+1);
        int temp = memberList[j];
        memberList[j] = memberList[i];
        memberList[i] = temp;
    }
    for(i = 0; i < numNodes; i++){
        int u = memberList[i];
	P->where[u] = P->n;
        SetAdd(P->moved, u);
    }
    for(i = 0; i < numNodes; ++i){
	int u = memberList[i];
	CommunityAddNode(newCom, u);
        CommunityDelNode(oldCom, u);
    }
    CommunityUpdate(oldCom, newCom, P->moved);
    Free(memberList);
    PartitionAddCommunity(P, newCom);
    _oldCom = c_id; 
    _newCom = P->n-1;
}


/*
    Measures
    TODO: Make it easier to swap out scoring functions, remove COMMUNITY * C requirement.
*/

double IntraEdgeDensity(COMMUNITY *C, int inEdges){
   return inEdges/(C->n *(C->n - 1) / 2.0);
}

double InterEdgeDensity(COMMUNITY *C, int edgesOut){
    int tot = C->G->n - C->n;
    //printf("tot = %d, edgesOut = %d\n", tot, edgesOut);
    return (double)edgesOut/(C->n * tot);
}

double NewmanAndGirvan(int edgesOut, int edgesIn, int gDegree){
    // Eq 15 on Pg 16 on the pdf viewer
    int inEdges = edgesIn;
    int cDegree = inEdges + edgesOut;
    return inEdges/gDegree - (cDegree/(2.0*gDegree) * cDegree/(2.0*gDegree));
}

double Conductance(int edgesIn, int edgesOut){
     return (double)(edgesOut/(edgesOut + edgesIn));
}

double HayesScore(COMMUNITY *C, int inEdges){ 
    //printf("inEdges = %d, C->n = %d\n", inEdges, C->n);
    if(C->n < 2)
        return 0;
    double eps = inEdges / ((C->n * (C->n-1))/2.0); 
    if(eps <= TARGET_EDGE_DENSITY)
	return 0;
    else
	return inEdges * (TARGET_EDGE_DENSITY); // inEdges*eps * (target/eps), to down-weight if eps is above the target
}


double PartitionAllScores(PARTITION *P) {
    int i, gDegree = 0;
    double total = 0.0;
    GRAPH *G = P->C[0]->G;
    for(i = 0; i < G->n; i++){
        gDegree += G->degree[i];
    }
    printf("G degree = %d\n", gDegree);
    for(i = 0; i < P->n; i++){
        // Run the expensive calculations once per community
	double score = 0.0;
        COMMUNITY *Com = P->C[i];
        //printf("Size of Community %d = %d\n", i, Com->n);
        if(Com->n > 0){
        int interEdges = CommunityEdgeOutwards(Com);
        int intraEdges = CommunityEdgeCount(Com);
        score += IntraEdgeDensity(Com, intraEdges);
        printf("Score of IaED is %g\n", score);
	double IeEd = InterEdgeDensity(Com, interEdges);
	printf("Score of IeED is %g\n", IeEd);
        score -= IeEd;
        double ng = NewmanAndGirvan(intraEdges, interEdges, gDegree);
        printf("Score of N&G is %g\n", ng);
        score += ng;
        double con = Conductance(intraEdges, interEdges);
        printf("Score of Conductance is %g\n", con);
        score -= con;
        double hayes = HayesScore(Com, intraEdges);
        printf("Hayes score = %g\n", hayes);
        score += hayes;
	printf("\nScore of community %d is %g\n\n", i, score);
        total += score;
        }
    }
    return total;
}


#define DEBUG 0
double ScorePartition(Boolean global, foint f){
    PARTITION *P = (PARTITION*) f.v;
#if 0   
    for(int i = 0; i < P->n; ++i){
	printf("Com %d, %d\n", i, P->C[i]->n);
    }
#endif
#if VERBOSE > 0
    printf("o = %d, n = %d, total n = %d\n", _oldCom, _newCom, P->n);
#endif
    if(_oldCom == -1 && _newCom == -1){
#if VERBOSE > 0
	printf("No moves have been made yet\n");	    
#endif
    }
    else{
	// In case of a split/merge reject P->C[n] will no longer exist
	if(_newCom == P->n){
	    int swap = _newCom;
	    _newCom = _oldCom;
	    _oldCom = swap;
	}
	int in, out, fail = 1;
	if(_oldCom != P->n){ 
	    COMMUNITY * old = P->C[_oldCom];
	    double oldBefore = old->score;
	    old->score = HayesScore(old, old->edgesIn);
	    P->total += old->score - oldBefore;
#if VERBOSE > 0	    
	    printf("os = %f, ob = %f change = %f ", old->score, oldBefore, old->score - oldBefore);
#endif
	}
	COMMUNITY * new = P->C[_newCom];
	double newBefore = newBefore = new->score;
	new->score = HayesScore(new, new->edgesIn);
	P->total += new->score - newBefore;
#if DEBUG
	in = CommunityEdgeCount(new);
	out = CommunityEdgeOutwards(new);
	printf("ns = %f, nb = %f change = %f\n", new->score, newBefore, new->score - newBefore);
	printf("FROM GROUND, In %d Out %d ns %f\n", CommunityEdgeCount(new), CommunityEdgeOutwards(new), HayesScore(new, CommunityEdgeCount(new)));
	if(in != new->edgesIn){
	    fail = 0;
	    printf("ERROR: New in %d vs %d\n", in, new->edgesIn);
	}
	if(out != new->edgesOut){
	    fail = 0;
	    printf("ERROR: New out %d vs %d\n", out, new->edgesOut);
	}
	assert(fail);
#endif
    }
#if VERBOSE > 0
    printf("Updated total = %f\n\n", P->total);
#endif
#if DEBUG   
    // This test basically defeats the purpose of incremental updating
    // Only use when ensuring robustness
    #endif
    return P->total;
}

// returns the CHANGE in score due to the move
double PerturbPartition(foint f) {
    //printf("Perturb\n");
    PARTITION *P = (PARTITION*) f.v;
    double before = P->total;
    int choice = drand48() * 3;
    //printf("Choice = %d\n", choice);
    
    if(choice == 0 && P->n > 1){
	MoveRandomNode(P);
	_moveOption = 0; 
    }
    else if(choice == 1 && P->n > 1){
	int c1, c2;
	do{
	    c1 = (int)(drand48() * P->n);
	    c2 = (int)(drand48() * P->n);
	}
	while(c1 == c2);
	MergeCommunities(P, c1, c2);
	_moveOption = 1; 
    }
    else{
	// Assumes that there are is at least 1 community with n > 1 	
	int c, numNodes;
	do{
	    c = drand48() * P->n;
	}
	while(P->C[c]->n == 1);
	numNodes = (drand48() * (P->C[c]->n - 1)) + 1;	
	SplitCommunity(P, c, numNodes);
	_moveOption = 2;
    }
    double after = ScorePartition(true, f);
#if VERBOSE > 0    
    printf("Before = %f, After = %f\n", before, after); 
#endif    
    return after-before;
}

Boolean MaybeAcceptPerturb(Boolean accept, foint f) {
#if VERBOSE > 0    
    printf("Accept? %d\n", accept);
#endif    
    PARTITION * P = (PARTITION *) f.v;
    double before = P->total;
    //printf("P->n = %d, _newCom = %d, _oldCom = %d\n", P->n, _newCom, _oldCom);

    if(accept) ; // do nothing?
    else {
	int * memberList = Calloc(sizeof(int), SetCardinality(P->moved));
	int i, j, n = SetToArray(memberList, P->moved);
	if(_moveOption == 1 || _moveDel){
	    PartitionAddCommunity(P, CommunityAlloc(P->G));
	}
	COMMUNITY * newCom = P->C[_newCom];
	COMMUNITY * oldCom = P->C[_oldCom];
	CommunityUpdate(newCom, oldCom, P->moved);
	for(i = 0; i < n; ++i){
	    int u = memberList[i];
	    if(!SetIn(P->common, u)) // Only delete nodes that weren't shared before the change -> (Only for overlapping communities)
		CommunityDelNode(newCom, u);
	    CommunityAddNode(oldCom, u);
	    P->where[u] = _oldCom;
	}

	Free(memberList);
	if(_moveOption == 2)
	    PartitionDelCommunity(P, _newCom);

	double after = ScorePartition(true, f);
	//printf("Maybe Before = %f, After = %f\n\n", before, after);
	if(before > after)
	    fprintf(stderr, "ERROR: Rejection failed, before > after\n");
#if 1	    
	if(_moveDel){
	    for(int i = 0; i < P->n; ++i)
		assert(P->C[i]->n > 0);
	}
#endif
    }
    SetEmpty(P->moved);
    SetEmpty(P->common);
    _moveDel = 0;
#if VERBOSE > 0   
    printf("Current total = %f\n\n", P->total);
#endif
    return accept;
}

void HillClimbing(PARTITION *P, int tries){
    int i;
    foint f;
    f.v = (void*)P;

    printf("Beginning score %g\n", P->total);
    for(i = 0; i < tries; i++) {
        //for(int l = 0; l < P->n; l++) printf("%d, %p, size = %d\n", l, P->C[l], P->C[l]->n);
	double delta = PerturbPartition(f);
	printf("delta = %g...", delta);

        if(delta > 0) {
            printf("accept\n");
	    MaybeAcceptPerturb(true, f);
        }
        else{
            printf("reject\n");
	    MaybeAcceptPerturb(false, f);
        }
    }
    printf("Final score %g\n", P->total);
}

// EXTRA_ASSERTS will significantly degrade performance
#define EXTRA_ASSERTS 0
void SAR(int iters, foint f){

    PARTITION * P = f.v;
    int fail = 1, in, out;
    double ground = 0, withInfo, stored = 0;
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * com = P->C[i];
	in = CommunityEdgeCount(com);
	out = CommunityEdgeOutwards(com);
#if VERBOSE > 2
	printf("\nCom %d FROM GROUND, In %d Out %d\n", i, in, out);
#endif
	if(in != com->edgesIn){
	    fail = 0;
	    printf("Com %d ERROR: Old from ground in %d vs stored in %d\n", i, in, com->edgesIn);
	}
	if(out != com->edgesOut){
	    fail = 0;
	    printf("Com %d ERROR: Old from ground out %d vs stored out %d\n", i, out, com->edgesOut);
	}
#if EXTRA_ASSERTS
	stored += com->score;
	ground += HayesScore(com, CommunityEdgeCount(com));
	withInfo += HayesScore(com, com->edgesIn);
	assert(P->total - stored < 0.001 && stored - P->total < 0.001);
	assert(P->total - ground < 0.001 && ground - P->total < 0.001);
	assert(P->total - withInfo < 0.001 && withInfo - P->total < 0.001);
#endif
    }
    assert(fail);

}


#define RANDOM_START 1
#define CHECK_OVERLAP 1
int main(int argc, char *argv[])
{

    int i, j;
    srand48(GetFancySeed(false));

    Boolean sparse=maybe, supportNames = true;
    FILE *fp = Fopen(argv[1], "r"); // edge list file is given on the command line
    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNames, false);
    printf("G has %d nodes\n", G->n);

    PARTITION *P = PartitionAlloc(G);
#if RANDOM_START
    int numCommunities = 8; // communities numbered 0 through numCommunities-1 inclusive
    printf("Starting with %d random communities\n", numCommunities);
    for(i=0; i<numCommunities; i++) PartitionAddCommunity(P, CommunityAlloc(G));

    for(i=0; i<G->n; i++) {
	int which = (int)(drand48() * numCommunities);
	CommunityAddNode(P->C[which], i);
        P->where[i] = which;
	//printf("%d->%d, ", i, which);
    }
#else
	COMMUNITY *C = CommunityAlloc(G);
	for(i=0; i<n; i++) if(!SetIn(nodesUsed,nodeArray[i]) || drand48() > 0.5) {
	    SetAdd(nodesUsed,nodeArray[i]); CommunityAddNode(C,nodeArray[i]); ++numAdded;
	}
	//printf("%d ", numAdded); fflush(stdout);
	assert(C->n >0 && C->n < G->n && C->n==numAdded);
	PartitionAddCommunity(P, C);
        //printf("Size of Community = %d\n", C->n);
	
    }
    //printf("\n%d communities, score = %g\n", P->n, ScorePartition(P));
    SetFree(nodesUsed);
#endif
    printf("%d\n", P->n);
    
    SET * found = SetAlloc(G->n);    
    for(int k = 0; k < P->n; ++k){
	if(P->C[k]->n == 0)
	    printf("ERROR: Com %d has 0 nodes\n", k);
#if CHECK_OVERLAP 
	int * memberList = Calloc(sizeof(int), P->C[k]->n);
	int l, m, n = SetToArray(memberList, P->C[k]->nodeSet);
	
	for(l = 0; l < n; ++l){
	    if(SetIn(found, memberList[l]))
		printf("ERROR: %d is in an overlapping community\n", memberList[l]);
	    SetAdd(found, memberList[l]);
	}
	Free(memberList);
#endif
    }
    SetFree(found);

    for(int i = 0; i < P->n; ++i){
	COMMUNITY * C = P->C[i];
	double s = HayesScore(C, CommunityEdgeCount(C)); 
	C->score = s;
	C->edgesIn = CommunityEdgeCount(C);
	C->edgesOut = CommunityEdgeOutwards(C);
	P->total += s;
	printf("%d has in %d, out %d, n %d, score %f\n", i, C->edgesIn, C->edgesOut, C->n, C->score);
    }


#if 0
    HillClimbing(P, 200);
#else
    foint f;
    f.v = P;
    SIM_ANNEAL *sa = SimAnnealAlloc(1, f, PerturbPartition, ScorePartition, MaybeAcceptPerturb, G->n * G->numEdges, 0, 0, SAR);
    //SimAnnealSetSchedule(sa, 50, 10);
    SimAnnealAutoSchedule(sa); // to automatically create schedule
    //sa->tInitial = sa->tDecay = sa->temperature = 0; // equivalent to hill climbing
    SimAnnealRun(sa); // returns >0 if success, 0 if not done and can continue, <0 if error
    // foint SimAnnealSol(SIM_ANNEAL *sa))
    SimAnnealFree(sa);
#endif
    int nodes = 0, biggest = 0, num, which=-1;
    for(int j = 0; j < P->n; ++j){
	num = P->C[j]->n;
	int inEdges = CommunityEdgeCount(P->C[j]);
	printf("Com %d has %d nodes, %d edges, with edge density %g\n", j, num, inEdges, inEdges/(num*(num-1)/2.0));
	nodes += num;
	if(HayesScore(P->C[j], inEdges) && num > biggest) {
	    biggest = num;
	    which = j;
	}
	
    }
    assert(nodes == G->n);
    printf("Final score = %f\n", P->total);
    printf("Biggest community is %d, with %d nodes:\n", which, biggest);
    PrintCommunity(P->C[which]);
    printf("Attempting Partition Free\n");
    PartitionFree(P);
    printf("Partition Free completed\n");
    GraphFree(G);
    fclose(fp);
    return 0;
}
