#include <stdio.h>
#include "misc.h"
#include "rand48.h"
#include "graph.h"
#include "sets.h"
#include "sim_anneal.h"


#define TARGET_EDGE_DENSITY 0.5
#define VERBOSE 0 // 0 = no noisy outpt, 3 = lots, 1..2 is intermediate
#define DEBUG 0


/************************** Community routines *******************/
typedef struct _community {
    int id, n;
    int * nodeSet;
    GRAPH *G; // the graph we came from
    double score;
    int edgesIn, edgesOut;
} COMMUNITY;

double (*pCommunityScore)(COMMUNITY * C) = NULL;  

/******************** Sets of non-overlapping Communities (partition) ***********/
typedef struct _communitySet {
    unsigned n; // current number of non-empty communities
    GRAPH *G; // the graph we came from
    COMMUNITY **C; // array of pointers to COMMUNITY
    int *whichCommunity; // Tells where each node belongs to which community
    int *whichMember; // Within the community, tells at which index the node is located at 
    //SET *common; // In case merge of 2 communities have overlap, record them (Only for useful for overlapping communities)
    double total; // Cumulative score of partition
    int * visited; // A bool vector for community update (Moved here to stop allocating and freeing repeadetly) 
    int * marked; // Marking which ones will be moved (Essentially SET * moved but in int * format)  
    int * toMove; // Marking, but only holds which ones to move instead of marking them to move 
    int numMoved; // Number of nodes that will be moved
} PARTITION;

COMMUNITY *CommunityAlloc(GRAPH *G, int id){
    COMMUNITY *C = Calloc(sizeof(COMMUNITY), 1);
    C->id = id; 
    C->score = 0;
    C->n = 0;
    C->G = G;
    C->edgesIn = 0;
    C->edgesOut = 0;
    C->nodeSet = Calloc(sizeof(int), G->n);
    return C;
}

void CommunityFree(COMMUNITY *C) {
    Free(C->nodeSet);
    Free(C);
}

// Make sure to do error checking outside of CommunityAddNode
// Remember to update PARTITION->whichCommunity and PARTITION->whichMember
COMMUNITY *CommunityAddNode(COMMUNITY *C, int * whichMember, int node) {
    C->nodeSet[C->n] = node;
    whichMember[node] = C->n;
    ++C->n;
    return C;
}

// Similarily, error check DelNode as well outside
// Shift around whichMember 
COMMUNITY *CommunityDelNode(COMMUNITY *C, int * whichMember, int node) {
    int index = whichMember[node];
        if(index != C->n - 1){
	int last = C->nodeSet[C->n-1];
	C->nodeSet[index] = last;
	whichMember[last] = index;
    } 
    C->n--; 
    return C;
}

// Potential move one node function that I should have probably thought of 6 months ago
void MoveOneNode(PARTITION * P, int node, int dest){ 
     
    COMMUNITY * oldCom = P->C[P->whichCommunity[node]];
    COMMUNITY * newCom = P->C[dest];
    GRAPH * G = P->G; 
#if VERBOSE > 2
    printf("node %d, org %d, dest %d\n", node, P->whichCommunity[node], dest);
    printf("BEFORE\noc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut); 
#endif
    int old = 0, new = 0;
    for(int i = 0; i < G->degree[node]; ++i){
	int neighbor = G->neighbor[node][i];
    #if DEBUG
	printf("%d ", neighbor);
    #endif
	if(neighbor == node){
	    // Self loops
	    // This works but I'm not entirely sure why
	    continue;
	}
	if(!P->visited[neighbor]){
	    if(P->marked[neighbor]){
	    #if DEBUG	
		printf("W ");
	    #endif
		--oldCom->edgesIn;
		++newCom->edgesIn;
	    }
	    else if(P->whichCommunity[neighbor] == oldCom->id){
	    #if DEBUG
		printf("O ");
	    #endif
		--oldCom->edgesIn;
		++old;
	    }
	    else if(P->whichCommunity[neighbor] == newCom->id){
	    #if DEBUG
		printf("N ");
	    #endif
		++newCom->edgesIn;
		++new;
	    }
	    else{
	    #if DEBUG
		printf("- ");
	    #endif
		++newCom->edgesOut;
		--oldCom->edgesOut;
	    }
	}
#if DEBUG
	else
	    printf("visited");
#endif
    }
    int diff = old - new;
    oldCom->edgesOut += diff;
    newCom->edgesOut += diff;
    P->visited[node] = 1;
    CommunityDelNode(oldCom, P->whichMember, node);
    CommunityAddNode(newCom, P->whichMember, node);
    P->whichCommunity[node] = dest;
#if VERBOSE > 2 
    printf("AFTER\noc size = %d, nc size = %d\noc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->n, newCom->n, oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut);
#endif


}

void PrintCommunity(COMMUNITY * C){
    char space = '\0';
    for(int i = 0; i < C->n; ++i) {
	if(space) putchar(space);
	printf("%s", C->G->name[C->nodeSet[i]]);
	space=' ';
    }
    printf("\n");
}

int NodeInDegree(PARTITION *P, COMMUNITY *C, int node){
    int i, inDeg = 0;
    for(i=0;i<C->G->degree[node];i++){
	if(P->whichCommunity[C->G->neighbor[node][i]] == C->id) ++inDeg; 
    }
    return inDeg;
}

int CommunityEdgeCount(COMMUNITY *C) {
    int numEdges=0;
    for(int i=0; i<C->n; i++) {
	for(int j=i+1; j<C->n; j++) {
	    int u = C->nodeSet[i], v = C->nodeSet[j];
	    assert(u < C->G->n && v < C->G->n);
	    if(GraphAreConnected(C->G,u,v)) ++numEdges;
	}
    }
    return numEdges;
}

int CommunityEdgeOutwards(PARTITION * P, COMMUNITY *C){
    int out = 0;
    for(int i = 0; i < C->n; ++i){
	int u = C->nodeSet[i];
	out += C->G->degree[u] - NodeInDegree(P, C, u);
    }
    return out;
}

// Returns an empty partition containing no communities (they're not even allocated)
PARTITION *PartitionAlloc(GRAPH *G) {
    PARTITION *P = Calloc(sizeof(PARTITION), 1);
    P->G=G;
    P->C = Calloc(sizeof(COMMUNITY**), G->n);
    P->whichCommunity = Calloc(sizeof(int), G->n);
    P->whichMember = Calloc(sizeof(int), G->n);
    P->visited = Calloc(sizeof(int), G->n);
    P->marked = Calloc(sizeof(int), G->n);
//  P->common = SetAlloc(G->n);
    P->toMove = Calloc(sizeof(int), G->n);
    P->numMoved = 0;
    return P;  
}

PARTITION *PartitionAddCommunity(PARTITION *P, COMMUNITY *C) {
    assert(P->n < P->G->n);
    int i;
    for(i=0; i<P->G->n; i++){
	if(!P->C[i]) {
            P->C[i] = C;
            P->n++;
	    int index = 0;
	    for(int j = 0; j < C->n; ++j){
		int u = C->nodeSet[j];
		P->whichCommunity[u] = i;
		P->whichMember[u] = index;
		++index;
	    }
	    return P;// find empty community slot
        }
    }
    printf("Updated P->n %d\n", P->n);
    return NULL;
}

PARTITION *PartitionDelCommunity(PARTITION *P, int c){
#if VERBOSE > 2 
    printf("Deleting com %d, P->n = %d\n", c, P->n);
#endif
    COMMUNITY * C = P->C[c];
    if(C){
	P->total -= C->score;
        P->n--;
	// if C is not the last community, write over the deleted community's info
        if(c != P->n){
	    int index = 0;
	    COMMUNITY * lastCom = P->C[P->n];
	    for(int i = 0; i < lastCom->n; i++){
		P->whichCommunity[lastCom->nodeSet[i]] = c;
		P->whichMember[lastCom->nodeSet[i]] = index;
		++index;
	    }
	    CommunityFree(C);
	    P->C[c] = P->C[P->n];
	    P->C[c]->edgesIn = lastCom->edgesIn;
	    P->C[c]->edgesOut = lastCom->edgesOut;
	    P->C[c]->id = c;
	}
	else
	    CommunityFree(C);

	P->C[P->n] = NULL;
    }

    return P;
}

void PartitionFree(PARTITION *P) {
    int i;
    for(i=0; i<P->n; i++) if(P->C[i]) {CommunityFree(P->C[i]);}
    Free(P->C);
    Free(P->whichCommunity);
    Free(P->whichMember);
    //SetFree(P->common);
    Free(P->visited);
    Free(P->marked);
    Free(P->toMove);
    Free(P);
}

// oldCom is where to move the nodes back in case of reject
// newCom is where to get the nodes from
// moveDel is if MoveRandomNode deletes a community if it moves a node where it is the only node in the community
static int _moveOption = -1, _oldCom = -1, _newCom = -1, _moveDel = 0;

// THESE MOVE OPTIONS ASSUME NO OVERLAPPING COMMUNITIES

void MoveRandomNode(PARTITION *P){  
    int u = P->G->n * drand48();
    int oldCom = P->whichCommunity[u];
    COMMUNITY * oc = P->C[oldCom]; 
    _oldCom = oldCom;
    int newCom;
    do{newCom = (int)(P->n * drand48());}
    while(newCom == P->whichCommunity[u]); 
    COMMUNITY * nc = P->C[newCom];	
    P->marked[u] = 1;
    P->toMove[P->numMoved++] = u;
    MoveOneNode(P, u, newCom);     
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
    int index = 0;
    for(int i = 0; i < C2->n; i++){
	int u = C2->nodeSet[i]; 
	//printf("Marked %d ", u);
	P->marked[u] = 1;
	P->toMove[P->numMoved++] = u;
    }
    int iters = C2->n;
    for(int i = 0; i < iters; ++i){
	MoveOneNode(P, C2->nodeSet[i], c1);
    }
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
    COMMUNITY *newCom = CommunityAlloc(P->G, P->n);     
    PartitionAddCommunity(P, newCom);

    int u;
    for(int i = 0; i < numNodes; i++){	
	do{
	    int index = (int)(drand48() * oldCom->n);
	    u = oldCom->nodeSet[index];
	    //printf("u %d, i %d, located %d\n", u, index, P->whichCommunity[u]);
	}
	while(P->marked[u]);
	P->marked[u] = 1;
	P->toMove[P->numMoved++] = u;
    } 
 
    for(int i = 0; i < P->numMoved; ++i){
	MoveOneNode(P, P->toMove[i], P->n - 1);
    }

    _oldCom = c_id; 
    _newCom = P->n-1;
    assert(P->numMoved == numNodes);
}


/*
    Measures

    TODO: Make it easier to swap out scoring functions
	  Prototypes
*/


double IntraEdgeDensity(COMMUNITY *C){
   return C->edgesIn / (C->n *(C->n - 1) / 2.0);
}

double InterEdgeDensity(COMMUNITY *C){
    int tot = C->G->n - C->n;
    //printf("tot = %d, edgesOut = %d\n", tot, edgesOut);
    return (double)C->edgesOut/(C->n * tot);
}

double NewmanAndGirvan(COMMUNITY * C){
    // Eq 15 on Pg 16 on the pdf viewer
   
    // FIXME: Double check what exactly is needed
    return 0;
    //return C->inEdges/C->gDegree - (cDegree/(2.0*gDegree) * cDegree/(2.0*gDegree));
}

double Conductance(COMMUNITY * C){
     return (double)(C->edgesOut/(C->edgesOut + C->edgesIn));
}

double HayesScore(COMMUNITY *C){ 
    //printf("inEdges = %d, C->n = %d\n", inEdges, C->n);
    if(C->n < 2)
        return 0;
    double eps = C->edgesIn / ((C->n * (C->n-1))/2.0); 
    if(eps <= TARGET_EDGE_DENSITY)
	return 0;
    else
	return C->edgesIn * (TARGET_EDGE_DENSITY); // inEdges*eps * (target/eps), to down-weight if eps is above the target
}


double ScorePartition(Boolean global, foint f){
    PARTITION *P = (PARTITION*) f.v;
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
	    old->score = pCommunityScore(old);
	    P->total += old->score - oldBefore;
#if VERBOSE > 0	    
	    printf("os = %f, ob = %f change = %f ", old->score, oldBefore, old->score - oldBefore);
#endif 

	}
	COMMUNITY * new = P->C[_newCom];
	double newBefore = new->score;
	new->score = pCommunityScore(new);
	P->total += new->score - newBefore;
    }
#if VERBOSE > 0
    printf("Updated total = %f\n\n", P->total);
#endif
    return P->total;
}

// returns the CHANGE in score due to the move
double PerturbPartition(foint f) {
    //printf("Perturb\n");
    PARTITION *P = (PARTITION*) f.v;
    double before = P->total;
    
    int choice = drand48() * 3;
    P->numMoved = 0;    
    for(int i = 0; i < P->G->n; ++i){
	P->visited[i] = 0;
	P->marked[i] = 0;
    }
#if DEBUG
    for(int i = 0; i < P->n; ++i){ 
    	COMMUNITY * temp = P->C[i];
	printf("Com %d, has %d nodes, in %d, out %d\n", i, temp->n, temp->edgesIn, temp->edgesOut);
    }
#endif
    int comNums = 1;
    for(int i = 0; i < P->G->n; ++i){	    
	if(P->whichCommunity[i] >= P->n){
	    comNums = 0;
	    printf("ERROR: Node %d is in Com %d > %d\n", i, P->whichCommunity[i], P->n);
	}	    
    }
    assert(comNums);

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
	int c, numNodes, tries = 0;
	do{
	    c = drand48() * P->n;
	    ++tries; 
	}
	while(P->C[c]->n == 1 && tries < 500);

	if(tries >= 500){
	    // If after a while, it can't find a community with n > 1,
	    // just move
	    MoveRandomNode(P);
	    _moveOption = 0;
	}
	else{
	    numNodes = (drand48() * (P->C[c]->n - 1)) + 1;	
	    SplitCommunity(P, c, numNodes);
	    _moveOption = 2;
	}
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
	if(_moveOption == 1 || _moveDel){
	    PartitionAddCommunity(P, CommunityAlloc(P->G, P->n));
	}
	COMMUNITY * newCom = P->C[_newCom];
	COMMUNITY * oldCom = P->C[_oldCom];

	
	// In case of reject -> Look at marked, MoveOneNode for each one 
	for(int i = 0; i < P->G->n; ++i){
	    P->visited[i] = 0; // Make visited a set?
	}

	for(int i = 0; i < P->numMoved; ++i){
	#if DEBUG
	    printf("%d->%d\n", i, _oldCom);
	#endif
	    MoveOneNode(P, P->toMove[i], _oldCom);
	}
	if(_moveOption == 2)
	    PartitionDelCommunity(P, _newCom);

	double after = ScorePartition(true, f);
	//printf("Maybe Before = %f, After = %f\n\n", before, after);
	if(before > after)
	    fprintf(stderr, "ERROR: Rejection failed, before > after\n");    
	if(_moveDel){
	    for(int i = 0; i < P->n; ++i){
		//printf("Com %d has %d nodes\n", i, P->C[i]->n);
		assert(P->C[i]->n > 0);
	    }
	}

    }
    //SetEmpty(P->common);
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
    int fail = 1, in, out, best_id = -1;
    double ground = 0, withInfo, stored = 0, best = -1;
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * com = P->C[i];
	in = CommunityEdgeCount(com);
	out = CommunityEdgeOutwards(P, com);
	if(com->score > best){
	    best = com->score;
	    best_id = com->id;
	}
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
	ground += pCommunityScore(com);
	withInfo += pCommunityScore(com);
	assert(P->total - stored < 0.001 && stored - P->total < 0.001);
	assert(P->total - ground < 0.001 && ground - P->total < 0.001);
	assert(P->total - withInfo < 0.001 && withInfo - P->total < 0.001);
#endif
    }
    assert(fail);
    printf("\tBest: Com %d, with n %d, score %g  \tP->n = %d, Total score = %g", best_id, P->C[best_id]->n, best, P->n, P->total);
}


#define RANDOM_START 1
#define CHECK_OVERLAP 1
int main(int argc, char *argv[])
{

    pCommunityScore = HayesScore;
    int i, j;
    srand48(GetFancySeed(false));

   
    Boolean sparse=maybe, supportNames = true;
    FILE *fp = Fopen(argv[1], "r"); // edge list file is given on the command line
    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNames, false);
    printf("G has %d nodes\n", G->n);

    PARTITION *P = PartitionAlloc(G);
#if RANDOM_START
    int numCommunities = 2; // communities numbered 0 through numCommunities-1 inclusive
    printf("Starting with %d random communities\n", numCommunities);
    for(i=0; i<numCommunities; i++) PartitionAddCommunity(P, CommunityAlloc(G, i));

    for(i=0; i<G->n; i++) {
	int which = (int)(drand48() * numCommunities);
	//printf("%d->%d, ", i, which);
	CommunityAddNode(P->C[which], P->whichMember, i);
        P->whichCommunity[i] = which;
    }
    
#else
    for(int i = 0; i < ){
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
	int n = P->C[k]->n;
	for(int l = 0; l < n; ++l){
	    int node = P->C[k]->nodeSet[l];
	    if(SetIn(found, node))
		printf("ERROR: %d is in an overlapping community\n", node);
	    SetAdd(found, node);
	}
	
#endif
    }
    SetFree(found);

    for(int i = 0; i < P->n; ++i){
	COMMUNITY * C = P->C[i];
	double s = pCommunityScore(C); 
	C->score = s;
	C->edgesIn = CommunityEdgeCount(C);
	C->edgesOut = CommunityEdgeOutwards(P, C);
	P->total += s;
	printf("%d has in %d, out %d, n %d, score %f\n", i, C->edgesIn, C->edgesOut, C->n, C->score);
    }

    	
#if 0
    HillClimbing(P, 200);
#else
    foint f;
    f.v = P;
    	
    SIM_ANNEAL *sa = SimAnnealAlloc(1, f, PerturbPartition, ScorePartition, MaybeAcceptPerturb, 10*G->n*G->numEdges,0,0,SAR);
    SimAnnealSetSchedule(sa, 64, 16);
    //SimAnnealAutoSchedule(sa); // to automatically create schedule
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
	if(HayesScore(P->C[j]) && num > biggest) {
	    biggest = num;
	    which = j;
	}
	
    }
    assert(nodes == G->n);
    printf("Final score = %f\n", P->total);
    printf("Biggest community is #%d, with %d nodes:\n", which, biggest);
    PrintCommunity(P->C[which]);
    printf("Attempting Partition Free\n");
    PartitionFree(P);
    printf("Partition Free completed\n");
    GraphFree(G);
    fclose(fp);
    return 0;
}
