// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.


// Warning: setting PA = 0 and ND = 1 turns of ALL checking... but could make this up to 10x faster
#define PARANOID_ASSERTS 0
#define NDEBUG 1

#include <stdio.h>
#include "misc.h"
#include "rand48.h"
#include "graph.h"
#include "sets.h"
#include "sim_anneal.h"

#define TARGET_EDGE_DENSITY 0.5
#define VERBOSE 0 // 0 = no noisy outpt, 3 = lots, 1..2 is intermediate
#define DEBUG 0
#define MOVE_ONLY 0
#define PRINT_ALL_COMMUNITIES 0

/************************** Community routines *******************/
typedef struct _community {
    int id, n;
    int * nodeSet;
    GRAPH *G; // the graph we came from
    double score;
    int edgesIn, edgesOut;
} COMMUNITY;

double (*pCommunityScore)(COMMUNITY * C, int fakeN) = NULL;  

/******************** Sets of non-overlapping Communities (partition) ***********/
typedef struct _communitySet {
    unsigned n; // current number of non-empty communities
    GRAPH *G; // the graph we came from
    COMMUNITY **C; // array of pointers to COMMUNITY
    int *whichCommunity; // Tells where each node belongs to which community
    int *whichMember; // Within the community, tells at which index the node is located at 
    //SET *common; // In case merge of 2 communities have overlap, record them (Only for useful for overlapping communities)
    double total; // Cumulative score of partition
    SET *visited; // A bool vector for community update (Moved here to stop allocating and freeing repeadetly) 
    SET *marked; // Marking which ones will be moved (Essentially SET * moved but in int * format)  
    int * toMove; // Marking, but only holds which ones to move instead of marking which nodes in the graph to move 
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
    C->G = NULL;
#if DEBUG
    printf("Freeing Com %d, %p\n", C->id, C);
#endif
    Free(C);
}

void PrintCommunity(COMMUNITY * C){
    /*char space = '\0';
    for(int i = 0; i < C->n; ++i) {
	if(space) putchar(space);
	printf("%s", C->G->name[C->nodeSet[i]]);
	space=' ';
    }
    printf("\n");*/

    for(int i = 0; i < C->n; ++i){
	printf("%d ", C->nodeSet[i]);
    }
    printf("\n");
}


// Make sure to do error checking outside of CommunityAddNode
COMMUNITY *CommunityAddNode(COMMUNITY *C, PARTITION * P, int node) {
#if DEBUG
    printf("Add C %d, node %d, before size = %d\n", C->id, node, C->n);
#endif 
    C->nodeSet[C->n] = node;
    P->whichMember[node] = C->n;
    P->whichCommunity[node] = C->id;
    ++C->n;

    return C;
}

// Similarily, error check DelNode as well outside 
COMMUNITY *CommunityDelNode(COMMUNITY *C, PARTITION * P, int node) {
#if DEBUG
    printf("Del C %d, node %d, In com %d index %d\n", C->id, node, P->whichCommunity[node], P->whichMember[node]);
#endif
    int index = P->whichMember[node];
    if(index != C->n - 1){
	int last = C->nodeSet[C->n-1];
	C->nodeSet[index] = last;
	P->whichMember[last] = index;
    }
    --C->n;  
    return C;
}

static void MoveOneNode(PARTITION * P, int node, int dest){  
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
	if(!SetIn(P->visited,neighbor)){
	    if(SetIn(P->marked, neighbor)){
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
    SetAdd(P->visited, node);
    
#if VERBOSE > 2 
    printf("AFTER\noc size = %d, nc size = %d\noc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->n, newCom->n, oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut);
#endif


}


int NodeInDegree(PARTITION *P, COMMUNITY *C, int node){
    int i, inDeg = 0;
    for(i=0;i<C->G->degree[node];i++){
	int neigh = C->G->neighbor[node][i];
    
	//printf("N %d in com %d\n",neigh, P->whichCommunity[neigh]);
    
	if(P->whichCommunity[C->G->neighbor[node][i]] == C->id) ++inDeg; 
    }
    //printf("Node %d in Com %d inDeg %d\n", node, P->whichCommunity[node], inDeg);

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
    //printf("Com %d size %d, In %d\n", C->id, C->n, numEdges);
    return numEdges;
}

int CommunityEdgeOutwards(PARTITION * P, COMMUNITY *C){
    int out = 0;
    for(int i = 0; i < C->n; ++i){
	int u = C->nodeSet[i];
	int degree = C->G->degree[u];
	int inDeg = NodeInDegree(P, C, u);
	out += degree - inDeg;
    
	//printf("node %d degree %d, in %d\n", u, degree, inDeg);
    
    }

//    printf("Com %d size %d, out %d\n", C->id, C->n, out);

    return out;
}

// Returns an empty partition containing no communities (they're not even allocated)
PARTITION *PartitionAlloc(GRAPH *G) {
    PARTITION *P = Calloc(sizeof(PARTITION), 1);
    P->G=G;
    P->C = Calloc(sizeof(COMMUNITY**), G->n);
    P->whichCommunity = Calloc(sizeof(int), G->n);
    P->whichMember = Calloc(sizeof(int), G->n);
    P->visited = SetAlloc(G->n);
    P->marked = SetAlloc(G->n);
//  P->common = SetAlloc(G->n);
    P->toMove = Calloc(sizeof(int), G->n);
    P->numMoved = 0;
    return P;  
}

PARTITION *PartitionAddCommunity(PARTITION *P, COMMUNITY *C) {
#if VERBOSE > 2 
    printf("Add com %d, Before update, P->n = %d\n", C->id, P->n);
#endif
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
    return NULL;
}

PARTITION *PartitionDelCommunity(PARTITION *P, int c){
    COMMUNITY * C = P->C[c];
#if VERBOSE > 2 
    printf("Deleting com %d size %d, P->n = %d\n", c, P->C[c]->n, P->n);
    printf("C = %p\n", C);
#endif
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
    #if DEBUG
	printf("Freed Community %d\n", c);
    #endif
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
    SetFree(P->visited);
    SetFree(P->marked);
    Free(P->toMove);
    Free(P);
}

// oldCom is where to move the nodes back in case of reject
// newCom is where to get the nodes from
// moveDel is if MoveRandomNode deletes a community if it moves a node where it is the only node in the community
// oldComN and newComN are values of n before a move is considered
// oldComIn, oldComOut, newComIn, newComOut are the respective values before a move is considered
// oldScore and newScore serve the same purpose
// oldDiff and newDiff record the change in scores of the communities if the hypothetical move is done
static int _moveOption = -1, _oldCom = -1, _newCom = -1, _moveDel = 0, _oldComN, _newComN, _oldComIn, _oldComOut, _newComIn, _newComOut;
static double _oldScore, _newScore, _oldDiff = 0, _newDiff = 0;
// THESE MOVE OPTIONS ASSUME NO OVERLAPPING COMMUNITIES

static void SaveCommunityInfo(COMMUNITY * oldCom, COMMUNITY * newCom){
    _oldComN = oldCom->n; 
    _newComN = newCom->n;
    _newComIn = newCom->edgesIn;
    _newComOut = newCom->edgesOut;
    _oldComIn = oldCom->edgesIn; 
    _oldComOut = oldCom->edgesOut;
    _oldScore = oldCom->score;
    _newScore = newCom->score;
#if DEBUG
    printf("SAVE. Set_oldComN %d, _newComN %d\n", _oldComN, _newComN);
#endif
}

static void MoveRandomNode(PARTITION *P){  
    int u = P->G->n * drand48();
    int oldCom = P->whichCommunity[u];
    COMMUNITY * oc = P->C[oldCom]; 
    _oldCom = oldCom;
    int newCom;
    do{newCom = (int)(P->n * drand48());}
    while(newCom == P->whichCommunity[u]); 
    COMMUNITY * nc = P->C[newCom];	
    SetAdd(P->marked, u);
    P->toMove[P->numMoved++] = u;
    SaveCommunityInfo(oc, nc);    
    MoveOneNode(P, u, newCom);    
#if VERBOSE > 1
    printf("Mv(%d,%d->%d) ", u, oldCom, newCom);
#endif
    _newCom = newCom;

    if(oc->n - 1 == 0){
	// In case oc only had 1 node and we moved it away, delete oc
	_moveDel = 1;
    }
#if VERBOSE > 2
    printf("o %d, n %d, P->n %d\n", _oldCom, _newCom, P->n);
#endif
    assert(P->numMoved == 1);
}

// Merges C1 and C2 into just C1
static void MergeCommunities(PARTITION *P, int c1, int c2){   
    COMMUNITY *C1 = P->C[c1];
    COMMUNITY *C2 = P->C[c2];
#if VERBOSE > 1    
    printf("Mg(%d,%d) |%d,%d| ", c1, c2, C1->n, C2->n);
#endif 


    SaveCommunityInfo(C2, C1);
 
    int index = 0;
    for(int i = 0; i < C2->n; i++){
	int u = C2->nodeSet[i]; 
	//printf("Marked %d ", u);
	SetAdd(P->marked, u);
	P->toMove[P->numMoved++] = u;
    }
    int iters = C2->n;
    for(int i = 0; i < iters; ++i){
	MoveOneNode(P, C2->nodeSet[i], c1);
    }
    // Automatically takes care of P->n--;
    //PartitionDelCommunity(P, c2);
    if(c1 == P->n){
	// In the case C1 is the last community, to preseve the list of communities
	// C1 is placed where C2 was previosuly.
	_newCom = c2;
    }
    else{
	_newCom = c1;    
    }
    _oldCom = c2; 

   assert(P->numMoved == _oldComN); 
}

// Splits P->C[c_id] into a new community that has numNodes random nodes from the old community
static void SplitCommunity(PARTITION *P, int c_id, int numNodes){
    COMMUNITY *oldCom = P->C[c_id];
    assert(numNodes > 0);
    assert(numNodes < oldCom->n);
#if VERBOSE > 1
    printf("Sp(%d,|%d|->%d) ", c_id, oldCom->n, numNodes);
#endif
    COMMUNITY *newCom = CommunityAlloc(P->G, P->n);     
    _oldComN = oldCom->n; 
    SaveCommunityInfo(oldCom, newCom);
    PartitionAddCommunity(P, newCom);

    int u;
    for(int i = 0; i < numNodes; i++){	
	do{
	    int index = (int)(drand48() * oldCom->n);
	    u = oldCom->nodeSet[index];
	#if DEBUG
	    printf("u %d, i %d, located %d\n", u, index, P->whichCommunity[u]);
	#endif
	}
	while(SetIn(P->marked, u));
	SetAdd(P->marked, u);
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

*/

// fakeN is the number of nodes inside the community I TELL the measure

// Because of the nature of how I implemented faster rejects, most of the time 
// its not the actual true number of nodes in the community

double IntraEdgeDensity(COMMUNITY *C, int fakeN){
    if(C->n < 2)
	return 0;
    return C->edgesIn / (fakeN *(fakeN - 1) / 2.0);
}

double InterEdgeDensity(COMMUNITY *C, int fakeN){
    int tot = C->G->n - fakeN;
    //printf("tot = %d, edgesOut = %d\n", tot, edgesOut);
    return (double)C->edgesOut/(fakeN * tot);
}

double NewmanAndGirvan(COMMUNITY * C, int fakeN){
    // Eq 15 on Pg 16 on the pdf viewer
   
    // FIXME: Double check what exactly is needed
    return 0;
    //return C->inEdges/C->gDegree - (cDegree/(2.0*gDegree) * cDegree/(2.0*gDegree));
}

double Conductance(COMMUNITY * C, int fakeN){
     return (double)(C->edgesOut/(C->edgesOut + C->edgesIn));
}

double HayesScore(COMMUNITY *C, int fakeN){ 
#if DEBUG
    printf("Hayes Com %d, inEdges = %d, C->n = %d\n", C->id, C->edgesIn, fakeN);
#endif
    if(fakeN < 2){
	return 0;
    }
    double eps = C->edgesIn / ((fakeN * (fakeN-1))/2.0); 
    if(eps <= TARGET_EDGE_DENSITY){
    #if VERBOSE > 2
	printf("eps too low %f\n", eps);
    #endif
	return 0;
    }
    else{
    #if VERBOSE > 2
	printf("Result = %g\n", C->edgesIn * eps * TARGET_EDGE_DENSITY/eps);
    #endif
	return C->edgesIn*eps * (TARGET_EDGE_DENSITY/eps); // to down-weight if eps is above the target					  
    }
}

double ScorePartition(Boolean global, foint f){
    
    PARTITION *P = (PARTITION*) f.v;
#if VERBOSE > 0
    printf("SCORING PARTITION o = %d, n = %d, P->n = %d, P->numMoved = %d\n", _oldCom, _newCom, P->n, P->numMoved);
    printf("ncn %d, ocn %d\n", _newComN, _oldComN);
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
	int in, out;
	if(_oldCom != P->n){ 
	    COMMUNITY * old = P->C[_oldCom];
	    double oldBefore = old->score;
	    int testN = _oldComN - P->numMoved;
	    double oldScore = pCommunityScore(old, testN);
	    _oldDiff = oldScore - oldBefore;
	    P->total += _oldDiff;
#if VERBOSE > 1
	    printf("os = %g, ob = %g change = %g \n", oldScore, oldBefore, _oldDiff);
#endif 

	}
	COMMUNITY * new = P->C[_newCom];
	double newBefore = new->score;
	int testN = _newComN + P->numMoved; 
	double newScore = pCommunityScore(new, testN);
	_newDiff = newScore - newBefore;
	P->total += _newDiff;
#if VERBOSE > 1
	printf("ns = %g, nb = %g change = %g \n", newScore, newBefore, _newDiff);
#endif 
    }
#if VERBOSE > 0
    printf("Updated total = %g\n\n", P->total);
#endif
    return P->total;
}

// returns the CHANGE in score due to the move
double PerturbPartition(foint f) {
    //printf("Perturb\n");
    PARTITION *P = (PARTITION *) f.v;
    double before = P->total; 
#if MOVE_ONLY
    int choice = 0;
#else
    int choice = drand48() * 3;
#endif
    P->numMoved = 0;    

    SetReset(P->visited);
    SetReset(P->marked);
#if DEBUG
    printf("Before perturb\n");
    for(int i = 0; i < P->n; ++i){ 
    	COMMUNITY * temp = P->C[i];
	printf("Com %d, has %d nodes, score %g, in %d, out %d\n", i, temp->n, temp->score, temp->edgesIn, temp->edgesOut);
	PrintCommunity(temp);
    }
    printf("P->total score = %g\n", P->total);
    int fail = 1;
    for(int i = 0; i < P->G->n; ++i){
    #if DEBUG
	printf("Node %d in Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
    #endif
	if(P->C[P->whichCommunity[i]]->nodeSet[P->whichMember[i]] != i){
	    printf("ERROR: Mem mismatch. Node %d, Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
	    fail = 0;
	}	
	if(P->C[P->whichCommunity[i]]->id != P->whichCommunity[i]){
	    printf("ERROR: Com mismatch. Node %d, Com %d index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
	    fail = 0;
	}
    }
    assert(fail);

    int comNums = 1;
    for(int i = 0; i < P->G->n; ++i){	    
	if(P->whichCommunity[i] >= P->n){
	    comNums = 0;
	    printf("ERROR: Node %d is in Com %d > %d\n", i, P->whichCommunity[i], P->n);
	}	    
    }
    assert(comNums);
#endif
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
	if(P->C[c1]->n > P->C[c2]->n)
	    MergeCommunities(P, c1, c2);
	else
	    MergeCommunities(P, c2, c1);
	_moveOption = 1;
    }
    else{ 
	#define MAX_TRIES 20	
	int c, numNodes, tries = 0;
	do
	    c = drand48() * P->n;
	while(P->C[c]->n == 1 && tries++ < MAX_TRIES);

	if(tries >= MAX_TRIES){
	    // If after a while, it can't find a community with n > 1,
	    // just move
	    MoveRandomNode(P);
	    _moveOption = 0;
	}
	else{
	    numNodes = (int)(drand48() * (P->C[c]->n / 2)) + 1;	
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
    
    COMMUNITY * newCom = NULL;
    COMMUNITY * oldCom = NULL;

    if(P->C[_newCom])
	newCom = P->C[_newCom];
    if(P->C[_oldCom])
	oldCom = P->C[_oldCom];
     
#if DEBUG
    printf("P->n = %d ", P->n); 
    if(newCom)
	printf("_newCom = %d size %d ", newCom->id, newCom->n);
    else
	printf("newCom DNE ");
    if(oldCom)
	printf("_oldCom = %d size %d ", oldCom->id, oldCom->n);
    else
	printf("oldCom DNE ");
    printf("numMoved = %d\n", P->numMoved);
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * temp = P->C[i];
	printf("Com %d has %d nodes, %d in, %d out, score %g\n", i, temp->n, temp->edgesIn, temp->edgesOut, temp->score);
    }

    printf("_newComN %d, _oldComN %d, oc->in %d, oc->out %d, nc->in %d, nc->out %d\n", _newComN, _oldComN, _oldComIn, _oldComOut, _newComIn, _newComOut);

#endif
    
    if(accept){ 
	// Check if both they exist
	// Actually move the nodes if its accepted
    	oldCom->score += _oldDiff;
	newCom->score += _newDiff;
		
	for(int i = 0; i < P->numMoved; ++i){
	    CommunityDelNode(oldCom, P, P->toMove[i]);
	    CommunityAddNode(newCom, P, P->toMove[i]);
	}

	if(_moveDel || _moveOption == 1){
	    PartitionDelCommunity(P, _oldCom);
	    oldCom = NULL;
	} 
	
    }
    else {


	// In the case of reject, update the information to the old information
#if DEBUG
	printf("REJECTION\n");
	printf("Before P->total %g, Score old %g, new %g\n", P->total, oldCom->score, newCom->score);	
#endif	
	oldCom->edgesIn = _oldComIn;
	oldCom->edgesOut = _oldComOut;
	newCom->edgesIn = _newComIn;
	newCom->edgesOut = _newComOut;
	oldCom->score = _oldScore;
	newCom->score = _newScore; 
	P->total -= _oldDiff;
	P->total -= _newDiff;
#if DEBUG
	printf("After P->total %g, Score old %g, new %g\n", P->total, oldCom->score, newCom->score);
#endif
	if(_moveOption == 2){
	    PartitionDelCommunity(P, P->n - 1);
	    newCom = NULL;
	}
    }
    
    if(oldCom)
	_oldComN = oldCom->n;
    else
	_oldComN = 0;
    if(newCom)
	_newComN = newCom->n;
    else
	_newComN = 0;	
#if DEBUG
    printf("Now holding true values, _oldComN %d, _newComN %d\n", _oldComN, _newComN);
#endif


    P->numMoved = 0;
#if DEBUG
    printf("After A/R\n");
    for(int i = 0; i < P->n; ++i){ 
    	COMMUNITY * temp = P->C[i];
	printf("Com %d, has %d nodes, score %g, in %d, out %d\n", i, temp->n, temp->score, temp->edgesIn, temp->edgesOut);
	PrintCommunity(temp);
    }

    //for(int i = 0; i < P->G->n; ++i){
//	printf("Node %d in Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
  //  }
  
    int fail = 1; 
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * C = P->C[i];
	assert(C->n > 0);
	int groundOut = CommunityEdgeOutwards(P, C);
	int groundIn = CommunityEdgeCount(C);  
	if(groundOut != C->edgesOut){
	    printf("ERROR: Com %d GroundOut %d != stored %d\n", i, groundOut, C->edgesOut);
	    fail = 0;
	}
	if(groundIn != C->edgesIn){
	    printf("ERROR: Com %d GroundIn %d != stored %d\n", i, groundIn, C->edgesIn);
	    fail = 0;
	}
    }

    assert(fail);
#endif

    _moveDel = 0;
#if VERBOSE > 0   
    printf("Current total = %f\n\n", P->total);
#endif
    _oldCom = -1; // This ensures that I don't score the partition with the wrong information
    _newCom = -1; // The score is properly updated on MaybeAccept
    
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

void SAR(int iters, foint f){

    PARTITION * P = f.v;
    int fail = 1, in, out, best_id = -1;
    double ground = 0, withInfo, stored = 0, best = -1, totalGround = 0, totalStored = 0;
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
#if DEBUG
	ground = pCommunityScore(com, com->n);
	if(fabs(ground - com->score) > 0.001){
	    printf("Com %d ERROR: Ground score %g, stored score %g\n", i, ground, com->score);
	    fail = 0;
	}
	totalStored += com->score;
	totalGround += ground;
	#endif
    }
#if 1
    printf("P->Total %g, totalStored %g, totalGround %g\n", P->total, totalStored, totalGround);
#endif
#if DEBUG
    assert(fabs(P->total - totalStored) < 0.001);
    assert(fabs(P->total - totalGround) < 0.001);
    assert(P->total >= 0);

#endif
    assert(fail);
    if(fail == 0){
	printf("fail\n");
	exit(1);
    }
	
    printf("\tBest: Com %d, with n %d, score %g  \tP->n = %d, Total score = %g", best_id, P->C[best_id]->n, best, P->n, P->total);
}

// Assumes the P is already properly allocated
// Directly modifies P
void PartitionRead(FILE * fp, PARTITION * P){
    char line[BUFSIZ];
    int numCom = 0;
    SET * overlapCheck = SetAlloc(P->G->n);
    SET * checkAll = SetAlloc(P->G->n);
    while(fgets(line, sizeof(line), fp)){
	line[strcspn(line, "\n")] = '\0';	
	int len = strlen(line);	
	char *token = strtok(line, " ");
	COMMUNITY * C = CommunityAlloc(P->G, numCom++);
	
	while(token != NULL){
	    int node = atof(token);
	    if(SetIn(overlapCheck, node)){
		printf("ERROR: %d node is in an overlapping community\nSorry, haven't yet implemented that yet\n", node);
		exit(1);
	    }
	    SetAdd(overlapCheck, node);
	    CommunityAddNode(C, P, node);
	    token = strtok(NULL, " ");
	    SetAdd(checkAll, node);
	}
	PartitionAddCommunity(P, C);
    } 
    

    COMMUNITY * Extra = CommunityAlloc(P->G, numCom++);
    for(int i = 0; i < P->G->n; ++i){
	if(!SetIn(checkAll, i)){
	    printf("WARNING: Node %d is not in any community. Putting it into catch all\n", i);
	    CommunityAddNode(Extra, P, i);
	}
    }

    if(Extra->n == 0){
	CommunityFree(Extra);	
	printf("All nodes accounted for, deleting catch all community\n");
    }
    else{
	PartitionAddCommunity(P, Extra);
	printf("Catch all community added into partition\n");
    }
    
    
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * C = P->C[i];
	C->edgesIn = CommunityEdgeCount(C);
	C->edgesOut = CommunityEdgeOutwards(P, C);
	double s = pCommunityScore(C, C->n); 
	C->score = s;
	P->total += s;     
    }
    fclose(fp);
}


#define RANDOM_START 0
#define CHECK_OVERLAP 1
int main(int argc, char *argv[])
{
    printf("Running with PARANOID_ASSERTS=%d, NDEBUG=%d\n", PARANOID_ASSERTS, NDEBUG);
    // Set which measure to use here
    pCommunityScore = HayesScore;
    
    
    int i, j;
    srand48(GetFancySeed(false));

    Boolean sparse=maybe, supportNames = true;
    FILE *fp = Fopen(argv[1], "r"); // edge list file is given on the command line  

    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNames, false);
    
    printf("G has %d nodes, %d edges\n", G->n, G->numEdges);

    PARTITION *P = PartitionAlloc(G);

    if(argc > 2){
	printf("Reading partition %s\n", argv[2]);
	PartitionRead(Fopen(argv[2], "r"), P);
    }
    else{
#if RANDOM_START
    int numCommunities = 2; // communities numbered 0 through numCommunities-1 inclusive
    printf("Starting with %d random communities\n", numCommunities);
    for(i=0; i<numCommunities; i++) PartitionAddCommunity(P, CommunityAlloc(G, i));

    for(i=0; i<G->n; i++) {
	int which = (int)(drand48() * numCommunities);
	printf("%d->%d, ", i, which);
	CommunityAddNode(P->C[which], P, i);
    }

#else

    printf("BFS-based communities: \n");
    SET *nodesUsed=SetAlloc(G->n); // cumulative set of nodes that have been put into a partition

    int nodeArray[G->n], distArray[G->n], numCom = 0;
    while(SetCardinality(nodesUsed) < G->n) {
	int seed;
	do { seed = (int)(drand48() * G->n); }
	while(SetIn(nodesUsed, seed));
	int numAdded=0, distance = 4; // should be far enough
	int n=GraphBFS(G, seed, distance, nodeArray, distArray); // list of nodes within "distance" of seed
	//printf("BFS(%d[%d])=%d", seed, G->degree[seed], n);
	assert(n>0 && nodeArray[0]==seed && distArray[seed]==0);
	COMMUNITY *C = CommunityAlloc(G, numCom);
	for(i=0; i<n; i++) if(!SetIn(nodesUsed,nodeArray[i]) && drand48() > 0.5) {
	    SetAdd(nodesUsed,nodeArray[i]); CommunityAddNode(C, P, nodeArray[i]); ++numAdded;
	}
	//printf("%d ", numAdded); fflush(stdout);
	//assert(C->n >0 && C->n < G->n && C->n==numAdded);
	if(C->n > 0){
	    PartitionAddCommunity(P, C);
	    ++numCom;
	}
	else
	    CommunityFree(C);
	//printf("Size of Community = %d\n", C->n);
    }

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
	C->edgesIn = CommunityEdgeCount(C);
	C->edgesOut = CommunityEdgeOutwards(P, C);
	double s = pCommunityScore(C, C->n); 
	C->score = s;
	
	P->total += s; 
    }

    }
    fclose(fp);
	
#if 0
    HillClimbing(P, 200);
#else
    foint f;
    f.v = P;
	
    SIM_ANNEAL *sa = SimAnnealAlloc(1, f, PerturbPartition, ScorePartition, MaybeAcceptPerturb, 10*G->n*G->numEdges /*100*/,0,0,SAR);
    if(G->n==2390 && G->numEdges==16127) {
	printf("Hmm, this looks like yeast.el, using canned schedule\n");
	SimAnnealSetSchedule(sa, 6, 3);
    }
    else
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
	if(HayesScore(P->C[j], P->C[j]->n) && num > biggest) {
	    biggest = num;
	    which = j;
	}
	
    }
    assert(nodes == G->n);
    printf("Final score = %f\n", P->total);
    if(which == -1)
	which = 0;
    printf("Biggest community is #%d, with %d nodes:\n", which, biggest);
    PrintCommunity(P->C[which]);

#if PRINT_ALL_COMMUNITIES
    printf("\n Start ALL COMMUNITIES\n");
    for(int i = 0; i < P->n; ++i){
	PrintCommunity(P->C[i]);
    }
#endif
    printf("Attempting Partition Free\n");
    PartitionFree(P);
    printf("Partition Free completed\n");
    GraphFree(G);
    return 0;
}
