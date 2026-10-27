// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.


// Warning: setting PA = 0 and ND = 1 turns of ALL checking... but could make this up to 10x faster
#define PARANOID_ASSERTS 0
#define NDEBUG 1
#include "communities.h"
#include "structs.h"
#include "measures.h"
#include <stdio.h>
#include "misc.h"
#include "rand48.h"
#include "graph.h"
#include "sets.h"
#include "sim_anneal.h"
//#include "pthread.h"
#include "scoreopts.c"
#include "unistd.h"


#define VERBOSE 0 // 0 = no noisy outpt, 3 = lots, 1..2 is intermediate
#define DEBUG 0
// Which perturbations PerturbPartition uses (merge and split are planned to be removed):
#define PERTURB_ALL       0 // move, merge and split (the old MOVE_ONLY 0)
#define PERTURB_MOVE_ONLY 1 // only move (the old MOVE_ONLY 1)
#define PERTURB_MOVE_SWAP 2 // move or swap, 50/50
#define PERTURB_MODE PERTURB_MOVE_ONLY
#define PRINT_ALL_COMMUNITIES 1 // Prints out all the community info at the end of the program execution
#define PRINT_ITERS 1 // Prints out iters per second on terminal (Note: Is quite expensive)
#define PRINT_ALL_SCORES 0 // Prints out the scores of the partition based on all the functions I've implemented 

/************************** Community routines *******************/

double (*pCommunityScore)(PARTITION * P, COMMUNITY * C, int fakeN) = NULL;  
//static int num_threads; 

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

// All move operations are based off of this function
// Repeadetely call this to handle all edgesIn + out
static void MoveOneNode(PARTITION * P, int node, int dest){  
    COMMUNITY * oldCom = P->C[P->whichCommunity[node]];
    COMMUNITY * newCom = P->C[dest];
    GRAPH * G = P->G;

    int actualOldIn = CommunityEdgeCount(oldCom);
    int actualOldOut = CommunityEdgeOutwards(P, oldCom);

    int actualNewIn = CommunityEdgeCount(newCom);
    int actualNewOut = CommunityEdgeOutwards(P, newCom);

    int fail = 0; 
    if(oldCom->edgesIn != actualOldIn){
	printf("ERROR: oldCom %d IN stored %d vs actual %d\n", oldCom->id, oldCom->edgesIn, actualOldIn);
	fail = 1;
    }
    if(oldCom->edgesOut != actualOldOut){
	printf("ERROR: oldCom %d OUT stored %d vs actual %d\n", oldCom->id, oldCom->edgesOut, actualOldOut);
	fail = 1;
    }

    if(newCom->edgesIn != actualNewIn){
	printf("ERROR: newCom %d IN stored %d vs actual %d\n", newCom->id, newCom->edgesIn, actualNewIn);
	fail = 1;
    }
    if(newCom->edgesOut != actualNewOut){
	printf("ERROR: newCom %d OUT stored %d vs actual %d\n", newCom->id, newCom->edgesOut, actualNewOut);
	fail = 1;
    }

    
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
    printf("old = %d, new = %d, diff = %d\n", old, new, old-new);
    printf("AFTER\noc size = %d, nc size = %d\noc in %d, oc out %d, nc in %d, nc out %d\n", oldCom->n, newCom->n, oldCom->edgesIn, oldCom->edgesOut, newCom->edgesIn, newCom->edgesOut);
#endif

    if(fail){
	exit(1);
    }

}


int NodeInDegree(PARTITION *P, COMMUNITY *C, int node){
    int i, inDeg = 0;
        //printf("START: NodeInDegree for %d in Com %d\n", node, P->whichCommunity[node]);
    for(i=0;i<C->G->degree[node];i++){
	int neigh = C->G->neighbor[node][i];

	//printf("Neighbor %d in com %d\n",neigh, P->whichCommunity[neigh]);
	int belongsTo = P->whichCommunity[neigh]; 
	if(belongsTo == C->id){
	    ++inDeg;
	}
    }

    //printf("END: Node %d in Com %d inDeg %d\n", node, P->whichCommunity[node], inDeg);

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
//  P->common = SetAlloc(G->n);
P->visited = SetAlloc(G->n); // A set for community update to avoid double counting  
P->marked = SetAlloc(G->n); // Marking which ones will be moved (Essentially SET * toMove but in int * format)  
P->toMove = Calloc(sizeof(int), G->n);    
P->numMoved = 0; 
P->iters = 0;
P->clk = clock();
return P;  
}

PARTITION *PartitionAddCommunity(PARTITION *P, COMMUNITY *C) {
#if VERBOSE > 2 
printf("Add com %d, Before update, P->n = %d\n", C->id, P->n);
#endif
assert(P->n < P->G->n);
int i = P->n;
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
else{
    printf("MEM ERROR: Attempting to write to an existing community???\n");
    exit(1);
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
// For a swap, oldCom is the community of the first node picked and newCom is the community of the second;
// each node goes to the other's community, so neither community changes size (and none is ever deleted)
// _moveOption is 0 = move, 1 = merge, 2 = split, MOVE_OPTION_SWAP = swap (PerturbPartition's "choice" uses the same numbers)
#define MOVE_OPTION_SWAP 3
#define MAX_SWAP_TRIES 1000 // attempts at finding two nodes in different communities before giving up on a swap
int _moveOption = -1, _oldCom = -1, _newCom = -1, _moveDel = 0, _oldComN, _newComN, _oldComIn, _oldComOut, _newComIn, _newComOut;
double _oldScore, _newScore, _oldDiff = 0, _newDiff = 0;

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

// Choose one random node from any community, move it to a new community
static void MoveRandomNode(PARTITION *P){  
    int u = P->G->n * drand48();
    int oldCom = P->whichCommunity[u];
    COMMUNITY * oc = P->C[oldCom]; 
    _oldCom = oldCom;
    int newCom;
    // Note: MoveRandomNode is not called if there is only 1 community (Split is called instead)
    // This function assumes this and does not check
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

// Swap counterpart of MoveOneNode: update edgesIn/edgesOut of u's community A and v's community B
// as if u and v traded places (u goes to B, v goes to A). Like MoveOneNode, this only touches the
// edge counts; nodeSet/whichCommunity/whichMember are only changed in MaybeAcceptPerturb on accept.
// MoveOneNode can't just be called twice: the check at its top compares the stored counts with counts
// recomputed from nodeSet, and the first call has already changed the stored ones. Its "marked"
// shortcut also assumes all marked nodes go to the same community, but u and v go opposite ways.
//
// Let a_x and b_x be the number of neighbors x has in A and in B, not counting u, v or self loops. Then
//   edgesIn:  A += a_v - a_u,  B += b_u - b_v
// and since edgesOut = (sum of degrees) - 2*edgesIn, with degrees not counting self loops,
//   edgesOut: A += deg(v) - deg(u) + 2*(a_u - a_v),  B += deg(u) - deg(v) + 2*(b_v - b_u)
// An edge between u and v crosses between A and B both before and after, so it changes nothing.
static void SwapTwoNodes(PARTITION *P, int u, int v){
    GRAPH *G = P->G;
    COMMUNITY *A = P->C[P->whichCommunity[u]], *B = P->C[P->whichCommunity[v]];
    int a_u = 0, b_u = 0, a_v = 0, b_v = 0, degU = 0, degV = 0;
    for(int i = 0; i < G->degree[u]; ++i){
	int w = G->neighbor[u][i];
	if(w == u) continue; // self loop
	++degU;
	if(w == v) continue; // the u--v edge, see above
	if(P->whichCommunity[w] == A->id) ++a_u;
	else if(P->whichCommunity[w] == B->id) ++b_u;
    }
    for(int i = 0; i < G->degree[v]; ++i){
	int w = G->neighbor[v][i];
	if(w == v) continue; // self loop
	++degV;
	if(w == u) continue; // the u--v edge, see above
	if(P->whichCommunity[w] == A->id) ++a_v;
	else if(P->whichCommunity[w] == B->id) ++b_v;
    }
    A->edgesIn += a_v - a_u;
    B->edgesIn += b_u - b_v;
    A->edgesOut += degV - degU + 2*(a_u - a_v);
    B->edgesOut += degU - degV + 2*(b_v - b_u);
}

// Choose two random nodes that are in different communities and swap their communities
// Returns 0, having changed nothing, if it can't find such a pair, so the caller can do something else
// Note: Like MoveRandomNode, this assumes there is more than 1 community
static int SwapRandomNodes(PARTITION *P){
    int u, v, tries = 0;
    do{
	if(tries++ >= MAX_SWAP_TRIES)
	    return 0;
	u = P->G->n * drand48();
	v = P->G->n * drand48();
    }
    while(P->whichCommunity[u] == P->whichCommunity[v]);
    COMMUNITY * oc = P->C[P->whichCommunity[u]];
    COMMUNITY * nc = P->C[P->whichCommunity[v]];
    _oldCom = oc->id;
    _newCom = nc->id;
    P->toMove[P->numMoved++] = u;
    P->toMove[P->numMoved++] = v;
    SaveCommunityInfo(oc, nc);
    SwapTwoNodes(P, u, v);
#if VERBOSE > 1
    printf("Sw(%d,%d: %d<->%d) ", u, v, _oldCom, _newCom);
#endif
    assert(P->numMoved == 2);
    return 1;
}

// Merges C1 and C2 into just C1
static void MergeCommunities(PARTITION *P, int c1, int c2){   
    COMMUNITY *C1 = P->C[c1];
    COMMUNITY *C2 = P->C[c2];
#if VERBOSE > 1    
    printf("Mg(%d,%d) |%d,%d| ", c1, c2, C1->n, C2->n);
#endif 
    SaveCommunityInfo(C2, C1);
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
	    // A swap sends one node each way, so neither community changes size
	    int sizeChange = (_moveOption == MOVE_OPTION_SWAP) ? 0 : P->numMoved;
	    // In case of a split/merge reject P->C[n] will no longer exist
	    if(_newCom == P->n){
		int swap = _newCom;
		_newCom = _oldCom;
		_oldCom = swap;
	    }
	    if(_oldCom != P->n){
		// In case oldCom exists, record the change for oldCom
		// If it doesn't exist, the only changes of the score will be reflected in newCom
		COMMUNITY * old = P->C[_oldCom];
		double oldBefore = old->score;
		int testN = _oldComN - sizeChange;
		double oldScore = pCommunityScore(P, old, testN);
		_oldDiff = oldScore - oldBefore;
		P->total += _oldDiff;
#if VERBOSE > 1
		printf("os = %g, ob = %g change = %g \n", oldScore, oldBefore, _oldDiff);
#endif 

	    }
	    COMMUNITY * new = P->C[_newCom];
	    double newBefore = new->score;
	    int testN = _newComN + sizeChange;
	    double newScore = pCommunityScore(P, new, testN);
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
    PARTITION *P = (PARTITION *) f.v;
    double before = P->total; 
    P->numMoved = 0;    

#if PRINT_ITERS
    {
    static iter;
    double elapsed = (double)(clock() - P->clk) / CLOCKS_PER_SEC;
    double iters_per_sec = P->iters++ / elapsed; 
    if(iter++ % 1000==0) printf("\r%-0f\t",iters_per_sec);
    }
#endif
    
#if DEBUG
    printf("Before perturb\n");
    int fail = 1;
    for(int i = 0; i < P->n; ++i){ 
    	COMMUNITY * temp = P->C[i];
	printf("Com %d, has %d nodes, score %g, in %d, out %d\n", i, temp->n, temp->score, temp->edgesIn, temp->edgesOut);
	PrintCommunity(temp);
	if(temp->edgesIn < 0 || temp->edgesOut < 0){
	    printf("ERROR: Negative edges detected!\n");
	    fail = 0;
 	}
    }
    printf("P->total score = %g\n", P->total);
    for(int i = 0; i < P->G->n; ++i){
	// Below line will print out a lot of garbage, only use for extreme debugging
	//printf("Node %d in Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
	if(P->C[P->whichCommunity[i]]->nodeSet[P->whichMember[i]] != i){
	    printf("ERROR: Mem mismatch. Node %d, Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
	    fail = 0;
	}	
	if(P->C[P->whichCommunity[i]]->id != P->whichCommunity[i]){
	    printf("ERROR: Com mismatch. Node %d, Com %d index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
	    fail = 0;
	}
    }
    
    if(!fail){
	exit(1);
    }

    int comNums = 1;
    for(int i = 0; i < P->G->n; ++i){	    
	if(P->whichCommunity[i] >= P->n){
	    comNums = 0;
	    printf("ERROR: Node %d is in Com %d > %d\n", i, P->whichCommunity[i], P->n);
	}	    
    }
    assert(comNums);
#endif

    #if PERTURB_MODE == PERTURB_MOVE_ONLY
	int choice = 0;
    #elif PERTURB_MODE == PERTURB_MOVE_SWAP
	int choice = (drand48() < 0.5) ? 0 : MOVE_OPTION_SWAP; // move or swap, 50/50
    #else
	int choice = drand48() * 3;
    #endif
    
    if(choice == 0 && P->n > 1){
	// Move one node only
	MoveRandomNode(P);
	_moveOption = 0;	
    }
    else if(choice == 1 && P->n > 1){
	// Merge two communities
	int c1, c2;
	do{
	    c1 = (int)(drand48() * P->n);
	    c2 = (int)(drand48() * P->n);
	}
	while(c1 == c2);
	// Merge smaller community into bigger community
	if(P->C[c1]->n > P->C[c2]->n)
	    MergeCommunities(P, c1, c2);
	else
	    MergeCommunities(P, c2, c1);
	_moveOption = 1;
    }
    else if(choice == MOVE_OPTION_SWAP && P->n > 1){
	// Swap two nodes that are in different communities
	if(SwapRandomNodes(P))
	    _moveOption = MOVE_OPTION_SWAP;
	else{
	    // Couldn't find two such nodes, just move
	    MoveRandomNode(P);
	    _moveOption = 0;
	}
    }
    else{ 
	// Split a community into 2
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
    SetEmpty(P->visited);
    SetEmpty(P->marked);


    return after-before;
}

Boolean MaybeAcceptPerturb(Boolean accept, foint f) {
#if VERBOSE > 0    
    printf("Accept? %d\n", accept);
#endif

    PARTITION * P = (PARTITION *) f.v;
        
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
	// Actually move the nodes if its accepted
    	oldCom->score += _oldDiff;
	newCom->score += _newDiff;
		
	if(_moveOption == MOVE_OPTION_SWAP){
	    // Neither community changes size, so just exchange the two nodes' places
	    int u = P->toMove[0], v = P->toMove[1];
	    int iu = P->whichMember[u], iv = P->whichMember[v];
	    oldCom->nodeSet[iu] = v;
	    newCom->nodeSet[iv] = u;
	    P->whichMember[u] = iv;
	    P->whichMember[v] = iu;
	    P->whichCommunity[u] = newCom->id;
	    P->whichCommunity[v] = oldCom->id;
	}
	else{
	    for(int i = 0; i < P->numMoved; ++i){
		CommunityDelNode(oldCom, P, P->toMove[i]);
		CommunityAddNode(newCom, P, P->toMove[i]);
	    }
	}

	// If we merged or moved from a community of size 1, delete it
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
    _oldCom = -1; // This ensures that I don't score the partition with the wrong information
    _newCom = -1; // The score is properly updated on MaybeAccept

    
    P->numMoved = 0;
#if DEBUG
    printf("After A/R\n");
    for(int i = 0; i < P->n; ++i){ 
    	COMMUNITY * temp = P->C[i];
	printf("Com %d, has %d nodes, score %g, in %d, out %d\n", i, temp->n, temp->score, temp->edgesIn, temp->edgesOut);
	PrintCommunity(temp);
    }
    
    //For extreme debugging
    /*for(int i = 0; i < P->G->n; ++i){
    	printf("Node %d in Com %d, index %d\n", i, P->whichCommunity[i], P->whichMember[i]);
    }*/

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
#if DEBUG
    SAR(0, f);
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

// Short for Simmulated Annealing Report
// Called every 1% to report best community
// Also makes sure each of the communities hold correct information
void SAR(int iters, foint f){

    PARTITION * P = f.v;
    int fail = 1, in, out, best_id = -1;
    double ground = 0, stored = 0, best = -999999, totalGround = 0, totalStored = 0;
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
	    printf("Com %d edgesIn ERROR: ground %d vs stored %d\n", i, in, com->edgesIn);
	}
	if(out != com->edgesOut){
	    fail = 0;
	    printf("Com %d edgesOut ERROR: ground %d vs stored %d\n", i, out, com->edgesOut);
	}
#if 1
	ground = pCommunityScore(P, com, com->n);
	if(fabs(ground - com->score) > 0.001){
	    printf("Com %d ERROR: Ground score %g, stored score %g\n", i, ground, com->score);
	    fail = 0;
	}
	totalStored += com->score;
	totalGround += ground;
#endif
    }
#if 0
    printf("P->Total %g, totalStored %g, totalGround %g\n", P->total, totalStored, totalGround);
#endif
#if 1
    assert(fabs(P->total - totalStored) < 0.001);
    assert(fabs(P->total - totalGround) < 0.001);
    assert(P->total >= 0);

#endif
    assert(fail);
    if(fail == 0){
	printf("fail\n");
	exit(1);
    }
	
    printf("\tP->n = %d Best: Com %d, with n %d, score %g  \tTotal score = %g", P->n, best_id, P->C[best_id]->n, best, P->total);
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
		
	char *token = strtok(line, " ");
	COMMUNITY * C = CommunityAlloc(P->G, numCom++);
	
	while(token != NULL){
	    int node = atof(token);
	    if(SetIn(overlapCheck, node)){
		printf("ERROR: %d node is in an overlapping community\nSorry, haven't implemented that yet\n", node);
		exit(1);
	    }
	    SetAdd(overlapCheck, node);
	    CommunityAddNode(C, P, node);
	    token = strtok(NULL, " ");
	    SetAdd(checkAll, node);
	}
	PartitionAddCommunity(P, C);
    } 
    // Extra community to make sure that all nodes are accounted for
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
    
    if(!pCommunityScore){
	printf("Note: Default objective function is Hayes\n");
    }
    for(int i = 0; i < P->n; ++i){
	COMMUNITY * C = P->C[i];
	C->edgesIn = CommunityEdgeCount(C);
	C->edgesOut = CommunityEdgeOutwards(P, C);
	double s = pCommunityScore(P, C, C->n); 
	C->score = s;
	P->total += s;     
    }
    fclose(fp);
}

void PrintAllScores(PARTITION * P){
    double (*ObjectiveFuncs[])(PARTITION * P, COMMUNITY *, int) = {HayesScore, IntraEdgeDensity, InterEdgeDensity, NewmanAndGirvan, Conductance};
    const char * FuncNames[] = {"HayesScore", "IntraEdgeDensity", "InterEdgeDensity", "NewmanAndGirvan", "Conductance"};
    for(int i = 0; i < 5; ++i){
	double score = 0; 
	for(int j = 0; j < P->n; ++j){
	    COMMUNITY * C = P->C[j];
	    score += ObjectiveFuncs[i](P, C, C->n);
	}
	if(i != 3){
	    printf("Average %s, score = %g\n", FuncNames[i] , score / P->n);
	}
	else{
	    printf("%s, score = %g\n", FuncNames[i], score);
	}
    }
}

#define RANDOM_START 0
#define CHECK_OVERLAP 1
int main(int argc, char *argv[])
{
    //num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    printf("Running with PARANOID_ASSERTS=%d, NDEBUG=%d\n", PARANOID_ASSERTS, NDEBUG);
    
    pCommunityScore = HayesScore;
            
    int i;
    srand48(GetFancySeed(false));

    Boolean sparse=maybe, supportNames = true;
    FILE *fp = Fopen(argv[1], "r"); // edge list file is given on the command line  

    // we don't care about directed here, graph is assumed undirected
    GRAPH *G = GraphReadEdgeList(NULL, fp, false, supportNames, false);
    
    printf("G has %d nodes, %d edges\n", G->n, G->numEdges);

    PARTITION *P = PartitionAlloc(G);
    
    int opt, parRead = 0, dir = 1, eval = 0; 
    while((opt = getopt(argc, argv, "p:s:e:")) != -1){
	switch(opt){
	    case('p'):
		// -p means user is giving partition to start from
		// start hillclimbing immediately
		printf("Reading partition %s\n", argv[2]);
		PartitionRead(Fopen(optarg, "r"), P);
		parRead = 1;
		break;
	    case('s'):
		// -s means we are specifying a objective/scoring function
		printf("Setting scoring function to %s\n", optarg);
		for(int i = 0; i < numScoreOpts; i++){
		    // Hopefully find a matching scoring function name
		    if(strstr(optarg, scoreOpts[i].name) != NULL){
			pCommunityScore = scoreOpts[i].func;
			dir = scoreOpts[i].optimizingForMax;
			printf("NOTE: Using score %s, optimizing dir = %d\n", scoreOpts[i].name, dir);
			break;
		    }
		}
		if(pCommunityScore == NULL){
		    printf("ERROR: Provided scoring function name: %s does not exist.\n", optarg);
		    exit(1);
		}
		break; 
	    case('e'):
		// -e means eval mode 
		// Read input file similar to partition read option but eval scores
		PartitionRead(Fopen(optarg, "r"), P);
		eval = 1;
		break;
	}
    } 
    if(pCommunityScore == NULL){
	printf("NOTE: Using default HayesScore\n");
	pCommunityScore = HayesScore;
	dir = 1;
    }
    

    // If the user does not provide a Partition, create a new one
    if(!parRead && !eval){
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
	printf("Starting num coms = %d\n", P->n);

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
		
	// Initiate scores for all communities
	for(int i = 0; i < P->n; ++i){
	    COMMUNITY * C = P->C[i];
	    C->edgesIn = CommunityEdgeCount(C);
	    C->edgesOut = CommunityEdgeOutwards(P, C);
	    double s = pCommunityScore(P, C, C->n); 
	    //printf("Score = %f, edgesIn %d, edgesOut %d\n", s, C->edgesIn, C->edgesOut);
	    C->score = s;
	    P->total += s; 
	}

    }
    fclose(fp);
	
#if 0
    HillClimbing(P, 200);
#else
    if(!eval){
	foint f;
	f.v = P;
	
	SIM_ANNEAL *sa = SimAnnealAlloc(dir, f, PerturbPartition, ScorePartition, MaybeAcceptPerturb, 200*G->numEdges /*100*/, 0,0,SAR);
	if(parRead == 0){
	    // If user didn't provide partition run actual sim anneal
	    if(G->n==2390 && G->numEdges==16127) {
		printf("Hmm, this looks like yeast.el/communities.in, using canned schedule\n");
		SimAnnealSetSchedule(sa, 3, 8);
	    }
	    else
		SimAnnealAutoSchedule(sa); // to automatically create schedule
	}
	else{
	    // hill climbing
	    // TODO: Test if partition read actually hill climbs
	    SimAnnealSetSchedule(sa, 0, 0);
	}
	    
	SimAnnealRun(sa); // returns >0 if success, 0 if not done and can continue, <0 if error
	SimAnnealFree(sa);

    }
#endif

#if 0

    // Prints the biggest community, not really useful

    int nodes = 0, biggest = 0, num, which=-1;
    for(int j = 0; j < P->n; ++j){
	num = P->C[j]->n;
	//int inEdges = CommunityEdgeCount(P->C[j]);
	//printf("Com %d has %d nodes, %d edges, with edge density %g\n", j, num, inEdges, inEdges/(num*(num-1)/2.0));
	nodes += num;
	if(HayesScore(P, P->C[j], P->C[j]->n) && num > biggest) {
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
#endif

#if PRINT_ALL_COMMUNITIES
    printf("\n Start ALL COMMUNITIES\n");
    for(int i = 0; i < P->n; ++i){
	PrintCommunity(P->C[i]);
	//printf("\n");
    }
#endif
    printf("Total time elapsed = %.2f seconds\n", (double)(clock() - P->clk)/ CLOCKS_PER_SEC); 
#if PRINT_ALL_SCORES
    PrintAllScores(P);
#endif
    PartitionFree(P);
    GraphFree(G);
    return 0;
}
