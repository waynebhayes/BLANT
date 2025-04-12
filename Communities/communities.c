#include <stdio.h>
#include "misc.h"
#include "rand48.h"
#include "graph.h"
#include "sets.h"
#include "sim_anneal.h"

/************************** Community routines *******************/
typedef struct _community {
    int id, n;
    SET *nodeSet;
    GRAPH *G; // the graph we came from
    double score; 
} COMMUNITY;

COMMUNITY *CommunityAlloc(GRAPH *G){
    COMMUNITY *C = Calloc(sizeof(COMMUNITY), 1);
    C->id = 0; // not sure if this will be handy or not...
    C->score = 0;
    C->n = 0;
    C->G = G;
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
    else
	printf("Already in\n");
    return C;
}

int NodeInDegree(COMMUNITY *C, int node){
   // Loop through neighbor of node and use the SetIn(C->NodeSet)

   int i, inDeg = 0;
   for(i=0;i<C->G->degree[node];i++){
        if(SetIn(C->nodeSet, C->G->neighbor[node][i])) inDeg++;
   }
   return inDeg;
}

// It's not an error to delete a node that's already not there
COMMUNITY *CommunityDelNode(COMMUNITY *C, int i) {
    if(SetIn(C->nodeSet, i)) {
	assert(C->n > 0);
	SetDelete(C->nodeSet, i);
	C->n--;
    }
    else
	printf("Cannot find\n");
    return C;
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
    // Use total degree - InDegree of each node in the community
    int out = 0;
    GRAPH *G = C->G;
    int *memberList = Calloc(C->n, sizeof(int));
    int i, j, k = SetToArray(memberList, C->nodeSet);
    for(i = 0; i < C->n; i++){
        int node = memberList[i];
        out += G->degree[node] - NodeInDegree(C, node);
    }
    Free(memberList);
    return out;
}

void PrintCommunity(COMMUNITY * C){
    int * memberList = Calloc(sizeof(int), C->n);
    int i, j, n = SetToArray(memberList, C->nodeSet);
    for(i = 0; i < n; ++i)
	printf("%d ", memberList[i]);
    printf("\nEnd list\n");
    Free(memberList);
}

/******************** Sets of non-overlapping Communities (partition) ***********/


typedef struct _communitySet {
    unsigned n; // current number of non-empty communities
    GRAPH *G; // the graph we came from
    COMMUNITY **C; // array of pointers to COMMUNITY
    int *where; // Tells where each node belongs to which community
    SET *moved; // Records moved nodes in case of rejection 
    SET *common; // In case merge of 2 communities have overlap, record them
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
	//printf("Total = %f ", P->total);
	P->total -= P->C[c]->score;
	//printf(" Com %d score %f, after %f\n", c, P->C[c]->score, P->total);
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

			    // oldCom and newCom a bit misleading
			    // oldCom is where to move the nodes back in case of reject
			    // newCom is where to get the nodes from
static int _moveOption = -1, _oldCom = -1, _newCom = -1;

 
void MoveRandomNode(PARTITION *P){

    
    int u = P->G->n * drand48();
    int old = P->where[u];
    printf("Mv(%d,%d", u, old);
   
/* 
    for(int i = 0; i < P->G->n; ++i)
	printf("%d ", P->where[i]);
    printf("\n");
*/     
    CommunityDelNode(P->C[old], u);
    _oldCom = old;
    int newCom;
    do{newCom = (int)(P->n * drand48());} // Assumes there is at least one node that is not in ALL communities
    while(newCom == P->where[u] || SetIn(P->C[newCom]->nodeSet, u));  
    CommunityAddNode(P->C[newCom], u);
    P->where[u] = newCom;
    printf("->%d) ", newCom);
    SetAdd(P->moved, u);
    _newCom = newCom;
    //printf("Moved %d from %d to %d\n", u, old, newCom);
}

void MergeCommunities(PARTITION *P, int c1, int c2){
    // Merges C1 and C2 into just C1

    COMMUNITY *C1 = P->C[c1];
    COMMUNITY *C2 = P->C[c2];
    printf("Mg(%d,%d) |%d,%d| ", c1, c2, C1->n, C2->n);
    /*printf("C1\n");
    PrintCommunity(C1);
    printf("C2\n");
    PrintCommunity(C2);*/
    SetIntersect(P->common, C1->nodeSet, C2->nodeSet);
    SetUnion(C1->nodeSet, C2->nodeSet, C1->nodeSet);
    int *memberList = Calloc(C2->n, sizeof(int));
    int i, j, n = SetToArray(memberList, C2->nodeSet);
    for(i = 0; i < C2->n; i++){
	int u = memberList[i]; 
	P->where[u] = c1;
	SetAdd(P->moved, u);
    }
    Free(memberList);
    C1->n = SetCardinality(C1->nodeSet);
    SetEmpty(C2->nodeSet);
    // Automatically takes care of P->n--;
    PartitionDelCommunity(P, c2);
    if(c1 == P->n){
	// In the case C1 is the last community, to preseve the list of communities
	// C1 is placed where C2 was previosuly.
	_newCom = c2;
	//printf("Special case\n");
    }
    else{
	_newCom = c1;    
    }
    _oldCom = P->n; 
   
    //printf("Set o to %d, n to %d\n", _oldCom, _newCom);
}

void SplitCommunity(PARTITION *P, int c_id, int numNodes){
    // Splits P->C[c_id] into a new community that has numNodes random nodes from the old community
    COMMUNITY *oldCom = P->C[c_id];
    if(numNodes < 1 || numNodes > oldCom->n - 1){
	// Always reject if split will make communities of size >2 
        fprintf(stderr, "Invalid split size: %d, originally had %d\n", numNodes, P->C[c_id]->n);
        return;
    }
    printf("Sp(%d,|%d|->%d) ", c_id, oldCom->n, numNodes);

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
        CommunityAddNode(newCom, u);
        CommunityDelNode(oldCom, u);
	SetAdd(P->moved, u);
    }
    Free(memberList);
    PartitionAddCommunity(P, newCom);
    _oldCom = c_id; 
    _newCom = P->n-1;
}


void TestCommunityRoutines(PARTITION *P){
    printf("Calling TestCommunityRoutines\n");
    assert(P->n > 2);
    printf("\nDone with test\n");
}

/*
double ScorePartitionHayes(PARTITION *P) {
    int i;
    double score=0;
    for(i=0; i<P->n; i++) if(P->C[i]) {
	int n = P->C[i]->n;
	if(n>1) {
	    int m = CommunityEdgeCount(P->C[i]);
            printf("m = %d, n = %d\n", m, n);
	    score += m*m/(n*(n-1)/2.0);
	}
        printf("\nScore for community %d = %g\n", i, score);
    }
    return score;

}
*/

/*

    Measures

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
    return inEdges * inEdges / ((C->n * C->n-1)/2.0);
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



double ScorePartition(Boolean global, foint f) {
    
    PARTITION *P = (PARTITION*) f.v;
    //printf("o = %d, n = %d, total n = %d\n", _oldCom, _newCom, P->n);
    printf("\n");
    if(_oldCom == -1 && _newCom == -1){
	printf("No moves have been made yet\n");	    
    }
    else{
	// In case of a split/merge reject P->C[n] will no longer exist
	if(_newCom == P->n){
	    int swap = _newCom;
	    _newCom = _oldCom;
	    _oldCom = swap;
	    //printf("Swapped newCom\n");
	}
	if(_oldCom != P->n){ 
	    COMMUNITY * old = P->C[_oldCom];
	    //printf("o->n %d, ", old->n);
	    double oldBefore = old->score; 
	    old->score = HayesScore(old, CommunityEdgeCount(old));
	    P->total += old->score - oldBefore;
	    //printf("os = %f, ob = %f change = %f ", old->score, oldBefore, old->score - oldBefore);
	}
	COMMUNITY * new = P->C[_newCom];
	//printf("n->n %d\n", new->n);
	double newBefore = newBefore = new->score; 
	new->score = HayesScore(new, CommunityEdgeCount(new));
	P->total += new->score - newBefore;
	//printf("ns = %f, nb = %f change = %f\n", new->score, newBefore, new->score - newBefore);
    }
    //printf("Updated total = %f\n", P->total);
    return P->total;
}



// returns the CHANGE in score due to the move
double PerturbPartition(foint f) {
    //printf("Perturb\n");
    PARTITION *P = (PARTITION*) f.v;
    double before = P->total;
    int choice = drand48() * 3;
    //if(sa) printf("Iter %lu ", sa->iter);
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
	//printf("Merging C1 %d, C2 %d, P->n %d\n", c1, c2, P->n);
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
    //printf("Finding after in P\n"); 
    double after = ScorePartition(true, f);
    printf("Before = %f, After = %f\n", before, after); 
    return after-before;
}

Boolean MaybeAcceptPerturb(Boolean accept, foint f) {
    printf("Accept? %d\n", accept);
    PARTITION * P = (PARTITION *) f.v;
    double before = P->total;
    //printf("P->n = %d, _newCom = %d, _oldCom = %d\n", P->n, _newCom, _oldCom);

    if(accept) ; // do nothing?
    else {
	int * memberList = Calloc(sizeof(int), SetCardinality(P->moved));
	int i, j, n = SetToArray(memberList, P->moved);
	// Could abstract moving these nodes into a function
	if(_moveOption == 1){
	    PartitionAddCommunity(P, CommunityAlloc(P->G));	   
	}
	/*printf("OLDCOM before\n");
	PrintCommunity(P->C[_oldCom]);
	printf("NEWCOM before\n");
	PrintCommunity(P->C[_newCom]);*/
	for(i = 0; i < n; ++i){
	    int u = memberList[i];
	    if(!SetIn(P->common, u)) // Only delete nodes that weren't shared before the change
		CommunityDelNode(P->C[_newCom], u);
	    CommunityAddNode(P->C[_oldCom], u);
	    P->where[u] = _oldCom;
	    //printf("Moved %d from %d to %d\n", u, _newCom, _oldCom);
	}
	/*
	printf("OLDCOM after\n");
	PrintCommunity(P->C[_oldCom]);

	printf("NEWCOM after\n");
	PrintCommunity(P->C[_newCom]);
	*/
	
	Free(memberList);
	if(_moveOption == 2)
	    PartitionDelCommunity(P, _newCom);
    //printf("Finding after in Maybe\n");

    double after = ScorePartition(true, f);
    //printf("Maybe Before = %f, After = %f\n\n", before, after);
    if(before > after)
	fprintf(stderr, "ERROR: Rejection failed, before > after\n");

    }
    SetEmpty(P->moved);
    SetEmpty(P->common);
        
    //printf("Current total = %f\n", P->total); 
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

void simAnnealReport(unsigned long iters, foint currentSolution){

}

int main(int argc, char *argv[])
{

    int i, j;
    srand48(time(0));

    Boolean sparse=maybe, supportNames = true;
    FILE *fp = Fopen(argv[1], "r"); // edge list file is given on the command line
    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNames, false);
    printf("G has %d nodes\n", G->n);

    PARTITION *P = PartitionAlloc(G);
    /*
    int numCommunities = 10; // communities numbered 0 through numCommunities-1 inclusive
    printf("Creating %d random communities: \n", numCommunities);
    for(i=0; i<numCommunities; i++) PartitionAddCommunity(P, CommunityAlloc(G));

    for(i=0; i<G->n; i++) {
	int which = (int)(drand48() * numCommunities);
	CommunityAddNode(P->C[which], i);
        P->where[i] = which;
	//printf("%d->%d, ", i, which);
    }


    //double s = ScorePartitionHayes(P);
    //printf("Hayes score = %g\n", s);
    //double s = ScorePartition(P);
    //printf("\nTotal Partition Score = %g\n", s);



    */
    PartitionFree(P);

    
    printf("BFS-based communities: \n");
    P = PartitionAlloc(G);
    SET *nodesUsed=SetAlloc(G->n); // cumulative set of nodes that have been put into a partition

    int nodeArray[G->n], distArray[G->n];
    while(SetCardinality(nodesUsed) < G->n) {
	int seed;
	do { seed = (int)(drand48() * G->n); }
	    while(SetIn(nodesUsed, seed));
	int numAdded=0, distance = 4; // should be far enough
	int n=GraphBFS(G, seed, distance, nodeArray, distArray); // list of nodes within "distance" of seed
	//printf("BFS(%d[%d])=%d", seed, G->degree[seed], n);
	assert(n>0 && nodeArray[0]==seed && distArray[seed]==0);
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


    //for(i=0; i<numCommunities; i++) CommunityEdgeCount(P->C[i]);
    //puts("");
    //s = ScorePartition(P);
    //printf("Partition score is %g\n", s);
    //TestCommunityRoutines(P);

    printf("%d\n", P->n);
    
    SET * found = SetAlloc(G->n);    
    for(int k = 0; k < P->n; ++k){
	if(P->C[k]->n == 0)
	    printf("ERROR: Com %d has 0 nodes\n", k);
/*	
	int * memberList = Calloc(sizeof(int), P->C[k]->n);
	int l, m, n = SetToArray(memberList, P->C[k]->nodeSet);
	
	for(l = 0; l < n; ++l){
	    if(SetIn(found, memberList[l]))
		printf("ERROR: Overlapping communities\n");
	    SetAdd(found, memberList[l]);
	}
	Free(memberList);
*/
 
    }
    SetFree(found);

    for(int i = 0; i < P->n; ++i){
	double s = HayesScore(P->C[i], CommunityEdgeCount(P->C[i])); 
	P->C[i]->score = s; 
	P->total += s;
    }


#if 0
    HillClimbing(P, 200);
#else
    foint f;
    f.v = P;
    	

    SIM_ANNEAL *sa = SimAnnealAlloc(1, f, PerturbPartition, ScorePartition, MaybeAcceptPerturb, 100, 0, 0, NULL);
    // Boolean SimAnnealSetSchedule(SIM_ANNEAL *sa, double tInitial, double tDecay);
    // void SimAnnealAutoSchedule(SIM_ANNEAL *sa); // to automatically create schedule
    // The above two functions need the accept/reject function to work
    sa->tInitial = sa->tDecay = sa->temperature = 0; // equivalent to hill climbing
    SimAnnealRun(sa); // returns >0 if success, 0 if not done and can continue, <0 if error
    // foint SimAnnealSol(SIM_ANNEAL *sa))
    SimAnnealFree(sa);
#endif
    for(int k = 0; k < P->n; ++k){
	if(P->C[k]->n == 0)
	    printf("ERROR: Com %d has 0 nodes\n", k); 
    }
    printf("Final score = %f\n", P->total);
    printf("Attempting Partition Free\n");
    PartitionFree(P);
    printf("Partition Free completed\n");
    GraphFree(G);
    return 0;
}
