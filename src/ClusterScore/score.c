#include "misc.h"
#include "graph.h"
#include "sets.h"
#include "rand48.h"


typedef struct _clustering {
    GRAPH *G;
    SET **clusters;
    unsigned nC, *clusSize;
} CLUSTERING;

// Allocate a clustering and return it with the zero'th cluster assigned to all nodes in G.
CLUSTERING *ClusteringAlloc(GRAPH *G) {
    CLUSTERING *C = Calloc(sizeof(CLUSTERING), G->n);
    C->G = G;
    C->nC=0;
    C->clusters = Calloc(G->n, sizeof(SET*));
    C->clusters[0] = SetAlloc(G->n);
    for(int i=0;i < G->n; i++) SetAdd(C->clusters[0], i);
    C->clusSize[0] = G->n;
    ++C->nC;
    return C;
}

double ScoreOneCluster(CLUSTERING *C, int who) {
    unsigned nodes[C->G->n], n = SetToArray(nodes, C->clusters[who]), m=0;
    assert(n>0);
    assert(n == C->clusSize[who]);
    if(n==1) return 0;
    for(int i=0; i<n; i++) for(int j=i+1;j<n;j++) if(GraphAreConnected(C->G,nodes[i],nodes[j])) ++m;
    double ED= m/(n*(n-1)/2.0);
    return m*ED;
}

double ScoreClustering(CLUSTERING *C) {
    double score=0;
    for(int c=0; C->clusters[c]; c++) { double CS = ScoreOneCluster(C, c); printf(" %g", CS); score += CS; }
    return score;
}


// Remove some random set of elements from C[who], and return the set of removed elements
void SplitCluster(CLUSTERING *C, unsigned who) {
    GRAPH *G=C->G;
    unsigned nodes[G->n], n = SetToArray(nodes, C->clusters[who]);
    assert(n>1);
    assert(C->clusters[C->nC]==NULL);
    C->clusters[C->nC]=SetAlloc(G->n);
    int which, num2move = 1+drand48()*(n-1)/2; // remove at least 1, and at most half FIXME: might this introduce a bias?
    // int which, num2move = 1+drand48()*(n-1);   // remove at least 1, and at most (n-1) FIXME: might this introduce a bias?
    for(int i=0; i<num2move; i++) {
	until(SetIn(C->clusters[who], nodes[(which = drand48()*n)])) ;
	SetDelete(C->clusters[who], nodes[which]) ;
	SetAdd(C->clusters[C->nC], nodes[which]);
    }
    C->clusSize[C->nC]= num2move;
    C->clusSize[who] -= num2move;
    C->nC++;
}

Boolean TrySplit(CLUSTERING *C){
    GRAPH *G = C->G;
    assert(0 < C->nC && C->nC < G->n-1);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    int oldNC = C->nC, who, oldSize=0, tries = 0;
    while(oldSize < 2) {
	who = oldNC * drand48();
	oldSize = C->clusSize[who];
	assert(++tries < G->n);
    }
    double oldScore = ScoreOneCluster(C, who);
    SplitCluster(C, who);
    assert(C->nC == oldNC+1);
    double newScore = ScoreOneCluster(C, who) + ScoreOneCluster(C, oldNC);
    if(newScore > oldScore) return true;
    else {
	SetUnion(C->clusters[who], C->clusters[who], C->clusters[oldNC]);
	SetFree(C->clusters[oldNC]);
	C->clusters[oldNC] = NULL;
	C->clusSize[oldNC]= 0;
	C->clusSize[who]  = oldSize;
	C->nC = oldNC;
	return false;
    }
}

Boolean TryMove(CLUSTERING *C){
    assert(0 < C->nC && C->nC < C->G->n);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    if(C->nC < C->G->n-1) return TrySplit(C);
    return false;
}

int main(int argc, char *argv[])
{
    Boolean sparse=true, supportNodeNames=true, weighted=false;
    FILE *fp=Fopen(argv[1], "r");
    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNodeNames, weighted);

    CLUSTERING *C = ClusteringAlloc(G);

    // At this point we have one big cluster consisting of the entire graph, so start splitting

    srand48(time(0));
    while(true) {
	printf("Status: %d clusters: ", C->nC);
	double fullScore = ScoreClustering(C);
	printf(" full score %g\n", fullScore);
	TryMove(C);
    }
    return 0;
}

