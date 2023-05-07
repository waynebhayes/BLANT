#include "misc.h"
#include "sets.h"
#include "graph.h"
#include "hashmap.h"
#include "priority-queue.h"
#include "rand48.h"

unsigned _Gn; // number of nodes in the full input network

typedef struct _commType {
    unsigned n,m,nc2,k, index;
    double ED, score;
    SET *nodes;
} CLUSTER;

static int ClusterCompare(foint f1, foint f2)
{
    CLUSTER *c1 = (CLUSTER*)f1.v, *c2 = (CLUSTER*)f2.v;
    // abstract "smallest-at-top" wins, so reverse it for largest
    if(c1->score < c2->score) return +1;
    if(c1->score > c2->score) return -1;
    return 0;
}


// read the output of blant-clusters.sh, looks like this (after it's put through sed to remove non-digits)
// n  m   nc2 k ED%     list
// 19 154 171 5 90.0585 241 44 129 184 227 168 178 112 256 197 293 189 633 106 14 636 171 15 190
// returns number of nodes
CLUSTER *ReadCluster(FILE *fp)
{
    unsigned i,v;
    CLUSTER *c = (CLUSTER*) Calloc(1,sizeof(CLUSTER));
    c->nodes = SetAlloc(_Gn);
    fscanf(fp, "%u%u%u%u%lf", &c->n, &c->m, &c->nc2, &c->k, &c->ED);
    c->ED /= 100; // blant outputs a percentage, we want a fraction
    for(i=0; i<c->n;i++) {fscanf(fp, "%u", &v); SetAdd(c->nodes,v);}
    return c;
}

CLUSTER **_cluster; // gets dynamically bigger as necessary, starting at _Gn
unsigned _numClus, _clusSize, *_finalMemberships; // array [_Gn] counting number of clusters each node in INPUT graph belongs to
double _overlapThresh = 0.5;
GRAPH *_clusterSimGraph;
SET *_finalComm;
#define MAX_GRAPH 1000000

enum Measure { undef, MEASURE_OMOD, MEASURE_EDN };
static enum Measure measure = undef;
double _stopT = 0; // stop threshold

PRIORITY_QUEUE *_PQ;

double scoreOfNodeInCommunity(CLUSTER *c, unsigned u, unsigned membershipCount)
{
    assert(measure != undef);
    if(measure==MEASURE_OMOD) {Apology("boo hoo, no OMOD"); return 0;}
      // return ( (kin[c][u] - ( degree[u]-kin[c][u] ) ) / degree[u] ) * ( 1 / ms ) * edgeDensity[c] * ( 1 / nc[c] )
    if(measure == MEASURE_EDN) return ( 1.0 / membershipCount ) * c->ED;
    return 0;
}


unsigned *sList(SET *s, unsigned *space)
{
    if(s->list) return s->list;
    else SetToArray(space, s);
    return space;
}

double communityScore(CLUSTER *c)
{
    int i;
    double CS = 0;
    unsigned member[c->n], *list = sList(c->nodes, member);
    assert(c->n==SetCardinality(c->nodes));
    for(i=0;i<c->n;i++) CS += scoreOfNodeInCommunity(c, list[i], _finalMemberships[list[i]]+1);
    return CS;
}

double _currentScore;

double potentialScore(CLUSTER *c)
{
    double P = _currentScore;
    unsigned clusNodes[c->n], *clist = sList(c->nodes, clusNodes), i,j;
    unsigned currentClusters[SetCardinality(_finalComm)], *ccList = sList(_finalComm, currentClusters);
    for(i=0;i<c->n;i++) {
	unsigned u = clist[i];
	if(_finalMemberships[u]!=0) {
	    for(j=0; j<SetCardinality(_finalComm); j++) {
	    CLUSTER *c2 = _cluster[ccList[j]];
	    if(SetIn(c2->nodes, u))
		P=P-scoreOfNodeInCommunity(c2, u, _finalMemberships[u])+scoreOfNodeInCommunity(c2, u, _finalMemberships[u] + 1);
	    }
	}
    }
    P+=communityScore(c);
    return P;
}



void addToResult(CLUSTER *c)
{
    double P=potentialScore(c), diff=P-_currentScore;
    if (diff < _stopT) return;
    _currentScore=P;
    int array[c->n], *list = sList(c->nodes, array), i;
    for(i=0;i<c->n;i++) {
	unsigned u = list[i];
	_finalMemberships[u]++;
    }
    SetAdd(_finalComm,c->index);
}

SET *_finalClusVisited; 

void expand(CLUSTER *c)
{
    if(GraphDegree(_clusterSimGraph, c->index)) {
	int j; 
	for(j=0;j<_clusterSimGraph->degree[c->index];j++) {
	    unsigned n=_clusterSimGraph->neighbor[c->index][j];
	    if(!SetIn(_finalClusVisited, n)) {
		double P=potentialScore(_cluster[n]), diff=P-_currentScore;
		if ( diff > _stopT ) {
		    _cluster[n]->score = P;
		    PriorityQueueInsert(_PQ,(foint)(void*)_cluster[n]);
		}
		SetAdd(_finalClusVisited,n);
	    }
	}
    }
}



int main(int argc, char *argv[])
{
    int i,j;
    assert(argc==3);
    _clusSize = _Gn = atoi(argv[1]);
    _stopT = atof(argv[2]);
    _cluster = Malloc(_clusSize*sizeof(CLUSTER));
    SET *intersect = SetAlloc(_Gn);
    _finalMemberships = Calloc(_Gn, sizeof(_finalMemberships[0]));

    Boolean sparse = true, names=false;
    _clusterSimGraph = GraphAlloc(MAX_GRAPH, sparse, names);
    GraphMakeWeighted(_clusterSimGraph);
    unsigned sim[MAX_GRAPH];

    while(!feof(stdin)) {
	CLUSTER *c = ReadCluster(stdin);

	// Either remember, or forget, this cluster based on overlap with previous ones
	for(i=0;i<_numClus;i++)  {
	    SetReset(intersect);
	    SetIntersect(intersect, _cluster[i]->nodes, c->nodes);
	    assert(SetCardinality(_cluster[i]->nodes) >= SetCardinality(c->nodes));
	    sim[i]=SetCardinality(intersect);
	    if(sim[i] / (1.0*SetCardinality(_cluster[i]->nodes)) > _overlapThresh) {
		//printf("Skipping cluster %d\n", _numClus);
		SetFree(c->nodes);
		Free(c);
		break;
	    }
	}
	if(i<_numClus) continue; // loop above terminated early, meaning c had too much overlap with existing cluster

	for(i=0;i<_numClus;i++) if(sim[i]>0) {
	    //printf("sim %d %d = %d\n",i,_numClus,sim[i]);
	    GraphSetWeight(_clusterSimGraph,i,_numClus,sim[i]);
	}

	if(_numClus == _clusSize) {
	    printf("Upgrading _clusSize from %d to %d\n", _clusSize, 2*_clusSize);
	    _clusSize *=2;
	    _cluster = Realloc(_cluster, _clusSize*sizeof(CLUSTER));
	}
	c->index = _numClus;
	_cluster[_numClus++] = c;
	fscanf(stdin, " ");
    }
#if 0 // debugging just to check we read things correctly
    printf("Read %d clusters:\n", _numClus);
    for(i=0; i<_numClus; ++i) {
	CLUSTER *c=_cluster[i];
	printf("Cluster %d: n=%d m=%d nc2=%d k=%d ED=%g%%: ", i, c->n, c->m, c->nc2, c->k, c->ED);
	SetPrint(_cluster[i]->nodes);
	puts("");
    }
#endif

    Warning("May need to load original graph to compute overlapping modularity");

    _PQ = PriorityQueueAlloc(_numClus, ClusterCompare, NULL);
    measure = MEASURE_EDN;

    // For each connected component, mark it off and pick one node at random to put in the PQ
    int numVisited = 0, Varray[_numClus], rootNode=0;
    SET *visitedDFS = SetAlloc(_numClus);
    _finalClusVisited = SetAlloc(_numClus);
    for(rootNode=0;rootNode<_numClus;rootNode++) {
	if(!SetIn(visitedDFS, rootNode)) {
	    int prevNum = numVisited;
	    GraphVisitCC(_clusterSimGraph, rootNode, visitedDFS, Varray, &numVisited);
	    assert(numVisited == SetCardinality(visitedDFS));
	    // pick a node at random from the most recent CC
	    int CCsize = numVisited-prevNum, randOffset=CCsize*drand48(), ccRandNode=Varray[prevNum+randOffset];
	    _cluster[ccRandNode]->score = communityScore(_cluster[ccRandNode]);
	    PriorityQueueInsert(_PQ, (foint)(void*)_cluster[ccRandNode]);
	    SetAdd(_finalClusVisited, ccRandNode);
	}
    }
    printf("Cluster graph has %d nodes; PQ size is %d\n", _numClus, PriorityQueueSize(_PQ));
    _finalComm = SetAlloc(_numClus);

    while (PriorityQueueSize(_PQ)) {
	CLUSTER *c = PriorityQueueNext(_PQ).v;
	addToResult(c);
	expand(c);
    }

    assert(measure!=MEASURE_OMOD); // printf "Qov=%s",Q/length(finalComm)
    printf("EDN=%g\n",_currentScore);

    unsigned community[_numClus], *list = sList(_finalComm, community);
    for(i=0; i<SetCardinality(_finalComm);i++) {
	unsigned index = list[i];
	CLUSTER *c=_cluster[index];
	printf("%d nodes, %d of %d edges from k %d (%g%%):", c->n, c->m, c->nc2, c->k, 100*c->ED);
	unsigned nodes[c->n], *node = sList(c->nodes, nodes);
	for(j=0;j<c->n;j++) {
	    unsigned u = node[j];
	    printf(" %d", u);
	}
	puts("");
    }
    return 0;
}
