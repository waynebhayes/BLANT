#include "misc.h"
#include "sets.h"
#include "graph.h"
#include <stdio.h>

// This should be an exact substitute for a similarly named bash+awk script in the blant-clusters.sh suite.
// Input file must have lines in the following format:
//   n m w v1 v2 v3 ...
// where n=number of nodes in the cluster, m is the number of edges, w is mean edge weight, and v_i are the member nodes.
// G->n is the maximum expected cluster size (easiest if it's just the number of nodes in the input graph).
void RemoveOverlapClusters(double overlap, GRAPH *G, FILE *fp) {
    unsigned i, numClus=0, kk[G->n];
    SET *S = SetAlloc(G->n);
    SET **cluster = Calloc(G->n, sizeof(SET*)); // actual number of clusters can probably be bigger...
    unsigned *edges = Calloc(G->n, sizeof(unsigned));
    float *edgeSum = Calloc(G->n, sizeof(float));
    while(!feof(fp)) {
	double edgeWgts;
	unsigned numNodes, edgeHits, k;
	int numRead = fscanf(fp, "%u%u%lf%u", &numNodes, &edgeHits, &edgeWgts, &k);
	if(numRead!=4) Fatal("couldn't find n m w k, got only %d values: %u %u %g %u",numRead, numNodes,edgeHits,edgeWgts,k);
	SetEmpty(S);
	for(i=0;i<numNodes;i++) {
	    char name[BUFSIZ]; strcpy(name,"<undefined>");
	    if(fscanf(fp, "%s ", name) != 1) Fatal("couldn't get nodeName");
	    unsigned node = GraphNodeNameToInt(G,name);
	    SetAdd(S, node);
	}
	Boolean add=true;
	for(i=0;i<numClus;i++) if(SetIntersectCount(S, cluster[i]) >= overlap*SetCardinality(S)) {
	    add=false; break;
	}
	if(add) {
	    if(numClus >= G->n) {
		Warning("numClus has reached the maximum of %u, not reading any more clusters", numClus);
		break;
	    }
	    cluster[numClus] = SetCopy(NULL, S);
	    edges[numClus]=edgeHits; edgeSum[numClus]=edgeWgts;kk[numClus]=k;
	    ++numClus; 
	}
    }
    for(i=0;i<numClus;i++) {
	unsigned j, memberList[G->n];
	int n = SetCardinality(cluster[i]);
	unsigned maxEdges=n*(n-1)/2;
	printf("%d %d %g %d", n ,edges[i], edgeSum[i], kk[i]);
	SetToArray(memberList, cluster[i]);
	for(j=0;j<SetCardinality(cluster[i]);j++) printf(" %s", G->name[memberList[j]]);
	puts("");
    }
}

const char const * const USAGE = "expecting \"OVERLAP_BOUND edgeList.el\" on command line";

int main(int argc, char *argv[]) {
    if(argc!=3) Fatal(USAGE);
    double overlap=atof(argv[1]);
    if(overlap<0 || overlap>1) Fatal("overlap %g must be in [0,1]", overlap);
    FILE *graphFile = Fopen(argv[2],"r"); assert(graphFile);
    GRAPH *G = GraphReadEdgeList(graphFile, false, false, false);
    RemoveOverlapClusters(overlap, G, stdin);
}
