#include "misc.h"
#include "sets.h"
#include "graph.h"
#include <stdio.h>

// This should be an exact substitute for a similarly named bash+awk script in the blant-clusters.sh suite.
// Input file must have lines in the following format:
//   n m w v1 v2 v3 ...
// where n=number of nodes in the cluster, m is the number of edges, w is mean edge weight, and v_i are the member nodes.
// G->n is the maximum expected cluster size (easiest if it's just the number of nodes in the input graph).
int RemoveSubsetClusters(int k, GRAPH *G, FILE *fp) {
    unsigned i, numClus=0;
    SET *S = SetAlloc(G->n);
    SET **cluster = Calloc(G->n, sizeof(SET*)); // actual number of clusters can probably be bigger...
    unsigned *edges = Calloc(G->n, sizeof(unsigned));
    float *edgeSum = Calloc(G->n, sizeof(float));
    while(!feof(fp)) {
	unsigned numNodes, edgeHits; double edgeWgts;
	int numRead = fscanf(fp, "%u%u%lf", &numNodes, &edgeHits, &edgeWgts);
	if(numRead != 3) Fatal("couldn't find n m w, got only %d values: %u %u %g", numRead, numNodes,edgeHits,edgeWgts);
	SetEmpty(S);
	for(i=0;i<numNodes;i++) {
	    char name[BUFSIZ]; strcpy(name,"<undefined>");
	    if(fscanf(fp, "%s ", name) != 1) Fatal("couldn't get nodeName");
	    unsigned node = GraphNodeNameToInt(G,name);
	    SetAdd(S, node);
	}
	Boolean add=true;
	for(i=0;i<numClus;i++) if(SetSubsetEq(S, cluster[i])) { add=false; break;} // eliminate only EXACT subsets for now
	if(add) {
	    if(numClus >= G->n) {
		Warning("numClus has reached the maximum of %u, not reading any more clusters", numClus);
		break;
	    }
	    cluster[numClus] = SetCopy(NULL, S);
	    edges[numClus]=edgeHits; edgeSum[numClus]=edgeWgts;
	    ++numClus; 
	}
    }
    for(i=0;i<numClus;i++) {
	unsigned j, memberList[G->n];
	int n = SetCardinality(cluster[i]);
	unsigned maxEdges=n*(n-1)/2;
	printf("%d %d %g %d", n ,edges[i], edgeSum[i], k);
	SetToArray(memberList, cluster[i]);
	for(j=0;j<SetCardinality(cluster[i]);j++) printf(" %s", G->name[memberList[j]]);
	puts("");
    }
}

const char const * const USAGE = "expecting \"k edgeList.el\" on command line";

int main(int argc, char *argv[]) {
    if(argc!=3) Fatal(USAGE);
    unsigned k=atoi(argv[1]);
    if(k<3 || k>9) Fatal("k must be between 3 and 9");
    FILE *graphFile = Fopen(argv[2],"r"); assert(graphFile);
    GRAPH *G = GraphReadEdgeList(graphFile, false, false, false);
    RemoveSubsetClusters(k, G, stdin);
}
