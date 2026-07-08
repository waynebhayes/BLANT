// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <stdio.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

// A "sub-canononical" is the (k-1)-canonical graphlet you get when you delete one of the k nodes of a k-graphlet.
// This program creates the list of k sub-canonicals, each of size (k-1), that result from deleting each node
// of a k-graphlet, across all k-graphlet canonicals.

static Gordinal_type *K; // The (k-1)-graphlet canonmap.
static Gint_type canon_list[MAX_CANONICALS];
static char canon_num_edges[MAX_CANONICALS];
Boolean directed=false;

int main(int argc, char* argv[]) {
    int i, j, k;
    if(argc == 3) {
	if(strcmp(argv[2], "directed")!=0) {
	    fprintf(stderr, "if given two arguments, the second must be the word \"directed\"\n");
	    exit(1);
	}
	directed=true;
    } else assert(argc==2);
    k = atoi(argv[1]); assert(2<=k && k<=8);
    if(directed && k > 6) {
	fprintf(stderr, "Error: directed graphs are only supported for k<=6 (got k=%d)\n", k);
	exit(1);
    }

#if SELF_LOOPS
    if (k>7) Fatal("cannot make_subcanon_maps when SELF_LOOPS");
    TINY_GRAPH* G = TinyGraphSelfAlloc(k);
#else
    TINY_GRAPH* G = TinyGraphAlloc(k,false,false);
#endif
    //Create k canon list
    char BUF[BUFSIZ];
    SET *connectedCanonicals = canonListPopulate(BUF, canon_list, k, canon_num_edges, directed);
    int numCanon = connectedCanonicals->maxElem;
    SetFree(connectedCanonicals);

    //Create canon map for k-1
    K = mapCanonMap(BUF, K, k-1, directed);

    /*
        For every canonical in canonical_list
            create tinygraph
            print canonical
            For every induced subgraph of size k-1
                get canonical from k-1
                print decimal for canonical
            print newline

    */
    TINY_GRAPH *g = TinyGraphAlloc(k-1,G->selfLoops,false);
    TSET induceTSET = TSET_NULLSET;
    Gint_type Gint = 0;
    int tsetBit;
    for (i = 0; i < numCanon; i++) {
        int canonical = canon_list[i];
        Int2TinyGraph(G, canonical);
        printf("%d ", canonical);
        //Reset induceTSET
        for (j = 0; j < k; j++) {
            induceTSET = TSET_NULLSET;
            for (tsetBit = 0; tsetBit < k; tsetBit++) {
                induceTSET <<= 1;
                if (tsetBit != j)
                    induceTSET |= 1;
            }
            g = TinyGraphInduced(g, G, induceTSET);
            Gint = TinyGraph2Int(g, k-1);
            printf("%d ", K[Gint]);
        }
        printf("\n");
    }
    assert(g && G);
    TinyGraphFree(g);
    TinyGraphFree(G);
    return 0;
}
