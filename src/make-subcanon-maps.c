#include <stdio.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

// A "sub-canononical" is the (k-1)-graphlet you get when you delete one of the k nodes of a k-graphlet.
// This program creates the list of k sub-canonicals, each of size (k-1), that result from deleting each node
// of a k-graphlet, across all k-graphlet canonicals.

static short int *K; // The (k-1)-graphlet canonmap.
static Gint_type canon_list[MAX_CANONICALS];

int main(int argc, char* argv[]) {
    int i, j, k;
    k=atoi(argv[1]);
    TINY_GRAPH* G = TinyGraphAlloc(k);

    //Create k canon list
    char BUF[BUFSIZ];
    SET *connectedCanonicals = canonListPopulate(BUF, canon_list, k);
    int numCanon = connectedCanonicals->n;
    SetFree(connectedCanonicals);

    //Create canon map for k-1
    K = mapCanonMap(BUF, K, k-1);

    /*
        For every canonical in canonical_list
            create tinygraph
            print canonical
            For every induced subgraph of size k-1
                get canonical from k-1
                print decimal for canonical
            print newline

    */
    TINY_GRAPH *g = NULL;
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
