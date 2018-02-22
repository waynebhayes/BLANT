#include <stdio.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

// Here we are allocating 256MB x sizeof(short int) = 512MB for the canon map.
static short int K[maxBk] __attribute__ ((aligned (8192)));

int main(int argc, char* argv[]) {
    int i, j, k;
    k=atoi(argv[1]);
    TINY_GRAPH* G = TinyGraphAlloc(k);

    
    char BUF[BUFSIZ];
    //Create canon list for k
    sprintf(BUF, CANON_DIR "/canon_list%d.txt", k);
    FILE *fp_ord=fopen(BUF, "r");
    assert(fp_ord);
    int numCanon;
    fscanf(fp_ord, "%d",&numCanon);
    int canon_list[numCanon];
    for(i=0; i<numCanon; i++) fscanf(fp_ord, "%d", &canon_list[i]);
    fclose(fp_ord);

    //Create canon map for k-1
    short int *Kf = createCanonMap(BUF, K, k-1);

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
    int Gint = 0;
    int tsetBit;
    for (i = 0; i < numCanon; i++) {
        int canonical = canon_list[i];
        BuildGraph(G, canonical);
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
