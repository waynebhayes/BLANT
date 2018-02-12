#include <stdio.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

static unsigned int Bk;//Bk=actual number of entries in the canon_map for given k - 1.

// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int K[maxBk] __attribute__ ((aligned (8192)));

int main(int argc, char* argv[]) {
    int i, j, k;
    k=atoi(argv[1]);
    TINY_GRAPH* G = TinyGraphAlloc(k);

    Bk = (1 <<((k-1)*(k-2)/2));
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
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k-1);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    short int *Kf = Mmap(K, Bk*sizeof(K[0]), Kfd);
    assert(Kf == K);

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
        BuildGraph(canonical, G);
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
