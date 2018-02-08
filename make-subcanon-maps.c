#include <stdio.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"


//needed for blant stuff
#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map
#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/canon_maps" // if you happen to put it there...

static unsigned int Bk, _k; // _k is the global variable storing k; Bk=actual number of entries in the canon_map for given k.

// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int K[maxBk] __attribute__ ((aligned (8192)));

//Needed for BuildGraph
static TINY_GRAPH *G;
static int k;

#define MMAP 1 
// Try to mmap, and if it fails, just slurp in the file (sigh, Windoze)
void *Mmap(void *p, size_t n, int fd)
{
#if MMAP
    void *newPointer = mmap(p, n, PROT_READ, MAP_PRIVATE|MAP_FIXED, fd, 0);
    if(newPointer == MAP_FAILED)
#endif
    {
#if !__WIN32__ && !__CYGWIN__ // it will always fail on Windoze so don't bother with a warning
	Warning("mmap failed");
#endif
	if(read(fd, p, n) != n)
	    Fatal("cannot mmap, or cannot read the file, or both");
    }
    return p;
}

static void BuildGraph(int Gint)
{
    int i, j, bitPos=0;
    int Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
    TinyGraphEdgesAllDelete(G);
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
    {
	for(j=i-1;j>=0;j--)
#else	// UPPER_TRIANGLE
    for(i=k-2;i>=0;i--)
    {
	for(j=k-1;j>i;j--)
#endif
	{
	    if(!Gint2) break;
	    int bit = (1 << bitPos);
	    if(Gint & bit)
		TinyGraphConnect(G,i,j);
	    Gint2 &= ~bit;
	    bitPos++;
	}
	if(!Gint2) break;
    }
}

#define LOWER_TRIANGLE 1
// Given a TINY_GRAPH and k, return the integer ID cretaed from one triangle (upper or lower) of the adjacency matrix.
static int TinyGraph2Int(TINY_GRAPH *g, int numNodes)
{
    int i, j, bitPos=0, Gint = 0, bit;
    
#if LOWER_TRIANGLE	// Prefer lower triangle to be compatible with Ine Melckenbeeck's Jesse code.
    for(i=numNodes-1;i>0;i--)
    {
        for(j=i-1;j>=0;j--)
#else   // UPPER_TRIANGLE // this is what we used in the original faye code and paper with Adib Hasan and Po-Chien Chung.
    for(i=numNodes-2;i>=0;i--)
    {
        for(j=numNodes-1;j>i;j--)
#endif
        {
	    if(TinyGraphAreConnected(g,i,j))
	    {
		bit = (1 << bitPos);
		Gint |= bit;
	    }
            bitPos++;
        }
    }
    return Gint;
}

int main(int argc, char* argv[]) {
    int i, j;
    k=atoi(argv[1]);
    _k = k;
    G = TinyGraphAlloc(k);

    Bk = (1 <<(k*(k-1)/2));
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

    //Getting and filling canon maps
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
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
    TINY_GRAPH *g;
    TSET induceTSET = TSET_NULLSET;
    int Gint = 0;
    for (i = 0; i < numCanon; i++) {
        int canonical = canon_list[i];
        BuildGraph(canonical);
        printf("%d ", canonical); 

        //Reset induceTSET
        for (j = 0; j < k; j++) {
            induceTSET = TSET_NULLSET;
            for (int tsetBit = 0; tsetBit < k; tsetBit++) {
                induceTSET <<= 1;
                if (tsetBit != j)
                    induceTSET |= 1;
            }
            g = TinyGraphInduced(g, G, induceTSET);
            Gint = TinyGraph2Int(g, k);
            printf("%d ", K[Gint]);
        }
        printf("\n");
    }

    free(g);
    free(G);

    return 0;
}