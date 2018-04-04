#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "rand48.h"
#include "heap.h"
#include "blant.h"

#define CANON_DIR "canon_maps"

static int _numCanon, _canonList[MAX_CANONICALS];

static short int _K[maxBk] __attribute__ ((aligned (8192)));

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array. 
*/
void mapCanonMap(char* BUF, short int *K, int k) {
    int Bk = (1 <<(k*(k-1)/2));
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    short int *Kf = Mmap(K, Bk*sizeof(K[0]), Kfd);
    assert(Kf == K);
}

int main(int argc, char* argv[]) {
    char* BUF[99999];
    mapCanonMap(BUF, _K, 7);
    while (true) {
        scanf("%s", BUF);
        printf("%i\n", _K[atoi(BUF)]);
    }
}