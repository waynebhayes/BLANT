#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "blant.h"
#include "sets.h"

int _numCanon;
short int* _K;

void SetGlobalCanonMaps(int k){
    assert(3 <= k && k <= 8);
    unsigned int _Bk = (1 <<(k*(k-1)/2));
    char BUF[BUFSIZ];
    Gint_type _canonList[MAX_CANONICALS];
    SET *_connectedCanonicals = canonListPopulate(BUF, _canonList, k);
    _numCanon = _connectedCanonicals->n;
    _K = (short int*) mapCanonMap(BUF, _K, k);

    sprintf(BUF, CANON_DIR "/perm_map%d.bin", k);
}

int getHammingDistance(int a, int b){
    // return number of different corresponding bits in the binary representation of a and b
    int xor = a^b;
    int i;
    int ans = 0;

    for(i=0; i < (sizeof(int) * 8); i++){
    	if((xor%2) != 0)
    		ans += 1;
    	xor = xor>>1;
    }

    if(a==b)
    	assert(ans == 0);

    assert(ans >= 0);
    return ans;
}

int main(int argc, char *argv[]){
	// USAGE: makeEHD 3 > ehd3.txt
	// prints -
	//        g1 g2 ehd
	//        g1 g3 ehd
	//        ...
	//        ...

	if(argc != 2){
		fprintf(stderr, "specify a 'k' value. USAGE: ./makeEHD 5 > ehd5.txt\n");
		exit(1);
	}

	int k = atoi(argv[1]);
	SetGlobalCanonMaps(k);  // construct lookup table

	int EHD[_numCanon][_numCanon];
	int i, j, c1, c2, id1, id2;

	// initialize the EHD matrix
	for(i=0; i<_numCanon; i++)
		for(j=0; j<_numCanon; j++)
			EHD[i][j] = k*k + 1;   // some big number

	int perms = pow(2, (k*(k-1))/2);  // total number of permutations of all canonicals.
	// loop through this number -> thus creating all possible canonical permutations

	// given a canonical number, store any integer ID (this can also be done by reading the canon_maps files)
	int canon2int[_numCanon];
	for(i=0; i<perms; i++){
		assert((_K[i] >= 0) && (_K[i] < _numCanon));
		canon2int[_K[i]] = i;
	}

	// MAIN loop
	for(c1=0; c1<_numCanon; c1++){  // loop through each *canonical*, connected or disconnected
		id1 = canon2int[c1];  // given a canonical, get an integer ID
		assert((id1>=0) && (id1<perms));

		for(id2=0; id2<perms; id2++){
			c2 = _K[id2];
			EHD[c1][c2] = MIN(EHD[c1][c2], getHammingDistance(id1, id2));
		}
	}

	// print the result
	for(i=0; i<_numCanon; i++){
		assert(EHD[i][i] == 0);
		for(j=0; j<_numCanon; j++){
            assert(EHD[i][j] < ((k*k) + 1));
			assert(EHD[i][j] == EHD[j][i]);
			printf("%d %d %d\n", i, j, EHD[i][j]);
		}
	}

	return 0;
}
