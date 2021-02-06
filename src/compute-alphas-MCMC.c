#include "combin.h"
#include "tinygraph.h"
#include "blant.h"
#include "misc.h"
#include <stdio.h>

int _alphaList[MAX_CANONICALS];
Gint_type _canonList[MAX_CANONICALS];

static int L;
static COMBIN *_Lcombin;
static int *_Darray;
static int _alpha;

/* if all consecutive elements within s share d-1 nodes then alpha++
   Used as function pointer for CombinAllPermutations.
   Iterates through permutations of connected d graphlets from the k graphlet and checks if an ordering of them
   has all consecutive graphlets sharing d-1 vertices in common.
   If so, increments alpha.
   Conceptually, this represents if this permutations of d nodes is possible to walk during sampling.
*/
Boolean _permuteDgraphlets(int size, int* array) {
	int g1, i;
	Boolean validWalk = true;
	for (g1 = 0; g1 < L-1; g1++) { //Unsure if this combination strategy works for d != 2
		TSET tset1 = 0;
		TSET tset2 = 0;
		for (i = 0; i < mcmc_d; i++) {
			TSetAdd(tset1, _Darray[_Lcombin->array[array[g1]]*mcmc_d+i]);
			TSetAdd(tset2, _Darray[_Lcombin->array[array[g1+1]]*mcmc_d+i]);
		}

		if (TSetCardinality(TSetIntersect(tset1, tset2)) != mcmc_d-1) {
			validWalk = false;
			break;
		}
	}
	if (validWalk) _alpha += 1;
	return 0; //Tell permutation to continue
}

/*  This function computes the number of ways SampleGraphlet MCMC can walk over each graphlet to sample it
	For example a path has 2 ways to walk over it.
	The number of ways this is possible is a bit higher than the number of Hamiltonion Paths in our graphlet.
	Given preallocated k graphlet and d graphlet. Assumes Gk is connected
*/
int ComputeAlpha(TINY_GRAPH *Gk, TINY_GRAPH *Gd, unsigned* combinArrayD, unsigned* combinArrayL, int k, int L) {
	assert(k >= 3 && k <=8); //TSET used limits to 8 bits of set represntation.
	_alpha = 0;
	int numDGraphlets = CombinChoose(k, mcmc_d); //The number of possible d graphlets in our k graphlet
	int Darray[numDGraphlets * mcmc_d]; //Vertices in our d graphlets
	_Darray = Darray; //static global variable for permuatation function
    COMBIN * Dcombin = CombinZeroth(k, mcmc_d, combinArrayD); // generates all the edges of g if mcmc_d == 2
	//Fill the S array with all connected d graphlets from Gk
	int SSize = 0;
	int i, j;
	if (mcmc_d != 2) {
		do {
			TSET mask = 0;
			for (i = 0; i < mcmc_d; i++) {
				TSetAdd(mask, Dcombin->array[i]);
			}
			Gd = TinyGraphInduced(Gd, Gk, mask);
			if (TinyGraphDFSConnected(Gd, 0))
			{
				for (i = 0; i < mcmc_d; i++)
				{
					Darray[SSize*mcmc_d+i] = Dcombin->array[i];
				}
				SSize++;
			}
		} while (CombinNext(Dcombin));
	} else {
		do { //if there is an edge between any two vertices in the graphlet
			if (TinyGraphAreConnected(Gk, combinArrayD[0], combinArrayD[1]))
			{ //add it to the Darray
				Darray[SSize*mcmc_d] = Dcombin->array[0];
				Darray[SSize*mcmc_d+1] = Dcombin->array[1];
				SSize++;
			}
		} while (CombinNext(Dcombin));
	}

	//for s over all combinations of L elements in S
	COMBIN *Lcombin = CombinZeroth(SSize, L, combinArrayL);
	_Lcombin = Lcombin;
	do {
		//add vertices in combinations to set.
		TSET mask = 0;
		for (i = 0; i < L; i++) {
			for (j = 0; j < mcmc_d; j++) {
				TSetAdd(mask, Darray[Lcombin->array[i]*mcmc_d+j]);
			}
		}
		//if Size of nodeset of nodes in each element of s == k then
		if (TSetCardinality(mask) == k)
		{
			int permArray[L];
			memset(permArray, L, sizeof(int));
			//all permutations of elements within s do
			CombinAllPermutations(L, permArray, _permuteDgraphlets);
		}
	} while (CombinNext(Lcombin));

	CombinFree(Dcombin);
	CombinFree(Lcombin);
	return _alpha / 2;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "USAGE: %s k\n", argv[0]);
        exit(-1);
    }
    int k = atoi(argv[1]);
    char BUF[BUFSIZ];
    TINY_GRAPH *gk = TinyGraphAlloc(k);
    TINY_GRAPH *gd = TinyGraphAlloc(mcmc_d);
    SET* _connectedCanonicals = canonListPopulate(BUF, _canonList, k);
    int numCanon = _connectedCanonicals->n;

    L = k - mcmc_d  + 1;
    unsigned combinArrayD[mcmc_d]; //Used to hold combinations of d graphlets from our k graphlet
	unsigned combinArrayL[L]; //Used to hold combinations of L d graphlets from Darray
	// create the alpha list
	int i;
	for (i = 0; i < numCanon; i++) {
		Int2TinyGraph(gk, _canonList[i]);
		TinyGraphEdgesAllDelete(gd);
		if (SetIn(_connectedCanonicals, i)) {
			_alphaList[i] = ComputeAlpha(gk, gd, combinArrayD, combinArrayL, k, L);
		}
		else _alphaList[i] = 0; // set to 0 if unconnected graphlet
	}

    printf("%d\n", numCanon);
	for (i = 0; i < numCanon; i++) {
		printf("%d ", _alphaList[i]);
	}
	printf("\n");

	TinyGraphFree(gk);
	TinyGraphFree(gd);
    return 0;
}
