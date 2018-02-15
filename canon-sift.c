/*
** graphette-sift.c: given a list of exactly 2 files on the command line
** that contain parallel-generated graphettes that were adjacent during
** the parallel run, assume the first file contains the "real" canonicals.
** Then, when reading the second file, we need map the faux canonicals
** to the "real" ones and note any genuinely new ones we find, and also
** re-map repeated canonicals in the later files back to the "real" ones.
** Here's the logic:
**
** STEP 1: READ THE FIRST FILENAME ON THE COMMAND LINE
   while not EOF
	read line
	if $1 == $2 then
	    store _canonicalGraph[_numCanon] and _canonicalSig[_numCanon]
	    ++_numCanon
	end if
	echo line to stdout
    end while

    STEP 2: read the second file
    while not eof
	read line
	if $1==$2 then // "new" possibly faux canonical
	    if $1 is isomorphic to any existing true canonical
		_fauxCanonSig[_numFauxCanon] = $1
		_fauxCanonTrue[_numFauxCanon] = trueCanon
		_fauxCanonPerm[_numFauxCanon] = perm of the true canonical to the faux one found above
		output $1, trueCanon, perm
	    else // it's a new true canonical
		store _canonicalGraph[_numCanon] and _canonicalSig[_numCanon]
		++_numCanon
		output line verbatim
	    end if
	else
	    assert $2 is equal to an existing _fauxCanonSig[F]
	    output $1, _fauxCanonTrue[F], _fauxCanonPerm[F] composed with $3
	end if
    end while

**
*/


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

// We use TinySets since they can be stored in a char (8 bits), and thus the adjacency matrix for k<=8 fits into 8 chars.
static TINY_GRAPH *G, *_canonicalGraph[MAX_CANONICALS];
static int k, _numCanonicals, _canonicalSig[MAX_CANONICALS];

static int _numFauxCanon, _fauxCanonSig[MAX_CANONICALS], _fauxCanonTrueSig[MAX_CANONICALS];


/* The array perm[] will hold the permutation between the non-canonical and canonical.
** More formally, if perm[i]=j, it means that you should map node i from the non-canonical
** to node j of the canonical in order to prove the isomorphism. Then, when faye is doing
** samples, it means that if we print out the sampled nodes in the order listed in the perm
** array, then columns representing the same canonical present a perfect local alignment.
*/
static int perm[MAX_TSET], _fauxCanonPerm[MAX_CANONICALS][MAX_TSET];

int main(int argc, char *argv[])
{
    int i;
    if(argc != 4) {
	fprintf(stderr, "USAGE: %s k File1 File2\n", argv[0]);
	exit(1);
    }
    k = atoi(argv[1]);
    FILE *fp1 = fopen(argv[2], "r"), *fp2=fopen(argv[3], "r");
    G = TinyGraphAlloc(k);
    for(i=0; i<MAX_CANONICALS; i++)
	_canonicalGraph[i] = TinyGraphAlloc(k);
    int Gint, Gcanon; // maxGint = (1<<(k*(k-1)/2));
    int isGraphlet, numRead; // isGraphlet is ignored for now and just echo'd through.
    char line[BUFSIZ];

    // STEP 1: read first file
    
    while(!feof(fp1))
    {
	char cPerm[MAX_TSET+1];
	fgets(line, sizeof(line), fp1);
	numRead = sscanf(line, "%d%d%s%d", &Gint, &Gcanon, cPerm, &isGraphlet);
	assert(strlen(cPerm) <= MAX_TSET && strlen(cPerm)==k);
	if(Gint == Gcanon) // it's a true canonical
	{
	    _canonicalSig[_numCanonicals] = Gint;
	    BuildGraph(G, Gint);
	    TinyGraphCopy(_canonicalGraph[_numCanonicals], G);
	    ++_numCanonicals;
	}
	printf("%d\t%d\t%s",Gint, Gcanon, cPerm);
	if(numRead == 4) printf(" %d\n", isGraphlet);
	else puts("");
	fscanf(fp1, " ");
    }
    fclose(fp1);

    // STEP 2: read the second file
    while(!feof(fp2))
    {
	int j;
	char cPerm[MAX_TSET+1];
	fgets(line, sizeof(line), fp2);
	numRead = sscanf(line, "%d%d%s%d", &Gint, &Gcanon, cPerm, &isGraphlet);
	assert(strlen(cPerm) <= MAX_TSET && strlen(cPerm)==k);
	if(Gint == Gcanon) // "new" canonical, possibly faux
	{
	    BuildGraph(G, Gint);
	    for(i=0; i<_numCanonicals; i++)
	    {
#if PERMS_CAN2NON // the permutation provided is from canonical to non-canonical
		if(TinyGraphsIsomorphic(perm, _canonicalGraph[i], G))
#else // the permutation provided is from non-canonical to canonical
		if(TinyGraphsIsomorphic(perm, G, _canonicalGraph[i]))
#endif
		{
		    _fauxCanonSig[_numFauxCanon] = Gint;
		    _fauxCanonTrueSig[_numFauxCanon] = _canonicalSig[i];
		    printf("%d\t%d\t", Gint, _canonicalSig[i]);
		    for(j=0; j<k; j++)
		    {
			_fauxCanonPerm[_numFauxCanon][j] = perm[j];
			printf("%d", perm[j]);
		    }
		    if(numRead == 4) printf(" %d\n", isGraphlet);
		    else puts("");
		    ++_numFauxCanon;
		    break;
		}
	    }
	    if(i == _numCanonicals) // it's a real and true new canonical
	    {
		printf("%d\t%d\t", Gint, Gint);
		    for(j=0; j<k; j++) printf("%d", j);
		if(numRead == 4) printf(" %d\n", isGraphlet);
		else puts("");
		_canonicalSig[_numCanonicals] = Gint;
		TinyGraphCopy(_canonicalGraph[_numCanonicals], G);
		++_numCanonicals;
	    }
	}
	else // this is equivalent to an existing canonical, real or faux we're not sure yet.
	{
	    for(i=0; i<_numFauxCanon; i++)
	    {
		if(Gcanon == _fauxCanonSig[i])
		    break;
	    }
	    if(i < _numFauxCanon) // it's equivalent to a previously found faux canonical
	    {
		printf("%d\t%d\t", Gint, _fauxCanonTrueSig[i]);
		//for(j=0; j<k; j++) printf("%d", _fauxCanonPerm[i][cPerm[j]-'0']);
		for(j=0; j<k; j++) printf("%c", cPerm[_fauxCanonPerm[i][j]]);
		if(numRead == 4) printf(" %d\n", isGraphlet);
		else puts("");
	    }
	    else // it's equiv to the NEW canoical that was in the second file, so just echo
	    {
		printf("%d\t%d\t%s", Gint, Gcanon, cPerm);
		if(numRead == 4) printf(" %d\n", isGraphlet);
		else puts("");
	    }
	}
	fscanf(fp2, " ");
    }
    fclose(fp2);

    return 0;
}
