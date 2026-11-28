// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include "combin.h"
#include "tinygraph.h"
#include "blant.h"
#include "misc.h"
#include <stdio.h>

Gint_type *_canonList, *_alphaList;
char *_canonNumEdges;

// NBE counts each distinct frontier node once.
Gint_type ComputeAlphaNode(TINY_GRAPH* g, int k) {
    return ComputeExpansionAlpha(g, k, false);
}

SET* _connectedCanonicals;
int main(int argc, char* argv[]) {
    if (argc != 2 && argc != 3 && argc != 4) {
        fprintf(stderr, "USAGE: %s k\nOr,  %s k ID\nOr,  %s k start end [inclusive]\n", argv[0], argv[0], argv[0]);
        exit(-1);
    }
    int k = atoi(argv[1]);
    char BUF[BUFSIZ];
    _canonList = Malloc(MAX_CANONICALS * sizeof(Gint_type));
    _alphaList = Malloc(MAX_CANONICALS * sizeof(Gint_type));
    _canonNumEdges = Malloc(MAX_CANONICALS * sizeof(char));
    TINY_GRAPH *g = TinyGraphAlloc(k, SELF_LOOPS, false);
    _connectedCanonicals = canonListPopulate(BUF, _canonList, k, _canonNumEdges,false);
    int numCanon = _connectedCanonicals->maxElem;

	int start, end;
	if(argc == 2) {
		start = 0, end = numCanon - 1;
	} else if (argc == 3) {
		start = atoi(argv[2]), end = atoi(argv[2]);
		if(start < 0 || start >= numCanon) {
			fprintf(stderr, "Invalid ID\n");
			exit(-1);
		}
	} else {
		start = atoi(argv[2]), end = atoi(argv[3]);
		if(start < 0 || start >= numCanon || end < 0 || end >= numCanon || start > end) {
			fprintf(stderr, "Invalid range\n");
			exit(-1);
		}
	}

    int i;
    for (i = start; i <= end; i++) {
	    Int2TinyGraph(g, _canonList[i]);
    if (!SetIn(_connectedCanonicals, i))
	_alphaList[i] = 0;
	    else
	_alphaList[i] = ComputeAlphaNode(g, k);
    }

    printf(GORDINAL_FMT "\n", numCanon);
    for (i = 0; i < numCanon; i++) {
	    printf(GINT_FMT " ", _alphaList[i]);
    }
    printf("\n");

    TinyGraphFree(g);
    SetFree(_connectedCanonicals);
    Free(_canonList);
    Free(_alphaList);
    Free(_canonNumEdges);

    return 0;
}
