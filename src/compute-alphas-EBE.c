#include "combin.h"
#include "tinygraph.h"
#include "blant.h"
#include "misc.h"
#include <stdio.h>

Gint_type _canonList[MAX_CANONICALS], _alphaList[MAX_CANONICALS];
char _canonNumEdges[MAX_CANONICALS];

Gint_type CountPath(TINY_GRAPH* g, TSET seen, TSET candidates, int k) {
    seen = TSetUnion(seen, candidates);

    Gint_type outbound[((int)(k+1)/2)*((int)k/2)], outboundSize = 0;
    if (TSetCardinality(seen) == k) return 1;

    //Fill outbound for seen set
    int i, j;
    for (i = 0; i < k; i++) {
        if (!TSetIn(seen, i)) continue;

        for (j = 0; j < k; j++) {
            if (i != j && !TSetIn(seen, j) && TinyGraphAreConnected(g, i, j)) {
                outbound[outboundSize++] = j;
            }
        }
    }

    Gint_type calculatedValue[k];
    for (i = 0; i < k; i++) calculatedValue[i] = 0;

    //Recurse through neighbors in outbound
    Gint_type total = 0;
    for (i = 0; i < outboundSize; i++) {
        if (!calculatedValue[outbound[i]]) {
            TSetEmpty(candidates);
            TSetAdd(candidates, outbound[i]);
            calculatedValue[outbound[i]] = CountPath(g, seen, candidates, k);
        }
        total += calculatedValue[outbound[i]];
    }
    return total;
}

Gint_type ComputeAlphaNode(TINY_GRAPH* g, int k) {
    Gint_type total = 0, i, j;
    TSET seen = 0, candidates;
    for (i = 0; i < k; i++) {
        for (j = i+1; j < k; j++) {
            if (!TinyGraphAreConnected(g, i, j))
                continue;

            TSetEmpty(candidates);
            TSetAdd(candidates, i);
            TSetAdd(candidates, j);

            total += CountPath(g, seen, candidates, k);
            TSetEmpty(candidates);
        }
    }
    return total;
}

SET* _connectedCanonicals;
int main(int argc, char* argv[]) {
    if (argc != 2 && argc != 3 && argc != 4) {
        fprintf(stderr, "USAGE: %s k\nOr,  %s k ID\nOr,  %s k start end [inclusive]\n", argv[0], argv[0], argv[0]);
        exit(-1);
    }
    int k = atoi(argv[1]);
    char BUF[BUFSIZ];
    TINY_GRAPH *g = TinyGraphAlloc(k, false, false);
    _connectedCanonicals = canonListPopulate(BUF, _canonList, k, _canonNumEdges);
    Gordinal_type numCanon = _connectedCanonicals->maxElem;

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
    if (!SetIn(_connectedCanonicals, i)) _alphaList[i] = 0;
    else _alphaList[i] = ComputeAlphaNode(g, k);
    }

    printf(GORDINAL_FMT "\n", numCanon);
    for (i = 0; i < numCanon; i++) {
	printf(GINT_FMT " ", _alphaList[i]);
    }
    printf("\n");

    TinyGraphFree(g);
    SetFree(_connectedCanonicals);

    return 0;
}
