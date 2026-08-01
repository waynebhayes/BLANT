// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
//
// Debug harness: for each canonical gint given on the command line, enumerate all
// k! vertex relabelings of the graphlet, which yields every decimal in its
// isomorphism class, and run L_K_Func_Sort on each distinct decimal. This checks
// whether the number of permutations L_K_Func_Sort searches (the DEBUG_ATTEMPTS
// counter, read back via the global "attempts") depends on which decimal of the
// class it starts from, or only on the isomorphism class itself.
//
// Build (not part of the normal build; only on explicit request):
//     make l-k-perm-test
// Usage:
//     ./l-k-perm-test k canon_gint [canon_gint ...]
#include "blant.h"
#include <stdio.h>
#include <stdlib.h>

unsigned int _k;      // normally defined in blant.c, which this harness does not link
Boolean _directed;    // ditto

#if DEBUG_ATTEMPTS
extern unsigned long attempts;      // defined in blant-utils.c under DEBUG_ATTEMPTS
extern Boolean logAttemptsToFile;   // ditto
#else
#error l-k-perm-test requires -DDEBUG_ATTEMPTS
#endif

// L_K_Func_Sort is defined in blant-utils.c but declared in no header.
Gint_type L_K_Func_Sort(Gint_type Gint, unsigned char permOut[], unsigned short olist[], Boolean computeOrbits);

static Gint_type decs[1 << 22];          // distinct decimals found so far
static unsigned long decsAtt[1 << 22];   // L_K_Func_Sort attempt count for each
static int numDecs;
static TINY_GRAPH *gCanon;

// Build the relabeled graph (node i of the new graph = node perm[i] of gCanon)
// and return its decimal.
static Gint_type relabelDecimal(const int perm[]) {
    TINY_GRAPH g2;
    memset(&g2, 0, sizeof(g2));
    g2.n = _k;
    g2.directed = gCanon->directed;
    g2.selfLoops = gCanon->selfLoops;
    for (int i = 0; i < _k; i++)
        for (int j = i + 1; j < _k; j++)
            if (TinyGraphAreConnected(gCanon, perm[i], perm[j]))
                TinyGraphConnect(&g2, i, j);
    return TinyGraph2Int(&g2, _k);
}

static void genPerms(int pos, int *perm, Boolean *used) {
    if (pos == _k) {
        Gint_type d = relabelDecimal(perm);
        for (int i = 0; i < numDecs; i++)
            if (decs[i] == d) return;   // already tested this decimal
        decs[numDecs] = d;
        Gint_type canon = L_K_Func_Sort(d, NULL, NULL, 0);
        if (canon != decs[0])
            Fatal("decimal %llu canonicalizes to %llu, expected %llu",
                (unsigned long long)d, (unsigned long long)canon, (unsigned long long)decs[0]);
        decsAtt[numDecs] = attempts;
        numDecs++;
        return;
    }
    for (int v = 0; v < _k; v++) {
        if (used[v]) continue;
        perm[pos] = v;
        used[v] = true;
        genPerms(pos + 1, perm, used);
        used[v] = false;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "USAGE: %s k canon_gint [canon_gint ...]\n", argv[0]);
        return -1;
    }
    int k = atoi(argv[1]);
    if (k < 1 || k > 9) {   // k! gets out of hand beyond k=9
        fprintf(stderr, "k must be in [1,9]\n");
        return -1;
    }
    _k = k;
    _directed = false;
    logAttemptsToFile = false;   // we read "attempts" directly instead

    long long factorial = 1;
    for (int i = 2; i <= k; i++) factorial *= i;
    fprintf(stderr, "enumerating %lld relabelings per canonical\n", factorial);

    for (int a = 2; a < argc; a++) {
        Gint_type canonGint = (Gint_type)strtoull(argv[a], NULL, 10);
        TINY_GRAPH g;
        memset(&g, 0, sizeof(g));
        g.n = _k;
        g.directed = _directed;
        g.selfLoops = 0;
        Int2TinyGraph(&g, canonGint);
        gCanon = &g;
        numDecs = 0;
        decs[0] = L_K_Func_Sort(canonGint, NULL, NULL, 0);
        if (decs[0] != canonGint)
            Warning("input %llu is not canonical; its canonical is %llu",
                (unsigned long long)canonGint, (unsigned long long)decs[0]);
        decsAtt[0] = attempts;
        numDecs = 1;

        int perm[MAX_K];
        Boolean used[MAX_K] = { false };
        genPerms(0, perm, used);

        unsigned long minAtt = decsAtt[0], maxAtt = decsAtt[0];
        int nVarying = 0;
        for (int i = 1; i < numDecs; i++) {
            if (decsAtt[i] < minAtt) minAtt = decsAtt[i];
            if (decsAtt[i] > maxAtt) maxAtt = decsAtt[i];
            if (decsAtt[i] != decsAtt[0]) nVarying++;
        }
        printf("canonical %llu: %d distinct decimals in class | canonical=%lu attempts | "
            "min=%lu max=%lu | %d other decimals vary (%s)\n",
            (unsigned long long)decs[0], numDecs, decsAtt[0], minAtt, maxAtt,
            nVarying, nVarying ? "VARY" : "all equal");
        for (int i = 1; i < numDecs && i <= 10; i++)
            if (decsAtt[i] != decsAtt[0])
                printf("    decimal %llu -> %lu attempts\n",
                    (unsigned long long)decs[i], decsAtt[i]);
    }
    return 0;
}
