// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
//
// Debug harness: run L_K_Func_Sort over every canonical read from a canon list
// file and do nothing else. The canonical gints are taken from the file rather
// than rediscovered, so no scan over all graphlets is done. Each canonical is
// canonicalized by L_K_Func_Sort itself (so file entries from the old
// CANON_ASCENDING_NEIGHBORS=0 definition are converted to their current canonical
// decimal on the fly). When compiled with -DDEBUG_ATTEMPTS the per-call attempt
// count is appended to attlog.txt.
//
// Build (not part of the normal build; only on explicit request):
//     make l-k-func-sort-test
// Usage:
//     ./l-k-func-sort-test k [canon_list_file]
// If canon_list_file is omitted, canon_maps/canon_list{k}.txt is used. The file
// format is the BLANT canon_list format: first line is the number of canonicals,
// then one line per canonical with the gint as the first column.
#include "blant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

unsigned int _k;      // normally defined in blant.c, which this harness does not link
Boolean _directed;    // ditto

// L_K_Func_Sort is defined in blant-utils.c but declared in no header.
Gint_type L_K_Func_Sort(Gint_type Gint, unsigned char permOut[], unsigned short olist[], Boolean computeOrbits);

int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "USAGE: %s k [canon_list_file]\n", argv[0]);
        return -1;
    }
    int k = atoi(argv[1]);
    if (k < 1 || k > MAX_K) {
        fprintf(stderr, "k must be in [1,%d]\n", MAX_K);
        return -1;
    }
    _k = k;
    _directed = false;

    char fname[BUFSIZ];
    if (argc == 3) {
        snprintf(fname, sizeof(fname), "%s", argv[2]);
    } else {
        snprintf(fname, sizeof(fname), "%s/canon_list%d.txt", DEFAULT_CANON_DIR, k);
    }
    FILE *fp = fopen(fname, "r");
    if (!fp) {
        fprintf(stderr, "cannot open canon list file %s\n", fname);
        return -1;
    }

    char line[BUFSIZ];
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "failed to read numCanon from %s\n", fname);
        return -1;
    }
    unsigned long long numCanon = strtoull(line, NULL, 10);
    if (numCanon == 0) {
        fprintf(stderr, "failed to read numCanon from %s\n", fname);
        return -1;
    }
    Gint_type *canon = Malloc((Gordinal_type)numCanon * sizeof(Gint_type));
    for (Gordinal_type i = 0; i < numCanon; i++) {
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "unexpected EOF in %s\n", fname);
            return -1;
        }
        char *tok = strtok(line, " \t");
        if (!tok) {
            fprintf(stderr, "malformed line in %s: %s", fname, line);
            return -1;
        }
        canon[(size_t)i] = (Gint_type)strtoull(tok, NULL, 10);
    }
    fclose(fp);

    for (Gordinal_type i = 0; i < numCanon; i++) {
        (void)L_K_Func_Sort(canon[(size_t)i], NULL, NULL, 0);
        if ((i + 1) % 100 == 0 || i + 1 == numCanon)
            fprintf(stderr, "checked %llu / %llu canonicals\n",
                (unsigned long long)(i + 1), (unsigned long long)numCanon);
    }
    Free(canon);
    return 0;
}
