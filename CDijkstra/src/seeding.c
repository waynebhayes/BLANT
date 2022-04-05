#include "seeding.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void seed_from_file(seed_t* seed, char* fname, unsigned int line_number, GRAPH* g) {
    FILE* seed_file = fopen(fname, "r");
    if (seed_file == NULL) {
        return;
    }

    unsigned int lines_skipped = 0;

    while (lines_skipped < line_number - 1) {
        switch (fgetc(seed_file)) {
        case EOF:
            return;

        case '\n':
            lines_skipped++;
        }
    }

    char* buf = NULL;
    ssize_t bufsize;
    getline(&buf, &bufsize, seed_file);

    char* delim = strchr(buf, ' ');
    *delim = '\0';

    unsigned int kval = atoi(buf);
    unsigned int node = GraphNodeName2Int(g, delim + 1);

    seed->kval = kval;
    seed->node = node;
}
