#include "odv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __ODV_COL_DELIMSTR " "
#define __ODV_COL_DELIM ' '
#define __ODV_LINE_DELIM '\n'

static odvrow_t* odvdata;
static int odvdataCapacity;
static int nOdvRows;

void parseOdvRow(odvrow_t* row, char* line) {
    row->nodeName = strtok(line, __ODV_COL_DELIMSTR);
    
    for (int orbitNumber = 0; orbitNumber < __ODV_N_ORBITS; orbitNumber++) {
        char* columnStr = strtok(NULL, __ODV_COL_DELIMSTR);
        double odvVal = atoi(columnStr);

        row->odvValues[orbitNumber] = odvVal;
    }
}

void freeOdvRow(odvrow_t* row) {
    free(row->nodeName);
}

void parseOdvFromFile(char* fname) {
    // open file
    FILE* odvFile = fopen(fname, "r");
    // TODO: fopen error checking

    // init odvdata
    nOdvRows = 0;
    odvdataCapacity = 1000; // init buffer to 1000 rows
    odvdata = malloc(sizeof(odvrow_t) * odvdataCapacity);

    // TODO: malloc error checking

    // read one line at a time, parse into odvrow_t, stick in odvdata
    while (!feof(odvFile)) {
        // lazily expand odvdata if necessary
        if (nOdvRows == odvdataCapacity) {
            odvdataCapacity *= 2;
            odvdata = realloc(odvdata, sizeof(odvrow_t) * odvdataCapacity); // TODO: error checking
        }

        odvrow_t* row = odvdata + nOdvRows;

        int nodeNameCapacity = 16;
        int nodeNameSize = 0;
        char* nodeName = calloc(sizeof(char), nodeNameSize);
        // TODO: calloc error checking

        char c = fgetc(odvFile);

        while (c != __ODV_COL_DELIM && c != EOF) {
            if (nodeNameSize == nodeNameCapacity) {
                nodeNameCapacity *= 2;
                nodeName = realloc(nodeName, sizeof(char) * nodeNameCapacity);
                // TODO: alloc error checking
            }

            nodeName[nodeNameSize++] = c;
            c = fgetc(odvFile);
        }

        // do this to stick a null character at the end
        if (nodeNameSize == nodeNameCapacity) {
            nodeNameCapacity++;
            nodeName = realloc(nodeName, sizeof(char) * nodeNameCapacity);
            // TODO: alloc error checking
        }

        int i = 0;
        double currentOdvVal = 0.0;
        c = fgetc(odvFile);

        while (c != __ODV_LINE_DELIM && c != EOF) {
            if (c != __ODV_COL_DELIM) {
                int digit = c - '0'; // i assume that the file format will be correct here
                currentOdvVal = 10 * currentOdvVal + digit;
            } else {
                row->odvValues[i++] = currentOdvVal;
                currentOdvVal = 0.0;
            }
            
            c = fgetc(odvFile);
        }

        row->odvValues[i++] = currentOdvVal;
        row->nodeName = nodeName;
        nOdvRows++;
    }
}

void freeOdvData() {
    for (int i = 0; i < nOdvRows; i++) {
        freeOdvRow(odvdata + i);
    }

    free(odvdata);
}

void getOdvValues(double* heuristicVals, int orbitNumber, char** nodeNames) {

}

// TODO: remove this, just a simple function to prove that it works
int main(int argc, char** argv) {
    char* odvfile = argv[1];
    parseOdvFromFile(odvfile);

    for (int i = 0; i < nOdvRows; i++) {
        odvrow_t* row = odvdata + i;
        printf("%i: %s ", i, row->nodeName);

        for (int j = 0; j < __ODV_N_ORBITS; j++) {
            printf("%d ", (int)(row->odvValues[j]));
        }

        printf("\n");
    }

    freeOdvData();

    return 0;
}
