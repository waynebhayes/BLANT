#ifndef BLANT_ODV_HEADER
#define BLANT_ODV_HEADER

#define __ODV_N_ORBITS 15

typedef struct {
    char* nodeName;
    double odvValues[__ODV_N_ORBITS];
} odvrow_t;

void parseOdvFromFile(char* fname);
void freeOdvInfo();

void getOdvValues(double* heuristicVals, int orbitNumber, char** nodeNames);

#endif
