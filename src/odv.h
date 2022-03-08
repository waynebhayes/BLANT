#ifndef BLANT_ODV_HEADER
#define BLANT_ODV_HEADER

void parseOdvFromFile(char* fname);
void freeOdvInfo();

void getOdvValues(double* heuristicVals, int orbitNumber, char** nodeNames, int nodes);

#endif
