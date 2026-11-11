#ifndef MEASURES_H
#define MEASURES_H

#include "communities.h"
#include "structs.h"
double IntraEdgeDensity(PARTITION * P, COMMUNITY *C, int fakeN);
double InterEdgeDensity(PARTITION * P, COMMUNITY *C, int fakeN);
double NewmanAndGirvan(PARTITION * P, COMMUNITY *C, int fakeN);
double Conductance(PARTITION * P, COMMUNITY *C, int fakeN);
double HayesScore(PARTITION * P, COMMUNITY *C, int fakeN);
double Modularity(PARTITION * P, COMMUNITY *C, int fakeN);
double Expansion(PARTITION * P, COMMUNITY * C, int fakeN);
double NormalizedCut(PARTITION * P, COMMUNITY * C, int fakeN);
double CPM(PARTITION * P, COMMUNITY * C, int fakeN);
extern int _newCom, _oldCom;
#endif
