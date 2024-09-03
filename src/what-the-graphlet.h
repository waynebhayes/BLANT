#ifndef WHAT_THE_GRAPHLET_H
#define WHAT_THE_GRAPHLET_H

#include "blant.h"


#define MAXK 9

#if MAXK == 9
#define MAX_SIZE 3160576
#define MAX_CANON 274668
#elif MAXK == 8
#define MAX_SIZE 133632
#define MAX_CANON 12346
#else
#define MAX_SIZE 8988
#define MAX_CANON 1044
#endif


#define CANON_MAP_FOLDER "canon_maps/"


extern void smaller_canon_map(Gint_type num, int k, Gint_type* return_canon, unsigned char* return_permutation);
extern void read_maps(int max_k);
extern Gordinal_type canon_to_ordinal(Gint_type canon, int k);


extern Gint_type map_non_canon[MAXK+1][MAX_SIZE];


#endif