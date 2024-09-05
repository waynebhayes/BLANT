#ifndef WHAT_THE_GRAPHLET_H
#define WHAT_THE_GRAPHLET_H

#include "blant.h"



#if MAX_K == 9
#define SMALL_MAP_SIZE 3160576
#elif MAX_K == 8
#define SMALL_MAP_SIZE 133632
#else
#define SMALL_MAP_SIZE 8988
#endif


#define CANON_MAP_FOLDER "canon_maps/"


extern Gint_type smaller_canon_map(Gint_type num, int k, unsigned char* return_permutation);
extern void read_smaller_canon_maps(int max_k);
extern Gordinal_type canon_to_ordinal(Gint_type canon, int k);


extern Gint_type map_non_canon[MAX_K+1][SMALL_MAP_SIZE];


#endif