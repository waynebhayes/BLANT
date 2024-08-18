#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "what-the-graphlet.h"



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

unsigned long long _canonical_count[] =
	{0, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864};

unsigned long long map_ordinal[MAXK+1][MAX_SIZE];

unsigned long long map_canon[MAXK+1][MAX_SIZE];
char map_permutation[MAXK+1][MAX_SIZE][MAXK];





void read_maps(int max_k) {
    char filename[100];
    int cnt;
    char buf[BUFSIZ];

    for(int k = 3; k <= max_k; k++) {
        sprintf(filename, "%s%s%d.txt", CANON_MAP_FOLDER, "smaller_canon_map", k);
        FILE* file = fopen(filename, "r");

        cnt = 0;
        while(fgets(buf, sizeof(buf), file)) {
            sscanf(buf, "%llu\t%llu\t%s\t%*s", &map_ordinal[k][cnt], &map_canon[k][cnt], map_permutation[k][cnt]);
            cnt++;
        }
        fclose(file);
    }
}



int permute(int num, char* perm, int n) {
    int permuted = 0;
    for(int i = 0; i < n; i++) {
        permuted |= ((num >> (n - i - 1)) & 1) << (n - (perm[i] - '0') - 1);
    }
    return permuted;
}

int reverse_permute(int num, char* perm, int n) {
    int permuted = 0;
    for(int i = 0; i < n; i++) {
        permuted |= ((num >> (n - (perm[i] - '0') - 1)) & 1) << (n - i - 1);
    }
    return permuted;
}

void permutation_composition(char* perm1, char* perm2, int n, char* permuted) {
    for(int i = 0; i < n; i++) {
        permuted[i] = perm1[perm2[i] - '0'];
    }
    permuted[n] = '\0';
}






void smaller_canon_map(unsigned long long num, int k, unsigned long long* return_canon, char* return_permutation) {
    // The base case
    if(k == 3) {
        ///////// Binary search the index ////////
        int index = -1, low = 0, high = _canonical_count[2] << 2, mid;
        while (low <= high) {
            mid = (low + high) / 2;
            if (map_ordinal[3][mid] == num) {
                index = mid;
                break;
            } else if (map_ordinal[3][mid] < num) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        ///////// Binary search the index ////////

        *return_canon = map_canon[3][index];
        strcpy(return_permutation, map_permutation[3][index]);
        return;
    }





    // Recursion for > 3
    // Getting the canonical for first k-1 nodes and their permutation
    unsigned long long prev_canon;
    char prev_perm[MAXK];
    // #define prev_perm return_permutation // Reuse the memory, no need to allocate new memory
    smaller_canon_map(num >> (k-1), k-1, &prev_canon, prev_perm);


    // Preparing the last row
    unsigned int last_row = num & ((1ll << (k-1))-1);
    last_row = reverse_permute(last_row, prev_perm, k-1);



    // Permutation for the first transformation
    prev_perm[k - 1] = k - 1 + '0';
    prev_perm[k] = '\0';

    


    // Using the pre-calculated canonicals

    ///////// Binary search the index ////////
    int index = -1, low = 0, high = _canonical_count[k - 1] << (k - 1), mid;
    while (low <= high) {
        mid = (low + high) / 2;
        if (map_ordinal[k][mid] == ((prev_canon << (k - 1)) | last_row)) {
            index = mid;
            break;
        } else if (map_ordinal[k][mid] < ((prev_canon << (k - 1)) | last_row)) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    ///////// Binary search the index ////////
    


    *return_canon = map_canon[k][index];
    // Because we want return_permutation(X) == prev_perm(map_permutation[k][index](X))
    permutation_composition(prev_perm, map_permutation[k][index], k, return_permutation);

    return;
}




// int main(int argc, char* argv[]) {

//     int k = atoi(argv[1]);

//     read_maps(k);



//     for(int num = 0; num < (1 << (k * (k - 1) / 2)); num++) {
//         unsigned long long canon1;
//         char perm1[MAXK];
//         smaller_canon_map(num, k, &canon1, perm1);

//         printf("%d\t%llu\t%s\n", num, canon1, perm1);
//     }    

//     // unsigned long long canon1;
//     // char perm1[MAXK];
//     // int num = 33;

//     // smaller_canon_map(num, 4, &canon1, perm1);

//     // printf("%d\t%llu\t%s\n", num, canon1, perm1);



//     return 0;
// }