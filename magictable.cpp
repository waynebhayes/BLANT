#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "blant.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

#define CANON_DIR "canon_maps"

static int _numCanon, _canonList[MAX_CANONICALS];

static short int _K[maxBk] __attribute__ ((aligned (8192)));

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array. 
*/
void mapCanonMap(char* BUF, short int *K, int k) {
    int Bk = (1 <<(k*(k-1)/2));
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    short int *Kf = Mmap(K, Bk*sizeof(K[0]), Kfd);
    assert(Kf == K);
}

uint64_t Upper2Lower(uint64_t Gint, int k)
{
    auto matrix = vector<vector<uint64_t>>(k, vector<uint64_t>(k, 0));
    unsigned j = 1, i = 0;
    uint64_t mask = 1, result;

    for (int bit = 0; bit < k*(k-1)/2; bit++) {
        if (Gint & mask) {
            matrix[j][i] = 1;
        }
        mask <<= 1;
        j++;
		if (j == k) {
			i++;
			j = i + 1;
		}
    }
    mask = 1;
    j = 0;
    i = 1;
    for (int bit = 0; bit < k*(k-1)/2; bit++) {
        if (matrix[j][i]) {
            result |= mask;
        }

        mask <<= 1;
        j++;
        if (j == i) {
            i++;
            j = 0;
        }

    }
}

//first, check 4th column matches blant canonical
//2nd, fix 5th column blant ordinale
int main(int argc, char* argv[]) {
    for (int k = 3; k <= 7; k++) {
        auto table = vector<vector<int>>();
        ofstream outfile;
        outfile.open("Table" + to_string(k) + ".txt");
        ifstream infile;
        infile.open("UpperToLower" + to_string(k) + ".txt");
        while (infile) {
            vector<int> temp;
            temp.reserve(7);
            for (int i = 0; i < 7; i++) {
                infile >> temp[i];
            }
            table.push_back(temp);
        }

        outfile.close();
        infile.close();


        //Testing
        for (int i = 0; i < table.size(); i++) {
            assert(table[i][3] == Upper2Lower(table[i][2], k));
        }
    }

    
    return 0;
}