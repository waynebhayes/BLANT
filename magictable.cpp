#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <memory.h>
#include <sys/mman.h>
#include <sys/time.h>   /*getrusage doesn't seem to work.*/
#include <sys/resource.h>
#include <unistd.h>
#include <algorithm>
#include <sstream>
#include <unordered_map>


extern "C"{
    void mapCanonMap(char* BUF, short int *K, int k);
    void *Mmap(void *p, size_t n, int fd);
    int canonListPopulate(char *BUF, int *canon_list, int k);
}

using std::sort;
using std::vector;
using std::unordered_map;
using std::to_string;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::ostream;

#define maxK 8
#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map
#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8

static int _numCanon, _canonList[MAX_CANONICALS];
static short int _K[maxBk] __attribute__ ((aligned (8192)));
const auto umap3 = unordered_map<uint64_t,uint64_t>{{3, 1}, {7, 2}};
const auto umap4 = unordered_map<uint64_t, uint64_t>{{11, 4}, {13, 3}, {15, 6}, {30, 5}, {31, 7}, {63, 8}};
const auto umap5 = unordered_map<uint64_t, uint64_t>{{75, 11}, {77, 10}, {79, 14}, {86, 9}, {87, 12},
 {94, 16}, {95, 17}, {117, 13}, {119, 19}, {127, 23}, {222, 20}, {223, 22}, {235, 18}, {236, 15}, {237, 21},
  {239, 24}, {254, 25}, {255, 26}, {507, 27}, {511, 28}, {1023, 29}};

//assuming little endian? bad filling matrix
uint64_t Upper2Lower(uint64_t Gint, int k)
{
    auto matrix = vector<vector<bool>>(k, vector<bool>(k, false));
    int j = k-1, i = k-2;
    uint64_t mask = 1, result = 0;
    int numBits = k*(k-1)/2;
    for (int bit = 0; bit < numBits; bit++) {
        if (Gint & mask) {
            matrix[j][i] = true;
            matrix[i][j] = true;
        }
        mask <<= 1;
        j--;
		if (j == i) {
			i--;
			j = k-1;
		}
    }
    mask = 1;
    j = k-2;
    i = k-1;
    for (int bit = 0; bit < k*(k-1)/2; bit++) {
        if (matrix[j][i]) {
            result |= mask;
        }
        mask <<= 1;
        j--;
        if (j == -1) {
            i--;
            j = i
            -1;
        }
    }
    return result;
}

bool sortcol5(const vector<uint64_t>& v1, const vector<uint64_t>& v2 ) {
    return v1[4] < v2[4];
}

bool sortcol2(const vector<uint64_t>& v1, const vector<uint64_t>& v2 ) {
    return v1[1] < v2[1];
}

ostream& operator<<(ostream& os, const vector<vector<uint64_t>> table) {
    for (auto row : table) {
        for (auto num : row) {
            os << num << ' ';
        }
        os << '\n';
    }
    return os;
}

//first, check 4th column matches blant canonical
//2nd, fix 5th column blant ordinale
//3rd do k >= 5 jesse
//4th do k < 5 jesse

int main(int argc, char* argv[]) {
    int connectedCount = 1;
    for (int k = 3; k <= 7; k++) {
        //load into memory
        auto table = vector<vector<uint64_t>>();
        ofstream outfile;
        outfile.open("newTransformations/UpperToLower" + to_string(k) + ".txt");
        ifstream infile;
        infile.open("transformations/UpperToLower" + to_string(k) + ".txt");
        if (!infile) {
            cerr << "Failed to open file\n";
            exit(EXIT_FAILURE);
        }
        char BUF[BUFSIZ];
        _numCanon = canonListPopulate(BUF, _canonList, k);
        mapCanonMap(BUF, _K, k);

        int j = 0;
        while (infile) {
            uint64_t temp;
            table.push_back(vector<uint64_t>());
            for (int i = 0; i < 7; i++) {
                infile >> temp;
                table[j].push_back(temp);
            }
            table[j].push_back(0);
            j++;
        }
        infile.close();
        table.pop_back();

        sort(table.begin(), table.end(), sortcol5);
        //Process
        for (int i = 0; i < table.size(); i++) {
            uint64_t lower, lowerOrdinal, lowerCanonical;
            lower = Upper2Lower(table[i][2], k);
            lowerOrdinal = _K[lower];
            lowerCanonical = _canonList[lowerOrdinal];

            table[i][3] = lowerCanonical;
            table[i][4] = lowerOrdinal;

            table[i][7] = 0;
            if (table[i][0]) {
                if (k <= 5) {
                    if (k == 3) {
                        table[i][7] = umap3.at(table[i][2]);
                    } else if (k == 4) {
                        table[i][7] = umap4.at(table[i][2]);
                    } else {
                        table[i][7] = umap5.at(table[i][2]);
                    }
                } else {
                    table[i][7] = connectedCount;
                }
                connectedCount++;
            } else {
                table[i][7] = 0;
            }
        }   
        sort(table.begin(), table.end(), sortcol2);
        //output
        outfile << table;
        outfile.close();
    }
    return 0;
}