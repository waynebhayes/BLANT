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

// extern "C"{
// typedef struct _tinyGraph;
// int TinyGraphBFS(TINY_GRAPH *G, int seed, int distance, int *nodeArray, int *distArray);
// TINY_GRAPH *TinyGraphAlloc(unsigned int n);
// #define TinyGraphFree free
// void BuildGraph(TINY_GRAPH* G, int Gint);
// int TinyGraph2Int(TINY_GRAPH *g, int numNodes);
// void mapCanonMap(char* BUF, short int *K, int k);
// int canonListPopulate(char *BUF, int *canon_list, int k);
// char** convertToEL(char* file); // from convert.cpp
// }

using namespace std;

#define maxK 8
#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map
#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8
#define MAX_ORBITS	79264	// This is the number of orbits for k=8
#define CANON_DIR "canon_maps"
#define MMAP 1

static int _numCanon, _canonList[MAX_CANONICALS];
static short int _K[maxBk] __attribute__ ((aligned (8192)));
const auto umap3 = unordered_map<uint64_t,uint64_t>{{3, 1}, {7, 2}};
const auto umap4 = unordered_map<uint64_t, uint64_t>{{7, 4}, {13, 3}, {15, 6}, {30, 5}, {31, 7}, {63, 8}};
const auto umap5 = unordered_map<uint64_t, uint64_t>{{15, 11}, {29, 10}, {31, 14}, {58, 9}, {59, 12},
 {62, 16}, {63, 17}, {126, 20}, {127, 22}, {185, 13}, {187, 19}, {191, 23}, {207, 18}, {220, 15}, {223, 24},
  {254, 25}, {255, 26}, {495, 27}, {511, 28}, {1023, 29}};

// Try to mmap, and if it fails, just slurp in the file (sigh, Windoze)
void *Mmap(void *p, size_t n, int fd)
{
#if MMAP
    void *newPointer = mmap(p, n, PROT_READ, MAP_PRIVATE|MAP_FIXED, fd, 0);
    if(newPointer == MAP_FAILED)
#endif
    {
#if !__WIN32__ && !__CYGWIN__ // it will always fail on Windoze so don't bother with a warning
        perror("mmap failed");
#endif
	int status;
	size_t numRead = 0;
        while(numRead < n && (status = read(fd, (char*)p + numRead, n-numRead)) > 0)
	    numRead += status;
        if(numRead < n) perror("cannot mmap, or cannot read the file, or both");
    }
    return p;
}

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array. 
*/
void mapCanonMap(char* BUF, short int *K, int k) {
    int Bk = (1 <<(k*(k-1)/2));
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    short int *Kf = (short int*) Mmap(K, Bk*sizeof(K[0]), Kfd);
    assert(Kf == K);
}

int canonListPopulate(char *BUF, int *canon_list, int k) {
    sprintf(BUF, CANON_DIR "/canon_list%d.txt", k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) perror("cannot find\n");
    int numCanon, i;
    fscanf(fp_ord, "%d",&numCanon);
    for(i=0; i<numCanon; i++) fscanf(fp_ord, "%d", &canon_list[i]);
    fclose(fp_ord);
    return numCanon;
}


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

            //Check if k<= 5 jesse
            // if (k <= 5 && table[i][0]) {
            //     stringstream ss;
            //     ss << "$(../SanaGV/graphette2dot -u -k " << k << " -d " << table[i][2] << " -t k" << k << "d" << table[i][3] << " -o k" << k << "d" << table[i][3] << ")\n";
            //     system(ss.str().c_str());
            // }
        }   
        sort(table.begin(), table.end(), sortcol2);
        //output
        outfile << table;
        outfile.close();
    }
    return 0;
}