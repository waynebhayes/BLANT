#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
/*
    This file assumes LOWER_TRIANGLE is defined in blant.h
*/

//Functions from libwayne and libblant.c
extern "C" {
    struct TINY_GRAPH;
    void mapCanonMap(char* BUF, short int *K, int k);
    void *Mmap(void *p, size_t n, int fd);
    int canonListPopulate(char *BUF, int *canon_list, int k);
    void BuildGraph(TINY_GRAPH* G, int Gint);
    TINY_GRAPH *TinyGraphAlloc(unsigned int n);
    int TinyGraphBFS(TINY_GRAPH *G, int seed, int distance, int *nodeArray, int *distArray);
    typedef unsigned char Boolean;
    Boolean TinyGraphDFSConnected(TINY_GRAPH *G, int seed);
    void Fatal(const char *fmt, ...);
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
using std::unordered_set;

#define maxK 8
#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map
#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8

//Columns
const int CONNECTED = 0;
const int UPPER_ORDINAL = 1;
const int UPPER_DECIMAL = 2;
const int LOWER_DECIMAL = 3;
const int LOWER_ORDINAL = 4;
const int NUM_ORBITS = 5;
const int NUM_NODES_FIRST_ORBIT = 6;
const int ORCA = 7;
const int FIRST_ODV_ORBIT_CON = 8;
const int FIRST_ODV_ORBIT_ALL = 9;
const int FIRST_ODV_ORBIT_CON_FAYE = 10;
const int FIRST_ODV_ORBIT_ALL_FAYE = 11;

const int TABLE_WIDTH = 12;

const int MIN_K = 3;
const int MAX_K = 7;

static int _numCanon, _canonList[MAX_CANONICALS];
static int _numCanonU, _canonListU[MAX_CANONICALS];
static short int _K[maxBk] __attribute__ ((aligned (8192)));
static TINY_GRAPH* G;
const char* DIR = "orca_jesse_blant_table/";

//Manually generated mapping from upper(FAYE) binary representation of connected canonical graphlets to GRAAL/Przulj numbering
const auto umap3 = unordered_map<uint64_t,uint64_t>{{3, 1}, {7, 2}};
const auto umap4 = unordered_map<uint64_t, uint64_t>{{11, 4}, {13, 3}, {15, 6}, {30, 5}, {31, 7}, {63, 8}};
const auto umap5 = unordered_map<uint64_t, uint64_t>{{75, 11}, {77, 10}, {79, 14}, {86, 9}, {87, 12},
 {94, 16}, {95, 17}, {117, 13}, {119, 19}, {127, 23}, {222, 20}, {223, 22}, {235, 18}, {236, 15}, {237, 21},
  {239, 24}, {254, 25}, {255, 26}, {507, 27}, {511, 28}, {1023, 29}};

int canonListPopulate(char *BUF, int *canon_list, int k, char c) {
    stringstream ss;
    ss << DIR << "/canon_list" << c << k << ".txt";
    ifstream infile;
    infile.open(ss.str());
    if(!infile) {
        cerr << "Cannot open: " << ss.str() << "\n";
    }
    int numCanon, i = 0;
    infile >> numCanon;
    for (int i = 0; i < numCanon; i++) {
        infile >> canon_list[i];
    }
    infile.close();
    return numCanon;
}

//assuming little endian? bad filling matrix.
//TODO make better copying functions in make_orbit_maps
uint64_t Upper2Lower(uint64_t Gint, int k) {
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

//Sort table predicate on lower/BLANT ordinal
bool sortLower(const vector<uint64_t>& v1, const vector<uint64_t>& v2 ) {
    return v1[LOWER_ORDINAL] < v2[LOWER_ORDINAL];
}

//Sort table predicate on upper/FAYE ordinal
bool sortUpper(const vector<uint64_t>& v1, const vector<uint64_t>& v2 ) {
    return v1[UPPER_ORDINAL] < v2[UPPER_ORDINAL];
}

//Prints table to output stream
ostream& operator<<(ostream& os, const vector<vector<uint64_t>> table) {
    for (auto row : table) {
        for (auto num : row) {
            os << num << ' ';
        }
        os << '\n';
    }
    return os;
}

int main(int argc, char* argv[]) {
    //Connected count starts at 1 because k=2 has a connected node
    int connectedCount = 1;
    for (int k = MIN_K; k <= MAX_K; k++) {
        auto orbitTable = vector<vector<int>>();
        auto orbitTableUpper = vector<vector<int>>();

        stringstream ss;
        ss << DIR << "/UpperToLower" << k << ".txt";
        ofstream outfile;
        outfile.open(ss.str());
        if (!outfile) {
            cerr << "Failed to open outputfile: " << ss.str() << '\n';
        }
        ifstream orbitInfile;

        if (G) {
            free(G);
        }
        G = TinyGraphAlloc(k);

        //Load canon_list and canon_map
        char BUF[BUFSIZ];
        _numCanon = canonListPopulate(BUF, _canonList, k, 'l');
        _numCanonU = canonListPopulate(BUF, _canonListU, k, 'u');
        mapCanonMap(BUF, _K, k);
        if (_numCanon != _numCanonU) {
            perror("Num canons not equal\n");
            exit(-1);
        }

        //Create table
        auto table = vector<vector<uint64_t>>(_numCanon, vector<uint64_t>(TABLE_WIDTH, 0));
        int lowerDecimal;

        //Fill table
        for (int i = 0; i < table.size(); i++) {
            table[i][UPPER_ORDINAL] = i;
            table[i][UPPER_DECIMAL] = _canonListU[i];
            lowerDecimal = Upper2Lower(table[i][UPPER_DECIMAL], k);
            table[i][LOWER_ORDINAL] = _K[lowerDecimal];
            table[i][LOWER_DECIMAL] = _canonList[table[i][LOWER_ORDINAL]];

            BuildGraph(G, table[i][LOWER_DECIMAL]);
            if (TinyGraphDFSConnected(G, 0)) {
                table[i][CONNECTED] = 1;
            } else {
                table[i][CONNECTED] = 0;
            }
        }

        //Load num nodes first orbit information
        ss.str("");
        ss << DIR << "num_nodes_first_orbit" << k << ".txt";
        orbitInfile.open(ss.str());
        if (!orbitInfile) {
            cerr << "Failed to open: " << ss.str() << "\n";
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < table.size(); i++) {
            orbitInfile >> table[i][NUM_NODES_FIRST_ORBIT];
        }
        orbitInfile.close();

        //Load upper orbit information and fill out table with it
        unordered_set<int> orbits;
        ss.str("");
        ss << DIR << "orbit_mapu" << to_string(k) << ".txt";
        orbitInfile.open(ss.str());
        if (!orbitInfile) {
            cerr << "Failed to open: " << ss.str() << "\n";
            exit(EXIT_FAILURE);
        }
        int num, i = 0;
        orbitInfile >> num;
        while (orbitInfile) {
            orbitTableUpper.push_back(vector<int>(k, 0));
            for (int j = 0; j < k; j++) {
                orbitInfile >> orbitTableUpper[i][j];
            }
            i++;
        }
        orbitInfile.close();
        orbitTableUpper.pop_back();

        if (orbitTableUpper.size() != table.size()) {
            cerr << "k: " << k << " Orbit table upper size: " << orbitTableUpper.size() << " Table: " << table.size() << '\n';
            exit(-1);
        }

        //Calculate orbit information from upper
        int numConnectedOrbits = 0;
        int numTotalOrbits = 0;
        for (int i = 0; i < orbitTableUpper.size(); i++) {
            orbits.clear();
            for (int j = 0; j < k; j++) {
                orbits.insert(orbitTableUpper[i][j]);
            }
            table[i][NUM_ORBITS] = orbits.size();
            if (table[i][0]) {
                table[i][FIRST_ODV_ORBIT_CON_FAYE] = numConnectedOrbits;
                numConnectedOrbits += table[i][NUM_ORBITS];
            }
            table[i][FIRST_ODV_ORBIT_ALL_FAYE] = numTotalOrbits;
            numTotalOrbits += table[i][NUM_ORBITS];
        }

        //Process in Upper Sorting
        for (int i = 0; i < table.size(); i++) {
            //table[i][NUM_ORBITS] = orbitTable[i][0];
            orbits.clear();
            for (int j = 0; j < orbitTableUpper[i].size(); j++) {
                orbits.insert(orbitTableUpper[i][j]);
            }
            //std::cout << "NUM_ORBITS: " << orbits.size() << '\n';
            table[i][NUM_ORBITS] = orbits.size();
        }

        //Begin Lower sorting things
        sort(table.begin(), table.end(), sortLower);

        //Load lower orbit information and compare with upper
        numConnectedOrbits = 0;
        numTotalOrbits = 0;
        i = 0;
        ss.str("");
        ss << DIR << "orbit_map" << k << ".txt";
        orbitInfile.open(ss.str());
        if (!orbitInfile) {
            cerr << "Failed to open: " << ss.str() << "\n";
            exit(EXIT_FAILURE);
        }
        orbitInfile >> num;
        while (orbitInfile) {
            orbitTable.push_back(vector<int>(k, 0));
            for (int j = 0; j < k; j++) {
                orbitInfile >> orbitTable[i][j];
            }
            i++;
        }
        orbitTable.pop_back();
        orbitInfile.close();
        if (orbitTable.size() != table.size()) {
            cerr << "Orbit table size: " << orbitTable.size() << " Table: " << table.size() << '\n';
            exit(EXIT_FAILURE);
        }

        //Process in lower
        for (int i = 0; i < orbitTable.size(); i++) {
            orbits.clear();
            for (int j = 0; j < k; j++) {
                orbitInfile >> orbitTable[i][j];
                orbits.insert(orbitTable[i][j]);
            }
            if (orbits.size() != table[i][NUM_ORBITS]) {
                cerr << k << '\n';
            }

            if (table[i][0]) {
                table[i][FIRST_ODV_ORBIT_CON] = numConnectedOrbits;
                numConnectedOrbits += table[i][NUM_ORBITS];
            }
            table[i][FIRST_ODV_ORBIT_ALL] = numTotalOrbits;
            numTotalOrbits += table[i][NUM_ORBITS];

            if (table[i][0]) {
                if (k <= 5) {
                    if (k == 3) {
                        table[i][ORCA] = umap3.at(table[i][UPPER_DECIMAL]);
                    } else if (k == 4) {
                        table[i][ORCA] = umap4.at(table[i][UPPER_DECIMAL]);
                    } else {
                        table[i][ORCA] = umap5.at(table[i][UPPER_DECIMAL]);
                    }
                } else {
                    table[i][ORCA] = connectedCount;
                }
                connectedCount++;
            } else {
                table[i][ORCA] = 0;
            }
        }

        //output
        sort(table.begin(), table.end(), sortUpper);
        outfile << table;
        outfile.close();
    }
    return 0;
}