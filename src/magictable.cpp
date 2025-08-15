#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include "blant-fundamentals.h"

// We assume LOWER_TRIANGLE is defined in blant-fundamentals.h; otherwise, Int2TinyGraph won't build a correct graph.
// .... except the ORCA numbering assumes we're using UPPER!  (Thus the tables umap3, umap4, umap5). We need lmaps...

//Functions from libwayne and libblant.c
extern "C" {
    struct TINY_GRAPH;
    void mapCanonMap(char* BUF, short int *K, int k);
    void Int2TinyGraph(TINY_GRAPH* G, int Gint);
    TINY_GRAPH *TinyGraphAlloc(unsigned int n, bool selfLoops, bool directed);
    int TinyGraphBFS(TINY_GRAPH *G, int seed, int distance, int *nodeArray, int *distArray);
    typedef unsigned char Boolean;
    Boolean TinyGraphDFSConnected(TINY_GRAPH *G, int seed);
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

//Column contents information. Order can be changed easily. 
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

static int _numCanon, _canonList[MAX_CANONICALS]; //Lower Canon List
static int _numCanonU, _canonListU[MAX_CANONICALS]; //Upper Canon List
static short int _K[maxBk] __attribute__ ((aligned (8192))); //Holds canon_map.bin for lookup
static TINY_GRAPH* G;
const char* DIR = "orca_jesse_blant_table/";

#if LOWER_TRIANGLE
#if 0 // Nahian: here's where to add stuff. All (or most) the numbers in lmap[3-5] need to change!
//Manually generated mapping from LOWER binary representation of connected canonical graphlets to GRAAL/Przulj numbering
const auto lmap3 = unordered_map<uint64_t,uint64_t>{{3, 1}, {7, 2}};
const auto lmap4 = unordered_map<uint64_t, uint64_t>{{11, 4}, {13, 3}, {15, 6}, {30, 5}, {31, 7}, {63, 8}};
const auto lmap5 = unordered_map<uint64_t, uint64_t>{{75, 11}, {77, 10}, {79, 14}, {86, 9}, {87, 12},
 {94, 16}, {95, 17}, {117, 13}, {119, 19}, {127, 23}, {222, 20}, {223, 22}, {235, 18}, {236, 15}, {237, 21},
  {239, 24}, {254, 25}, {255, 26}, {507, 27}, {511, 28}, {1023, 29}};
#define MAP3 lmap3
#define MAP4 lmap4
#define MAP5 lmap5
#warning "sorry, we don't have lmaps yet!"
#else
#endif 
//Manually generated mapping from upper(FAYE) binary representation of connected canonical graphlets to GRAAL/Przulj numbering
const auto umap3 = unordered_map<uint64_t,uint64_t>{{3, 1}, {7, 2}};
const auto umap4 = unordered_map<uint64_t, uint64_t>{{11, 4}, {13, 3}, {15, 6}, {30, 5}, {31, 7}, {63, 8}};
const auto umap5 = unordered_map<uint64_t, uint64_t>{{75, 11}, {77, 10}, {79, 14}, {86, 9}, {87, 12},
 {94, 16}, {95, 17}, {117, 13}, {119, 19}, {127, 23}, {222, 20}, {223, 22}, {235, 18}, {236, 15}, {237, 21},
  {239, 24}, {254, 25}, {255, 26}, {507, 27}, {511, 28}, {1023, 29}};
#define MAP3 umap3
#define MAP4 umap4
#define MAP5 umap5
#endif

int canonListPopulate(char *BUF, int *canon_list, int k, char c) {
    stringstream ss;
    ss << DIR << "/canon_list" << c << k << ".txt";
    ifstream infile;
    infile.open(ss.str());
    if(!infile) {
        cerr << "Cannot open: " << ss.str() << "\n";
    }
    int numCanon;
    infile >> numCanon;
    for (int i = 0; i < numCanon; i++) {
        infile >> canon_list[i];
#if 0 // leave this out for now since magictable is using old canon_list files without the connected Boolean or numEdges.
	int connected, numEdges;
        infile >> connected;
        infile >> numEdges;
#endif
    }
    infile.close();
    return numCanon;
}

//Converts between upper and lower decimal forms of a graphette
uint64_t Upper2Lower(uint64_t Gint, int k) {
    auto matrix = vector<vector<bool>>(k, vector<bool>(k, false));

    int i, j, bitPos=0, bit;
    uint64_t Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
    for(i=k-2;i>=0;i--) {
        for(j=k-1;j>i;j--) {
            if(!Gint2) break;
            bit = (1 << bitPos);
            if(Gint & bit) {
                matrix[i][j] = true;
                matrix[j][i] = true;
            }
            Gint2 &= ~bit;
            bitPos++;
        }
	    if(!Gint2) break;
    }

    bitPos=0, Gint2 = 0;
    for(i=k-1;i>0;i--)  {
        for(j=i-1;j>=0;j--) {
            if(matrix[i][j]) {
                bit = (1 << bitPos);
                Gint2 |= bit;
	        }
            bitPos++;
        }
    }
    return Gint2;
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
    int max_k;
    //Connected count starts at 1 because k=2 has a connected node
    if (argc == 2)
        max_k = atoi(argv[1]);
    else {
        cerr << "Please enter one argument: maxK. minK is 3. argc: " << argc << std::endl;
        exit(EXIT_FAILURE);
    }
    int connectedCount = 1;
    for (int k = MIN_K; k <= max_k; k++) {
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
        G = TinyGraphAlloc(k,false,false);

        //Load canon_list and canon_map
        char BUF[BUFSIZ];
        _numCanon = canonListPopulate(BUF, _canonList, k, 'l');
        _numCanonU = canonListPopulate(BUF, _canonListU, k, 'u');
        mapCanonMap(BUF, _K, k);
        if (_numCanon != _numCanonU) { //Sanity Assertion
            perror("Num canons not equal\n");
            exit(-1);
        }

        //Create table
        auto table = vector<vector<uint64_t>>(_numCanon, vector<uint64_t>(TABLE_WIDTH, 0));
        int lowerDecimal;

        //Fill table with graphlet information
        for (size_t i = 0; i < table.size(); i++) {
            table[i][UPPER_ORDINAL] = i;
            table[i][UPPER_DECIMAL] = _canonListU[i];
            lowerDecimal = Upper2Lower(table[i][UPPER_DECIMAL], k);
            table[i][LOWER_ORDINAL] = _K[lowerDecimal];
            table[i][LOWER_DECIMAL] = _canonList[table[i][LOWER_ORDINAL]];

            Int2TinyGraph(G, table[i][LOWER_DECIMAL]);
            table[i][CONNECTED] = TinyGraphDFSConnected(G,0 > 0);
        }

        //Load num nodes first orbit information and put in table
        ss.str("");
        ss << DIR << "num_nodes_first_orbit" << k << ".txt";
        orbitInfile.open(ss.str());
        if (orbitInfile) {
            for (size_t i = 0; i < table.size(); i++) {
                orbitInfile >> table[i][NUM_NODES_FIRST_ORBIT];
            }
            orbitInfile.close();
        } else if (k < 8) {
            cerr << "WARNING: Expected to find num_nodes_first_orbit" << k << ".txt\n";
        }

        //Load upper orbit information
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
        orbitTableUpper.pop_back(); //EOF after space results in row of zeroes at end. Remove

        //Sanity assertion
        if (orbitTableUpper.size() != table.size()) {
            cerr << "k: " << k << " Orbit table upper size: " << orbitTableUpper.size() << " Table: " << table.size() << '\n';
            exit(-1);
        }

        //Calculate orbit information from upper and fill out table
        int numConnectedOrbits = 0;
        int numTotalOrbits = 0;
        for (size_t i = 0; i < orbitTableUpper.size(); i++) {
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


        //Begin processing that requires table to be sorted in lower order
        sort(table.begin(), table.end(), sortLower);

        //Load lower orbit information
        i = 0;
        ss.str("");
        ss << DIR << "orbit_map" << k << ".txt";
        orbitInfile.open(ss.str());
        if (!orbitInfile) {
            cerr << "WARNING: Failed to open: " << ss.str() << "\n";
        } else {
            orbitInfile >> num;
            while (orbitInfile) {
                orbitTable.push_back(vector<int>(k, 0));
                for (int j = 0; j < k; j++) {
                    orbitInfile >> orbitTable[i][j];
                }
                i++;
            }
            orbitTable.pop_back(); //Pop off empty vector
            orbitInfile.close();
        }

        if (orbitTable.size() != table.size()) { //Sanity Assertion
            cerr << "Orbit table size: " << orbitTable.size() << " Table: " << table.size() << '\n';
            exit(EXIT_FAILURE);
        }

        //Calculate lower orbit information and fill out table
        numConnectedOrbits = 0;
        numTotalOrbits = 0;
        for (size_t i = 0; i < orbitTable.size(); i++) {
            orbits.clear();
            for (int j = 0; j < k; j++) {
                orbitInfile >> orbitTable[i][j];
                orbits.insert(orbitTable[i][j]);
            }
            if (orbits.size() != table[i][NUM_ORBITS]) {
                static int prev_k;
                if(prev_k && k != prev_k){
                    cerr << "WARNING: (orbits.size() != table[i][NUM_ORBITS]) for k="<<k<<'\n';
                    prev_k = k;
                }
            }

            if (table[i][0]) {
                table[i][FIRST_ODV_ORBIT_CON] = numConnectedOrbits;
                numConnectedOrbits += table[i][NUM_ORBITS];
            }
            table[i][FIRST_ODV_ORBIT_ALL] = numTotalOrbits;
            numTotalOrbits += table[i][NUM_ORBITS];

        //Jesse/ORCA use Przulj numbering for lower k<=5. AFter they use lower ordering
            if (table[i][CONNECTED]) {
                if (k <= 5) {
                    if (k == 3) {
                        if (umap3.find(table[i][UPPER_DECIMAL]) == umap3.end()) //Map sanity check
                            cerr << "Cannot find: " << table[i][UPPER_DECIMAL] << " in umap" << k << "\n";
                        else
                            table[i][ORCA] = umap3.at(table[i][UPPER_DECIMAL]);
                    } else if (k == 4) {
                        if (umap4.find(table[i][UPPER_DECIMAL]) == umap4.end()) //Map sanity check
                            cerr << "Cannot find: " << table[i][UPPER_DECIMAL] << " in umap" << k << "\n";
                        else
                        table[i][ORCA] = umap4.at(table[i][UPPER_DECIMAL]);
                    } else {
                        if (umap5.find(table[i][UPPER_DECIMAL]) == umap5.end()) //Map sanity check
                            cerr << "Cannot find: " << table[i][UPPER_DECIMAL] << " in umap" << k << "\n";
                        else
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

        //output table sorted in upper order
        sort(table.begin(), table.end(), sortUpper);
        outfile << table;
        outfile.close();
    }
    return 0;
}
