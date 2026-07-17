// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
//Requires c++11 for stoull
#include <string>
#include <vector>
#include <cstring>
#include <cstdio>
#include <math.h>
#include <unistd.h>
#include "graphette2dotutils.h"
#include "Graphette.h"


using std::cout;
using std::cerr;
using std::string;
using std::stringstream;
using std::ofstream;
using std::ifstream;
using std::vector;

const int RADIUS_SCALING = 25;
const string USAGE = "USAGE: graphette2dot <-k number_of_nodes> <-b bitstring | -d decimal_representation | -i lower_ordinal> <-o output_filename> [-n input_filename] [-u | -l] [-t title] [-e edge width] [-a enable orbits] [-p disable circular] -h for verbose help\n";
const string GRAPH_ARGS = "";
const string NODE_ARGS = "shape = \"none\", fontsize = 12.0";
const string EDGE_ARGS = "";
const string TITLE_ARGS = "fontsize = 24.0";
const string DECIMAL_INPUT_WARNING = "Warning. Decimal input was used with k > 11. Edge information may have been lost.\n";
const double PI  = 3.141592653589793238463;
const vector<string> COLORS = {"black", "red", "lawngreen", "orange", "blue", "yellow", "indigo"};
const int DEFAULT_EDGE_WIDTH = 1;

enum class OutputMode { none, circular, planar };

typedef struct _Graphette2DotParams {
	int k = 0;
	vector<Graphette> graphettes;
	string outputFile, namesFile, graphTitle;
	int edgewidth = DEFAULT_EDGE_WIDTH;
	TriangularRepresentation triangularRepresentation = TriangularRepresentation::lower;
	OutputMode outputMode = OutputMode::circular;
	bool showOrbits = false;
	bool directed = false;
} Graphette2DotParams;

void parseInput(int argc, char* argv[], Graphette2DotParams& params);
void printUsage();
void printHelp();

void createDotfileFromBit(const Graphette2DotParams& params);
string getPos(int i, int k);
void writeEdges(ofstream& outfile, const vector<Graphette>& graphettes, int edgewidth, bool directed);
void printGraphConversionInstruction(const Graphette2DotParams& params);

//Functions from libwayne and libblant.c
extern "C" {
	#include <sets.h>
	SET *SetAlloc(unsigned int n);
	#if DYNAMIC_CANON_MAP
	#define maxK 16
	#else
    SET *canonListPopulate(char *BUF, int *canon_list, int k, char *canon_num_edges, bool directed);
	#define maxK 8
	#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8
	#define MAX_ORBITS	79264	// This is the number of orbits for k=8
	int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][maxK],  int orbit_canon_mapping[MAX_ORBITS], char orbit_canon_node_mapping[MAX_ORBITS], int numCanon, int k, bool directed);
	void mapCanonMap(char* BUF, short int *K, int k, bool directed);
	#endif
	#define maxBk 1 << (maxK*(maxK-1)/2) // maximum number of entries in the canon_map - need to change for k=6 directed graphs
}
#if !DYNAMIC_CANON_MAP
// Canon Maps Loading
static unsigned _Bk;
static int _numCanon;
static int _canonList[MAX_CANONICALS];
static SET *_connectedCanonicals;
static int _numOrbits;
static int _orbitList[MAX_CANONICALS][maxK];
static char _canonNumEdges[MAX_CANONICALS];
static char _orbitCanonNodeMapping[MAX_ORBITS];
static int _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)
static short int _K[maxBk] __attribute__ ((aligned (8192)));

void loadCanonMaps(int k, bool directed) {
	char BUF[BUFSIZ];
	_Bk = ( directed ? (1 <<(k*(k-1))) : (1 <<(k*(k-1)/2)) );
	_connectedCanonicals = canonListPopulate(BUF, _canonList, k, _canonNumEdges, directed);
	_numCanon = _connectedCanonicals->maxElem;
	_numOrbits = orbitListPopulate(BUF, _orbitList, _orbitCanonMapping, _orbitCanonNodeMapping, _numCanon, k, directed);
	mapCanonMap(BUF, _K, k, directed);
}
#endif
int main(int argc, char* argv[]) {
	//Parse input passing variables by reference.
	Graphette2DotParams params; 
	parseInput(argc, argv, params);
	createDotfileFromBit(params);
	printGraphConversionInstruction(params);
	return EXIT_SUCCESS;
}

unsigned __int128 stou128(const string& s) {
	unsigned __int128 result = 0;
	for (char c : s) {
		if (c < '0' || c > '9') break;
		result = result * 10 + (c - '0');
	}
	return result;
}

/**
 * Parses command line input.
 * Doesn't allow for repeated inputs.
 * Prints usage and exits if invalid input is passed.
 * */
void parseInput(int argc, char* argv[], Graphette2DotParams& params) {
	bool input = false, matrixType = false;
	unsigned __int128 inputDecimalNum = 0;
	int opt;
	vector<unsigned __int128> lowerOrdinals;
	vector<unsigned __int128> inputDecimals;
	vector<string> inputBitStrings;

	while((opt = getopt(argc, argv, "k:b:d:i:o:t:e:Dapnhul")) != -1)
    {
		switch(opt)
		{
		case 'k':
			if (params.k > 0) {
				cerr << "Only one k is allowed\n";
				printUsage();
			}
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			params.k = atoi(optarg);
			break;

		case 'b':
			input = true;
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			inputBitStrings.push_back(optarg);
			break;

		case 'd':
			input = true;
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			inputDecimalNum = stou128(optarg);
			inputDecimals.push_back(inputDecimalNum);
			break;

		case 'i':
			input = true;
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			lowerOrdinals.push_back(std::stoull(optarg));
			break;

		case 'o':
			if (!params.outputFile.empty()) {
				cerr << "Only one output file is allowed.\n";
				printUsage();
			}
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			params.outputFile = optarg;
			break;

		case 't':
			if (!params.graphTitle.empty()) {
				cerr << "Only one title allowed.\n";
				printUsage();
			}
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			params.graphTitle = optarg;
			break;

		case 'e':
			if (params.edgewidth != DEFAULT_EDGE_WIDTH) {
				cerr << "Only one edge width allowed.\n";
				printUsage();
			}
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			params.edgewidth = atoi(optarg);
			break;
		
		case 'D':
			params.directed = true;
			break;
		#if !DYNAMIC_CANON_MAP
		case 'a':
			params.showOrbits = true;
			break;
		#endif
		case 'p':
			params.outputMode = OutputMode::planar;
			break;

		case 'n':
			if (!params.namesFile.empty()) {
				cerr << "Only one names file is allowed.\n";
				printUsage();
			}
			if (!*(optarg)) {
				cerr << "No following argument to " << opt << '\n';
				printUsage();
			}
			params.namesFile = optarg;
			break;

		case 'h':
			printHelp();
			break;

		case 'u':
			if (params.triangularRepresentation != TriangularRepresentation::none) {
				cerr << "Only one matrix type allowed.\n";
				printUsage();
			}
			params.triangularRepresentation = TriangularRepresentation::upper;
			break;

		case 'l':
			if (params.triangularRepresentation != TriangularRepresentation::none) {
				cerr << "Only one matrix type allowed.\n";
				printUsage();
			}
			params.triangularRepresentation = TriangularRepresentation::lower;
			break;

		default:
			cerr << "Unrecognized argument: " << opt << "\n";
			printUsage();
		}
	}

	if (params.triangularRepresentation == TriangularRepresentation::none)
		params.triangularRepresentation = TriangularRepresentation::lower;

	int numGraphettes = inputBitStrings.size() + lowerOrdinals.size() + inputDecimals.size();

	if (numGraphettes < 1 || numGraphettes > COLORS.size()) {
		cerr << "Must specify between 1 and " << COLORS.size() << " graphlets as input.\n";
		printUsage();
	}

	if (params.k < 3 || params.outputFile.size() < 0) {
		cerr << "You must specify a graphette size k > 2\n";
		printUsage();
	}
	#if !DYNAMIC_CANON_MAP
	if (numGraphettes > 1 && params.showOrbits) {
		cerr << "Displaying orbits is currently disabled for multiple inputs.\n";
		printUsage();
	}

	if (params.showOrbits && params.triangularRepresentation == TriangularRepresentation::upper) {
		cerr << "Ordinal input is only allowed for lower triangular inputs\n";
	}

	if (params.k > 11 && inputDecimals.size() > 0)
		cerr << DECIMAL_INPUT_WARNING;
	loadCanonMaps(params.k, params.directed);
	#endif
	for (const auto& bitString : inputBitStrings) {
		#if !DYNAMIC_CANON_MAP
		params.graphettes.emplace_back(bitString, params.triangularRepresentation, params.k, _K, params.directed);
		#else
		params.graphettes.emplace_back(bitString, params.triangularRepresentation, params.k, params.directed);
		#endif
	}
	#if !DYNAMIC_CANON_MAP
	for (const auto& ordinal : lowerOrdinals) {
		if (ordinal < 0 || ordinal >= _numCanon) {
			cerr << "Ordinal for k: " << params.k << " must be between 0 and " << _numCanon << '\n';
			exit(EXIT_FAILURE);
		}
		params.graphettes.emplace_back(ordinal, params.triangularRepresentation, params.k, _canonList, params.directed);
	}

	for (const auto& decimal : inputDecimals) {
		params.graphettes.emplace_back(appendLeadingZeros(toBitString(decimal, params.k), params.k, params.directed), params.triangularRepresentation, params.k, _K, params.directed);
	}

	if (!lowerOrdinals.empty()) {
		if (params.k < 3 || params.k > 8) {
			cerr << "Ordinal input is only allowed for k between 3 and 8 (inclusive)\n";
			exit(EXIT_FAILURE);
		}
	}
	#else
	for (const auto& decimal : inputDecimals) {
		params.graphettes.emplace_back(appendLeadingZeros(toBitString(decimal, params.k), params.k, params.directed), params.triangularRepresentation, params.k, params.directed);
	}
	#endif
}

void printUsage() {
	cerr << USAGE;
	exit(EXIT_SUCCESS);
}

//Contains useful information about the assumptions made in the program.
void printHelp() {
	std::cout << USAGE 
			  << "You must specify the number of nodes with -k\n"
			  << "Currently number of nodes is limited to 16 if decimal input is used\n"
			  << "This is because k=17 requires 136 bits to store the decimal input\n"
			  << "You must specify at least one bitstring, decimal, or ordinal input with -b -d or -i respectively\n"
			  << "-u and -l specify if all input is upper or lower triangular row major. Lower is assumed for ordinal input.\n"
			  << "Lower triangular row major is assumed\n"
			  << "You must specify an output file name with -o\n"
			  << "Please do not include file extension for output. The program with generate a .dot\n"
			  << "If no names file is selected, nodes will be named 0, 1, ...., (k-1)\n"
			  << "Names file parsing assumes one name per line\n"
			  << "You may specify a title with -t\n"
			  << "You may specify an edge width with -e. 1 is default\n"
			  << "If less names than nodes, additional nodes will be labeled by their index #\n"
			  << "If more names than nodes, additional names will be ignored\n";
	exit(EXIT_SUCCESS);
}

/**Creates .dot file from given input.
 * First, node names are listed.
 * If a names file was specified, the nodes are labeled with their name.
 * If the names file has a different number of names than k,
 * extra names become isolated nodes and additional nodes aren't labeled.
 * Then the edges are written to the file.
 * */
void createDotfileFromBit(const Graphette2DotParams& params) {
	int finalBitstringSize = ( params.directed ? (params.k * (params.k - 1)) : (params.k * (params.k - 1)) / 2 );
	int size;

	for (const auto& graphette : params.graphettes) {
		size = graphette.bitstring.size();
		if (finalBitstringSize != size) {
			cerr << "Input size does not match number of nodes.\n"
				<< "Expected Bitstring Size given k = " << params.k << " is:  "
				<< finalBitstringSize << "\nInput Bitstring Size: " << size << "\n";
			exit(EXIT_FAILURE);
		}
	}

	ofstream outfile;
	stringstream ss;
	ss << params.outputFile << ".dot";
	outfile.open(ss.str());
	if (!outfile) {
		cerr << "Unable to create Dot File " << ss.str() << "\n";
		exit(EXIT_FAILURE);
	}
	
	outfile << (params.directed ? "digraph {\n" : "graph {\n") << GRAPH_ARGS;

	int i = 0;
	if (params.namesFile != "") {
		std::ifstream infile;
		infile.open(params.namesFile);
		if (infile) {
			string nodeName;
			while (std::getline(infile, nodeName) && i < params.k) {
				outfile << 'n' << i << " [label=<" << nodeName;
				#if !DYNAMIC_CANON_MAP
				if (params.showOrbits)
					outfile << "<BR />" << "<I>" << _orbitList[params.graphettes[0].lowerOrdinal][i] - _orbitList[params.graphettes[0].lowerOrdinal][0] << "</I>";
				#endif
				outfile << ">, ";
				if (params.outputMode == OutputMode::circular)
					outfile << ", pos=\"" << getPos(i, params.k) << "!\"";
				outfile << NODE_ARGS << ";]\n";
				i++;
			}
			if (i < params.k) {
				cerr << "Warning: Less nodes in names file than -k.\n"
					 << "Number of nodes: " << params.k << " Names file number of nodes: " << i << "\n";
			}

			while (std::getline(infile, nodeName))
				i++;
			if (i > params.k) {
				cerr << "Warning: More nodes in names file than -k.\n"
			         << "Number of nodes: " << params.k << " Names file number of nodes: " << i << "\n";								
			}

			infile.close();
		} else {
			cerr << "Could not open name file\n";
		}
	}
	while (i < params.k) {
		outfile << 'n' << i << "[label=<" << i;
		#if !DYNAMIC_CANON_MAP
		if (params.showOrbits)
			outfile << "<BR />" << "<I>" << _orbitList[params.graphettes[0].lowerOrdinal][i] - _orbitList[params.graphettes[0].lowerOrdinal][0] << "</I>";
		#endif
		outfile << ">, ";
		if (params.outputMode == OutputMode::circular)
			outfile << "pos=\"" << getPos(i, params.k) << "!\"";
		outfile  << NODE_ARGS << ", width=.25" << ";]\n";
		i++;
	}

	writeEdges(outfile, params.graphettes, params.edgewidth, params.directed);

	if (params.graphTitle != "") {
		outfile << "labelloc=\"b\";\n"
				<< "label=\"" << params.graphTitle << "\"\n"
				<< TITLE_ARGS << '\n';
	}
	outfile << "}";
	outfile.close();
}

string getPos(int i, int k) {
	stringstream ss;
	ss << RADIUS_SCALING * k * cos(PI /2 - (2 * PI / k * i)) << ", " << RADIUS_SCALING * k * sin(PI / 2 - (2 * PI / k * i));
	return ss.str();
}

void writeEdges(ofstream& outfile, const vector<Graphette>& graphettes, int edgewidth, bool directed) {
	string penwidth = "";
	if (edgewidth != 1) {
		penwidth = (", penwidth=" + std::to_string(edgewidth));
	}

	unsigned i, j;
	for (size_t c = 0; c < graphettes.size(); c++) {
		switch (graphettes[c].triangularRepresentation) {
			case TriangularRepresentation::lower:
				i = 1;
				j = 0;
				break;
			case TriangularRepresentation::upper:
				i = 0;
				j = 1;
				break;
			case TriangularRepresentation::none: exit(EXIT_FAILURE); break;
		}
		if(directed) {
			i = 0;
			j = 1;
		}
		for (size_t k = 0; k < graphettes[c].bitstring.size(); k++) {
			if (graphettes[c].bitstring[k] == '1')
				outfile << "n" << i << (directed ? " -> " : " -- ") << "n" << j << "[color=" << COLORS[c] << penwidth << ", " << EDGE_ARGS << "]" << "\n";
			else if (graphettes[c].bitstring[k] != '0')
				cerr << "Unknown input: " << graphettes[c].bitstring[k] << " in input bitstring.\n";
			switch (graphettes[c].triangularRepresentation) {
				case TriangularRepresentation::lower:
					j++;
					if(directed&&(j == i)) j++;
					if (j >= (directed ? graphettes[c].k : i)) {
						i++;
						j = 0;
					}
					break;
				case TriangularRepresentation::upper:
					j++;
					if(directed&&(j == i)) j++;
					if (j >= graphettes[c].k) {
						i++;
						j = (directed ? 0 : i + 1);
					}
					break;
				case TriangularRepresentation::none: exit(EXIT_FAILURE); break;
			}
		}
	}
}

void printGraphConversionInstruction(const Graphette2DotParams& params) {
	stringstream ss;
	ss << "neato ";
	if (params.outputMode == OutputMode::circular)
		ss <<"-n ";
	ss << "-Tpdf \"" << params.outputFile << ".dot\" -o \"" << params.outputFile << ".pdf\"";
	std::cout << ss.str() << std::endl;
}
