// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <string>
#include <stdint.h>
#include "graphette2dotutils.h"
#define DYNAMIC_CANON_MAP 1

enum class TriangularRepresentation { none, lower, upper };

class Graphette {
	public:
	#if DYNAMIC_CANON_MAP
	Graphette(const std::string& bitstring, TriangularRepresentation triangularRepresentation, int k, bool directed);
	#else
	Graphette(const std::string& bitstring, TriangularRepresentation triangularRepresentation, int k, short int* _K, bool directed);
	Graphette(unsigned __int128 lowerOrdinal, TriangularRepresentation triangularRepresentation, int k, int* _canonList, bool directed);
	unsigned __int128 lowerOrdinal;
	#endif
	std::string bitstring;
	unsigned __int128 decimal;
	TriangularRepresentation triangularRepresentation;
	int k;
};
