// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <string>
#include <stdint.h>
#include "graphette2dotutils.h"

enum class TriangularRepresentation { none, lower, upper };

class Graphette {
	public:
	Graphette(const std::string& bitstring, TriangularRepresentation triangularRepresentation, int k, short int* _K);

	Graphette(uint64_t lowerOrdinal, TriangularRepresentation triangularRepresentation, int k, int* _canonList);
 
	std::string bitstring;
	uint64_t lowerOrdinal;
	uint64_t decimal;
	TriangularRepresentation triangularRepresentation;
	int k;
};
