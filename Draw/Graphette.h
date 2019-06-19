#include <string>
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
