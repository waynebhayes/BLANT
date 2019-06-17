#include "Graphette.h"

Graphette::Graphette(const std::string& bitstring, TriangularRepresentation triangularRepresentation, int k, short int* _K) : 
    bitstring(bitstring), 
    triangularRepresentation(triangularRepresentation), 
    k(k) {
        this->decimal = toDecimal(bitstring, k);
        if (triangularRepresentation == TriangularRepresentation::lower)
            this->lowerOrdinal = _K[this->decimal];
        else
            this->lowerOrdinal = -1;
}

Graphette::Graphette(uint64_t lowerOrdinal, TriangularRepresentation triangularRepresentation, int k, int* _canonList) : 
    lowerOrdinal(lowerOrdinal), 
    triangularRepresentation(triangularRepresentation),
    k(k) {
    this->decimal = _canonList[lowerOrdinal];
    this->bitstring = appendLeadingZeros(toBitString(this->decimal, k), k);
}
