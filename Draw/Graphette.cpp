// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
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
