// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include "graphette2dotutils.h"
#include <sstream>
#include <iostream>
#include <algorithm>

using std::stringstream;
using std::string;

string appendLeadingZeros(const string& inputBitstring, int k) {
	size_t bitStringFinalSize = static_cast<size_t>(k * (k - 1) / 2);
	size_t inputSize = inputBitstring.size();
	if (bitStringFinalSize > inputSize) {
		size_t leadingZeros = bitStringFinalSize - inputSize;
		stringstream ss;
		for (size_t i = 0; i < leadingZeros; i++) {
			ss << '0';
		}
		ss << inputBitstring;
		return ss.str();
	}
	return inputBitstring;
}

//Converts 64 bit decimal input into a bit string
string toBitString(uint64_t inputDecimalNum, int k) {
    stringstream ss;

	//Convert input decimal to reversed bitstring
	while (inputDecimalNum) {
		ss << (inputDecimalNum & 1);
		inputDecimalNum /= 2;
	}

	//Reverse string to correct reverse order of bits
	string result = ss.str();
	std::reverse(result.begin(), result.end());
	return result;
}

uint64_t toDecimal(const std::string& inputBitstring, int k) {
    int _Bk = k * (k - 1) / 2;
    if (_Bk < inputBitstring.size()) {
        std::cerr << "Input Bitstring too long for k=" << k << ". Max Size: " << _Bk << " Found: " << inputBitstring.size() << "\n";
        exit(EXIT_FAILURE);
    }
    uint64_t ret = 0;
    for (auto i = inputBitstring.cbegin(); i != inputBitstring.cend(); i++) {
        if (*i != '0' && *i != '1') {
            std::cerr << "Invalid Bitstring: " << inputBitstring << "\n";
            exit(EXIT_FAILURE);
        }
        ret |= *i - 48;
        ret <<= 1;
    }
    return ret;
}
