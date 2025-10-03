// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <string>
#include <stdint.h>

std::string appendLeadingZeros(const std::string& inputBitstring, int k);

//Converts 64 bit decimal input into a bit string
std::string toBitString(uint64_t inputDecimalNum, int k);

uint64_t toDecimal(const std::string& inputBitstring, int k);
