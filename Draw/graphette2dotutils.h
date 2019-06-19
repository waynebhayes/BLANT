#include <string>

std::string appendLeadingZeros(const std::string& inputBitstring, int k);

//Converts 64 bit decimal input into a bit string
std::string toBitString(uint64_t inputDecimalNum, int k);

uint64_t toDecimal(const std::string& inputBitstring, int k);
