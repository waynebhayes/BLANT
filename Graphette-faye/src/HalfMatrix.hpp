// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef HALFMATRIX_HPP
#define HALFMATRIX_HPP

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

typedef long long unsigned int ullint;

class HalfMatrix{
public:
    HalfMatrix(){};
    HalfMatrix(ullint n, std::vector<bool>& bitVector);
    HalfMatrix(ullint n, ullint decimalNumber);
    HalfMatrix(ullint n);
    HalfMatrix(const HalfMatrix& m);   // Copy constructor
    ~HalfMatrix();

    ullint length();
    void print();
    bool& operator() (ullint row, ullint col);
    //bool operator() (ullint row, ullint col) const;
    HalfMatrix& operator= (const HalfMatrix& m);
    void clear();
private:
    ullint len_;
       bool* bitArray_ = NULL;
       void encodeBitArray(ullint decimalNumber);
};
#endif