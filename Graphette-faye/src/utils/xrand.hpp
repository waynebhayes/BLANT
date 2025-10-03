// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef XRAND_HPP
#define XRAND_HPP
#include <bits/stdc++.h>

typedef unsigned long long  ullint;

ullint xrand(ullint begin, ullint end); //the range is [begin, end)
void xshuffle(std::vector<ullint>& nodes, ullint len);
#endif