// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <iosfwd>

template <typename charT, typename traits>
int fileno(const std::basic_ios<charT, traits>& stream);
