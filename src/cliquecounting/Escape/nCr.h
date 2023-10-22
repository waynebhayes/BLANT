#ifndef ESCAPE_NCR_H_
#define ESCAPE_NCR_H_

#include <fstream>
#include <sstream>

using namespace Escape;

double nCr[1001][101];

void populate_nCr()
{
    std::ifstream file;
    std::string ifname = "../Escape/nCr.txt";
    std::cout << ifname << std::endl;
    file.open(ifname);
    if (!file.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        return;
    }
    std::cout << "In populate_nCr" << std::endl;


    for(int row = 0; row < 1001; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < 101; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;
            // std::cout << "val = " << val << std::endl;

            std::stringstream convertor(val);
            convertor >> nCr[row][col];
            // std::cout << nCr[row][col] << " ";
        }
        //std::cout << std::endl;
    }
}


#endif
