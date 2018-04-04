#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

uint64_t Upper2Lower(uint64_t Gint, int k)
{
    auto matrix = vector<vector<uint64_t>>(k, vector<uint64_t>(k, 0));
    unsigned j = 1, i = 0;
    uint64_t mask = 1, result;

    for (int bit = 0; bit < k*(k-1)/2; bit++) {
        if (Gint & mask) {
            matrix[j][i] = 1;
        }
        mask <<= 1;
        j++;
		if (j == k) {
			i++;
			j = i + 1;
		}
    }
    mask = 1;
    j = 0;
    i = 1;
    for (int bit = 0; bit < k*(k-1)/2; bit++) {
        if (matrix[j][i]) {
            result |= mask;
        }

        mask <<= 1;
        j++;
        if (j == i) {
            i++;
            j = 0;
        }

    }
}

//first, check 4th column matches blant canonical
//2nd, fix 5th column blant ordinale
int main(int argc, char* argv[]) {
    for (int k = 3; k <= 7; k++) {
        auto table = vector<vector<int>>();
        ofstream outfile;
        outfile.open("transformations/Table" + to_string(k) + ".txt");
        ifstream infile;
        infile.open("transformations/UpperToLower" + to_string(k) + ".txt");
        if (!infile) {
            cerr << "Failed to open file\n";
            exit(EXIT_FAILURE);
        }
        int j = 0;
        while (infile) {
            int temp;
            table.push_back(vector<int>());
            for (int i = 0; i < 7; i++) {
                infile >> temp;
                table[j].push_back(temp);
            }
            j++;
        }

        outfile.close();
        infile.close();


        //Testing
        for (int i = 0; i < table.size(); i++) {
            int lowerCanonical = 0;
            lowerCanonical = Upper2Lower(table[i][2], k);
            assert(table[i][3] == lowerCanonical);
        }
    }
    return 0;
}