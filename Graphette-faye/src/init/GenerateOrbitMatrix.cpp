// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include "GenerateOrbitMatrix.hpp"

using namespace std;

void generateOrbitMatrix(ullint n){
    string inname = "data/canon_list"+to_string(n)+".txt";
    string outname = "data/orbit_map"+to_string(n)+".txt";
    ifstream fin(inname);
    ofstream fout(outname);
    ullint num, cgraph, orbitId = 0;
    vector<ullint> canonical;

    fin >> num;
    while(fin >> cgraph){
        Graphette g = Graphette(n, cgraph);
        vector<ullint> idList(n);        
        vector<vector<ullint>> orbits = g.orbits();
        for(auto orbit : orbits){
            for(auto node : orbit){
                idList[node] = orbitId;
            }
            orbitId++;
        }
        for(ullint i=0; i < idList.size(); i++){
                fout << idList[i] << " ";
            }
        fout << endl;
    }
    fin.close();
    fout.close();
}
