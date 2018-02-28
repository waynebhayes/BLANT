#include <iostream>
#include <fstream>
#include <string.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;
typedef long long unsigned int ullint;
vector< vector<ullint> > orbitId_;

ullint strToUllint(string word1,int size){
    ullint x=0,m=1;
    for(int i=size-1;i>=0;i--)
    {
       x+=(word1[i]-'0')*m;
       m*=10;
    }
    return x;
}

void orbit(ullint l, string Permutation, int k, vector<ullint>& o){
     o.resize(k);
    for(ullint i = 0; i < k; i++){
       o[i]=orbitId_[l][Permutation[i]-'0'];
    }
}

void Fatal(string s)
{
    cerr << "ERROR: " << s << endl;
    cerr<<"USAGE: k old_canon_map orbit_map canon_list new_canon_map\n";
    exit(1);
}

int main(int argc,char* argv[]){
    if(argc!=6) Fatal("not enough arguments");
    bool check=1, same=1;
    int k=argv[1][0]-'0';
    ullint m, decimal,numOrbitId_, new_int, old_int;
    ullint e = k*(k-1)/2;
    ullint p = pow(2,e);
    ullint t = 0;
    ifstream fcanon_map(argv[2]), forbit_map(argv[3]);
    ifstream fcanon_list(argv[4]);
    if(fcanon_map.fail())Fatal("old canon_map file not found");
    if(forbit_map.fail())Fatal("orbit_map file not found");
    if(fcanon_list.fail())Fatal("canon_list file not found");
    string canonicalPermutation_new;
    string canonicalPermutation;

    vector<ullint> canonicalGraphette;
    fcanon_list >> m; //reading the number of canonical graphettes
    forbit_map >> numOrbitId_; //reading the number of orbit ids
    for(ullint i = 0; i < m; i++){
    fcanon_list >> decimal;
    canonicalGraphette.push_back(decimal);
    //reading orbit ids for the canonical graphette decimal
    vector<ullint> ids(k);
     for(short j = 0; j < k; j++){
    forbit_map >> ids[j];
     }
     orbitId_.push_back(ids);
     }
    fcanon_list.close();
    forbit_map.close();
    ifstream fcanon_map_new(argv[5]);
    if(fcanon_map_new.fail())Fatal("new canon_map file not found");
    string line, word1;
    ullint t1=0;
    vector<ullint> o1(k);
    vector<ullint> o2(k);
    int n=0,m12=0;
    bool b;
    while(t1<p){
    	fcanon_map >> old_int;
    	//fcanon_map_new >> new_int;
    	getline(fcanon_map_new,line);
        m12=line.find_first_of(" \t,:-;");
    	word1 = line.substr(0,m12);
        line = line.substr(m12+1);
   	canonicalPermutation_new = line.substr(0,line.find_first_of(" \t,:-"));
    	new_int=strToUllint(word1,word1.size());
    	check &= (new_int == old_int);
        fcanon_map >> canonicalPermutation;
    	ullint l = lower_bound(canonicalGraphette.begin(), canonicalGraphette.end(), new_int) - canonicalGraphette.begin();
    	orbit(l,canonicalPermutation_new,k,o1);
    	orbit(l,canonicalPermutation,k,o2);
    	for(int r=0; r<k; r++)
    		check &= (o1[r]==o2[r]);
    	if(check==0){cout << "line no: " << t1+1 << " old file: " << old_int << " " << canonicalPermutation << " new file: " << new_int << " " << canonicalPermutation_new << endl; check = 1;same = 0;}
    	t1++;}
    fcanon_map.close();
    fcanon_map_new.close();
    if(same==1) return 0;
    else return 1;
}
