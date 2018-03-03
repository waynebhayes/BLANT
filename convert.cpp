#include<cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include<cstdio>
#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;


extern "C"
{
    int BlantAddEdge(int a, int b);
    char** convertToEL(char* file);
}

int numofEdges = 0;
int find_type(string input)
{
	regex graphml(".*(\\.xml)");
	regex gml(".*(\\.gml)");
	regex leda(".*(\\.leda)");
	regex csv(".*(\\.csv)");
	regex lgf(".*(\\.lgf)");
	regex edgelist(".*(\\.el)");
	if(regex_match(input,lgf))
	return 6;
	if(regex_match(input,graphml))
	return 1;
	if(regex_match(input,gml))
	return 2;
	if(regex_match(input,edgelist))
	return 3;
	if(regex_match(input,csv))
	return 4;
	if(regex_match(input,leda))
	return 5;

    return 0;
}

void convert_leda( string filename){
	ifstream ifs;
	ifs.open(filename,ifstream::in);
	ofstream ofs;
	ofs.open("/tmp/example.txt");
	
	
	 string line;
	 regex pattern("([a-zA-Z\\d]*)(\\s)([a-zA-Z\\d]*)(\\s)([a-zA-Z\\d]*)(\\s)\\|\\{(.*)\\}\\|(.*)");
	while(!ifs.eof()){
		getline(ifs,line);
		if( regex_match(line,pattern))
		{
			int check = 0;
			int target_start =0;
			 string start,target;
			for(unsigned int i = 0; i < line.size(); i++){
				if(isspace(line[i]) && check == 0){
					check = 1;
					start = line.substr(0,i);
					ofs<<start<<" ";
				}
				if(!isspace(line[i]) && check == 1){
					check = 2;
					target_start = i;
				}
				if(isspace(line[i]) && check == 2){
					target = line.substr(target_start,i-target_start);
					ofs<<target<< endl;
					break;
				}
			}
		}
	}
	ifs.close();
	ofs.close();
	
}





void convert_lgf( string filename){
	ifstream ifs;
	ifs.open(filename,ifstream::in);
	ofstream ofs;
	ofs.open("/tmp/example.txt");
	
	string line;
	regex pattern("^([\\da-zA-Z]+)(\\s+)([\\da-zA-Z]+)(\\s+)(\\d+)(\\s+)(\\d+)(\\s*)$");

		
	while(!ifs.eof()){
		getline(ifs,line);
		if(regex_match(line,pattern)){
			int check = 0;
			int target_start =0;
			unsigned int i = 0;
			while(i < line.length())
			{
				if(isspace(line[i]) && check == 0){
					check = 1;
					string start =  line.substr(0,i);
					ofs<<start<<" ";
				}
				if(!isspace(line[i]) && check == 1){
					check = 2;
					target_start = i;
				}
				if(isspace(line[i]) && check == 2){
					check = 0;
					string target = line.substr(target_start,i-target_start);
					ofs<<target<< endl;
						break;
				}
				i++;
			}
			getline(ifs,line);
		}
	}
	
	ifs.close();
	ofs.close();
}



void convert_grapgml(string inp)
{

    regex source("(.*)(source)(.*)");
    regex target("(.*)(target)(.*)");
    ifstream inFile;
    ofstream myfile;
    inFile.open(inp);
    myfile.open ("/tmp/example.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    string line;
    string sub;
    while (inFile.good()) {

	getline (inFile,line);
	if(regex_match(line,source))
{
		size_t  s_place = line.find("source=");
		string re = line.substr(s_place);
		int n = 7;
		if(re[7] == '"')
			n++;
		string sub = re.substr(n);
		char f = ' ';
                if(n == 8)
			f = '"';
		string final_source;
		for(string::iterator it=sub.begin(); it!=sub.end(); ++it)
		{
			if(*it == f)
				break;
			final_source.push_back(*it);
			
		}
		size_t  t_place = re.find("target=");
		string t_re = re.substr(t_place);
		int x = 7;
		if(t_re[7] == '"')
			x++;
		string t_sub = t_re.substr(n);
		string final_target;
		for(string::iterator tit=t_sub.begin(); tit!=t_sub.end(); ++tit)
		{
			if(*tit == f || *tit == '>'||*tit == '/' )
				break;
			final_target.push_back(*tit);
			
		}
		myfile << final_source << "  " << final_target << "\n";
    }
}
    myfile.close();
    inFile.close();

}
void convert_gml(string inp)
{    
    regex source("(.*)(source)(.*)");
    regex target("(.*)(target)(.*)");
    ifstream inFile;
    ofstream myfile;
    inFile.open(inp);
    myfile.open ("/tmp/example.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    string line;
    string sub;
    while (inFile.good()) {

	getline (inFile,line);
	if(regex_match(line,source))
		
		myfile << line.substr(9) << ' ';
	if(regex_match(line,target))
		myfile << line.substr(9) << "\n";
continue;
	


    }
    myfile.close();
    inFile.close();

}




void convert_csv( string filename){
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	ofstream ofs;
	ofs.open("/tmp/example.txt");
	
	string line;
	regex pattern("(.*);(.*)");
	
	while(!ifs.eof()){
		getline(ifs,line);
		string start,target;
		if(regex_match(line,pattern)){
			size_t semi = line.find(";");
			target = line.substr(semi);
			start = line.substr(0,line.length()-target.length());
			target = target.substr(1,target.length()-1);
			ofs<<start<<" ";
			ofs<<target<< endl;
		}
	}
	ifs.close();
	ofs.close();
}

map<int,string> toString(string filename){

	ifstream ifs;
	ifs.open(filename,ifstream::in);
	ofstream ofs;
	ofs.open("/tmp/hashed.txt");
	
	vector<int> Nodes;
	vector<int>::iterator it;
	map<string,int> ref;
	map<int,string> result;
	int start = 0;
	
	string line;
	regex pattern(".*");
	
	while(!ifs.eof()){
		getline(ifs,line);
		if(regex_match(line,pattern)){
			size_t pos = line.find(" ");

			if(pos != string::npos)
			{
				string firstNode = line.substr(0,pos);
				string secondNode = line.substr(pos+1);
				
				if(ref.find(firstNode) == ref.end()){
					Nodes.push_back(start);
					ref[firstNode] = start;
					result[start] = firstNode;
					start++;
					numofEdges ++;
				}
				if(ref.find(secondNode) == ref.end()){
					Nodes.push_back(start);
					ref[secondNode] = start;
					result[start] = secondNode;
					start++;
					numofEdges ++;
				}
				ofs<<ref[firstNode]<<" "<<ref[secondNode]<<endl;
				
			}
		}
		
	}
	ifs.close();
	ofs.close();
	return result;
}



string  convert(string filename)
 {

    string result = "example.txt";
    switch(find_type(filename))
{
	
	case 0: 
	{cout << "Wrong" << endl;
		exit(1);}
	case 1: 
	{ 	convert_grapgml(filename);
		break;}
	case 2: 
	{	convert_gml(filename);
		break;}
	case 3:
	{
		return filename;
	}
	case 4: 
	{	convert_csv(filename);
		break;}
	case 5: 
	{	convert_leda(filename);
		break;}
	case 6: 
	{	convert_lgf(filename);
		break;}

}

    return result;
}

void addAllEdges(string filename)
{
	
	regex pattern("(.)+(\\s)(.)+");
	string line;
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	while(!ifs.eof()){
	getline(ifs,line);
	if(regex_match(line,pattern)){
		size_t pos = line.find(" ");

		if(pos !=  string::npos)
		{
			 string firstNode = line.substr(0,pos);
			 string secondNode = line.substr(pos+1);
		//	 cout << stoi(firstNode) << "   " << stoi(secondNode) << endl;
			 BlantAddEdge(stoi(firstNode),stoi(secondNode));

		}
	}
}

}


char** getStringNodes(string filename)
{
map<int,string> strnodes = toString(filename);
char ** result = new char * [numofEdges];
for(int i = 0;i< numofEdges;i++)
{
result[i] = new char[strnodes[i].size()+1];
strcpy(result[i],strnodes[i].c_str());
}
return result;

}


void clean(char** result)
{
remove("/tmp/example.txt");
remove("/tmp/hashed.txt");

}

char** convertToEL(char* file)
{
string filename = file;
string converted  = convert(filename);
char ** temp = getStringNodes(converted);
/*
for(int i = 0;i< numofEdges;i++)
{
cout << temp[i]<< endl;
}*/
addAllEdges("/tmp/hashed.txt");
return temp;
}



#if 0
int main()
{
string converted  = convert("test.el");
char ** temp = getStringNodes(converted);

for(int i = 0;i< numofEdges;i++)
{
cout << temp[i]<< endl;
}
addAllEdges("hashed.txt");
RunBlantEdgesFinished(temp);
clean(temp);
return 0;

}

#endif
