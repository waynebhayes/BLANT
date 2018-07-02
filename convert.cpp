#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <string.h>
using namespace std;

extern "C"
{
    int BlantAddEdge(int a, int b);
    char** convertToEL(char* file);
}
int numofEdges = 0;
int numofNodes = 0;
vector <string> stringEdge;
vector <int> intEdge;
map<string,int> nodeIndex;
vector <string> stringNode;


int find_type(string input)
{
	regex graphml(".*(\\.xml)");
	regex gml(".*(\\.gml)");
	regex leda(".*(\\.leda)");
	regex csv(".*(\\.csv)");
	regex lgf(".*(\\.lgf)");
	regex edgelist(".*(\\.el)");

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
	if(regex_match(input,lgf))
		return 6;
    return 0;
}



void convert_leda(string filename){
	ifstream ifs;
	ifs.open(filename,ifstream::in);
	string line;

	regex pattern("(.*)(\\s)(.*)(\\s)(.*)(\\s)\\|\\{(.*)\\}\\|(.*)");
	while(!ifs.eof()){
		getline(ifs,line);
		if( regex_match(line,pattern))
		{
			int check = 0;
			int target_start =0;
			string start,target;
			for(unsigned int i = 0; i < line.size(); i++)
				{
				if(isspace(line[i]) && check == 0){
					check = 1;
					start = line.substr(0,i);
				}
				if(!isspace(line[i]) && check == 1){
					check = 2;
					target_start = i;
				}
				if(isspace(line[i]) && check == 2){
					target = line.substr(target_start,i-target_start);
		stringEdge.push_back(start);
		stringEdge.push_back(target);
		numofEdges ++;

					break;
				}

			}
		}
	}
	ifs.close();	
}


void convert_lgf( string filename){
	ifstream ifs;
	ifs.open(filename,ifstream::in);
	string line;
	regex pattern("^(.+)(\\s+)(.+)(\\s+)(\\d+)(\\s+)(\\d+)(\\s*)$");

		
	while(!ifs.eof()){
		getline(ifs,line);
		if(regex_match(line,pattern)){
			int check = 0;
			int target_start =0;
			unsigned int i = 0;
			string start,target;
			while(i < line.length())
			{
				if(isspace(line[i]) && check == 0){
					check = 1;
					start =  line.substr(0,i);
				}
				if(!isspace(line[i]) && check == 1){
					check = 2;
					target_start = i;
				}
				if(isspace(line[i]) && check == 2){
					check = 0;
					target = line.substr(target_start,i-target_start);
		stringEdge.push_back(start);
		stringEdge.push_back(target);
		numofEdges ++;

						break;
				}
				i++;

			}
			getline(ifs,line);
		}
	}
	ifs.close();
}



void convert_grapgml(string inp)
{

    regex source("(.*)(source)(.*)");
    regex target("(.*)(target)(.*)");
    ifstream inFile;
    inFile.open(inp);
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
		stringEdge.push_back(final_source);
		stringEdge.push_back(final_target);
		numofEdges ++;

    }
}
    inFile.close();
}
void convert_gml(string inp)
{    
    regex source("(.*)(source)(.*)");
//  regex target("(.*)(target)(.*)");
    ifstream inFile;
    inFile.open(inp);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    string line,s,t;
    while (inFile.good()) {
	getline (inFile,line);
	if(regex_match(line,source))
		{
//		cout << line << " ";
		s = line.substr(9);
		getline (inFile,line);
//		cout << line << endl;
		t = line.substr(9);
		stringEdge.push_back(s);
		stringEdge.push_back(t);
		numofEdges ++;
		}
   }
    inFile.close();
}

void convert_el(string filename){
	regex pattern("(.)+(\\s)(.)+");
	string line;
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	while(!ifs.eof()){
	getline(ifs,line);
	
	if(regex_match(line,pattern)){
		size_t pos = line.find("	");
		if(pos !=  string::npos)
		{
			 string firstNode = line.substr(0,pos);
			 string secondNode = line.substr(pos+1);
			 stringEdge.push_back(firstNode);
			 stringEdge.push_back(secondNode);
			 numofEdges++;

		}
		pos = line.find(" ");
		if(pos !=  string::npos)
		{
			 string firstNode = line.substr(0,pos);
			 string secondNode = line.substr(pos+1);
			 stringEdge.push_back(firstNode);
			 stringEdge.push_back(secondNode);
			 numofEdges++;

		}
	}
}

}


void convert_csv( string filename){
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	string line;
	regex pattern("(.*);(.*)");
	string source,target;
	while(!ifs.eof()){
		getline(ifs,line);
		string start,target;
		if(regex_match(line,pattern)){
			size_t semi = line.find(";");
			target = line.substr(semi);
			start = line.substr(0,line.length()-target.length());
			target = target.substr(1,target.length()-1);
			stringEdge.push_back(start);
			stringEdge.push_back(target);
			numofEdges++;
		}
	}
	ifs.close();
}

void indexStringNodes(){

 	 for (std::vector<string>::const_iterator i = stringEdge.begin(); i != stringEdge.end(); ++i)
	{
		string s = *i;
		if(nodeIndex.find(s) == nodeIndex.end())
		{
			nodeIndex[s] = numofNodes;	
			numofNodes++;
			stringNode.push_back(s);
		}		
	}

}

void hashEdges()
{
	 for (std::vector<string>::const_iterator i = stringEdge.begin(); i != stringEdge.end(); ++i)
	{
		intEdge.push_back(nodeIndex[*i]);
			}

}

void addAllEdges()
{
	int s,t;
	int c = 0;
	 for (std::vector<int>::const_iterator i = intEdge.begin(); i != intEdge.end(); ++i)
	{
		c++;
		if(c==1)
			s = *i;
		if(c==2)
		{
			t = *i;
			BlantAddEdge(s,t);
			c = 0;
		}
	}
}

void convert(string filename)
 {
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
		convert_el(filename);
		break;
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
}


char** getStringNodes()
{

char ** result = new char * [numofNodes];
int index = 0;
 for (std::vector<string>::const_iterator i = stringNode.begin(); i != stringNode.end(); ++i)
{
result[index] = new char[(*i).size()+1];
strcpy(result[index],(*i).c_str());
index++;
}
return result;
}



char** convertToEL(char* file)
{
string filename = file;
convert(filename);
indexStringNodes();
hashEdges();
char ** temp =  getStringNodes();
addAllEdges();
return temp;
}
