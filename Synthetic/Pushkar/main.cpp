#include "../../libwayne/C++/SanaGraphBasis/Graph.hpp"
#include <filesystem>
#include <dlfcn.h>
#include <vector>
#include <queue>
#include <unordered_map>

constexpr std::filesystem::path script_path = std::filesystem::current_path();
constexpr std::filesystem::path library_path = script_path/"Snap-6.0/";

std::pair<int, int> dobfs(
	const std::unordered_map<int, std::vector<int>>& graph,
    	int x,
    	std::unordered_map<int, int>& components,
    	int cnum){
    int nodes = 0, edges = 0;
    std::queue<int> queue; //bfs queue
    queue.push(x);
    components[x] = cnum;
    while(!queue.empty()){
    	int src = queue.front();
    	queue.pop();
    	++nodes;
    	edges += graph.at(src).size();
    	for(int n : graph.at(src)){
    	    if(components.find(n) == components.end()){
    	    	queue.push(n);
    		components[n] = cnum;
    	    }
    	}
    }	
    return {nodes, edges/2};
}

int getcomponents(
	const std::unordered_map<int, std::vector<int>>& graph,
	const std::vector<std::pair<int, int> csize){
    int n = graph.size();
    std::vector<int> components(n, -1);
    int cnum = 0;
    for(int i = 0; i < n; ++i){
	if(components[i] == -1){
	    std::pair<int, int> numnodenumedge = dobfs(graph, i, components, cnum);
	    csize[cnum] = numnodenumedge;
	    ++cnum;
	}
    }
    return cnum;
}

int main(){
    void* handle = dlopen(library_path, RTLD_LAZY);
    if(!handle){	
	std::cerr << "Failed to load library: " << dlerror() << std::endl;
    	dlclose(handle);
    	return 1;
    }
    dlerror();
    //insert_code_here	
        
    dlclose(handle);	
    return 0;
}
