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

void getprops(const string &input_file_name, int evalues = -1) {
    string output_file_name = "/tmp/" + fs::path(input_file_name).filename().string() + to_string(getpid());

    // Load graph dynamically
    PUNGraph snap_graph = LoadEdgeList(input_file_name.c_str(), 0, 1);
    if (!snap_graph) {
        cerr << "Failed to load graph from " << input_file_name << endl;
        return;
    }

    int nodes = snap_graph->GetNodes();
    int edges = snap_graph->GetEdges();
    
    vector<set<int>> adj_list_graph(nodes);
    
    ifstream infile(input_file_name);
    int a, b;
    while (infile >> a >> b) {
        if (adj_list_graph[a].find(b) == adj_list_graph[a].end()) {
            edges++;
        }
        adj_list_graph[a].insert(b);
        adj_list_graph[b].insert(a);
    }
    infile.close();

    // Clustering Coefficient & Eccentricity
    double cc_sum = 0;
    int diameter = -1;
    vector<double> ev;
    map<int, int> khop;
    map<int, int> degree_dist;

    for (int nodeid = 0; nodeid < nodes; nodeid++) {
        double cc = GetNodeClustCf(snap_graph, nodeid);
        cc_sum += cc;
        int ecc = GetNodeEcc(snap_graph, nodeid, false);
        diameter = max(diameter, ecc);
        
        TIntPr64V nodevec;
        GetNodesAtHops(snap_graph, nodeid, nodevec, false);
        for (auto item : nodevec) {
            int k = item.GetVal1();
            int nodes_at_k = item.GetVal2();
            khop[k] += nodes_at_k;
        }
    }

    // Degree Distribution
    TIntPr64V degvec;
    TSnap::GetDegCnt(snap_graph, degvec);
    for (auto item : degvec) {
        degree_dist[item.GetVal1()] = item.GetVal2();
    }

    // Eigenvalues
    TFltV peigv;
    int eig_count = (evalues == -1) ? nodes / 2 : min(evalues / 2, nodes / 2);
    if (eig_count > 0) {
        GetEigVals(snap_graph, eig_count, peigv);
        for (int i = 0; i < peigv.Len(); i++) {
            ev.push_back(peigv[i]);
        }
        sort(ev.rbegin(), ev.rend());
    }

    // Betweenness Centrality
    TIntFlt64H nbw;
    TIntPrFlt64H ebw;
    GetBetweennessCentr(snap_graph, nbw, ebw, 1.0, false);

    // Output results
    cout << "\n########################################################### Global" << endl;
    cout << "Eigenvalues " << ev.size() << ": ";
    for (double v : ev) {
        cout << fixed << setprecision(4) << v << " ";
    }
    cout << endl;
    
    cout << "Nodes: " << nodes << " Edges: " << edges << endl;
    cout << "Avg. Clustering Coefficient: " << cc_sum / nodes << endl;
    cout << "Diameter: " << diameter << endl;

    cout << "K-hop distribution: ";
    for (auto &[k, v] : khop) {
        cout << v << " ";
    }
    cout << endl;

    cout << "Degree Distribution: ";
    for (int i = 0; i <= degree_dist.rbegin()->first; i++) {
        cout << degree_dist[i] << " ";
    }
    cout << endl;

    cout << "nodeName clusCoff eccentricity node_betweenness" << endl;
    for (auto &[nodeid, val] : nbw) {
        cout << nodeid << " " << nbw[nodeid] << endl;
    }

    cout << "node1 node2 edge_betweenness" << endl;
    for (auto &[edge, bw] : ebw) {
        cout << edge.GetVal1() << " " << edge.GetVal2() << " " << bw << endl;
    }

    cout << "########################################################### End-of-output\n";
}

int main(int argc, char* argv[]) {
    if (argc == 2) {
        getprops(argv[1]);
    } else if (argc == 4 && string(argv[1]) == "-e") {
        int evalues = stoi(argv[2]);
        getprops(argv[3], evalues);
    } else {
        cerr << "Input error! Usage: ./main [-e eigValsToCompute] inputFile" << endl;
        return 1;
    }

    dlclose(handle); // Close the library when done
    return 0;
}
