#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdlib>    // For atoi
#include <algorithm>
#include <iomanip>
#include <unistd.h>   // For getpid
#include "Snap-6.0/snap-core/Snap.h"

using namespace std;
using namespace TSnap;

// Define a type alias for the graph structure using std::map (instead of unordered_map)
typedef map<int, vector<int> > GraphMap;

// A helper function to extract the filename from a path.
string get_filename(const string &path) {
    size_t pos = path.find_last_of("/\\");
    if (pos == string::npos)
        return path;
    else
        return path.substr(pos + 1);
}

// Modified BFS that works with a GraphMap and a components vector.
// The membership test is replaced by checking if components[n] == -1.
pair<int, int> dobfs(const GraphMap &graph, int x, vector<int> &components, int cnum) {
    int nodes = 0, edges = 0;
    queue<int> q; // BFS queue
    q.push(x);
    components[x] = cnum;
    while (!q.empty()) {
        int src = q.front();
        q.pop();
        ++nodes;
        // Use find() because std::map in C++98 does not have at()
        typename GraphMap::const_iterator it_src = graph.find(src);
        if (it_src != graph.end()) {
            edges += it_src->second.size();
            const vector<int>& neighbors = it_src->second;
            for (size_t i = 0; i < neighbors.size(); i++) {
                int n = neighbors[i];
                if (components[n] == -1) { // If not visited
                    q.push(n);
                    components[n] = cnum;
                }
            }
        }
    }
    return make_pair(nodes, edges / 2);
}

// Compute connected components.
// Here we assume that graph keys range from 0 to n-1.
int getcomponents(const GraphMap &graph, vector< pair<int, int> > &csize) {
    int n = graph.size();
    vector<int> components(n, -1);
    int cnum = 0;
    for (int i = 0; i < n; ++i) {
        if (components[i] == -1) {
            pair<int,int> numnodenumedge = dobfs(graph, i, components, cnum);
            if (csize.size() <= static_cast<unsigned int>(cnum)) {
                csize.push_back(numnodenumedge);
            } else {
                csize[cnum] = numnodenumedge;
            }
            ++cnum;
        }
    }
    return cnum;
}

void getprops(const string &input_file_name, int evalues = -1) {
    // Create an output filename using a helper function and stringstream conversion instead of to_string.
    ostringstream oss;
    oss << getpid();
    string output_file_name = "/tmp/" + get_filename(input_file_name) + oss.str();

    // Load graph dynamically using TSnap.
#if OLD_SNAP_CODE_FOR_INT_ONLY
    PUNGraph snap_graph = LoadEdgeList<PUNGraph>(TStr(input_file_name.c_str()), 0, 1);
    if (snap_graph.Empty()) {
        cerr << "Failed to load graph from " << input_file_name << endl;
        return;
    }
#endif
    PUNGraph snap_graph = TUNGraph::New();
    map<string,int> idmap;
    vector<string> id2name;
    int nextId = 0;

    ifstream gfile(input_file_name.c_str());
    string su, sv;
    while (gfile >> su >> sv) {
        // assign each unique node-string an integer ID
        if (!idmap.count(su)) {
            idmap[su] = nextId;
	    id2name.push_back(su);
            snap_graph->AddNode(nextId++);
        }
        if (!idmap.count(sv)) {
            idmap[sv] = nextId;
	    id2name.push_back(sv);
            snap_graph->AddNode(nextId++);
        }
        // add the undirected edge
        snap_graph->AddEdge(idmap[su], idmap[sv]);
    }
    gfile.close();

    // Now it’s safe to call TSnap routines on snap_graph:
    int nodes = snap_graph->GetNodes();
    int edges = snap_graph->GetEdges();
    
    // Build an adjacency list graph. Note: we use vector<set<int> > for storing unique neighbors.
    vector< set<int> > adj_list_graph(nodes);
    ifstream infile(input_file_name.c_str());
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
        
        TIntPrV nodevec;
        GetNodesAtHops(snap_graph, nodeid, nodevec, false);
        for (int i = 0; i < nodevec.Len(); ++i) {
            int k = nodevec[i].GetVal1();
            int nodes_at_k = nodevec[i].GetVal2();
            khop[k] += nodes_at_k;
        }
    }

    // Degree Distribution
    TIntPrV degvec;
    TSnap::GetDegCnt(snap_graph, degvec);
    for (int i = 0; i < degvec.Len(); ++i) {
        degree_dist[degvec[i].GetVal1()] = degvec[i].GetVal2();
    }

    // Eigenvalues
    TFltV peigv;
    int eig_count = (evalues == -1) ? nodes / 2 : min(evalues / 2, nodes / 2);
    if (eig_count > 0) {
        GetEigVals(snap_graph, eig_count, peigv);
        for (int i = 0; i < peigv.Len(); i++) {
            ev.push_back(peigv[i]);
        }
        // Sort eigenvalues in descending order.
        sort(ev.rbegin(), ev.rend());
    }

    // Betweenness Centrality
    TIntFltH nbw;
    TIntPrFltH ebw;
    GetBetweennessCentr(snap_graph, nbw, ebw, 1.0, false);

    // Output results.
    cout << "\n########################################################### Global" << endl;
    cout << "Eigenvalues " << ev.size() << ": ";
    for (size_t i = 0; i < ev.size(); i++) {
        cout << fixed << setprecision(4) << ev[i] << " ";
    }
    cout << endl;
    
    cout << "Nodes: " << nodes << " Edges: " << edges << endl;
    cout << "Avg. Clustering Coefficient: " << cc_sum / nodes << endl;
    cout << "Diameter: " << diameter << endl;

    cout << "K-hop distribution: ";
    for (map<int, int>::iterator it = khop.begin(); it != khop.end(); ++it) {
        cout << it->second << " ";
    }
    cout << endl;

    cout << "Degree Distribution: ";
    if (!degree_dist.empty()) {
        int max_deg = degree_dist.rbegin()->first;
        for (int i = 0; i <= max_deg; i++) {
            cout << degree_dist[i] << " ";
        }
    }
    cout << endl;

    cout << "nodeName clusCoff eccentricity node_betweenness" << endl;
    for (TIntFltH::TIter it = nbw.BegI(); it != nbw.EndI(); it++) {
        TInt nodeid = it.GetKey();
        TFlt val = it.GetDat();
        cout << id2name[nodeid] << " " << val << endl;
    }

    cout << "node1 node2 edge_betweenness" << endl;
    for (TIntPrFltH::TIter it = ebw.BegI(); it != ebw.EndI(); it++) {
        TIntPr nodepid = it.GetKey();
        TFlt val = it.GetDat();
        cout << id2name[nodepid.GetVal1()] << " : " << id2name[nodepid.GetVal2()] << " " << val << " " << endl;
    }

    cout << "########################################################### End-of-output\n";
}

int main(int argc, char* argv[]) {
    if (argc == 2) {
        getprops(argv[1]);
    } else if (argc == 4 && string(argv[1]) == "-e") {
        int evalues = atoi(argv[2]);
        getprops(argv[3], evalues);
    } else {
        cerr << "Input error! Usage: ./main [-e eigValsToCompute] inputFile" << endl;
        return 1;
    }
    return 0;
}
