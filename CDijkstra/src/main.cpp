
#include "main.h"
#include "skiplist0.h"

using namespace std;

// Function to create adjacency matrix and return the node-to-index mapping
vector<pair<string, int>> createAdjacencyMatrix(const string& filename, vector<vector<int>>& adjMatrix) {
    // Set to store unique nodes (ensuring no duplicates)
    set<string> nodes_set;

    // Open the file for reading
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        exit(1);  // Exit if file cannot be opened
    }

    string line;
    vector<pair<string, string>> edges;

    // Reading the file and extracting edges and nodes
    while (getline(file, line)) {
        stringstream ss(line);
        string start_node, end_node;

        ss >> start_node >> end_node;

        // Insert nodes into set (this ensures only unique nodes are stored)
        nodes_set.insert(start_node);
        nodes_set.insert(end_node);

        edges.push_back({ start_node, end_node });
    }

    file.close();

    // Convert the set of nodes to a vector for indexing
    vector<string> nodes(nodes_set.begin(), nodes_set.end());

    // Map nodes to indices (order of nodes is important)
    unordered_map<string, int> node_index;
    for (int i = 0; i < nodes.size(); ++i) {
        node_index[nodes[i]] = i;
    }

    // Store the node name to index mapping in a vector of pairs

    vector<pair<string, int>> nodeIndexMapping;
    for (const auto& pair : node_index) {
        nodeIndexMapping.push_back(pair);  // Add each pair (name, index) to the vector
    }

    // Create the adjacency matrix
    int numNodes = nodes.size();
    adjMatrix = vector<vector<int>>(numNodes, vector<int>(numNodes, 0));

    // Fill the adjacency matrix based on edges
    for (const auto& edge : edges) {
        int start_index = node_index[edge.first];  // Get index for start node
        int end_index = node_index[edge.second];   // Get index for end node

        // Set the corresponding cells in the adjacency matrix for undirected graph
        adjMatrix[start_index][end_index] = 1;  // Edge from start_node to end_node
        adjMatrix[end_index][start_index] = 1;  // Edge from end_node to start_node (symmetry)
    }

    // Print the adjacency matrix
    cout << "\nAdjacency Matrix:" << filename << endl;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }

    return nodeIndexMapping;  // Return the node-to-index mapping for later use
}

// Function to get the connected nodes of a specific node in graph1
vector<int> getConnectedNodes(const vector<int>& nodeIndices, const vector<vector<int>>& adjMatrix) {
    vector<int> connectedNodes;
    for (int nodeIndex : nodeIndices) {
        // Check the row corresponding to the node index in graph1
        for (int i = 0; i < adjMatrix[nodeIndex].size(); ++i) {
            // If there's a connection (1) between nodeIndex and node i, add i to the connected nodes list
            if (adjMatrix[nodeIndex][i] == 1) {
                connectedNodes.push_back(i);
            }
        }
    }
    return connectedNodes;
}
vector<int> getConnectedNodes(int nodeIndex, const vector<vector<int>>& adjMatrix) {
    vector<int> connectedNodes;

    // Check the row corresponding to the node index in graph1
    for (int i = 0; i < adjMatrix[nodeIndex].size(); ++i) {
        // If there's a connection (1) between nodeIndex and node i, add i to the connected nodes list
        if (adjMatrix[nodeIndex][i] == 1) {
            connectedNodes.push_back(i);
        }
    }
    return connectedNodes;
}


// Function to manually search for a node name in the vector and return its index, or -1 if not found
int findIndex(const vector<pair<string, int>>& nodeIndexMapping, const string& nodeName) {
    for (int i = 0; i < nodeIndexMapping.size(); ++i) {
        if (nodeIndexMapping[i].first == nodeName) {
            return nodeIndexMapping[i].second;  // Return the index if the node is found
        }
    }
    return -1;  // Return -1 if the node name is not found
}
void mapNamesToIndices(const string& filename, unordered_map<string, int>& nodeIndexMapping, vector<pair<int, int>>& alignmentList) {
    // Open the seed file for reading
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        exit(1);  // Exit if file cannot be opened
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string node1_name, node2_name;

        // Read two node names (assuming each line has two names separated by a space)
        ss >> node1_name >> node2_name;

        // Check if the node names exist in the mapping
        if (nodeIndexMapping.find(node1_name) != nodeIndexMapping.end() && nodeIndexMapping.find(node2_name) != nodeIndexMapping.end()) {
            int node1_index = nodeIndexMapping[node1_name];  // Get index for the first node
            int node2_index = nodeIndexMapping[node2_name];  // Get index for the second node

            // Add the pair (node1_index, node2_index) to the alignment list
            alignmentList.push_back({ node1_index, node2_index });
        }
        else {
            cout << "Node names " << node1_name << " or " << node2_name << " not found in the index mapping." << endl;
        }
    }

    file.close();
}

void mapNamesToIndicesFromFile(const string& filename,
    const vector<pair<string, int>>& nodeIndexMapping1,
    const vector<pair<string, int>>& nodeIndexMapping2,
    vector<pair<int, int>>& alignmentList) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        return;  // Exit if file cannot be opened
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string node1_name, node2_name;
        ss >> node1_name >> node2_name;

        // Map the node names to indices
        int node1_index = findIndex(nodeIndexMapping1, node1_name);
        int node2_index = findIndex(nodeIndexMapping2, node2_name);

        // If both nodes are found in the respective mappings, store the pair in alignmentList
        if (node1_index != -1 && node2_index != -1) {
            alignmentList.push_back({ node1_index, node2_index });
        }
        else {
            cout << "Node names not found in the index mapping: "
                << node1_name << ", " << node2_name << endl;
        }
    }

    file.close();
}

ostream& operator<<(ostream& os, const vector<int>& vec) {
    for (const int& element : vec) {
        os << element << " ";  // Print each element of the vector followed by a space
    }
    return os;
}

ostream& operator<<(ostream& os, const pair<int, int>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

// Overload the << operator for vector<pair<int, int>>
ostream& operator<<(ostream& os, const vector<pair<int, int>>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1) {  // Add a comma between elements, but not after the last one
            os << ", ";
        }
    }
    os << "]";
    return os;
}
void DisplayNodetoIndex(vector<pair<string, int>> nodeIndexMapping, const string& filename) {
    cout << "\nNode Name to Index Mapping :" << filename << endl;
    for (const auto& pair : nodeIndexMapping) {
        cout << "Node " << pair.first << " has Index: " << pair.second << endl;
    }
}

void displayConnectedNodes(int SeedNodeGraph1, const vector<int>& connectedNodes1,
    int SeedNodeGraph2, const vector<int>& connectedNodes2) {
    // Display the connected nodes for Graph 1
    cout << "\nConnected Nodes of Node " << SeedNodeGraph1 << " (Graph1):" << endl;
    for (int connectedNode : connectedNodes1) {
        cout << connectedNode << " ";  // Output the index of the connected node
    }
    cout << endl;

    // Display the connected nodes for Graph 2
    cout << "\nConnected Nodes of Node " << SeedNodeGraph2 << " (Graph2):" << endl;
    for (int connectedNode : connectedNodes2) {
        cout << connectedNode << " ";  // Output the index of the connected node
    }
    cout << endl;
}
void displayConnectedNodes(const vector<int>& connectedNodes) {
    for (int connectedNode : connectedNodes) {
        cout << connectedNode << " ";  // Output the index of the connected node
    }
    cout << endl;
}

vector<tuple<int, int, double>>ReadSimFIle(const string& filename, const vector<pair<string, int>>& nodeIndexMapping1, const vector<pair<string, int>>& nodeIndexMapping2) {
    //Step1:read sim file
    ifstream similarityFile(filename);
    if (!similarityFile) {
        cerr << "Error opening similarity file!" << endl;
        exit(1);  // Exit if the file cannot be opened
    }
    //Step2: convert the content in sim file to vector<tuple<int, int, double>> similarityPairs
    string node1_name, node2_name;
    double similarity;
    //create a vector of tuple to store the Pairs in Sim file for comparision
    vector<tuple<int, int, double>> similarityPairs;
    while (similarityFile >> node1_name >> node2_name >> similarity) {
        // Convert the node names to their corresponding indices using the map
         // Use the findIndex function to get the node indices from the mappings
        int node1_index = findIndex(nodeIndexMapping1, node1_name);
        int node2_index = findIndex(nodeIndexMapping2, node2_name);

        // If both nodes exist in their respective mappings
        if (node1_index != -1 && node2_index != -1) {
            // Add the pair and its similarity value to the similarityPairs vector
            similarityPairs.push_back(make_tuple(node1_index, node2_index, similarity));
        }
    }
    similarityFile.close();
    cout << "\nNode Pairs with Similarity Values:" << endl;
    for (const auto& simPair : similarityPairs) {
        cout << "(" << get<0>(simPair) << ", " << get<1>(simPair) << ") has similarity " << get<2>(simPair) << endl;
    }
    return similarityPairs;
}

int main() {
    //input file names that will be used: include g1, g2, seed and sim
    string filenameGraph1 = "graph1.txt";  // read from file(may need use .el file)
    string filenameGraph2 = "graph2.txt";
    string filenameSeed = "seed.txt";  // Replace with your actual file path
    string filenameSim = "sim.txt";
    

    //store the g1 and g2 into adjMatrix
    vector<vector<int>> adjMatrix1;  // Adjacency matrix will be filled inside createAdjacencyMatrix
    vector<pair<string, int>> nodeIndexMapping1 = createAdjacencyMatrix(filenameGraph1, adjMatrix1);
    vector<vector<int>> adjMatrix2;  // Adjacency matrix will be filled inside createAdjacencyMatrix
    vector<pair<string, int>> nodeIndexMapping2 = createAdjacencyMatrix(filenameGraph2, adjMatrix2);


    // Access the node-to-index mapping
    DisplayNodetoIndex(nodeIndexMapping1,filenameGraph1);
    DisplayNodetoIndex(nodeIndexMapping2,filenameGraph2);
    

    // Create a new alignment list
    vector<pair<int, int>> alignmentList;
    //vector stores the integer representation of seed node pairs
    mapNamesToIndicesFromFile(filenameSeed, nodeIndexMapping1, nodeIndexMapping2, alignmentList);
    
    //vector stores the integer representation of seedNode of g1 and g2
    vector<int>SeedNodeGraph1;
    vector<int>SeedNodeGraph2;

    // Output the alignmentList
    cout << "Alignment List: " << endl;
    //store the integer representation of seedNode of g1 and g2 
    // for computing the Natural Product of possible new node pairs
    for (const auto& pair : alignmentList) {
        SeedNodeGraph1.push_back(pair.first);
        SeedNodeGraph2.push_back(pair.second);
        cout << "(" << pair.first << ", " << pair.second << ")" << endl;
    }
    

    //computing the possible new node pairs(Natural Product)
    // Step1
    //SeedNodeGraph1&2 reference to the adjacency matrix to search for neighbours of aligned nodes
    vector<int> connectedNodes1 = getConnectedNodes(SeedNodeGraph1, adjMatrix1);
    vector<int> connectedNodes2 = getConnectedNodes(SeedNodeGraph2, adjMatrix2);

    // Display the connected nodes
    cout << "\nConnected Nodes of Seed " << SeedNodeGraph1 << " (Graph1):" << endl;
    displayConnectedNodes(connectedNodes1);
    cout << endl;
    cout << "\nConnected Nodes of Seed " << SeedNodeGraph2 << " (Graph2):" << endl;
    displayConnectedNodes(connectedNodes2);
    cout << endl;
    //Step2 Natural product
    vector<pair<int, int>> nodePairs;
    for (int node1 : connectedNodes1) {
        for (int node2 : connectedNodes2) {
            nodePairs.push_back({ node1, node2 });  // Add pair (node1, node2) to the list
        }
    }
    // Display the node pairs
    cout << "\nNode Pairs (Natural Product of Connected Nodes):" << endl;
    for (const auto& pair : nodePairs) {
        cout << "(" << pair.first << ", " << pair.second << ")" << endl;
    }

    //Check if the Natural Product node pairs exist in the sim file(only care about those have a sim)
    //Step1:read sim file and convert it to interger representation 
    vector<tuple<int, int, double>> similarityPairs = ReadSimFIle(filenameSim, nodeIndexMapping1, nodeIndexMapping2);
    //Step2: check for the matched node pairs
    vector<tuple<int, int, double>> matchedTuples;
    cout << "\nSimilarity Pairs found in Node Pairs:" << endl;
    // Iterate over each similarity pair: MIGHT BE OPTIMIZADE
    //All pairs find in the similarityPairs should be stored in the skiplist
    //some error happened when passing it to the skiplist so i just skip that part
    //need modification; passing tuple should be fine, unencapusulate it in the skiplist
    for (const auto& nodePair : nodePairs) {
        for (const auto& simPair : similarityPairs) {
            // Check if the pair matches the one in the nodePairs vector
            if (get<0>(simPair) == nodePair.first && get<1>(simPair) == nodePair.second) {
                matchedTuples.push_back(make_tuple(nodePair.first, nodePair.second, get<2>(simPair)));
                cout << "(" << nodePair.first << ", " << nodePair.second << ") with similarity " << get<2>(simPair) << endl;
            }
        }
    }
    //add pair: extract candidate pair from skiplist and add it to the alignmentList
    //newAlignPair=skipList.node//for now i will use the one i find in the similarity file directly
    //read from vector of the tuple and for each element run recursively

    //add one candidate from the skiplist to the alignmentList
    for (const auto& matchedTuple : matchedTuples) {
        int first = get<0>(matchedTuple);
        int second = get<1>(matchedTuple);
        alignmentList.push_back(make_pair(first, second));
        cout << alignmentList;
        //calculate new neighbours: some might be in the aligned nodes already. Be careful for them
        vector<int> newConnectedNodes1 = getConnectedNodes(first, adjMatrix1);
        vector<int> newConnectedNodes2 = getConnectedNodes(second, adjMatrix2);

        // Display the new connected nodes
        cout << "\nConnected Nodes of Node " << SeedNodeGraph1 << " (Graph1):" << endl;
        displayConnectedNodes(newConnectedNodes1);
        cout << endl;
        cout << "\nConnected Nodes of Node " << SeedNodeGraph2 << " (Graph2):" << endl;
        displayConnectedNodes(newConnectedNodes2);
        cout << endl;

        //calculate the natural product new connectedNodes
        vector<pair<int, int>> newNodePairs;
        for (int node1 : newConnectedNodes1) {
            for (int node2 : newConnectedNodes2) {
                newNodePairs.push_back({ node1, node2 });  
            }
        }
        //output the new connectedNodes
        cout << "\nNode Pairs (Natural Product of Connected Nodes):" << endl;
        for (const auto& pair : newNodePairs) {
            cout << "(" << pair.first << ", " << pair.second << ")" << endl;
        }
    }


//run over the process again: stoppoint?

    return 0;
}



