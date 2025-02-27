#include "graph.h"
#include <fstream>
#include <sstream>
#include <set>
#include <cstdlib>

// Definition: createAdjacencyMatrix
std::vector<std::pair<std::string, int>>
createAdjacencyMatrix(const std::string& filename,
                      std::vector<std::vector<int>>& adjMatrix)
{
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: file not found: " << filename << std::endl;
        exit(1); 
    }

    std::set<std::string> nodes_set;
    std::vector<std::pair<std::string, std::string>> edges;

    // Read edges
    {
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string start_node, end_node;
            ss >> start_node >> end_node;
            nodes_set.insert(start_node);
            nodes_set.insert(end_node);
            edges.push_back({start_node, end_node});
        }
    }
    file.close();

    // Build node-index map
    std::vector<std::string> nodes(nodes_set.begin(), nodes_set.end());
    std::unordered_map<std::string, int> node_index;
    node_index.reserve(nodes.size());
    for (int i = 0; i < (int)nodes.size(); ++i) {
        node_index[nodes[i]] = i;
    }

    // Create the adjacency matrix
    int numNodes = (int)nodes.size();
    adjMatrix.assign(numNodes, std::vector<int>(numNodes, 0));

    // Fill adjacency
    for (auto& edge : edges) {
        int start_i = node_index[edge.first];
        int end_i   = node_index[edge.second];
        adjMatrix[start_i][end_i] = 1;
        adjMatrix[end_i][start_i] = 1;
    }

    // Debug print
    /*
    std::cout << "\nAdjacency Matrix: " << filename << std::endl;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            std::cout << adjMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    

    // Convert node_index to a vector of (name, index) pairs
    std::vector<std::pair<std::string, int>> nodeIndexMapping;
    nodeIndexMapping.reserve(node_index.size());
    for (auto& it : node_index) {
        nodeIndexMapping.push_back({ it.first, it.second });
    }

    return nodeIndexMapping;
}

// Overload: getConnectedNodes(vector<int>, matrix)
std::vector<int> getConnectedNodes(const std::vector<int>& nodeIndices,
                                   const std::vector<std::vector<int>>& adjMatrix)
{
        std::set<int> uniqueConnectedNodes;  // Use a set to ensure uniqueness

        for (int node : nodeIndices) {
            for (size_t i = 0; i < adjMatrix[node].size(); ++i) {
                if (adjMatrix[node][i] != 0) {  // Check if there's an edge
                    uniqueConnectedNodes.insert(i);  // Insert ensures no duplicates
                }
            }
        }

    // Convert set back to vector before returning
    return std::vector<int>(uniqueConnectedNodes.begin(), uniqueConnectedNodes.end());
}
std::vector<int> getConnectedNodes(int nodeIndex,
    const std::vector<std::vector<int>>& adjMatrix)
{
std::set<int> uniqueConnectedNodes;  // Use a set to ensure uniqueness

// Traverse the adjacency row for the single node
for (size_t i = 0; i < adjMatrix[nodeIndex].size(); ++i) {
if (adjMatrix[nodeIndex][i] != 0) {  // Check if there's an edge
uniqueConnectedNodes.insert(i);  // Insert ensures no duplicates
}
}

// Convert set back to vector before returning
return std::vector<int>(uniqueConnectedNodes.begin(), uniqueConnectedNodes.end());
}



// findIndex
int findIndex(const std::vector<std::pair<std::string, int>>& nodeIndexMapping,
              const std::string& nodeName)
{
    for (int i = 0; i < (int)nodeIndexMapping.size(); ++i) {
        if (nodeIndexMapping[i].first == nodeName) {
            return nodeIndexMapping[i].second;
        }
    }
    return -1;
}

void ReadSeed(const std::string& seedFilename,
    const std::vector<std::pair<std::string,int>>& nodeIndexMapping1,
    const std::vector<std::pair<std::string,int>>& nodeIndexMapping2,
    std::vector<int>& SeedNodeGraph1,
    std::vector<int>& SeedNodeGraph2)
{
std::ifstream file(seedFilename);
if (!file) {
std::cerr << "Error opening file: " << seedFilename << std::endl;
return;  // or exit(1), depending on your preference
}

std::string line;
while (std::getline(file, line)) {
if (line.empty()) {
  continue;  // skip empty lines
}
std::stringstream ss(line);
std::string node1_name, node2_name;
ss >> node1_name >> node2_name;

// Find each name's index using the vectors
int idx1 = findIndex(nodeIndexMapping1, node1_name);
int idx2 = findIndex(nodeIndexMapping2, node2_name);

if (idx1 != -1 && idx2 != -1) {
  // Valid indices, so append them to the seed vectors
  SeedNodeGraph1.push_back(idx1);
  SeedNodeGraph2.push_back(idx2);
} else {
  // If either name is not found, warn and skip
  std::cerr << "Warning: node name not found in mapping(s): ("
            << node1_name << ", " << node2_name << ")\n";
}
}

file.close();
}


// DisplayNodetoIndex
void DisplayNodetoIndex(std::vector<std::pair<std::string, int>> nodeIndexMapping,
                        const std::string& filename)
{
    std::cout << "\nNode Name->Index Mapping: " << filename << std::endl;
    for (auto& p : nodeIndexMapping) {
        std::cout << "Node " << p.first << " has Index " << p.second << std::endl;
    }
}

// displayConnectedNodes (two sets)
void displayConnectedNodes(int SeedNodeGraph1, const std::vector<int>& connectedNodes1,
                           int SeedNodeGraph2, const std::vector<int>& connectedNodes2)
{
    std::cout << "\nConnected Nodes of Node " << SeedNodeGraph1 << ":\n";
    for (auto cn : connectedNodes1) {
        std::cout << cn << " ";
    }
    std::cout << "\n\nConnected Nodes of Node " << SeedNodeGraph2 << ":\n";
    for (auto cn : connectedNodes2) {
        std::cout << cn << " ";
    }
    std::cout << std::endl;
}

// displayConnectedNodes (single set)
void displayConnectedNodes(const std::vector<int>& connectedNodes)
{
    for (auto cn : connectedNodes) {
        std::cout << cn << " ";
    }
    std::cout << std::endl;
}

// ReadSimFIle
std::vector<std::vector<double>> ReadSimFile(
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping1,
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping2,
    const std::string& filename)
{
    // Initialize matrix with default value -1 (indicating no similarity)
    size_t size1 = nodeIndexMapping1.size();
    size_t size2 = nodeIndexMapping2.size();
    std::vector<std::vector<double>> similarityMatrix(size1, std::vector<double>(size2, -1.0));

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening similarity file: " << filename << std::endl;
        return similarityMatrix;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string node1_name, node2_name;
        double similarity;
        ss >> node1_name >> node2_name >> similarity;

        int idx1 = findIndex(nodeIndexMapping1, node1_name);
        int idx2 = findIndex(nodeIndexMapping2, node2_name);

        if (idx1 != -1 && idx2 != -1) {
            similarityMatrix[idx1][idx2] = similarity;
        } else {
            std::cerr << "Warning: Node names not found for similarity: " 
                      << node1_name << ", " << node2_name << std::endl;
        }
    }

    file.close();
    return similarityMatrix;
}



// Overloads for printing
std::ostream& operator<<(std::ostream& os, const std::vector<int>& vec)
{
    for (auto el : vec) {
        os << el << " ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::pair<int,int>& p)
{
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os,
                         const std::vector<std::pair<int,int>>& v)
{
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << "(" << v[i].first << ", " << v[i].second << ")";
        if (i + 1 < v.size()) os << ", ";
    }
    os << "]";
    return os;
}
