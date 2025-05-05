#include "graph.h"
#include <sstream>
#include <set>
#include <cstdlib>

bool buildNameMappings(
    const std::string& filename,
    std::unordered_map<std::string, int>& nameToIndex,
    std::vector<std::string>& indexToName,
    int estimatedNodeCount
) {
    // Reserve space if you know roughly how many nodes to expect
    nameToIndex.reserve(estimatedNodeCount);
    indexToName.reserve(estimatedNodeCount);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file: " << filename << "\n";
        return false;
    }

    int currentIndex = 0;
    std::string nodeName;
    while (file >> nodeName) {
        // If nodeName isn't in the map, add it
        if (nameToIndex.find(nodeName) == nameToIndex.end()) {
            nameToIndex[nodeName] = currentIndex;
            indexToName.push_back(nodeName);
            currentIndex++;
        }
    }
    file.close();

    return true;
}

bool AdjMatrix(const std::string& filename,
    const std::unordered_map<std::string,int>& nameToIndex,
    std::vector<std::vector<bool>>& adjMatrix)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
    std::cerr << "Error: Could not open file: " << filename << "\n";
    return false;
    }

    std::string node1, node2;
    while (file >> node1 >> node2) {
    // Convert each name to its index
    auto it1 = nameToIndex.find(node1);
    auto it2 = nameToIndex.find(node2);
    if (it1 == nameToIndex.end() || it2 == nameToIndex.end()) {
    // This would be unusual if first pass was done correctly
    // handle errors or skip
    std::cerr << "Warning: Node not found in map: " 
       << node1 << " or " << node2 << "\n";
    continue;
    }

    int idx1 = it1->second;
    int idx2 = it2->second;

    // Mark adjacency as true 
    adjMatrix[idx1][idx2] = true;
    adjMatrix[idx2][idx1] = true;
    }
    file.close();
    return true;
}


// Overload: getConnectedNodes(vector<int>, matrix)
std::vector<int> getConnectedNodes(const std::vector<int>& nodeIndices,
                                   const std::vector<std::vector<bool>>& adjMatrix)
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
    const std::vector<std::vector<bool>>& adjMatrix)
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


void mapNamesToIndicesFromFile(
    const std::string& filenameSeed,
    std::unordered_map<std::string, int>& nodeIndexMapping1,
    std::unordered_map<std::string, int>& nodeIndexMapping2,
    std::vector<std::pair<int, int>>& alignmentList
) {
    std::ifstream file(filenameSeed);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filenameSeed << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string name1, name2;

        if (!(iss >> name1 >> name2)) {
            std::cerr << "Error: Invalid line format in " << filenameSeed << ": " << line << std::endl;
            continue;
        }

        // Ensure both names have an index
        if (nodeIndexMapping1.find(name1) == nodeIndexMapping1.end()) {
            int newIndex = nodeIndexMapping1.size();
            nodeIndexMapping1[name1] = newIndex;
        }

        if (nodeIndexMapping2.find(name2) == nodeIndexMapping2.end()) {
            int newIndex = nodeIndexMapping2.size();
            nodeIndexMapping2[name2] = newIndex;
        }

        // Store the alignment as index pairs
        alignmentList.emplace_back(nodeIndexMapping1[name1], nodeIndexMapping2[name2]);
    }

    file.close();
}
// findIndex
int findIndex(const std::unordered_map<std::string,int>& nodeIndexMapping,
    const std::string& nodeName)
{
// Use .find() to locate the key
auto it = nodeIndexMapping.find(nodeName);
if (it == nodeIndexMapping.end()) {
return -1; // or handle "not found"
}
return it->second; // return the integer index
}


void ReadSeed(const std::string& seedFilename,
    const std::unordered_map<std::string,int>& nodeIndexMapping1,
    const std::unordered_map<std::string,int>& nodeIndexMapping2,
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
std::vector<std::vector<float>> ReadSimFile(
    const std::unordered_map<std::string, int>& nodeIndexMapping1,
    const std::unordered_map<std::string, int>& nodeIndexMapping2,
    const std::string& filename)
{
    // Initialize matrix with default value -1 (indicating no similarity)
    size_t size1 = nodeIndexMapping1.size();
    size_t size2 = nodeIndexMapping2.size();
    std::vector<std::vector<float>> similarityMatrix(size1, std::vector<float>(size2, -1.0));

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening similarity file: " << filename << std::endl;
        return similarityMatrix;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string node1_name, node2_name;
        float similarity;
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


void PrintNameToIndex(const std::unordered_map<std::string, int>& nameToIndex)
{
    std::cout << "\n=== nameToIndex Map ===\n";
    for (const auto& pair : nameToIndex) {
        // pair.first is the node name (std::string)
        // pair.second is the index (int)
        std::cout << "Node \"" << pair.first 
                  << "\" => index " << pair.second << '\n';
    }
}

void PrintIndexToName(const std::vector<std::string>& indexToName)
{
    std::cout << "\n=== indexToName Vector ===\n";
    for (int i = 0; i < (int)indexToName.size(); ++i) {
        std::cout << "Index " << i 
                  << " => Node \"" << indexToName[i] << "\"\n";
    }
}

void PrintAdjMatrix(const std::vector<std::vector<bool>>& adjMatrix)
{
    // If you want to show column indices:
    std::cout << "    ";
    for (int col = 0; col < (int)adjMatrix.size(); ++col) {
        std::cout << col << " ";
    }
    std::cout << "\n";

    // Print each row
    for (int row = 0; row < (int)adjMatrix.size(); ++row) {
        // Print row index
        std::cout << row << ":  ";

        for (int col = 0; col < (int)adjMatrix[row].size(); ++col) {
            // Convert the bool to 1 or 0
            std::cout << (adjMatrix[row][col] ? 1 : 0) << " ";
        }
        std::cout << "\n";
    }
}