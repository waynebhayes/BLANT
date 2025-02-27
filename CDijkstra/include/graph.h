#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <tuple>

// Declarations of all functions previously in your single .cpp (except main)

// Build adjacency matrix
std::vector<std::pair<std::string, int>>
createAdjacencyMatrix(const std::string& filename,
                      std::vector<std::vector<int>>& adjMatrix);

// Get connected nodes (overload for a list of nodes)
std::vector<int> getConnectedNodes(const std::vector<int>& nodeIndices,
                                   const std::vector<std::vector<int>>& adjMatrix);

// Get connected nodes (overload for a single node)


// Find a node index by name
int findIndex(const std::vector<std::pair<std::string, int>>& nodeIndexMapping,
              const std::string& nodeName);

// Map names to indices from a file
void ReadSeed(const std::string& seedFilename,
    const std::vector<std::pair<std::string,int>>& nodeIndexMapping1,
    const std::vector<std::pair<std::string,int>>& nodeIndexMapping2,
    std::vector<int>& SeedNodeGraph1,
    std::vector<int>& SeedNodeGraph2);

void mapNamesToIndicesFromFile(const std::string& filename,
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping1,
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping2,
    std::vector<std::pair<int, int>>& alignmentList);

// Display node→index mapping
void DisplayNodetoIndex(std::vector<std::pair<std::string, int>> nodeIndexMapping,
                        const std::string& filename);

// Display connected nodes (two sets at once)
void displayConnectedNodes(int SeedNodeGraph1, const std::vector<int>& connectedNodes1,
                           int SeedNodeGraph2, const std::vector<int>& connectedNodes2);

// Display connected nodes (single set)
void displayConnectedNodes(const std::vector<int>& connectedNodes);

// Read similarity file and return vector of (idx1, idx2, similarity)
std::vector<std::vector<double>> ReadSimFile(
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping1,
    const std::vector<std::pair<std::string, int>>& nodeIndexMapping2,
    const std::string& filename);
std::vector<int> getConnectedNodes(int nodeIndex,
    const std::vector<std::vector<int>>& adjMatrix);
// Overloaded operators for printing
std::ostream& operator<<(std::ostream& os, const std::vector<int>& vec);
std::ostream& operator<<(std::ostream& os, const std::pair<int, int>& p);
std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<int, int>>& v);


