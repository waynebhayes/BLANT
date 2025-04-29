#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <fstream>



bool buildNameMappings(
    const std::string& filename,
    std::unordered_map<std::string, int>& nameToIndex,
    std::vector<std::string>& indexToName,
    int estimatedNodeCount
);

bool AdjMatrix(const std::string& filename,
    const std::unordered_map<std::string,int>& nameToIndex,
    std::vector<std::vector<bool>>& adjMatrix);

// Get connected nodes (overload for a list of nodes)
std::vector<int> getConnectedNodes(const std::vector<int>& nodeIndices,
                                   const std::vector<std::vector<bool>>& adjMatrix);

// Get connected nodes (overload for a single node)


// Find a node index by name
int findIndex(const std::unordered_map<std::string,int>& nodeIndexMapping,
              const std::string& nodeName);

// Map names to indices from a file
void ReadSeed(const std::string& seedFilename,
    const std::unordered_map<std::string,int>& nodeIndexMapping1,
    const std::unordered_map<std::string,int>& nodeIndexMapping2,
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
    const std::unordered_map<std::string, int>& nodeIndexMapping1,
    const std::unordered_map<std::string, int>& nodeIndexMapping2,
    const std::string& filename);
std::vector<int> getConnectedNodes(int nodeIndex,
    const std::vector<std::vector<bool>>& adjMatrix);

// Overloaded operators for printing
std::ostream& operator<<(std::ostream& os, const std::vector<int>& vec);
std::ostream& operator<<(std::ostream& os, const std::pair<int, int>& p);
std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<int, int>>& v);

void PrintNameToIndex(const std::unordered_map<std::string, int>& nameToIndex);
void PrintIndexToName(const std::vector<std::string>& indexToName);
void PrintAdjMatrix(const std::vector<std::vector<bool>>& adjMatrix);

std::vector<std::vector<double>> ReadSimFile(
    const std::unordered_map<std::string, int>& nodeIndexMapping1,
    const std::unordered_map<std::string, int>& nodeIndexMapping2,
    const std::string& filename);