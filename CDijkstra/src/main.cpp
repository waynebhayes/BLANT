#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>  
#include <chrono>
#include "graph.h"
#include "skiplist.h" 

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // Input file names
    std::string Graph1 = "test/graph1.el";
    std::string Graph2 = "test/graph2.el";
    std::string Seed   = "test/seed.txt";
    std::string Sim    = "test/sim.txt";
    int numNodes1=10000;//number of nodes in graph 1
    int numNodes2=10000;//number of nodes in graph 2

    // 1. read graph and map string to int
    std::unordered_map<std::string, int> nameToIndex1;
    std::vector<std::string> indexToName1;
    if (!buildNameMappings(Graph1, nameToIndex1, indexToName1, numNodes1)) {
        std::cerr << "Failed building node mapping.\n";
        return 1;
    }

    std::unordered_map<std::string, int> nameToIndex2;
    std::vector<std::string> indexToName2;
    if (!buildNameMappings(Graph2, nameToIndex2, indexToName2, numNodes2)) {
        std::cerr << "Failed building node mapping.\n";
        return 1;
    }

    //2. read graph again to build adjacency matrix
    numNodes1 = nameToIndex1.size();
    std::vector<std::vector<bool>> adjMatrix1(
        numNodes1,
        std::vector<bool>(numNodes1, false)
    );
    if (!AdjMatrix(Graph1, nameToIndex1, adjMatrix1)) {
        std::cerr << "Failed to fill adjacency matrix.\n";
        return 1;
    }

    numNodes2 = nameToIndex2.size();
    std::vector<std::vector<bool>> adjMatrix2(
        numNodes2,
        std::vector<bool>(numNodes2, false)
    );
    if (!AdjMatrix(Graph2, nameToIndex2, adjMatrix2)) {
        std::cerr << "Failed to fill adjacency matrix.\n";
        return 1;
    }


    // 2. Display node->index mappings
    PrintNameToIndex(nameToIndex1);
    PrintNameToIndex(nameToIndex2);
    //    Display index->node mappings
    PrintIndexToName(indexToName1);
    PrintIndexToName(indexToName2);
    //    Display the matrix
    PrintAdjMatrix(adjMatrix1);
    PrintAdjMatrix(adjMatrix2);

    // 3. Read seed alignments
    std::vector<std::pair<int,int>> alignmentList;
    std::cout << "\nAlignment List: \n";
    for (auto& p : alignmentList) {
        std::cout << "(" << p.first << ", " << p.second << ")\n";
    }

    // 4. Separate aligned nodes for adjacency lookups
    std::vector<int> SeedNodeGraph1, SeedNodeGraph2;
    for (auto& p : alignmentList) {
        SeedNodeGraph1.push_back(p.first);
        SeedNodeGraph2.push_back(p.second);
    }

    // 5. Initialize alignment process
    bool alignmentInProgress = true;
    SkipList skiplist(20, 0.5);
    int iterationCount = 0;
    const int maxIterations = 10;
    // Main loop: continue until no candidates meet the threshold
    while (alignmentInProgress&& iterationCount < maxIterations) {
        iterationCount++;
        // 6. Get connected neighbors for seeds
        std::vector<int> connectedNodes1 = getConnectedNodes(SeedNodeGraph1, adjMatrix1);
        std::vector<int> connectedNodes2 = getConnectedNodes(SeedNodeGraph2, adjMatrix2);

        std::cout << "\nConnected Nodes of Seed Graph 1:\n";
        displayConnectedNodes(connectedNodes1);

        std::cout << "\nConnected Nodes of Seed Graph 2:\n";
        displayConnectedNodes(connectedNodes2);

        // 7. Natural product of neighbor sets
        std::vector<std::pair<int, int>> candidatePairs;
        for (int n1 : connectedNodes1) {
            for (int n2 : connectedNodes2) {
                candidatePairs.push_back({n1, n2});
            }
        }

        // Step 2: Filter out pairs that are already in the alignment list
        candidatePairs.erase(
            std::remove_if(candidatePairs.begin(), candidatePairs.end(),
                [&](const std::pair<int, int>& pair) {
                    return std::find(alignmentList.begin(), alignmentList.end(), pair) != alignmentList.end();
                }),
            candidatePairs.end()
        );

        std::cout << "\nFiltered Node Pairs (Natural Product):\n";
        for (auto& np : candidatePairs) {
            std::cout << "(" << np.first << ", " << np.second << ")\n";
        }

        // 8. Read similarity file
        std::vector<std::vector<double>> similarityMatrix = ReadSimFile(nameToIndex1, nameToIndex2, Sim);


        // 9. Add matching pairs to skip list
        for (auto& np : candidatePairs) {
        double similarity = similarityMatrix[np.first][np.second];

        if (similarity > 0) {  // Only insert if similarity exists and is positive
            std::cout << "Inserting: (" << np.first << ", " << np.second << ") with similarity " << similarity << "\n";
            skiplist.insertElement(similarity, np.first, np.second);
        }
        }
        // 10. Display skip list
        skiplist.displayList();

        // 11. Pop one candidate from skip list
        auto tup = skiplist.pop(0.1);
        double key;
        int first, second;
        std::tie(key, first, second) = tup;
        // Check if a candidate was found
        if (key != -1.0) {
            std::cout << "Popped node => key: " << key
                      << ", VertexA: " << first
                      << ", VertexB: " << second << "\n";

            // Add to alignment list
            alignmentList.push_back({ first, second });
            std::cout << "[DEBUG] Using maxNode: (" << first << ", " << second
                      << ") with value=" << key << "\n";
            // Update seeds for the next iteration
            SeedNodeGraph1.push_back(first);
            SeedNodeGraph2.push_back(second);
        } else {
            // No more candidates, stop the loop
            std::cout << "[DEBUG] Skiplist is empty, no candidate pairs.\n";
            alignmentInProgress = false;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\nExecution Time: " << elapsed.count() << " seconds\n";

    // Final output
    std::cout << "\nFinal Alignment List:\n";
    for (const auto& pair : alignmentList) {
        std::cout << "(" << pair.first << ", " << pair.second << ")\n";
    }

    return 0;
}