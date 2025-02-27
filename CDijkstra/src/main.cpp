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
    std::string filenameGraph1 = "rsrc/graph1.el";
    std::string filenameGraph2 = "rsrc/graph2.el";
    std::string filenameSeed   = "rsrc/seed.txt";
    std::string filenameSim    = "rsrc/sim.txt";

    // 1. Build adjacency matrices
    std::vector<std::vector<int>> adjMatrix1;
    std::vector<std::pair<std::string,int>> nodeIndexMapping1 =
        createAdjacencyMatrix(filenameGraph1, adjMatrix1);

    std::vector<std::vector<int>> adjMatrix2;
    std::vector<std::pair<std::string,int>> nodeIndexMapping2 =
        createAdjacencyMatrix(filenameGraph2, adjMatrix2);

    // 2. Display node→index mappings
    DisplayNodetoIndex(nodeIndexMapping1, filenameGraph1);
    DisplayNodetoIndex(nodeIndexMapping2, filenameGraph2);

    // 3. Read seed alignments
    std::vector<std::pair<int,int>> alignmentList;
    mapNamesToIndicesFromFile(filenameSeed, nodeIndexMapping1, nodeIndexMapping2, alignmentList);

    std::cout << "\nAlignment List: \n";.…
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
        std::vector<std::vector<double>> similarityMatrix = ReadSimFile(nodeIndexMapping1, nodeIndexMapping2, filenameSim);


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