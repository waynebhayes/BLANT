#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>  
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <unordered_map>
#include <unordered_set>
#include "graph.h"
#include "skiplist.h" 

int main()
{
    auto start_total = std::chrono::high_resolution_clock::now();
    // Input file names
    std::string Graph1 = "test/RNorvegicus.el";
    std::string Graph2 = "test/SPombe.el";
    std::string Seed   = "test/RNorvegicus-SPombe-Seed.txt";
    std::string Sim    = "test/RNorvegicus-SPombe.sim";
    int numNodes1=2330;//number of nodes in graph 1
    int numNodes2=4711;//number of nodes in graph 2

    auto start_mapping = std::chrono::high_resolution_clock::now();
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
    auto end_mapping = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_mapping = end_mapping - start_mapping;
    std::cout << "[TIME] Graph mapping: " << elapsed_mapping.count() << " seconds\n";

    auto start_adj_matrix = std::chrono::high_resolution_clock::now();
    //2. read graph again to build adjacency matrix
    numNodes1 = nameToIndex1.size();
    std::vector<std::vector<bool>> adjMatrix1(numNodes1,std::vector<bool>(numNodes1, false));
    if (!AdjMatrix(Graph1, nameToIndex1, adjMatrix1)) {
        std::cerr << "Failed to fill adjacency matrix.\n";
        return 1;
    }

    numNodes2 = nameToIndex2.size();
    std::vector<std::vector<bool>> adjMatrix2(numNodes2,std::vector<bool>(numNodes2, false));
    if (!AdjMatrix(Graph2, nameToIndex2, adjMatrix2)) {
        std::cerr << "Failed to fill adjacency matrix.\n";
        return 1;
    }

    auto end_adj_matrix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_adj_matrix = end_adj_matrix - start_adj_matrix;
    std::cout << "[TIME] Adjacency matrix construction: " << elapsed_adj_matrix.count() << " seconds\n";

    // 2. Display node->index mappings
    //PrintNameToIndex(nameToIndex1);
    //PrintNameToIndex(nameToIndex2);
    //    Display index->node mappings
    //PrintIndexToName(indexToName1);
    //PrintIndexToName(indexToName2);
    //    Display the matrix
    //PrintAdjMatrix(adjMatrix1);
    //PrintAdjMatrix(adjMatrix2);
    
    // 3. read seed
    std::vector<int> SeedNodeGraph1;
    std::vector<int> SeedNodeGraph2;
    ReadSeed(Seed,nameToIndex1,nameToIndex2, SeedNodeGraph1,SeedNodeGraph2);
    std::cout<<SeedNodeGraph1;
    std::cout<<SeedNodeGraph2;
    //// 4. Read similarity file
    auto start_sim = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<float>> similarityMatrix = ReadSimFile(nameToIndex1, nameToIndex2, Sim);
    auto end_sim = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sim = end_sim - start_sim;
    std::cout << "[TIME] Reading similarity file: " << elapsed_sim.count() << " seconds\n";
    
    std::ofstream outFile("output.txt",std::ios::app);
    if (!outFile.is_open()) {
        std::cerr << "Error opening output.txt\n";
        return 1;
    }
    outFile << "===== file reading time =====\n";
    outFile << "Graph mapping time: " << elapsed_mapping.count() << " seconds\n";
    outFile << "Adjacency matrix construction time: " << elapsed_adj_matrix.count() << " seconds\n";
    outFile << "Similarity file reading time: " << elapsed_sim.count() << " seconds\n\n";


    // 5. Initialize alignment process
    auto start_alignment = std::chrono::high_resolution_clock::now();
    bool alignmentInProgress = true;
    SkipList skiplist(20, 0.5);
    int iterationCount = 0;
    //const int maxIterations = 10;
    std::vector<int> connectedNodes1;
    std::vector<int> connectedNodes2;
    // Main loop: continue until no candidates meet the threshold
    while (alignmentInProgress) {
        if(iterationCount == 0){
            connectedNodes1 = getConnectedNodes(SeedNodeGraph1, adjMatrix1);
            connectedNodes2 = getConnectedNodes(SeedNodeGraph2, adjMatrix2);
        }

        //std::cout << "\nConnected Nodes of Seed Graph 1:\n";
        //displayConnectedNodes(connectedNodes1);
        //std::cout << "\nConnected Nodes of Seed Graph 2:\n";
        //displayConnectedNodes(connectedNodes2);

        //7. Calculating natural product and insert into skiplist
        //   throw away already aligned nodes when doing so
        for (int n1 : connectedNodes1)
        {
        // Check if n1 is already in the seed set
        if (std::find(SeedNodeGraph1.begin(), SeedNodeGraph1.end(), n1) != SeedNodeGraph1.end())
            continue;

            for (int n2 : connectedNodes2)
            {
            // Check if n1 is already in the seed set
            if (std::find(SeedNodeGraph1.begin(), SeedNodeGraph1.end(), n1) != SeedNodeGraph1.end()) 
                continue;
            double sim = similarityMatrix[n1][n2];
            //std::cout << "Candidate pair (" << n1 << ", " << n2 << ") has similarity " << sim << "\n";
            if (sim >0){
                //std::cout << "Inserting: (" << n1 << ", " << n2 << ") with similarity " << sim << "\n";
                skiplist.insertElement(sim, n1, n2);
                }
            }
        }

<<<<<<< HEAD
        // 8. Display skip list
        //skiplist.displayList();

        // 9. Pop one candidate from skip list
        auto tup = skiplist.pop(0.1);  // returns (key, first, second)
        double key;
        int first, second;
        std::tie(key, first, second) = tup;


        if (key != -1.0) {
        //std::cout << "Popped node => key: " << key
              //<< ", VertexA: " << first
              //<< ", VertexB: " << second << "\n";

        //10. Add to aligned list
        SeedNodeGraph1.push_back(first);
        SeedNodeGraph2.push_back(second);
        //std::cout << "[DEBUG] Using maxNode: (" << first << ", " << second
        //      << ") with value=" << key << "\n";
        connectedNodes1 = getConnectedNodes(first, adjMatrix1);
        connectedNodes2 = getConnectedNodes(second, adjMatrix2);
    } else {
    // No more candidates, stop the loop
    std::cout << "[DEBUG] Skiplist is empty, no candidate pairs.\n";
    alignmentInProgress = false;
    }
    iterationCount++;
}

auto end_alignment = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed_alignment = end_alignment - start_alignment;
std::cout << "[TIME] Alignment process: " << elapsed_alignment.count() << " seconds\n";


    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;
    //std::cout << "\nExecution Time: " << elapsed.count() << " seconds\n";

    std::cout << "\nFinal Alignment (name->name):\n";
    outFile << "Final Alignment (name -> name):\n";
    for (size_t i = 0; i < SeedNodeGraph1.size(); ++i) {
        // Indices in each graph
    int idx1 = SeedNodeGraph1[i]; 
    int idx2 = SeedNodeGraph2[i];
    // Convert to node names
    std::string node1Name = indexToName1[idx1];
    std::string node2Name = indexToName2[idx2];
    // Print the pair
    std::cout << "(" << node1Name << " => " << node2Name << ")\n";
    outFile << "(" << node1Name << " => " << node2Name << ")\n";
    }
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_total = end_total - start_total;
    std::cout << "[TIME] Total execution time: " << elapsed_total.count() << " seconds\n";

    outFile << "[TIME] Alignment process: " << elapsed_alignment.count() << " seconds\n";
    outFile << "===== End of Run ===== ";
    outFile.close();
    return 0;
}