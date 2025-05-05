#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>  
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cstdlib>   
#include "rand48.h"
#include "graph.h"
#include "skiplist.h" 
#include "ECScore.h"

#define SKIPLIST_MAX_LEVEL 20 //define the max level of skip list
#define SKIPLIST_PROB 0.5 //define the probability of promote to upper level for a node
#define MAX_DISCARD_TIMES 10 //define the how many times we can discard bad pairs if they have bad EC scores until we stop

int main(int argc, char* argv[])
{   
    //command line input
    float delta = 0.1;
    double ECthresh = 0.5;
    std::string Sim = "";
    std::string Seed = "";
    std::string Graph1 = "";
    std::string Graph2 = "";

    if(argc!= 7){
        std::cerr<<"try again: delta ECthresh simFile seedFile g1.el g2.el\n";
        return 1;
    }
    else{
        delta = std::stod(argv[1]);
        ECthresh = std::stod(argv[2]);
        Sim = argv[3];
        Seed = argv[4];
        Graph1 = argv[5];
        Graph2 = argv[6];
    }
    srand48(time(NULL));
    auto start_total = std::chrono::high_resolution_clock::now();
    // Input file names
    //./build/cdijkstra 0.1 0.5 ../../SANA/sequence/graphlet+seq.1/RNorvegicus-SPombe.sim test/RNorvegicus-SPombe-Seed.txt ../../SANA/networks/RNorvegicus.el ../../SANA/networks/SPombe.el
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
    //Newly added: initially calculate EC scores
    int globalEA = computeEA(SeedNodeGraph1, SeedNodeGraph2, adjMatrix1, adjMatrix2);
    int globalE1 = computeE(SeedNodeGraph1, adjMatrix1);
    int globalE2 = computeE(SeedNodeGraph2, adjMatrix2);
    double GlobalEC1 = computeEC(globalEA, globalE1);
    double GlobalEC2 = computeEC(globalEA, globalE2);
    double GlobalS3 = computeS3(globalEA, globalE1, globalE2);
    
    std::cout << "EA: " << globalEA << "\n";
    std::cout << "E1: " << globalE1 << "\n";
    std::cout << "E2: " << globalE2 << "\n";
    std::cout <<"The initial globalEC1: "<<GlobalEC1<< "\n";
    std::cout <<"The initial globalEC2: "<<GlobalEC2<< "\n";
    std::cout <<"The initial globalS3: "<<GlobalS3<< "\n";
    std::cout << "Seed size: " << SeedNodeGraph1.size() << "\n";

    // 5. Initialize alignment process
    auto start_alignment = std::chrono::high_resolution_clock::now();
    bool alignmentInProgress = true;
    SkipList skiplist(SKIPLIST_MAX_LEVEL, SKIPLIST_PROB);
    int iterationCount = 0;
    //const int maxIterations = 10;
    int discardCount = 0;
    int alignCount = 0;
    std::vector<int> connectedNodes1;
    std::vector<int> connectedNodes2;
    std::vector<int> discardedNodes1;
    std::vector<int> discardedNodes2;
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
        //auto insert_start = std::chrono::high_resolution_clock::now();
        for (int n1 : connectedNodes1)
        {
        // Check if n1 is already in the seed set
        if (std::find(SeedNodeGraph1.begin(), SeedNodeGraph1.end(), n1) != SeedNodeGraph1.end())
            continue;
        if (std::find(discardedNodes1.begin(), discardedNodes1.end(), n1) != discardedNodes1.end())
            continue;

            for (int n2 : connectedNodes2)
            {
            // Check if n2 is already in the seed set
            if (std::find(SeedNodeGraph2.begin(), SeedNodeGraph2.end(), n2) != SeedNodeGraph2.end()) 
                continue;
            if (std::find(discardedNodes2.begin(), discardedNodes2.end(), n2) != discardedNodes2.end()) 
                continue;
            float sim = similarityMatrix[n1][n2];
            //std::cout << "Candidate pair (" << n1 << ", " << n2 << ") has similarity " << sim << "\n";
            if (sim >0){
                //std::cout << "Inserting: (" << n1 << ", " << n2 << ") with similarity " << sim << "\n";
                skiplist.insertElement(sim, n1, n2);
                //for debug randomness
                //skiplist.displayList();

                }
            }
        }
        //auto insert_end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed_insertion = insert_end - insert_start;
        //std::cout << "Time taken for insertion: " << elapsed_insertion.count() << " seconds\n";
        // 8. Display skip list
        //skiplist.displayList();

        // 9. Pop one candidate from skip list
        auto tup = skiplist.pop(delta);  // returns (key, first, second) // call this constant DELTA above
        float key;
        int first, second;
        std::tie(key, first, second) = tup;
        if (key != -1.0) {
            
            int localEA = incrementalEA(first, second, SeedNodeGraph1, SeedNodeGraph2, adjMatrix1, adjMatrix2);
            int localE1 = incrementalE(first, SeedNodeGraph1, adjMatrix1);
            int localE2 = incrementalE(second, SeedNodeGraph2, adjMatrix2);
            //output local integer
            std::cout<<"The number of node1 add back to the graph1 is"<<localE1<<"\n";
            std::cout<<"The number of node2 add back to the graph2 is"<<localE2<<"\n";
            std::cout<<"The number of aligned edges is"<<localEA<<"\n";
            std::cout<<"localEC1:"<<double(localEA)/localE1<<"; localEC2:"<<double(localEA)/localE2<<"\n";
            int current_globalEA = localEA+globalEA;
            int current_globalE1 = localE1+globalE1;
            int current_globalE2 = localE2+globalE2;
            double current_GlobalEC1 = computeEC(current_globalEA,current_globalE1);
            double current_GlobalEC2 = computeEC(current_globalEA, current_globalE2);
            double current_GlobalS3 = computeS3(current_globalEA, current_globalE1, current_globalE2);
            std::cout<<"current_EC1: "<<current_GlobalEC1<< "\n";
            std::cout<<"current_EC2: "<<current_GlobalEC2<< "\n";
            std::cout<<"current_S3: "<<current_GlobalS3<< "\n";

        //std::cout << "Popped node => key: " << key
              //<< ", VertexA: " << first
              //<< ", VertexB: " << second << "\n";
            if(current_GlobalEC1 >= ECthresh && current_GlobalEC2 >= ECthresh){
                //10. Add to aligned list
            SeedNodeGraph1.push_back(first);
            SeedNodeGraph2.push_back(second);
            //std::cout << "[DEBUG] Using maxNode: (" << first << ", " << second
            //      << ") with value=" << key << "\n";
            connectedNodes1 = getConnectedNodes(first, adjMatrix1);
            connectedNodes2 = getConnectedNodes(second, adjMatrix2);
            globalEA = current_globalEA;
            globalE1 = current_globalE1;
            globalE2 = current_globalE2;
            GlobalEC1 = current_GlobalEC1;
            GlobalEC2 = current_GlobalEC2;
            GlobalS3 = current_GlobalS3;
            alignCount++;
            discardCount = 0;
            }
            else{
                std::cout << "This pair does not have good EC.\n";
                discardedNodes1.push_back(first);
                discardedNodes2.push_back(second);
                discardCount++;
                alignCount = 0;
                if(iterationCount > 500 && discardCount >=MAX_DISCARD_TIMES && alignCount == 0){
                    alignmentInProgress = false;
                }
            }
            
    } else {
    // No more candidates, stop the loop
    if(key == -1.0){
         std::cout << "[DEBUG] Skiplist is empty, no candidate pairs.\n";
    alignmentInProgress = false;
    }
    else{
        std::cout <<"This pair does not have good ec.\n";
    }
    }
    iterationCount++;
    std::cout<<"EC1: "<<GlobalEC1<< "\n";
    std::cout<<"EC2: "<<GlobalEC2<< "\n";
    std::cout<<"S3: "<<GlobalS3<< "\n";
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
    std::cout << "(" << node1Name << "	 " << node2Name << ")\n";
    outFile << "(" << node1Name << "	 " << node2Name << ")\n";
    }
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_total = end_total - start_total;
    std::cout << "[TIME] Total execution time: " << elapsed_total.count() << " seconds\n";
    std::cout<<"total number of edges used in the alignment EA E1 E2:"<<globalEA<<" "<<globalE1<<" "<<globalE2;
    outFile << "[TIME] Alignment process: " << elapsed_alignment.count() << " seconds\n";
    outFile << "===== End of Run ===== ";
    outFile.close();
    return 0;
}
