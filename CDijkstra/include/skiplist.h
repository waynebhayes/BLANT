#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <string.h>
#include <cstring>
#include <cstdint>  // For int8_t
#include <cstdlib>
#include <algorithm>
#include <cmath>

// Define constants.
const int ROWS = 41456; // human protein
const int COLS = 16127; // yeast protein

// Declare candidateMatrix externally.
extern std::vector<bool> candidateMatrix;

// Class declaration for Node.
class Node {
public:
    double key;
    int VertexA;
    int VertexB; // node pair info
    Node **forward; // pointer array for forward links

    // Constructor and destructor.
    Node(double key, int level, int vA, int vB);
    ~Node();
};

// Class declaration for SkipList.
class SkipList {
private:
    int MAXLVL;   // maximum level for this skip list
    float P;      // fraction of nodes with level i pointers that also have level i+1 pointers
    int level;    // current level of skip list
    Node *header; // pointer to header node

public:
    // Constructor and destructor.
    SkipList(int MAXLVL, float P);
    ~SkipList();

    // Member functions.
    int randomLevel();
    Node* createNode(double key, int level, int vA, int vB);
    void insertElement(double key, int vA, int vB);
    void deleteElement(double key, int vA, int vB);
    Node* searchElement(double key);
    Node* randomSelect(Node* start);
    void displayList();
    double topValue();
    int currentLevel();
    std::tuple<double, int, int> pop(double delta);
};


