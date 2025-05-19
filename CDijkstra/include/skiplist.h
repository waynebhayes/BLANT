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
#include "rand48.h"


// Class declaration for Node.
class Node {
public:
    double key;
    int VertexA;
    int VertexB; // node pair info
    int mylevel;
    Node **forward; // pointer array for forward links

    // Constructor and destructor.
    Node(float key, int level, int vA, int vB);
    ~Node();
};

// Class declaration for SkipList.
class SkipList {
private:
    int MAXLVL;   // maximum level for this skip list
    float P;      // fraction of nodes with level i pointers that also have level i+1 pointers
    int level;    // current level of skip list
    Node *header; // pointer to header node
    std::vector<bool> candidateMatrix;
    int ROWS, COLS;
public:
    // Constructor and destructor.
    SkipList(int maxlvl, float prob, int matrix_rows, int matrix_cols);
    ~SkipList();

    // Member functions.
    int randomLevel();
    Node* createNode(float key, int level, int vA, int vB);
    void insertElement(float key, int vA, int vB);
    void deleteElement(float key, int vA, int vB);
    Node* searchElement(float key);
    Node* randomSelect(Node* start);
    void displayList();
    float topValue();
    int currentLevel();
    void deleteVector(const std::vector<std::tuple<float, int, int>>& deletions);
    std::tuple<float, int, int> pop(float delta);
};


