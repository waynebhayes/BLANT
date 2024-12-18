#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <set>
#include <fstream>

using namespace std;



class Node {
public:
    int height;//level
    double value;//sim
    pair<int, int> nodePair;;
    vector<Node*> next;
    Node(int nodeA, int nodeB, double value, int level)
        : nodePair(nodeA, nodeB), value(value), next(level + 1, nullptr) {}
};

class skipList {
private:
    Node* head;
    int Level;
    unordered_map<int, set<int>> adjacencyList;  // Adjacency list for storing neighbors

public:
    skipList();
    int getCurrentLevel();
    void insert(int nodeA, int nodeB, double value);  // Insert value
    bool search(int nodeA, int nodeB);  // Search for a value
    void display();         // Display the skip list
    void addEdge(int node1, int node2);  // Add edge to adjacency list
    vector<int> getNeighbors(int node);  // Get neighbors of a node
    void extendAlignment(int node1, int node2);  // Extend alignment using neighbors
};

// Constructor for skip list
skipList::skipList() {
    head = new Node(0,0,0,0);  // Initialize the skip list with the 0 levels
    Level = 0;                             // Starting at level 0
}

int skipList::getCurrentLevel() {
    return Level;
}

void skipList::insert(int nodeA, int nodeB, double value) {
    // Unpack the tuple
    int newLevel = 0;

    // Determine the level of the new node
    while (newLevel < getCurrentLevel() && (rand() % 2) == 1) {
        newLevel++;
    }

    // Update the list structure if the new level exceeds the current level
    if (Level < newLevel) {
        head->next.resize(newLevel + 1, nullptr);
        Level = newLevel;
    }

    Node* current = head;
    vector<Node*> Update(Level + 1, nullptr);

    // Traverse the list to find the correct position for insertion
    for (int i = Level; i >= 0; i--) {
        while (current->next[i] && current->next[i]->value <value) { // Compare using 'value'
            current = current->next[i];
        }
        Update[i] = current;
    }

    // Check if the value already exists
    current = current->next[0];
    if (current == nullptr || current->nodePair != make_pair(nodeA, nodeB)) { // Compare using 'value'
        Node* newNode = new Node(nodeA, nodeB, value, Level);
        // Insert the new node at all levels up to the new level
        for (int i = 0; i <= newLevel; i++) {
            newNode->next[i] = Update[i]->next[i];
            Update[i]->next[i] = newNode;
        }
        cout << "Node pair (" << nodeA << ", " << nodeB << ") inserted successfully.\n";
    }
    else {
        cout << "Node pair (" << nodeA << ", " << nodeB << ") already exists.\n";
    }
}




vector<int> skipList::getNeighbors(int node) {
    vector<int> neighbors;
    if (adjacencyList.find(node) != adjacencyList.end()) {
        for (int neighbor : adjacencyList[node]) {
            neighbors.push_back(neighbor);
        }
    }
    return neighbors;
}


bool skipList::search(int nodeA, int nodeB) {
    Node* current = head;

    // Traverse through levels, comparing pairs
    for (int i = Level; i >= 0; i--) {
        while (current->next[i] &&
            (std::get<0>(current->next[i]->nodePair) < nodeA ||
                (std::get<0>(current->next[i]->nodePair) == nodeA && std::get<1>(current->next[i]->nodePair) < nodeB))) {
            current = current->next[i];
        }
    }

    // Move to the lowest level and compare the pair
    current = current->next[0];

    // Check if the node pair matches the search pair
    if (current != nullptr &&
        std::get<0>(current->nodePair) == nodeA && std::get<1>(current->nodePair) == nodeB) {
        cout << "Element (" << nodeA << ", " << nodeB << ") found.\n";
        return true;
    }
    else {
        cout << "Element (" << nodeA << ", " << nodeB << ") not found.\n";
        return false;
    }
}


void skipList::display() {
    cout << "Skip List (Graph-Graph Relationship):" << endl;
    for (int i = Level; i >= 0; i--) {
        Node* current = head->next[i];
        cout << "Level " << i << ": ";
        while (current != nullptr) {
            // Display the pair (nodeA, nodeB) at each level
            cout << "(" << std::get<0>(current->nodePair) << ", " << std::get<1>(current->nodePair) << ") ";
            current = current->next[i];
        }
        cout << endl;
    }
}




