#include "skiplist.h"
#include <chrono>
// Define the global candidateMatrix.
std::vector<bool> candidateMatrix(ROWS * COLS, false);

// ------------------ Node Member Functions ------------------

// Constructor: Allocate forward array and initialize.
Node::Node(double key, int level, int vA, int vB)
{
    this->key = key;
    this->VertexA=vA;
    this->VertexB=vB;
    // Allocate memory to forward 
    forward = new Node*[level+1];
    // Fill forward array with 0(NULL)
    std::memset(forward, 0, sizeof(Node*)*(level+1));   
}
// Destructor: Free allocated forward array.
Node::~Node(){
    delete[] forward;
}

// ------------------ SkipList Member Functions ------------------

// Constructor: Initialize skip list with header node.
SkipList::SkipList(int MAXLVL, float P)
{
    this->MAXLVL = MAXLVL;
    this->P = P;
    level = 0;
    // create header node and initialize key to -1
    header = new Node(-1, MAXLVL, -1, -1);
}
// Destructor: Delete all nodes.
SkipList::~SkipList(){
    Node * current = header->forward[0];
    while(current != nullptr){
        Node *next = current->forward[0];
        delete current;  // This calls Node's destructor to free its forward array.
        current = next;
    }
    delete header;
}

// Generate a random level for node.
int SkipList::randomLevel()
{
    float r = drand48();
    int lvl = 0;
    while(r < P && lvl < MAXLVL)
    {
        lvl++;
        r = drand48();
    }
    return lvl;
}
// create new node
Node* SkipList::createNode(double key, int level, int vA, int vB)
{
    Node *n = new Node(key, level, vA, vB);
    return n;
}

// Insert given key in skip list
void SkipList::insertElement(double key, int vA, int vB)
{
    //auto start = std::chrono::high_resolution_clock::now();

    int index = vA*COLS + vB;
    if (candidateMatrix[index]){
        return;
    }
    candidateMatrix[index]=true;
    Node *current = header;

    // create update array and initialize it
    Node *update[MAXLVL+1]={nullptr};
    memset(update, 0, sizeof(Node*)*(MAXLVL+1));

    /*    start from highest level of skip list
        move the current pointer forward while key 
        is greater than key of node next to current
        Otherwise inserted current in update and 
        move one level down and continue search
    */
    for(int i = level; i >= 0; i--)
    {
        while(current->forward[i] != NULL &&
              current->forward[i]->key <= key)
            current = current->forward[i];
        update[i] = current;
    }

    /* reached level 0 and forward pointer to 
       right, which is desired position to 
       insert key. 
    */
    current = current->forward[0];

    /* if current is NULL that means we have reached
       to end of the level or current's key is not equal
       to key to insert that means we have to insert
       node between update[0] and current node */
    if (current == NULL || current->key > key)
    {
        // Generate a random level for node
        int rlevel = randomLevel();

        /* If random level is greater than list's current
           level (node with highest level inserted in 
           list so far), initialize update value with pointer
           to header for further use */
        if(rlevel > level)
        {
            for(int i=level+1;i<rlevel+1;i++)
                update[i] = header;

            // Update the list current level
            level = rlevel;
        }

        // create new node with random level generated
        Node* n = createNode(key, rlevel, vA, vB);

        // insert node by rearranging pointers 
        for(int i=0;i<=rlevel;i++)
        {
            n->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = n;
        }
        //cout<<"Successfully Inserted key "<<key<<"\n";
    }
    else//insert same value behind the original one
    {
        // Pick a random level for the *new* duplicate node
        int rlevel = randomLevel();

        // If the new level is higher than the skip list's current level, update the list level
        if (rlevel > level) {
            for (int i = level + 1; i <= rlevel; i++) {
                update[i] = header;  // so we can link from the header if needed
            }
            level = rlevel;
        }

        // Create the new node
        Node* n = createNode(key, rlevel, vA, vB);

        // Splice the new node *behind* 'current' at each level up to rlevel
        for (int i = 0; i <= rlevel; i++) {

        // If there are multiple duplicates, 'update[i]' might not be
        // the node whose .forward[i] == current, so we walk forward if needed.
        while (update[i]->forward[i] != current) {
            update[i] = update[i]->forward[i];
        }

        // Now link in the new node *after* 'current'
        n->forward[i] = current->forward[i];
        current->forward[i] = n;
        }
    }
    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;

    //std::cout << "Time taken: " << elapsed.count() << " seconds\n";
}
// Delete an element from the skip list.
void SkipList::deleteElement(double key, int vA, int vB)
{
    //auto start = std::chrono::high_resolution_clock::now();
    Node *current = header;

    // create update array and initialize it
    Node *update[MAXLVL+1];
    memset(update, 0, sizeof(Node*)*(MAXLVL+1));

    /*    start from highest level of skip list
        move the current pointer forward while key 
        is greater than key of node next to current
        Otherwise inserted current in update and 
        move one level down and continue search
    */
    for(int i = level; i >= 0; i--)
    {
        while(current->forward[i] != NULL  &&
              current->forward[i]->key < key)
            current = current->forward[i];
            update[i] = current;
    }

    /* reached level 0 and forward pointer to 
       right, which is possibly our desired node.*/
    current = current->forward[0];
    while (current != NULL and current->key){
        if(current->VertexA == vA and current->VertexB == vB)
            break;
        current = current ->forward[0];
    }
    // If current node is target node
    if(current != NULL and current->key == key and current->VertexA == vA and current->VertexB == vB)
    {
        /* start from lowest level and rearrange
           pointers just like we do in singly linked list
           to remove target node */

        //compute the update
           for (int i = level; i >= 0; i--) {
            Node* cur = update[i];
            while (cur->forward[i] != NULL && cur->forward[i] != current) {
                cur = cur->forward[i];
            }
            update[i] = cur;
        }
    
        for(int i=0;i<=level;i++)
        {
            /* If at level i, next node is not target 
               node, break the loop, no need to move 
              further level */
            if(update[i]->forward[i] != current)
                break;

            update[i]->forward[i] = current->forward[i];
        }

        // Remove levels having no elements 
        while(level>0 &&
              header->forward[level] == 0)
            level--;
         //cout<<"Successfully deleted key "<<key<<"\n";
        delete current;
        int index = vA * COLS + vB;
        if (index >= 0 && index < (int)candidateMatrix.size())
            candidateMatrix[index] = false;
    }
    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;

    //std::cout << "Time taken for deletion: " << elapsed.count() << " seconds\n";
}

// Search for element in skip list
Node* SkipList::searchElement(double key)
{
    Node *current = header;

    /*    start from highest level of skip list
        move the current pointer forward while key 
        is greater than key of node next to current
        Otherwise inserted current in update and 
        move one level down and continue search
    */
    for(int i = level; i >= 0; i--)
    {
        while(current->forward[i] &&
               current->forward[i]->key < key)
            current = current->forward[i];

    }

    /* reached level 0 and advance pointer to 
       right, which is possibly our desired node*/
    Node *candidate = current->forward[0];

    // If current node have key equal to
    // search key, we have found our target node
    if (candidate && std::fabs(candidate->key - key) < 1e-5)
    {
        //randomly select among all candiate with the same keys during traversal
        std::cout << "Found exact key: " << key << "\n";
        return randomSelect(candidate);
    }
    else
    {
        // No exact match: current now is the greatest node with key < search key (unless it's header)
        // candidate (if not null) is the smallest node with key > search key.
        // We choose the "closest" node based on absolute difference.
        
        // If there is no valid predecessor (current is header), return candidate.
        if (current == header)
        {
            if (candidate)
            {
                std::cout << "Exact key not found. Closest key is: " 
                          << candidate->key << "\n";
                return randomSelect(candidate);
            }
            else
            {
                std::cout << "Exact key not found. No candidate available.\n";
                return candidate; // candidate is null
            }
        }
        // If there is no candidate on the right, return the predecessor.
        if (!candidate)
        {
            std::cout << "Exact key not found. Closest key is: " << current->key << "\n";
            return current;
        }
        
        // Both predecessor and candidate exist. Compare differences.
        //double diffPredecessor = key - current->key;
        //double diffCandidate = candidate->key - key;
        /*if (diffPredecessor <= diffCandidate)
        {
            std::cout << "Exact key not found. Closest key is: " << current->key << "\n";
            return current;
        }
        else
        {
            std::cout << "Exact key not found. Closest key is: " << candidate->key << "\n";
            return candidate;
        }
        */
        
        //always use candidate because it has higher sim
        std::cout << "Exact key not found. Closest key is: " << candidate->key << "\n";
            return randomSelect(candidate);
    }

}

// Display skip list level wise
void SkipList::displayList()
{
    std::cout<<"\n*****Skip List*****"<<"\n";
    for(int i=0;i<=level;i++)
    {
        Node *node = header->forward[i];
        std::cout<<"Level "<<i<<": ";
        while(node != NULL)
        {
            std::cout<<node->key<<" ";
            node = node->forward[i];
        }
        std::cout<<"\n";
    }
}

int SkipList::currentLevel(){
    return level;
}
double SkipList::topValue(){
    Node *current = header;
    for(int i = level; i >= 0; i--){
        while(current->forward[i] != NULL)
            current = current->forward[i];
    }
    return(current->key);
}

std::tuple<double, int, int> SkipList::pop(double delta)
{
double threshold = topValue() - delta;
double randomKey = threshold + ((topValue() - threshold) * drand48());
std::cout<<randomKey;
Node* chosen = searchElement(randomKey);
if (!chosen) {
    // If no node exists, return a sentinel tuple.
    return std::make_tuple(-1.0, -1, -1);
}

// Record the vertex pair of the chosen node.
double chosenKey = chosen->key;
int chosenA = chosen->VertexA;
int chosenB = chosen->VertexB;

std::vector<std::tuple<double, int, int>> nodesToDelete;
Node* current = header->forward[0];
while (current != nullptr) {
    if (current->VertexA == chosenA || current->VertexB == chosenB) {
        //std::cout << "[DEBUG] Processing node key: " << current->key << "\n";
        //std::cout << "[DEBUG] Forward[0] is at address: " << current->forward[0] << "\n";
        nodesToDelete.push_back(std::make_tuple(current->key, current->VertexA, current->VertexB));
    }
    current = current->forward[0];
}
for (const auto& nodeInfo : nodesToDelete) {
    double key = std::get<0>(nodeInfo);
    int vA = std::get<1>(nodeInfo);
    int vB = std::get<2>(nodeInfo);
    deleteElement(key, vA, vB);
}

/*// Traverse level 0 (the bottom level) of the skip list.
Node* current = header->forward[0];
while (current != nullptr) {
    // Save next pointer, since deletion will affect current->forward[0].
    Node* next = current->forward[0];
    std::cout << "[DEBUG] Processing node key: " << current->key << "\n";
    std::cout << "[DEBUG] Forward[0] is at address: " << current->forward[0] << "\n";

    // If this node involves either vertex of the chosen node...
    if (current->VertexA == chosenA || current->VertexB == chosenB) {
        // Update the existence matrix accordingly.
        int index = current->VertexA * COLS + current->VertexB;
        candidateMatrix[index] = 0;  // Mark as removed.
        // Remove this node from the skip list.
        deleteElement(current->key, current->VertexA, current->VertexB);
    }
    current = next;
}
*/




    return std::make_tuple(chosenKey, chosenA, chosenB);
}


Node* SkipList::randomSelect(Node* start)
{
    // If start is null, return null.
    if (!start)
        return nullptr;

    Node* selected = start;
    int count = 1;

    // Traverse the level-0 chain while the key remains equal to start->key.
    Node* temp = start->forward[0];
    while (temp && fabs(temp->key - start->key) < 1e-5)
    {
        count++;
        // With probability 1/count, select the current node.
        if (drand48() < 1.0/count)
            selected = temp;
        temp = temp->forward[0];
    }
    return selected;
}

