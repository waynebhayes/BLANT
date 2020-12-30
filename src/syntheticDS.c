#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include "syntheticDS.h"
#include "graph.h"
#include "sets.h"

// DICTIONARY
// wrapper around 'uthash'
// credits - https://raw.githubusercontent.com/troydhanson/uthash/master/src/uthash.h

int dictionary_create(Dictionary* dictionary){
    dictionary->hashTable = NULL;
}

int dictionary_get(Dictionary* dictionary, int key, int default_value){
    KeyValue* s;
    HASH_FIND_INT(dictionary->hashTable, &key, s);
    if (s == NULL){
        return default_value;
    }else{
        return s->value;
    }
}

void dictionary_set(Dictionary* dictionary, int key, int value){
    KeyValue* s;
    HASH_FIND_INT(dictionary->hashTable, &key, s);
    if (s == NULL){
        // no such key exists
        s = (KeyValue*) malloc(sizeof(KeyValue));
        s->key = key;
        s->value = value;
        HASH_ADD_INT(dictionary->hashTable, key, s);
    }else{
        // overwrite
        s->value = value;
    }
}

// iterator
KeyValue* getIterator(Dictionary* dictionary){
    KeyValue* iterator = dictionary->hashTable;
    return iterator;
}

int getNext(KeyValue** iterator, int* key, int* value){
    if((iterator==NULL) || (*iterator == NULL))
        return -1;
    memcpy(key, &((*iterator)->key), sizeof(int));
    memcpy(value, &((*iterator)->value), sizeof(int));
    *iterator = (*iterator)->hh.next;
    return 0;
}


// STACK
int create_stack(RevertStack* stack, int size){
    stack->tos = -1;
    stack->size = size;
    stack->space = (Change*) malloc(size * sizeof(Change));
    if (stack->space != NULL)
        return 0;
    else
        return -1;
}

int init_stack(RevertStack* stack){
    stack->tos = -1;
    if (stack->space != NULL)
        return 0;
    else
        return -1;
}

int push(RevertStack* stack, Change elt){
    if (stack->tos == (stack->size-1))
        return -1;
    (stack->tos) += 1;
    memcpy(&(stack->space)[stack->tos], &elt, sizeof(Change));
    return 0;
}

int pop(RevertStack* stack, Change* elt){
    if (stack->tos == -1)
        return -1;
    memcpy(elt, &(stack->space)[stack->tos], sizeof(Change));
    stack->tos -= 1;
    return 0;
}


// sort
// comparison function definition (to be used by qsort)
int compare_ints(const void* a, const void* b){
    // this function is from -> https://en.cppreference.com/w/c/algorithm/qsort
    int arg1 = *(const int*)a;
    int arg2 = *(const int*)b;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int compare_doubles(const void* a, const void* b){
    // this function is from -> https://en.cppreference.com/w/c/algorithm/qsort
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int compare_int_tuples(const void* a, const void* b){
    // sort on the 0th key
    int arg1 = (*(const IntTuple*)a).nums[0];
    int arg2 = (*(const IntTuple*)b).nums[0];

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

// median
int getIntMedian(int* nums, int start, int end){
    // array must have MAX(start+1,end+1) elements
    int n = end-start+1;
    if (n%2 == 0){
        return (nums[start + (n/2)] + nums[start + (n/2) - 1]) / 2.0;
    }else{
        return nums[start + (n/2)];
    }
}

double getDoubleMedian(double* nums, int start, int end){
    // array must have MAX(start+1,end+1) elements
    int n = end-start+1;
    if (n%2 == 0){
        return ((double) (nums[start + (n/2)] + nums[start + (n/2) - 1])) / 2.0;
    }else{
        return nums[start + (n/2)];
    }
}

// poisson distribution
double PoissonDistribution(double l, int k){
    //->   (e^(-l) x l^k) / k!
    double r = exp(-l);
    int i;
    for(i=k; i>0; i--) // divide by k!
    r *= l/i;
    return r;
}


double getDoubleBinSize(int n, double localClustCoff[n], double* scratchspace){
    // sort
    // getIQR
    // return (2IQR) / (n^1/3)
    //assert(n >= 5);

    double* sorted = scratchspace;
    memcpy(sorted, localClustCoff, n * sizeof(double));
    qsort(sorted, n, sizeof(double), compare_doubles);

    double q1, q3;

    if (n%2 == 0){
        q1 = getDoubleMedian(sorted, 0, (n/2)-1);
        q3 = getDoubleMedian(sorted, n/2, n-1);
    }else{
        q1 = getDoubleMedian(sorted, 0, n/2);
        q3 = getDoubleMedian(sorted, n/2, n-1);
    }

    double obs = ((double) (2 * (q3-q1)) / cbrt((double) n));
    if (obs < 0)
        return obs;
    if (obs < 0.0001)
        return 0.05;

    return obs;
}

int getIntegerBinSize(int n, int GDVcolumn[n], int* scratchspace){
    // sort
    // getIQR
    // return  (2IQR)/(n^1/3)
    //assert(n>=5);

    int* sorted = scratchspace;
    memcpy(sorted, GDVcolumn, n * sizeof(int));
    qsort(sorted, n, sizeof(int), compare_ints);

    int q1,q3;

    if (n%2 == 0){
        q1 = getIntMedian(sorted, 0, (n/2)-1);
        q3 = getIntMedian(sorted, n/2, n-1);
    }else{
        q1 = getIntMedian(sorted, 0, n/2);
        q3 = getIntMedian(sorted, n/2, n-1);
    }

    //assert(q3 >= q1);
    double obs = ((double) (2 * (q3-q1)) / cbrt((double) n));
    int returnVal = (int) ceil(obs);

    if (returnVal < 0)
        return returnVal;

    if (returnVal < 1)
        return 1;

    return returnVal;
}

// returns a random node which is 'd' hops away from node 'src' - by doing a bfs traversal (till d levels)
// incase there are no nodes at 'd' hops from 'src', then return the farthest random node
int getRandomNodeAtHops(GRAPH* G, int src, int d){

    int i,j,k,start,end,maxd;
    SET* visited = SetAlloc(G->n);
    SetAdd(visited, src);
    i = -1;
    j = 0;

    int* queue = (int*) malloc(sizeof(int) * G->n);
    int* distance = (int*) malloc(sizeof(int) * G->n);

    queue[j] = src;
    distance[j] = 0;

    maxd = 0;
    start = end = 0;

    while(i<j){
        i += 1;
        for(k=0; k < G->degree[queue[i]]; k++){
            if(!SetIn(visited, (G->neighbor[queue[i]])[k])){
                j += 1;
                queue[j] = (G->neighbor[queue[i]])[k];
                distance[j] = distance[i] + 1;
                SetAdd(visited, (G->neighbor[queue[i]])[k]);

                if(distance[j] > d)
                    break;

                if(distance[j] > maxd){
                    maxd = distance[j];
                    start = end = j;
                }

                if(distance[j] == maxd)
                    end = j;

            }
        }
    }

    // return random vertex b/w start&end
    int ans = start + (int)((end-start+1) * drand48());
    ans = queue[ans];
    free(queue);
    free(distance);
    return ans;
}

// returns a random neighbor
int getRandomConnectedNode(GRAPH* G, int src){
    int degree = G->degree[src];
    if (degree == 0)
        return -1;
    int random = drand48() * degree;
    return G->neighbor[src][random];
}

// khop distribution from a sample; also computes nodes sorted by num. of shortest paths
void sampleKHop(GRAPH* G, Dictionary* khop, double quality, int nodesBySp[G->n]){
    /*
    The exact khop distribution for any graph can be computed by (this is one possible algorithm) -

    do a BFS starting from every node
        -during this BFS, keep track of how many nodes are a distance 1/2/3/... away from the current source
        -keep adding the number of nodes found at distance 1/2/3/... ; into a global array or dictionary

    But, we only use the median of the khop distribution (in the main calling code); so it's okay to take
    a sample of the nodes (every node is selected with probability p=quality); and the above algorithm is executed.

    Now, this function also returns an array of nodes, sorted by the number of shortest paths going through them,
    which can be computed in the same bfs traversal.

    Once you do a bfs traversal and reach the leaf nodes (you now have a tree),
    you can backtrack to the source; while aggregating the shortest paths, using the following recursive algorithm

    do a BFS starting from every node (SOURCE)
        - Now, for any node, the number of shortest paths going through it is: 1 + the sum of SPs going through it's child nodes (and just 1 for leaves)
          (every node that is discovered is actually a terminating point of a SP from SOURCE)

    Note: If there's an SP b/w nodes u & v, then the count for u and v (and all nodes b/w them in the path) is incremented TWICE (once during bfs from `u` and once from `v`)
    Again, the exact count of SPs doesn't matter, as we're interested in the relative node counts.

    Complexity: k * (2 * (n+m))   ,where k is the number of nodes a bfs is done for, n = nodes, m = edges
    */

    int n, i, j, k, d, p, q, this, neigh;
    double random;
    if (quality>1.0) quality = 0.999;
    if (quality<0.0) quality = -1;

    for(i=0; i<G->n; i++)
        nodesBySp[i] = 0;  // *TEMPORARILY*, before line 380, nodesBySp[i] stores num of shortest paths going through `i`th vertex
                            // later, it will be sorted where nodes with lesser num of shortest paths going through have lower indexes

    dictionary_create(khop);
    int* queue = (int*) malloc(sizeof(int) * G->n);
    int* parent = (int*) malloc(sizeof(int) * G->n);
    int* distance = (int*) malloc(sizeof(int) * G->n);

    for(n=0; n < (G->n); n++){ // outer loop (select source vertex)
        random = drand48();
        if (random > quality)
            continue;

        // do bfs on the src vertex 'n'
        SET* visited = SetAlloc(G->n);
        int nodesBySp_local[G->n];
        for(i=0; i<G->n; i++) nodesBySp_local[i] = 0;

        SetAdd(visited, n);
        nodesBySp_local[n] = 1;
        i=-1;
        j=0;
        queue[j] = n;
        parent[j] = -1;
        distance[j] = 0;

        while(i<j){
            i += 1;
            this = queue[i];
            for(k=0; k < G->degree[this]; k++){
                neigh = (G->neighbor[queue[i]])[k];
                if(!SetIn(visited, neigh)){
                    j += 1;
                    queue[j] = neigh;
                    parent[j] = i;
                    distance[j] = distance[i] + 1;
                    SetAdd(visited, neigh);
                    assert(nodesBySp_local[neigh] == 0);
                    nodesBySp_local[neigh] = 1;
                }
            }
        }

        SetFree(visited);

        // aggregate distances into khop dictioanry
        p=0;
        while(p<=j){
            d = distance[p];
            q = p;
            while((q<=j) && (distance[q] == d))
                q += 1;
            dictionary_set(khop, d, dictionary_get(khop, d, 0) + (q-p+1));
            p = q + 1;
        }

        // backtrack from outermost nodes towards the source, using parent pointers (to accumulate SPs going through each node)
        assert(i==j);
        while(j>=0){
            assert(nodesBySp_local[queue[j]] >= 0);
            nodesBySp[queue[j]] = nodesBySp_local[queue[j]];
            if(parent[j]>=0){
                nodesBySp_local[queue[parent[j]]] += nodesBySp_local[queue[j]];
                assert(nodesBySp_local[queue[parent[j]]] >= 0);
            }
            j -= 1;
        }
    }

    // sort nodesBySp
    IntTuple* scratch = (IntTuple*) malloc(sizeof(IntTuple) * G->n);
    for(i=0; i<G->n; i++){
        scratch[i].nums[0] = nodesBySp[i];
        scratch[i].nums[1] = i;
    }
    qsort(scratch, G->n, sizeof(IntTuple), compare_int_tuples);
    for(i=0; i<G->n; i++){
        nodesBySp[i] = (scratch[i]).nums[1];  // nodes (sorted by SPs going through them)
        assert(nodesBySp[i]<G->n);
    }

    free(scratch);
    free(queue);
    free(distance);
}

// NORMALIZED by sum-of-values (let's say x nodes at 1 hop, y nodes at 2 hops, ...; then normalization constant = x+y+...)
void print_khop_sample(Dictionary* khop){
    // line 403 does the normalization
    KeyValue* iter = getIterator(khop);
    int k,v;
    int maxk = -1;
    int valsum = 0;

    while(getNext(&iter, &k, &v) == 0){
        if(k > maxk)
            maxk = k;
        valsum += v;
    }

    if (valsum == 0){
        fprintf(stderr, "0 0 0 0 0...0\n");
        return;
    }

    int i;
    for(i=0; i<=k; i++){
        fprintf(stderr, "%g ", (double) dictionary_get(khop, i, 0)/ (double) valsum);
    }
    fprintf(stderr, "\n");
}

// returns -1 if medians are equal, 0 if khop[0] has bigger median, 1 if khop[1] has bigger median
int compareKHopByMedian(Dictionary* khop[2], int medians[2], int MAX_Keys[2]){
    KeyValue* iter;
    int count, k, v, i, l, mi;

    MAX_Keys[0] = MAX_Keys[1] = -1;

    for(i=0; i<2; i++){
        count = 0;

        iter = getIterator(&(khop[i]));
        while((getNext(&iter, &k, &v)) == 0){
            MAX_Keys[i] = MAX(MAX_Keys[i], k);
            count += v;
        }

        mi = (int) count/2;
        count = 0;
        for(l=0; l<=MAX_Keys[i]; l++){
            count += dictionary_get(&(khop[i]), l, 0);
            if (count >= mi){
                medians[i] = l;
                break;
            }
        }
    }

    if(medians[0] < medians[1]){
        return 1;
    }else if(medians[1] < medians[0]){
        return 0;
    }else{
        return -1;
    }
}



