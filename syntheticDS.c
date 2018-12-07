#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"
#include "syntheticDS.h"
#include "graph.h"
#include "sets.h"

// DICTIONARY
// wrapper around 'uthash'
// https://raw.githubusercontent.com/troydhanson/uthash/master/src/uthash.h

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

int getRandomNode(GRAPH* G, int src, int d){

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

void sampleKHop(GRAPH* G, Dictionary* khop, double quality){

    int n, i, j, k, d, p, q;
    double random;
    if (quality>1.0) quality = 0.999;
    if (quality<0.0) quality = -1;

    dictionary_create(khop);
    int* queue = (int*) malloc(sizeof(int) * G->n);
    int* distance = (int*) malloc(sizeof(int) * G->n);

    for(n=0; n < (G->n); n++){ // outer loop (select source vertex)
        random = drand48();
        if (random > quality)
            continue;

        // do bfs on 'n'
        SET* visited = SetAlloc(G->n);
        SetAdd(visited, n);
        i=-1;
        j=0;
        queue[j] = n;
        distance[j] = 0;

        while(i<j){
            i += 1;
            for(k=0; k < G->degree[queue[i]]; k++){
                if(!SetIn(visited, (G->neighbor[queue[i]])[k])){
                    j += 1;
                    queue[j] = (G->neighbor[queue[i]])[k];
                    distance[j] = distance[i] + 1;
                    SetAdd(visited, (G->neighbor[queue[i]])[k]);
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
    }

    free(queue);
    free(distance);
}

void print_khop(Dictionary* khop){
    /*NORMALIZED*/
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


// old implementation of Clustering Coeffienct Objective
/*
double AdjustClustCoff(const int x, const int y, const int connected, GRAPH* G, int localConnections[2][_maxNodes], double avg_cc[2]){

    // Should be called AFTER calling GraphConnect or GraphDisconnect
    // right now localConnections and avg_cc are in sync.

    assert(abs(connected) == 1);
    int i, j, node, nodecount, c, nc2;

    double sumchange = 0.0;
    double oldcc, newcc;

    SET* xn = SetAlloc(G->n);
    for(i=0; i < G->degree[x]; i++)
        if (G->neighbor[x][i] != y)
            SetAdd(xn, G->neighbor[x][i]);
    
    nodecount = 0;  // nodes which are connected both to x & y
    for(i=0; i < G->degree[y]; i++){
        node = G->neighbor[y][i];
        if ((node!=x) && (SetIn(xn, node))){
            nodecount += 1;
            assert(G->degree[node] >= 2);
            nc2 = (G->degree[node] * (G->degree[node] - 1))/2;
            c = localConnections[1][node];  // the original num of connections
            sumchange += (connected * (1/nc2));
            localConnections[1][node] = c + connected;  // the new num of connections
        }
    }

    SetFree(xn);

    // for x
    if ((G->degree[x] - connected) < 2){ // original degree
        c = localConnections[1][x];  // the old connections
        assert(c == 0);
        oldcc = 0;
    }else{
        nc2 = ((G->degree[x] - connected) * ((G->degree[x] - connected) - 1))/2;
        c = localConnections[1][x];  // the old connections
        oldcc = ((double) c) / ((double) nc2);
    }

    if(G->degree[x] < 2){  // current degree
        c += (connected * nodecount);
        assert(c == 0);
        newcc = 0;
    }else{
        nc2 = (G->degree[x] * (G->degree[x] - 1))/2;
        c += (connected * nodecount);
        newcc = ((double) c) / ((double) nc2);
    }
    localConnections[1][x] = c;
    sumchange += (newcc-oldcc);


    // for y
    if ((G->degree[y] - connected) < 2){ // original degree
        c = localConnections[1][y];  // the old connections
        assert(c == 0);
        oldcc = 0;
    }else{
        nc2 = ((G->degree[y] - connected) * ((G->degree[y] - connected) - 1))/2;
        c = localConnections[1][y];  // the old connections
        oldcc = ((double) c) / ((double) nc2);
    }

    if(G->degree[y] < 2){  // current degree
        c += (connected * nodecount);
        assert(c == 0);
        newcc = 0;
    }else{
        nc2 = (G->degree[y] * (G->degree[y] - 1))/2;
        c += (connected * nodecount);
        newcc = ((double) c) / ((double) nc2);
    }
    localConnections[1][y] = c;
    sumchange += (newcc-oldcc);

    
    // sanity check
    int ideal_con[G->n];
    getConnections(G, ideal_con); // slow!
    for(i=0; i<G->n; i++)
        assert(ideal_con[i] == localConnections[1][i]);

    avg_cc[1] += (sumchange / G->n);
    assert((avg_cc[1] >= 0) && (avg_cc[1] <= 1));
    return fabs(avg_cc[1] - avg_cc[0]);
}
*/

/*
double ClustCoffObjective(GRAPH* G[2], int localConnections[2][_maxNodes], double avg_cc[2]){

    // old implementation : compute fabs() b/w Average clust coffs of target & synthetic

    double t;
    double cc_sum[2] = {0.0, 0.0};
    int i, j, degree, nc2;

    for(i=0; i<2; i++){
        for(j=0; j < G[i]->n; j++){
            degree = G[i]->degree[j];
            nc2 = (degree * (degree-1)) / 2;
            if (degree > 1)
                t = ((double) localConnections[i][j]) / ((double) nc2);  // local clustering coefficient
            else
                t = 0.0;
            
            assert((t >= 0) && (t <= 1));
            cc_sum[i] += t;
        }
    }

    avg_cc[0] = cc_sum[0] / (G[0]->n);
    avg_cc[1] = cc_sum[1] / (G[1]->n);

    fprintf(stderr, "clustering coefficients: target=%g synthetic=%g\n", avg_cc[0], avg_cc[1]);

    double returnVal = (double) fabs(avg_cc[0] - avg_cc[1]);
    return returnVal;
}*/
