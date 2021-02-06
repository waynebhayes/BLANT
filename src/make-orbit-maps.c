#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "blant.h"

static int k;
static Gint_type canon_list[MAX_CANONICALS];

Boolean suitablePerm(int permutation[], int adj[k][k+1]);
void permuteNodes(int permutation[], int adj[k][k+1], long *D);
void makeOrbit(int permutation[], long orbit[]);
void getCycle(int permutation[], int cycle[], int seed, int current, Boolean visited[]);

Boolean nextPermutation(int permutation[])
{
    int t, i;
    for(i=k-1;i>0;i--)
    {
	if(permutation[i]>permutation[i-1])
	{
	    int j;
	    for(j=k-1;j>i-1;j--){
		if(permutation[i-1]<permutation[j]){
		    int t=permutation[i-1];
		    permutation[i-1]=permutation[j];
		    permutation[j]=t;
		    break;
		}
	    }
	    int l=i;
	    for(j=k-1;j>l;j--)
	    {
		if(i<j){
		    t=permutation[i];
		    permutation[i]=permutation[j];
		    permutation[j]=t;
		    i++;
		}
	    }
	    return 1;
	}
    }
    return 0;
}

void getDecimal(int adj[k][k+1],long *D){

    int i, j, bitPos=0;
    *D=0;
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
       for(j=i-1;j>=0;j--)
#else	// UPPER_TRIANGLE
    for(i=k-2;i>=0;i--)
	for(j=k-1;j>i;j--)
#endif
	{
	    if(adj[i][j]==1) *D+= (1 << bitPos);
		bitPos++;
	}

}

void getGraph(long int decimal, int adj[k][k+1]){
    int i,j;
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
       for(j=i-1;j>=0;j--)
#else	// UPPER_TRIANGLE
    for(i=k-2;i>=0;i--)
	for(j=k-1;j>i;j--)
#endif
	{
	    adj[i][j]=decimal%2;
            adj[j][i]=adj[i][j];
            decimal/=2;
   	}

     for(i=0; i<k; i++)
        for(j=0; j<k; j++)
            adj[i][k]+=adj[i][j];
}

void orbits(long int Decimal,long orbit[]){
    int permutation[k], i, j;

    //initially,the permutation is (0, 1, ..., numNodes_-1)
    //and each node is in its own orbit. we'll merge orbits later
    for(i = 0; i < k; i++){
        permutation[i]=i;
        orbit[i]=i;
    }
    int adj[k][k+1];
    for(i=0; i<k; i++)
     for(j=0; j<k+1; j++)
      adj[i][j]=0;
    getGraph(Decimal,adj);
    long d=0;

    while( nextPermutation(permutation) ){
        d=0;
        //ruling out searching into unnecessary permutations
        if(!suitablePerm(permutation,adj)) continue;
        permuteNodes(permutation,adj,&d);
        //Check for automorphism
        if(Decimal == d){

              makeOrbit(permutation, orbit);
        }
    }
}

Boolean suitablePerm(int permutation[], int adj[k][k+1]){
    int node;
    for(node = 0; node < k; node++){
	if(adj[node][k] != adj[permutation[node]][k]){
	    return false;
	}
    }
    return true;
}

void permuteNodes(int permutation[], int adj[k][k+1], long* D){
    //To store adjMatrix_ after permutation
    int pAdj[k][k+1], i, j;
    for(i=0; i<k; i++)
     for(j=0; j<k+1; j++)
        pAdj[i][j]=0;
    //Apply the permutation
    for(i = 0; i < k; i++){
        for(j =i+1; j < k; j++){
            int pi = permutation[i], pj = permutation[j];
	    pAdj[pi][pj] = adj[i][j];
	    pAdj[pj][pi] = adj[i][j];
        }
    }

    getDecimal(pAdj, D);

}
void makeOrbit(int permutation[], long orbit[]){
    int i, j;
    Boolean visited[k];
    for(i=0;i<k;i++)
        visited[i]=0;
    for(i = 0; i < k; i++){
        if(!visited[i]){
            //finding out each cycle at a time
            int cycle[k+1];
            cycle[0]=0;
	    int minOrbit=i;
            getCycle(permutation, cycle, i, i, visited);

            for(j=1; j<=cycle[0]; j++)
                minOrbit = MIN(orbit[cycle[j]], minOrbit);
            for(j=1; j<=cycle[0]; j++)
                orbit[cycle[j]] = minOrbit;
        }
    }
}

void getCycle(int permutation[], int cycle[], int seed, int current, Boolean visited[]){
    cycle[++cycle[0]]=current;
    visited[current] = true;
    if(permutation[current] != seed)
        getCycle(permutation, cycle, seed, permutation[current], visited);
}

int main(int argc, char* argv[]){
    if(argc != 2) {
	fprintf(stderr, "USAGE: %s k\n", argv[0]);
	exit(1);
    }
    k=atoi(argv[1]);
    assert(k > 2 && k <= 8);
    //reading data from canon_list file
    char BUF[BUFSIZ];
    SET *connectedCanonicals = canonListPopulate(BUF, canon_list, k);
    int numCanon = connectedCanonicals->n;
    SetFree(connectedCanonicals);
    int numOrbits = 0, i, j;
    long orbit[numCanon][k];
    for(i=0; i<numCanon; i++){
    	orbits(canon_list[i],orbit[i]);
        for(j=0;j<k;j++){
	    if(orbit[i][j]==j)
		orbit[i][j]=numOrbits++;
	    else
                orbit[i][j]=orbit[i][orbit[i][j]];
        }
    }

    //output
    printf("%d\n", numOrbits);
    for(i=0; i<numCanon; i++){
	for(j=0;j<k;j++) {
	    printf("%ld ", orbit[i][j]);
	}
	printf("\n");
    }
    return 0;
}
