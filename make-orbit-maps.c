#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "blant.h"


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static int k;


bool suitablePerm(int permutation[], int adj[k][k+1]);
void permuteNodes(int permutation[], int adj[k][k+1], long *D);
void makeOrbit(int permutation[], long orbit[]);
void getCycle(int permutation[], int cycle[], int seed, int current, bool visited[]);

bool nextPermutation(int permutation[])
{
   int t;
   for(int i=k-1;i>0;i--)
   {
       if(permutation[i]>permutation[i-1])
          {
          	for(int j=k-1;j>i-1;j--){
           		if(permutation[i-1]<permutation[j]){ 
             			int t=permutation[i-1];
             			permutation[i-1]=permutation[j];
             			permutation[j]=t;
             			break;
               		}	
             	}              
           	int l=i;
           	for(int j=k-1;j>l;j--)
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
    int permutation[k];

    //initially,the permutation is (0, 1, ..., numNodes_-1)
    //and each node is in its own orbit. we'll merge orbits later
    for(int i = 0; i < k; i++){
        permutation[i]=i;
        orbit[i]=i;
    } 
    int adj[k][k+1];
    for(int i=0; i<k; i++)
     for(int j=0; j<k+1; j++)
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

bool suitablePerm(int permutation[], int adj[k][k+1]){
      for(int node = 0; node < k; node++){
        if(adj[node][k] != adj[permutation[node]][k]){
            return false;
        }
   }
   return true;
}

void permuteNodes(int permutation[], int adj[k][k+1], long* D){
 
 //To store adjMatrix_ after permutation
    int pAdj[k][k+1];
    for(int i=0; i<k; i++)
     for(int j=0; j<k+1; j++)
        pAdj[i][j]=0;
    //Apply the permutation
    for(int i = 0; i < k; i++){
        for(int j =i+1; j < k; j++){
            int pi = permutation[i], pj = permutation[j]; 
                pAdj[pi][pj] = adj[i][j];
                pAdj[pj][pi] = adj[i][j];
        }
  }

   getDecimal(pAdj, D);
   
}
void makeOrbit(int permutation[], long orbit[]){

    bool visited[k];
    for(int i=0;i<k;i++)
        visited[i]=0;
    for(int i = 0; i < k; i++){
        if(!visited[i]){
            //finding out each cycle at a time
            int cycle[k+1];
            cycle[0]=0;
	    int minOrbit=i;
            getCycle(permutation, cycle, i, i, visited);
            
            for(int j=1; j<=cycle[0]; j++)
                minOrbit = MIN(orbit[cycle[j]], minOrbit);
            for(int j=1; j<=cycle[0]; j++)
                orbit[cycle[j]] = minOrbit;
        }
    }
}

void getCycle(int permutation[], int cycle[], int seed, int current, bool visited[]){
    cycle[++cycle[0]]=current;
    visited[current] = true;
    if(permutation[current] != seed)
        getCycle(permutation, cycle, seed, permutation[current], visited);
}

void main(int argc,char* argv[]){
    if(argc != 4) {
	fprintf(stderr, "USAGE: %s k canon_list orbit_maps\n", argv[0]);
	exit(1);
    }
    k=atoi(argv[1]);
    assert(k > 2 && k <= 8);
    //reading data from canon_list file
    FILE *fp_ord=fopen(argv[2],"r");
    assert(fp_ord);
    int numCanon;
    fscanf(fp_ord, "%d",&numCanon);
    int canon_list[numCanon];
    for(int i=0; i<numCanon; i++) 
	fscanf(fp_ord, "%d", &canon_list[i]);
    fclose(fp_ord);
    //making orbit_map file
    int x = 0;
    FILE *fp;
    fp=fopen(argv[3], "w");
    for(int i=0; i<numCanon; i++){
    	long orbit[k];
    	orbits(canon_list[i],orbit);
        for(int i=0;i<k;i++){
                if(orbit[i]==i){
                orbit[i]=x++;
             }
        else{
                orbit[i]=orbit[orbit[i]];
            } 
        } 
    for(int j=0;j<k;j++)
	fprintf(fp,"%d ", orbit[j]);
    fprintf(fp,"\n");}

    fclose(fp);
    //writing the value of total number of orbits in the beginning of orbit_map file
    char txt1[]="sed -i '1i'";
    char* txt2=malloc(1);
    int co=0;
    while(x!=0){
    	txt2[co]=(x%10+'0');
	txt2=(char*)realloc(txt2,(co+2));
	x/=10;co++;
    }
    char* txt;
    char tmp;
    txt = malloc(strlen(txt1)+strlen(txt2)+strlen(argv[3])+1); 
    strcpy(txt, txt1);
    for(int i=co-1;i>co/2;i--){
	tmp=txt2[i];
	txt2[i]=txt2[co-i-1];
	txt2[co-i-1]=tmp;
    }
    strcat(txt, txt2);
    strcat(txt, " ");
    strcat(txt, argv[3]);
    system(txt);
}
