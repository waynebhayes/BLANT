#include <stdio.h>
#include <sys/file.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "blant.h"

static int k; // number of bits required to store ints up to 2^(k^2)
static long numBitValues; // the actual value 2^(k^2)



unsigned long bitArrayToDecimal(int bitMatrix[k][k], char Permutations[], int numBits){
    unsigned long num=0;
    int lf=0;
    for(int i = 0; i < k; i++)
    for(int j = 0; j < k; j++){
	if(i==j) continue;
	num+=(((unsigned long)bitMatrix[(int)Permutations[i]][(int)Permutations[j]]) << (numBits-1-lf)); 
	lf++;
    }
    return num;
}

void decimalToBitArray(int bitMatrix[k][k], unsigned long D){
    for(int i=k-1; i>=0; i--)
	for(int j=k-1; j>=0; j--){
	    if(i==j) continue;
	    bitMatrix[i][j] = D%2;
	    D = D/2;
	}
}




typedef unsigned char xChar[5];//40 bits for saving index of canonical decimal and permutation

static xChar* data;
static bool* done;
static unsigned long canonicalDecimal[1000000];//placeholder

unsigned long power(int x, int y){
    assert(x>0 && y>=0);
    if(y==0)return 1;
    return (unsigned long)x*power(x,y-1);
}

void encodeChar(xChar ch, long indexD, long indexP){ 
    unsigned long x=(unsigned long)indexD+(unsigned long)indexP*(1<<20);//2^20 > 1e6 (placeholder values, enough for k=5)
    for(int i=4; i>=0; i--){
	ch[i]=(char)(x%(1<<8));
	x>>=8;
    }
}

void decodeChar(xChar ch, long* indexD, long* indexP){

    unsigned long x=0,m;
    int y=0,w;

    for(int i=4; i>=0; i--){
	w=(int)ch[i];
	m=(1<<y);
	x+=w*m;
	y+=8;
    }
    unsigned long z=(1<<20);
    *indexD=x%z;
    *indexP=x/z;
}

long factorial(int n) {
    assert(n>=0);
    if(n==0)return 1L;
    return (long)n*factorial(n-1);
}

bool nextPermutation(int permutation[]) {
    for(int i=k-1;i>0;i--) {
	if(permutation[i]>permutation[i-1]) { 
	    for(int j=k-1;j>i-1;j--){
		if(permutation[i-1]<permutation[j]){ 
		    int t=permutation[i-1];
		    permutation[i-1]=permutation[j];
		    permutation[j]=t;
		    break;
		}
	    }
	    int l=i;
	    for(int j=k-1;j>l;j--) {
		if(i<j){
		    int t=permutation[i];
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

void canon_map(bool directed){
    FILE *fcanon = stdout;

    int numBits = k*(k-1);
    int bitMatrix[k][k];

    for(int i=0; i<numBitValues; i++) assert(i>=0), done[i]=0;
    canonicalDecimal[0]=0;
    int f=factorial(k);
    char Permutations[f][k];
    int tmpPerm[k];
    for(int i=0;i<k;i++)tmpPerm[i]=i;

    //saving all permutations
    for(int i=0;i<f;i++){
	assert(i>=0);
	for(int j=0; j<k; j++)
	    Permutations[i][j]=tmpPerm[j];
	nextPermutation(tmpPerm);
    }
    done[0]=1;
    encodeChar(data[0],0,0);
    long num_canon=0;

    //finding canonical forms of all graphettes
    for(int t=1; t<numBitValues; t++){
	assert(t>=0);
	if(done[t]) continue;
	done[t]=1; // this is a new canonical, and it the lowest by construction
	encodeChar(data[t],++num_canon,0);
	canonicalDecimal[num_canon]=t;

	int num = 0;
	decimalToBitArray(bitMatrix, t);
	for(int nP=1; nP<f; nP++) // now go through all the permutations to compute the non-canonicals of t.
	{
	    assert(nP>0);
	    num=bitArrayToDecimal(bitMatrix, Permutations[nP], numBits);
	    assert(num>=0);
	    if(!done[num]){
		done[num]=true;
		encodeChar(data[num],num_canon,nP);
	    }
	}

    }

    // It's easiest to generate the permumations from the non-canonical to the canonical,
    // but the inverse is far more useful. Here we invert all the permutations so that
    // they map from canonical nodes to the non-canonical. That way, we can spend most of
    // our time thinking in "canonical space", and if we want to know where to find canonical
    // node j in a particular non-canonical, we use perm[j].
    if(PERMS_CAN2NON){
	int tmp[k];
	for(int i=0; i<f; i++){
	    for(int j=0; j<k; j++)
	    tmp[(int)Permutations[i][j]]=j;
	    for(int j=0; j<k; j++)
	    Permutations[i][j]=tmp[j];
	}
    }
    fprintf(stderr, "Finished computing... now writing out canon_map file\n"); fflush(stderr);

    //saving canonical decimal and permutation in the file
    long canonDec, canonPerm;
    TINY_GRAPH *G = TinyGraphAlloc(k,0,directed);
    for(unsigned long i=0; i<numBitValues; i++){
	char printPerm[k+1];
	printPerm[k]='\0';
	canonDec=0;canonPerm=0;
	decodeChar(data[i],&canonDec,&canonPerm);
	assert(canonDec >= 0);
	assert(canonPerm >= 0);
	for(int p=0;p<k;p++) printPerm[p] = '0' + Permutations[canonPerm][p];
	fprintf(fcanon,"%lu\t%s", canonicalDecimal[canonDec], printPerm);
	if(canonPerm == 0) {
	    TinyGraphEdgesAllDelete(G);
	    Int2TinyGraph(G, i);
	    int nodeArray[k], distArray[k];
	    fprintf(fcanon, "\t%c", TinyGraphNumEdges(G));
	    int u,v,sep='\t';
	    for(u=0;u<k;u++)for(v=(1-directed)*u;v<k;v++) if(TinyGraphAreConnected(G,u,v)) {
		fprintf(fcanon, "%c%d,%d",sep,u,v); sep=' ';
	    }
	}
	putc('\n', fcanon);
    }
    //fprintf(fcanon,"%ld",canonicalDecimal[0]);
}

static char USAGE[] = "USAGE: $0 k directed [1 or 0]";

int main(int argc, char* argv[]){
    if(argc != 3){fprintf(stderr, "expecting exactly two arguments, which are k and directed\n%s\n",USAGE); exit(1);}
    k = atoi(argv[1]);
    bool directed = atoi(argv[2]);
    numBitValues = (1UL << (k*(k-1)));
    assert(numBitValues>0);
    data = malloc(sizeof(xChar)*numBitValues);
    done = malloc(sizeof(bool)*numBitValues);
    canon_map(directed);

    return 0;
}
