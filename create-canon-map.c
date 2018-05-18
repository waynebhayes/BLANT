#include <stdio.h>
#include <sys/file.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "blant.h"

#if __WORDSIZE < 64
    typedef graphletInt uint64;
    typedef unsigned long uint32;
#else
    typedef unsigned long uint64;
    typedef unsigned int uint32;
#endif

static int k;
#if maxK > 8
typedef uint64 graphletInt;
#else
typedef uint32 graphletInt;
#endif
static graphletInt q;

static TINY_GRAPH *G;
#if LOWER_TRIANGLE

graphletInt bitArrayToDecimal(int bitarray[k][k], char Permutations[], int bitVectorSize){
        graphletInt num=0;
        int lf=0;
        for(int i = 1; i < k; i++)
                for(int j=0; j < i; j++){
                        num+=(((graphletInt)bitarray[(int)Permutations[i]][(int)Permutations[j]]) << (bitVectorSize-1-lf));
                        lf++;
                }
        return num;
}

void decimalToBitArray(int bitarray[k][k], graphletInt D){
        for(int i=k-1; i>=1; i--)
                for(int j=i-1; j>=0; j--){
                        bitarray[i][j] = D%2;
                        bitarray[j][i]=bitarray[i][j];
                        D = D/2;
                }
}


#else

graphletInt bitArrayToDecimal(int bitarray[k][k], char Permutations[], int bitVectorSize){
	graphletInt num=0;
	int lf=0;
	for(int i = 0; i < k-1; i++)
        	for(int j=i+1; j < k; j++){
                	num+=(((graphletInt)bitarray[(int)Permutations[i]][(int)Permutations[j]]) << (bitVectorSize-1-lf));
                	lf++;
		}
	return num;
}

void decimalToBitArray(int bitarray[k][k], graphletInt D){
	for(int i=k-2; i>=0; i--)
        	for(int j=k-1; j>i; j--){
                	bitarray[i][j] = D%2;
                	bitarray[j][i]=bitarray[i][j];
                	D = D/2;
                }
}

#endif



typedef unsigned char xChar[5];//40 bits for saving index of canonical decimal and permutation

static xChar* data;
static bool* check;
static graphletInt canonicalDecimal[274668];//274668 canonical graphettes for k=9

graphletInt power(int x, int y){
	if(y==0)return 1;
	return (graphletInt)x*power(x,y-1);
}

void encodeChar(xChar ch, long indexD, long indexP){

	graphletInt x=(graphletInt)indexD+(graphletInt)indexP*power(2,19);//19 bits for canonical decimal index
	graphletInt z=power(2,8);
	for(int i=4; i>=0; i--){
		ch[i]=(char)(x%z);
		x/=z;
	}
}

void decodeChar(xChar ch, long* indexD, long* indexP){

	graphletInt x=0,m;
	int y=0,w;

	for(int i=4; i>=0; i--){
		w=(int)ch[i];
		m=power(2,y);
		x+=w*m;
		y+=8;
	}
	graphletInt z=power(2,19);
	*indexD=x%z;
	*indexP=x/z; 
}

long factorial(int n)
{
	if(n==0)return 1;
	return (long)n*factorial(n-1);
}

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

void canon_map(void){
	FILE *fcanon = stdout;

	int bitVectorSize = (k*(k-1))/2;
	graphletInt D;
	int bitarray[k][k];

	for(graphletInt i=0; i<q; i++)check[i]=0;
	canonicalDecimal[0]=0;
	long f=factorial(k);
	char Permutations[f][k];
	//char Permutations2[f][k];
	int tmpPerm[k];
	for(int i=0;i<k;i++)tmpPerm[i]=i;

	//saving all permutations
	for(long i=0;i<f;i++){
		for(int j=0; j<k; j++)
		{
			Permutations[i][j]=tmpPerm[j];
		}
		nextPermutation(tmpPerm);
	}
	check[0]=1;
	encodeChar(data[0],0,0);
	long num_canon=0;

	//finding canonical forms of all graphettes
	for(graphletInt t=1; t<q; t++){
		if(check[t]) continue;
		check[t]=1;
		encodeChar(data[t],++num_canon,0);
		canonicalDecimal[num_canon]=t;

		graphletInt num = 0;
		decimalToBitArray(bitarray, t);
		for(long nP=1; nP<f; nP++)//while( nextPermutation(permutation) )
		{
			num=bitArrayToDecimal(bitarray, Permutations[nP], bitVectorSize);
		
			if(!check[num]){
				check[num]=true;
				encodeChar(data[num],num_canon,nP);
			}
		}

	}

	//saving canonical decimal and permutation in the file
	long canonDec, canonPerm;
	if(PERMS_CAN2NON){
		int tmp[k];
		for(int i=0; i<f; i++){
			for(int j=0; j<k; j++) 
				tmp[(int)Permutations[i][j]]=j;
			for(int j=0; j<k; j++)
				Permutations[i][j]=tmp[j]; 
		}
	}

	for(graphletInt i=0; i<q; i++){
		canonDec=0;canonPerm=0;
		decodeChar(data[i],&canonDec,&canonPerm);
		assert(canonDec >= 0);
		assert(canonPerm >= 0);
		fprintf(fcanon,"%llu\t%llu\t", i,canonicalDecimal[canonDec]);
		for(int p=0;p<k;p++)
			fprintf(fcanon,"%d", Permutations[canonPerm][p]);
		if(canonPerm == 0) {
			G = TinyGraphAlloc(k);
			BuildGraph(G, i);
			int nodeArray[k], distArray[k];
			if(TinyGraphBFS(G, 0, k, nodeArray, distArray) == k)
				fprintf(fcanon, " 1"); // connected, thus graphlet
			else
	    		fprintf(fcanon, " 0"); // disconnected, thus not graphlet
		}
		fprintf(fcanon,"\n");
	}
	//fprintf(fcanon,"%lld",canonicalDecimal[0]);
}

static char USAGE[] = "USAGE: $0 k";

int main(int argc, char* argv[]){
	if(argc != 2){fprintf(stderr, "expecting exactly one argument, which is k\n%s\n",USAGE); exit(1);}
	k = atoi(argv[1]);
	assert(k <= maxK);
	if(k<=8) {
		q = 1 << k*(k-1)/2;
	} else {
		assert(sizeof(uint64) >= 8);
		q = 1LL << k*(k-1LL)/2;
	}
	data = malloc(sizeof(xChar)*q);
	check = malloc(sizeof(bool)*q);
	canon_map();
	return 0;
}

