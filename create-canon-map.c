#include <stdio.h>
#include <sys/file.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>


#ifndef k
#define k 7
#endif
#if k <= 8
#define q (1 << k*(k-1)/2)
#else
#define q (1LL << k*(k-1LL)/2)
#endif

typedef unsigned char xChar[5];//40 bits for saving index of canonical decimal and permutation

static xChar data[q];
static long int canonicalDecimal[274668LL];//274668 canonical graphettes for k=9

long power(int x, int y){
        if(y==0)return 1;
        return x*power(x,y-1);
}

void encodeChar(xChar ch, long indexD, long indexP){

        long int x=indexD+indexP*power(2,19);//19 bits for canonical decimal index
        long z=power(2,8);
        for(int i=4; i>=0; i--){
                ch[i]=(char)(x%z);
                x/=power(2,8);
        }
}

void decodeChar(xChar ch, long* indexD, long* indexP){

        long int x=0,m;
        int y=0,w;

        for(int i=4; i>=0; i--){
                w=(int)ch[i];
                m=power(2,y);
                x+=w*m;
                y+=8;
        }
        long z=power(2,19);
       *indexD=x%z;
       *indexP=x/z;                  
}

long factorial(int n)
{
  	if(n==0)return 1;
  	return n*factorial(n-1);
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
    bool check[q];

    FILE *fcanon = stdout;
  
    int bitVectorSize = (k*(k-1))/2;
    long int D;
    int bitarray[k][k];
    
    for(long int i=0; i<q; i++)check[i]=0;
    canonicalDecimal[0]=0;
    long f=factorial(k);
    char Permutations[f][k];
    int tmpPerm[k];
    for(int i=0;i<k;i++)tmpPerm[i]=i;
   
    //saving all permutations
    for(long i=0;i<f;i++){
       	for(int j=0; j<k; j++)
		Permutations[i][j]=tmpPerm[j]+'0';
    	nextPermutation(tmpPerm);
    }

    encodeChar(data[0],0,0);
    unsigned long num_canon=0;

    //finding canonical forms of all graphettes
    for(unsigned long int t=1; t<q; t++){
        if(check[t]) continue;
	check[t]=1;
	encodeChar(data[t],++num_canon,0);
	canonicalDecimal[num_canon]=t;
             
        unsigned long int num = 0;
        D=t;
        for(int i=k-2; i>=0; i--){
            for(int j=k-1; j>i; j--){
                bitarray[i][j] = D%2;
                bitarray[j][i]=bitarray[i][j];
                D = D/2;
            }
        }

	int permutation[k];
	for(int i = 0; i < k; i++){
            permutation[i]=i;
        }
        long num_perm=0;

	//finding all isomorphic forms of graphette with decimal t
        while( nextPermutation(permutation) ){
            num=0;
            int f=0;
            for(int i = 0; i < k-1; i++)
                for(int j=i+1; j < k; j++){
		   num+=(((long)bitarray[permutation[i]][permutation[j]]) << (bitVectorSize-1-f));
                   f++;
                }
	    num_perm++;
            if(!check[num]){
		check[num]=true;
		encodeChar(data[num],num_canon,num_perm);
	    }
        }

    }

	//saving canonical decimal and permutation in the file
        long canonDec, canonPerm;

        for(int i=0; i<q; i++){
	    decodeChar(data[i],&canonDec,&canonPerm);
            fprintf(fcanon,"%lld ", canonicalDecimal[canonDec]);
            for(int p=0;p<k;p++)
                fprintf(fcanon,"%c", Permutations[canonPerm][p]);
            fprintf(fcanon,"\n");
        }
}

int main(int argc, char* argv[]){
   canon_map();
   return 0;
}
