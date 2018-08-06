#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>

#define sl3 8
#define VP3 6
#define sl4 64
#define VP4 24
#define sl5 1024
#define VP5 120
#define sl6 32768
#define VP6 720
#define sl7 2097152
#define VP7 5040


typedef int permDimen3[VP3+1]; 
typedef permDimen3 decDimen3[sl3]; 
typedef int permDimen4[VP4+1]; 
typedef permDimen4 decDimen4[sl4]; 
typedef int permDimen5[VP5+1]; 
typedef permDimen5 decDimen5[sl5]; 
typedef int permDimen6[VP6+1]; 
typedef permDimen6 decDimen6[sl6]; 
typedef int* permDimen7; 
typedef permDimen7 decDimen7[sl7]; 

static int can_dec[sl7];

static decDimen3 can_perm3;
static decDimen4 can_perm4;
static decDimen5 can_perm5;
static decDimen6 can_perm6;
static decDimen7 can_perm7;

typedef unsigned long long int ullint;

int** canonical(ullint x, int k1, ullint* can_dec, int* num_canPerm3,int permutation[5040][7], int c);

ullint factorial(int n)
{
        if(n==0)return 1;
        return (ullint)n*factorial(n-1);
}

bool nextPermutation(int permutation[], int s)
{
        int t, i;
        for(i=s-1;i>0;i--)
        {
                if(permutation[i]>permutation[i-1])
                {
			int j;
                        for(j=s-1;j>i-1;j--){
                                if(permutation[i-1]<permutation[j]){
                                        int t=permutation[i-1];
                                        permutation[i-1]=permutation[j];
                                        permutation[j]=t;
                                        break;
                                }
                        }
                        int l=i;
                        for(j=s-1;j>l;j--)
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

ullint bitArrayToDecimal(int k, int** bitarray, int Permutations[], int bitVectorSize){
        ullint num=0;
        int lf=0,i,j;
        for(i = 0; i < k-1; i++)
                for(j=i+1; j < k; j++){
                        num+=(((ullint)bitarray[Permutations[i]][Permutations[j]]) << (bitVectorSize-1-lf));
                        lf++;
                }
        return num;
}

int **permTable(ullint Rows, int Cols);

int** decimalToBitArray(int k, ullint D){
        int i,j;
	int** bitarray=permTable(k,k);
	for(i=k-2; i>=0; i--)
                for(j=k-1; j>=i+1; j--){
                        bitarray[i][j] = D%2;
                        bitarray[j][i]=bitarray[i][j];
                        D = D/2;
                }
	for(i=0;i<=k-1;i++)bitarray[i][i]=0;
	return bitarray;
}



int **permTable(ullint Rows, int Cols)
{
    	int **pTable = (int **)malloc(Rows * sizeof(int *));
    	ullint row;
    	for (row = 0; row < Rows; row++) {
    		pTable[row] = (int *)malloc(Cols * sizeof(int));
    	}
    	return pTable;
}

int *createColumns(int cols)
{
    return (int*)malloc(cols * sizeof(int));
}

void RepermTable(int **pTable, ullint Rows, int Cols)
{
     	pTable = (int **)realloc(pTable, Rows * sizeof(int *));
     	ullint row;
     	for (row = 0; row < Rows; row++) {
     		pTable[row] = (int *)malloc(Cols * sizeof(int));
     	}
}

void free_table(int **table, ullint Rows)
{
    	ullint row;
    	for (row = 0; row < Rows; row++) {
    		free(table[row]);
    	}
    	free(table);
}

//given some groups of nodes, it creates all possible permutations such that no nodes from one group is swapped with any other node from a different group
int** selectedPerms(int grp[], int s, int *pSize)
{
	int gValue=grp[0], grps=1, grpm[s], i, j;
    	for(i=0; i<s; i++) grpm[i]=0;
    	int permSize=1;
    	for(i=0; i<s; i++)
    	{
    	    	if(gValue==grp[i]) grpm[grps-1]++;
            	else {
                    	gValue=grp[i];
                    	permSize*=factorial(grpm[grps-1]);
                    	grps++;
                    	grpm[grps-1]++;
            	}
    	}
     	permSize*=factorial(grpm[grps-1]);
     	*pSize=permSize;
     	int **sPerms=permTable(permSize, s);
     	int ns=s;
     	int z, w, a;
     	int temp=permSize, temp1;

     	for(z=grps-1; z>=0; z--){
        	int m=grpm[z];
        	int gArr[m];
        	a=0;
        	temp1=permSize/temp;
        	temp/=factorial(m);
        	for(w=0; w<temp1; w++){
        	        for(i=0; i<m; i++) gArr[i]=ns-m+i;
                	do{
                	        for(i=a*temp; i<(a+1)*temp; i++)
                	                for(j=0; j<m; j++) sPerms[i][ns-m+j] = gArr[j];
                	        a++;
                	}while(nextPermutation(gArr, m));
        	}
        	ns-=m;
     	}
     	return sPerms;
}

//finds the minimum degree and the nodes that has the minimum degree
void degree(int k1, int** bitArray, int* minDegree, int minDegNodes[], int *num_minDegNodes)
{
	int minDeg=k1;
	int p, i, j, s=0, x;
	int tmp=0;
	for(i=0; i<k1; i++){
		s=0;
                for(j=0; j<k1; j++){
			s+=bitArray[i][j];
		} 
		if(minDeg>s){
			minDeg=s;
			minDegNodes[0]=i;
			p=1;
			tmp=1;
		}
		else if(minDeg==s){
			minDegNodes[p++]=i;
			tmp++;
		}
	}
	*minDegree=minDeg;
	*num_minDegNodes=tmp;
}

//makes bitArray for the subgraph with nodes 1,2,..,k-minDeg-1 after the first column is minimized by "minimize_firstCol"
void make_subgraph(int k1, int** bitArray1, int minDeg, int** bitArray2, int perm[])
{
	int i, j;
	for(i=1; i<k1-1-minDeg; i++){
		for(j=0; j<i; j++){
			bitArray2[i][j]=bitArray1[perm [i+1]][perm[j+1]];
			bitArray2[j][i]=bitArray2[i][j];
		}
	}
}

void modify_perm(int can_perm[], int perm[], int start, int end, int Perm[])
{
	int u=end-start+1;
	int m, i, tmp_perm[u];
	for(i=0;i<u;i++)tmp_perm[i]=Perm[i+start];
	for(i=start; i<=end; i++){
		perm[i]=tmp_perm[can_perm[i-start]];
	}
}	

void swap(int perm[], int i, int j){
	int u=perm[i];
	perm[i]=perm[j];
	perm[j]=u;
}

//minimizes the first column 
void minimize_firstCol(int k1, int minDegNode, int bitArr[], int perm[], int minDeg)
{
	int i, j;
	swap(perm, minDegNode, 0);
	j=k1-1;
	if(minDeg!=0){
	for(i=k1-1; i>0; i--)
		if(bitArr[perm[i]]==1)swap(perm, i, j--);}
}

//After the first column and the subgraph with nodes 1,2,...,k-minDeg-1 is minimized, this function minimizes the rest part of the graphette.
void part_of_dec(int k1, int perm[], int start, int end, int** bitArray, int c, int grp[],int u)
{
	if(c<u&&start!=end){
        	int i, j=end;
        	for(i=end; i>=start; i--)
                	if(bitArray[perm[i]][perm[c]]==1)swap(perm, i, j--);
		if(j<=end-2){
			part_of_dec(k1,perm, j+1, end, bitArray, c+1, grp, u);
			for(i=j+1; i<=end; i++) grp[i-start]=j+1;
		}
		else if(j==end-1) grp[end]=end;
		if(j>=start+1){
			part_of_dec(k1,perm, start, j, bitArray, c+1, grp, u);
		}
	}
}

//reads the files
void getCanons(int i){
	int u=(1<<(i*(i-1)/2));
	char buffer[256]={0};
        snprintf(buffer, 255, "data%d.txt", i);
        FILE *file=fopen(buffer, "r");
        if(!file) {
                printf("%s not found\n", buffer);
                return;
        }
        int l=0;
        char ch;
        ullint temp;
        int tempI;
	int j=0;
	int tempIndex[5040], inx;
        while(l<u)
        {
   
        	j=0;
                temp=0;
                while(1){
                	ch=fgetc(file);
                	if(ch==' ' || ch=='\n') break;
                	temp=10*temp+(ch-'0');
                }
                if(i==7) can_dec[l]=temp;
                while(1){
                	tempI=0;
                	while(1){
                        	ch=fgetc(file);
                        	if(ch==' '|| ch=='\n') break;
                        	tempI=10*tempI+(ch-'0');
                	}
                	if(ch!='\n'){
				switch(i){
                        	        case 3: can_perm3[l][j]=tempI; break;
                                	case 4: can_perm4[l][j]=tempI; break;
                                	case 5: can_perm5[l][j]=tempI; break;
                                	case 6: can_perm6[l][j]=tempI; break;
                                	case 7: tempIndex[j]=tempI; break;
                        	}
			}
                	j++;
                	if(ch=='\n') break;
            	}
	    	switch(i){
            		case 3: can_perm3[l][VP3]=j-1;break;
            		case 4: can_perm4[l][VP4]=j-1;break;
            		case 5: can_perm5[l][VP5]=j-1;break;
            		case 6: can_perm6[l][VP6]=j-1;break;
            		case 7: can_perm7[l]=createColumns(j);can_perm7[l][0]=j-1; for(inx=1; inx<j; inx++) can_perm7[l][inx]=tempIndex[inx-1];break; 
            	}
	    	l++;
        }
        fclose(file);
}

int** create(int k1, int minDegNode,int minDeg, int** bitArray,ullint* can_dec, ullint x, int* num_canPerm, int permutation[5040][7], int c)
{
	ullint d;
	int i,i1,j;
	int perm[k1];
	 int bitvectorsize1=k1*(k1-1)/2;
	for(i=0; i<k1; i++)perm[i]=i;
	if(k1-minDeg-1>0){minimize_firstCol(k1,minDegNode,bitArray[minDegNode],perm,minDeg);
	if(k1-minDeg-1>2){int s=k1-1-minDeg;
	int **bitArray2=permTable((k1-minDeg-1),(k1-minDeg-1));
	make_subgraph(k1,bitArray,minDeg,bitArray2,perm);
	int perm1[k1-minDeg-1];
	for(i=0; i<k1-minDeg-1; i++)perm1[i]=i;
	int bitvectorsize2=(k1-1-minDeg)*(k1-minDeg-2)/2;
	d=bitArrayToDecimal(k1-1-minDeg,bitArray2,perm1,bitvectorsize2);
	free_table(bitArray2,k1-minDeg-1);}}
	int** sub_perm;
	int num_canPerm2;
	int t;
	ullint cd;
	if(k1-1-minDeg>7) sub_perm=canonical(d,k1-1-minDeg,&cd,&num_canPerm2,permutation,c+1);
	else{
		switch(k1-1-minDeg){
			case 0:
				num_canPerm2=1;
				break;
			case 1: 
				sub_perm=permTable(1,1);
				sub_perm[0][0]=0;
				num_canPerm2=1;
				break;
			case 2: 
				sub_perm=permTable(2,2);
				sub_perm[0][0]=0;sub_perm[0][1]=1;
				sub_perm[1][0]=1;sub_perm[1][1]=0;
				num_canPerm2=2;
				break;
			case 3: 
				sub_perm=permTable(can_perm3[d][6],3);
				for(t=0; t<can_perm3[d][6];t++)
					for(j=0; j<3; j++)
						sub_perm[t][j] = permutation[can_perm3[d][t]][j+4]-4;
				num_canPerm2=t;
				break;
			case 4: 
				sub_perm=permTable(can_perm4[d][24],4);
                                for(t=0; t<can_perm4[d][24];t++)
                                        for(j=0; j<4; j++)
                                                sub_perm[t][j] = permutation[can_perm4[d][t]][j+3]-3;
				num_canPerm2=t;
				break;
			case 5: 
				sub_perm=permTable(can_perm5[d][120],5);
                                for(t=0; t<can_perm5[d][120];t++)
                                        for(j=0; j<5; j++)
                                                sub_perm[t][j] = permutation[can_perm5[d][t]][j+2]-2;
				num_canPerm2=t;

				break;
			case 6: 
				sub_perm=permTable(can_perm6[d][720],6);
                                for(t=0; t<can_perm6[d][720];t++)
                                        for(j=0; j<6; j++)
                                                sub_perm[t][j] = permutation[can_perm6[d][t]][j+1]-1;
				num_canPerm2=t;
				break;
			case 7: 
				sub_perm=permTable(can_perm7[d][0],7);
                                for(t=0; t<can_perm7[d][0];t++)
                                        for(j=0; j<7; j++)
                                                sub_perm[t][j] = permutation[can_perm7[d][t+1]][j];
				num_canPerm2=t;
				break;
		}
	}

	int grp[minDeg];
	ullint fac=factorial(k1-1);
	int u=k1-minDeg; 
	int pSize;
	int l,j1=0,i2;
	ullint y,z=x+1;
	int **Perm1;
	Perm1=permTable(fac, k1);
	for(i=0; i<minDeg; i++) grp[i]=u;
	int** sPerms;
	int m;
	j1=0;
	int Perm[k1], perm2[k1];
	for(i=0;i<k1;i++) Perm[i]=perm[i];
	if(k1-minDeg-1>0){
	for(i=0;i<num_canPerm2;i++){

		modify_perm(sub_perm[i],perm,1,k1-minDeg-1,Perm);
		if(minDeg>1){
			for(i=0; i<minDeg; i++) grp[i]=u;
			part_of_dec(k1,perm,k1-minDeg,k1-1,bitArray,1,grp,u);
			sPerms = selectedPerms(grp, minDeg, &pSize);
			for(i=0;i<k1;i++) perm2[i]=perm[i];}
		else
			pSize=1;

		for( i1=0;i1<pSize;i1++){
			if(minDeg>1)modify_perm(sPerms[i1],perm2,k1-minDeg,k1-1,perm);
			y=bitArrayToDecimal(k1,bitArray,perm,bitvectorsize1);
			if(i==0&&i1==0)z=y;
			if(z>y){
				z=y;
				j1=0;
                                for(l=0;l<k1;l++)
                                        Perm1[j1][l]=perm[l];
				j1++;
			}
			else if(z==y){
				for(l=0;l<k1;l++)
                                        Perm1[j1][l]=perm[l];
				j1++;
			}
		}
		if(minDeg>1)free_table(sPerms,pSize);
		if(minDeg==0&&c==0)break;
	}
	}
	else if(c>0&&k1-minDeg-1==0){
		int tmpPerm[k1];
        	for(i=0;i<k1;i++)tmpPerm[i]=i;
		for(i=0;i<fac;i++)
        	{
                	for(j=0; j<k1; j++)
                	{
                        	Perm1[i][tmpPerm[j]]=j;
                	}
                nextPermutation(tmpPerm,k1);
        	}
		z=x;
	}
	else{
		for(i=0;i<k1;i++)Perm1[0][i]=i;
		z=x;
	}
	if(k1-1-minDeg>7) free_table(sub_perm,factorial(k1-1-minDeg));
        else{
                switch(k1-1-minDeg){
			case 1:
				free_table(sub_perm,1);break;
			case 2:
				free_table(sub_perm,2);break;
                        case 3:
                                free_table(sub_perm,can_perm3[d][6]);break;

                        case 4:
                                free_table(sub_perm,can_perm4[d][24]);break;

                        case 5:
                                free_table(sub_perm,can_perm5[d][120]);break;

                        case 6:
                                free_table(sub_perm,can_perm6[d][720]);break;

                        case 7:
                                free_table(sub_perm,can_perm7[d][0]);break;

                }
	}
	*can_dec=z;
	*num_canPerm=j1;
	return Perm1; 	
}

int** canonical(ullint x, int k1, ullint* can_dec, int* num_canPerm3, int permutation[5040][7], int c)
{
	int **Perm;
        int **bitArray;
	int i,j;
        bitArray=decimalToBitArray(k1,x);
	int minDegree=k1;
	int minDegNodes[k1];
	int num_minDegNodes;
	int** canPerm;
	int num_canPerm;	
        degree(k1,bitArray,&minDegree,minDegNodes,&num_minDegNodes);
	j=0;
	int l,j1;
	int m;
	ullint z;
	ullint w=x+1;
	ullint fac=factorial(k1);
	Perm=permTable(fac, k1);
	ullint fac1=fac/k1;
	for(i=0;i<num_minDegNodes;i++){
                canPerm=create(k1,minDegNodes[i],minDegree,bitArray,&z,x,&num_canPerm,permutation,c);

                if(i==0)w=z; 
		if(w>z){
			w=z;
			j=0;
			for(j=0;j<num_canPerm;j++)
			for(l=0;l<k1;l++)
				Perm[j][l]=canPerm[j][l];
		}
		else if(w==z){
			for(j1=j;j1<j+num_canPerm;j1++)
                        for(l=0;l<k1;l++)
                                Perm[j1][l]=canPerm[j1-j][l];
                 	j=j1;       
        	}
		*num_canPerm3=j;
		free_table(canPerm,fac1);
		if(c==0)if(minDegree==0||minDegree==k1-1)break;
	}
	free_table(bitArray,k1);
	*can_dec=w;
	return Perm;
}

int main(int argc, char *argv[])
{
	if(argc < 3){fprintf(stderr, "expecting at least two arguments, which are value of k, and decimal numbers"); exit(1);}
	int i,j;
        int tmpPerm[7];
        for(i=0;i<7;i++)tmpPerm[i]=i;
	int permutation[5040][7];

        //saving all permutations
        for(i=0;i<5040;i++)
	{
         	for(j=0; j<7; j++)
         	{
         		permutation[i][tmpPerm[j]]=j;
                }
        	nextPermutation(tmpPerm,7);
        }

	//reading the canon_map files
	for(i=0;i<5;i++)
  	  	getCanons(i+3);

	int k=atoi(argv[1]);
	ullint num;
	ullint f=factorial(k);
	for(i=2;i<argc;i++){
		num=atoi(argv[i]);
		int numPerm;
		ullint d;
		int **perm;
		perm=canonical(num,k,&d,&numPerm,permutation,0);
		printf("d=%d ",d);
		for(j=0;j<k;j++)printf("%d",perm[0][j]);
		printf("\n");
		free_table(perm,f);
	}
	
return 0;
}


