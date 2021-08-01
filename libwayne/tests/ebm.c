// Empirical Brown's Method, adapted from ebm.py (Python version)
// see Poole (2016) 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "misc.h"
#include "stats.h"

#define MAX_M 20000 // rows: bound on number of distinct GO terms we expect in one alignment
#define MAX_N 20000 // columns: bound on number of proteins in the smaller species

static float data[MAX_M+1][MAX_N+1]; // waste row 0 and column 0 to be consistent with AWK version.
double pVal[MAX_M+1];
static int NR, NF; // just as in AWK, indexed from *1* (not 0 as in normal C)

// No line number, always uses current line NR
static void TransformData(STAT *normSamples, int n) {
    int j;
    for(j=1;j<=n;j++) data[NR][j] = -2*log(StatECDF(normSamples, data[NR][j]));
}


void ReadLines(FILE *fp) {
    static int i;
    static char buf[30*MAX_N], *startWord, *tab; //estimate maximum 30 characters per number (indexed from 0 as normal)
    STAT *samples = StatAlloc(0,0,0,0,true), *normSamples = StatAlloc(0,0,0,0,true);

    fgets(buf, sizeof(buf), fp); // read and discard header line--it's NR==0 because we discard it before AWK.

    NR=0;
    while(fgets(buf, sizeof(buf), fp)) {
	++NR; // only gets incremented if there's actually a line there
	assert(NR<=MAX_M);
	int col=0, n=0, len; // len = strlen(word)
	startWord=buf;
	do {
	    assert(!isspace(*startWord));
	    ++col;
	    tab=strchr(startWord,'\t');
	    if(tab) {assert(*tab=='\t'); *tab='\0'; len=strlen(startWord);}
	    else {len=strlen(startWord); assert(startWord[--len]=='\n'); startWord[len]='\0';}
	    if(col==1) ; // do nothing, it's the name of the line
	    else if(col==2) {//fprintf(stderr,"line %d col %d word %s ", NR,col, startWord);
		pVal[NR]=atof(startWord); // fprintf(stderr,"%g\n",pVal[NR]);
	    }
	    else {assert(col-2>=1 && col-2<=MAX_N); StatAddSample(samples, (data[NR][++n]=atof(startWord)));}
	    if(tab) {*tab='\t'; startWord=tab+1;} // put the tab back, but not the ending newline
	} while(tab);
	if(!NF) {assert(NR==1); NF=col;}
	else assert(col == NF);
	assert(n==NF-2);
	assert(samples->n == n);
	//fprintf(stderr, "DATA[%d]",NR);for(i=1;i<=n;i++)fprintf(stderr, " %g",data[NR][i]);fprintf(stderr,"\n");

	for(i=1;i<=n;i++){
	    data[NR][i] = (data[NR][i] - StatMean(samples))/StatStdDev(samples);
	    StatAddSample(normSamples, data[NR][i]);
	}
	//fprintf(stderr, "NORM[%d]",NR);for(i=1;i<=n;i++)fprintf(stderr," %g",data[NR][i]);fprintf(stderr, "\n");
	//fprintf(stderr, "ECDF[%d]",NR);for(i=1;i<=n;i++)fprintf(stderr, " [%d][%d][%g]=%g",NR,i,data[NR][i],StatECDF(samples,data[NR][i]));fprintf(stderr,"\n");
	assert(normSamples->n == n);
	TransformData(normSamples, n);
	//fprintf(stderr, "-LOG[%d]",NR); for(i=1;i<=n;i++)fprintf(stderr, " %g",data[NR][i]); fprintf(stderr, "\n");
	StatReset(samples); StatReset(normSamples);
    }
    //fprintf(stderr, "last line has pVal[%d]=%g\n", NR, pVal[NR]);
}

int main(int argc, char *argv[])
    {
	ReadLines(stdin);
	int m=NR, n=NF-2, i, j, k;
	if(m<1) Fatal("need at least one variable");
	assert(m<=MAX_M);
	COVAR *C[MAX_M+1][MAX_M+1];

	// Compute covariances across all samples of all input variables.
	for(i=1;i<=m;i++) for(j=i+1;j<=m;j++){
	    C[i][j]=CovarAlloc();
	    for(k=1;k<=n;k++) CovarAddSample(C[i][j], data[i][k], data[j][k]);
	}
	// Now perform Empirical Browns Method
	double df_fisher, Expected, covar, cov_sum, df_brown, Var, c, x, pProd, p_brown, log_p_brown, log_pProd;
	df_fisher = Expected = 2.0*m;
	cov_sum = 0;
	for(i=1;i<=m;i++) { 
	    //fprintf(stderr, "row %4d", i);
	    for(j=i+1;j<=m;j++) {
		covar=Covariance(C[i][j]);
		//fprintf(stderr, " %.5f", covar);
		cov_sum += covar;
	    }
	    //fprintf(stderr, "\n");
	}
	//printf("cov_sum %g\n",cov_sum);
	Var = 4.0*m+2*cov_sum;
	c = Var/(2.0*Expected);
	df_brown = 2.0*SQR(Expected)/Var;
	if(df_brown > df_fisher) {
	    df_brown = df_fisher;
	    c = 1.0;
	}
	//printf("c = %g\n", c);
	pProd=1; log_pProd=0;
	x=0; // twice the sum of logs of p-values
	for(i=1;i<=m;i++)
	    if(pVal[i]>0) {
		x += -log(pVal[i]); log_pProd+=log(pVal[i]); pProd *= pVal[i];
	    }
	    else fprintf(stderr, "skipping pVal[%d]=%g\n",i, pVal[i]);
	x *= 2;
	//printf("x = %g\n", x);
	log_p_brown = logChi2_pair((int)(df_brown+0.5), 1.0*x/c);
	p_brown = Exp(log_p_brown);

	if(p_brown < pProd)
	    Fatal("Oops, something wrong: p_brown should be < product(pVals), but p_brown = %g while product = %g", p_brown, pProd);
	printf("%g = %g bits; %g = %g bits\n", p_brown, -log_p_brown/log(2.0), pProd, -log_pProd/log(2));
    }
