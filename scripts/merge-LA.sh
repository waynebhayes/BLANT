#!/bin/sh
USAGE="$0 G1.el G2.el Patrick-Alig-File"
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
s3=${s3:?"environment variable s3 should be set to lower bound of desired S3"}

NL='
'

G1="$1"
G2="$2"
A="$3"
#if [ "`tail -c1 "$A"`" != "$NL" ]; then
#    #echo "Appending newline to Patrick alig file" >&2
#    (cat "$A"; echo) > $TMPDIR/A
#    A=$TMPDIR/A
#fi

hawk '
    function S3(A,  u,v,sum,numer,denom) {
	if(isarray(A)) {
	    for(u in A)for(v in A) if(u<v) {
		sum=edge[1][u][v] + edge[2][A[u]][A[v]];
		if(sum>0) denom++;
		if(sum==2) numer++;
	    }
	}
	if(denom==0) return 0;
	return numer/denom;
    }
    function S3pair(U,i,j, v) { delete U;
	for(v in LA[i]) {
	    if(v in LA[j]) {if(LA[i][v]==LA[j][v]) U[v]=LA[i][v];} # the pair is mutual
	    else # v is NOT in the other alignment, and if its pairing is also not in the other alignment, add it to U
		if(!(j in LApair[v","LA[i][v]])) U[v]=LA[i][v];
	}
	for(v in LA[j]) {
	    if(v in LA[i]) {if(LA[j][v]==LA[i][v]) U[v]=LA[j][v];} # the pair is mutual
	    else # v is NOT in the other alignment, and if its pairing is also not in the other alignment, add it to U
		if(!(i in LApair[v","LA[j][v]])) U[v]=LA[j][v];
	}
	if(!isarray(U)) {U[0]=0; delete U[0]; return 0;}
	return S3(U);
    }
    ARGIND<=2{edge[ARGIND][$1][$2]=edge[ARGIND][$2][$1]=1; ++D[ARGIND][$1];++D[ARGIND][$2]}
    ARGIND==3{
	if(NF!=3)exit;
	n1=split($2,V1,","); n2=split($3,V2,",");
	ASSERT(n1==n2,"size mismatch");
	for(i=1;i<=n1;i++) {LApair[V1[i]","V2[i]][FNR]=1; LA[FNR][V1[i]]=V2[i]}
	ASSERT(S3(LA[FNR])>='"$s3"', "S3 too low on line "FNR);

	j=FNR; for(i=1;i<FNR;i++) {
	    li=length(LA[i]); lj=length(LA[j]); s3=S3pair(P,i,j); lp=length(P);
	    if(lp>=MIN(1.5*MAX(li,lj)+1, biggest*0.9) && s3>='"$s3"') { printf "%d,%g\t",lp,s3; comma="";
		for(v in P){printf "%s%s",comma,v;comma=","}
		printf "\t"; comma="";
		for(v in P){printf "%s%s",comma,P[v];comma=","}
		print ""
	    }
	    if(lp>biggest)biggest=lp;
	}
    }
    ' "$G1" "$G2" "$A"
