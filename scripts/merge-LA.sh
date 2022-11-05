#!/bin/bash
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
USAGE="$0 G1.el G2.el Patrick-Alig-File"
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

s3=${s3:?"environment variable s3 should be set to lower bound of seeds on both input and merged output"}

G1="$1"
G2="$2"
A="$3"

tail -c1 "$A" | cmp - <(echo) || die "last line is missing trailing newline character; append one using:$NL${TAB}echo '' >> '$A'"

hawk '
    # NOTE: awk does not allow declarations of local variables, but it DOES allow a function to be called with fewer
    # arguments than parameters, and all parameters are local to the function. Thus, to "declare" variables that are
    # local to a function, simply list their names as extra parameters in the function declaration. By convention,
    # these extra parameters are listed after a "wide" section of whitespace in the parameter list. So for example
    # the S3 function below really only has one parameter: the alignment A; all other "parameters" are actually
    # local variables. The programmer is responsible for ensuring all such local variables are listed as extra
    # parameters; any variables that are used but not so declared will overwrite same-named global variables.
    # Also, in awk, FNR is a built-in variable storing the current line number of the current file (lines numbered from 1).
    # (There is also a built-in variable NR, which stores the *cumulative* number of lines seen so far in all input files.)
    #
    # Compute the s3 score of an alignment of a set of node pairs {(v1,v2),...} stored as an array with elements A[v1]=v2
    function S3(A,  u,v,sum,numer,denom) {
	if(isarray(A)) { # the "for" loop below will fail if A is empty.
	    for(u in A)for(v in A) if(u<v) {
		sum=edge[1][u][v] + edge[2][A[u]][A[v]];
		if(sum){ denom++; if(sum==2) numer++;}
	    }
	}
	if(denom==0) return 0;
	return numer/denom;
    }
    # Compute the S3 score resulting from merging LAs from input lines i and j; put the merged LA in U and return its s3 score.
    function S3pair(U,i,j, v) {
	delete U;
	for(v in LA[i]) {
	    if(v in LA[j]) {if(LA[i][v]==LA[j][v]) U[v]=LA[i][v];} # the aligned node pair is in both seeds
	    else # v is NOT in the other alignment, and if its pairing is also not in the other alignment, add it to U
		if(!(j in LApair[v","LA[i][v]])) U[v]=LA[i][v]; # neither node exists in LA[j] so can be added without conflict
	}
	# Now, exactly the same "for" loop as above but with i <-> j swapped.
	for(v in LA[j]) {
	    if(v in LA[i]) {if(LA[j][v]==LA[i][v]) U[v]=LA[j][v];}
	    else
		if(!(i in LApair[v","LA[j][v]])) U[v]=LA[j][v];
	}
	# if nothing was added to U, make it an explicit empty set so that length(U) works later and returns 0
	# (otherwise length(U) would return an error, "length() can only be applied to array tyes")
	if(!isarray(U)) {U[0]=0; delete U[0]; return 0;}
	return S3(U);
    }
    ARGIND<=2{edge[ARGIND][$1][$2]=edge[ARGIND][$2][$1]=1; ++D[ARGIND][$1];++D[ARGIND][$2]} # read in two edge lists
    ARGIND==3{ # start reading Patrick-formatted seed lines
	if(NF!=3)exit;
	# split 2nd and 3rd columns by commas into vertex arrays V1 and V2, with n1=|V1| and n2=|V2|.
	n1=split($2,V1,","); n2=split($3,V2,",");
	ASSERT(n1==n2, "size mismatch");
	for(i=1;i<=n1;i++) { # store the alignments in two ways:
	    # First, LApairs is a "set" indexed (in an associative array) by aligned node pairs, with a second index listing
	    # the set of lines that node pair appeared on.  eg. the aligned node pair (ABC,DEF) appearing in the seed on line
	    # 17 is encoded by setting LApairs["ABC,DEF"][17]=1.
	    LApair[V1[i]","V2[i]][FNR]=1;
	    # Second, we store the alignment vector similar to how SANA does it, with the full alignment on line 17 stored
	    # in LA[17] with a second index listing all the nodes in V1 and the value being the aligned node in V2. eg in
	    # the same case as above with "ABC" in V1 and "DEF" in V2, we would have L[17]["ABC"]="DEF"
	    LA[FNR][V1[i]]=V2[i]
	}
	ASSERT(S3(LA[FNR])>='"$s3"', "S3 too low on line "FNR);

	j=FNR; # try matching this latest line with all the previously seen ones
	for(i=1;i<FNR;i++) {
	    li=length(LA[i]); lj=length(LA[j]); s3=S3pair(P,i,j); lp=length(P);
	    # HEURISTIC ALERT: we try to return only "relevant, reasonably large" merges; the current heuristic that
	    # seems to work is: either 90% the size of the biggest so far, otherwise at least 50% bigger than the input sets
	    if(lp>=MIN(biggest*0.9, 1.5*MAX(li,lj)+1) && s3>='"$s3"') { printf "%d,%g\t",lp,s3;
		comma=""; for(v in P){printf "%s%s",comma,  v; comma=","} printf "\t";
		comma=""; for(v in P){printf "%s%s",comma,P[v];comma=","} printf "\n";
	    }
	    if(lp>biggest)biggest=lp;
	}
    } ' "$G1" "$G2" "$A"
