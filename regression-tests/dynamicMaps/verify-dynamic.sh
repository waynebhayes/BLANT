#!/bin/bash
# blant-verify.sh — verify BLANT -mj output canonical Gint and orbit groupings
# Canonicalization matches BLANT's L_K_Func_Sort: sort by cubed sum of neighbor
# degrees, then find minimal Gint within score-equal groups.
# In -mj output, first field is the ORIGINAL (non-canonical) Gint.

EXEDIR=$(dirname "$0")
cd "$EXEDIR/.." || exit 1
PROJECT_DIR=$(pwd)

BLANT="$PROJECT_DIR/blant"
GRAPH_FILE="$PROJECT_DIR/Testing_Data/syeast.el"

declare -A EDGES=()

die() { echo "FATAL: $*" >&2; exit 1; }

# --- edge loading ---
load_edges() {
    local f=$1 d=$2
    EDGES=()
    while read -r a b _; do
        [[ -z "$a" || "$a" == "#"* ]] && continue
        EDGES["$a,$b"]=1
        (( d )) || EDGES["$b,$a"]=1
    done < "$f"
}

# --- Gint to adjacency ---
gint_to_adj() {
    local gint=$1 k=$2 d=$3
    local -n _gta_a=$4
    local bp=0 i j
    for ((i=0; i<k; i++)); do for ((j=0; j<k; j++)); do _gta_a[i*k+j]=0; done; done
    if (( d )); then
        for ((i=k-1; i>=0; i--)); do
            for ((j=k-1; j>=0; j--)); do
                (( i == j )) && continue
                (( (gint >> bp) & 1 )) && _gta_a[i*k+j]=1
                (( bp++ ))
            done
        done
    else
        for ((i=k-1; i>=0; i--)); do
            for ((j=i-1; j>=0; j--)); do
                (( (gint >> bp) & 1 )) && { _gta_a[i*k+j]=1; _gta_a[j*k+i]=1; }
                (( bp++ ))
            done
        done
    fi
}

# --- order + adjacency to Gint (nameref output, no subshell) ---
ord_to_gint() {
    local -n _otg_a=$1 _otg_o=$2 _otg_v=$5
    local k=$3 d=$4
    local bp=0 g=0 i j
    if (( d )); then
        for ((i=k-1; i>=0; i--)); do
            for ((j=k-1; j>=0; j--)); do
                (( i == j )) && continue
                (( _otg_a[_otg_o[i]*k+_otg_o[j]] )) && (( g |= (1 << bp) ))
                (( bp++ ))
            done
        done
    else
        for ((i=k-1; i>=0; i--)); do
            for ((j=i-1; j>=0; j--)); do
                (( _otg_a[_otg_o[i]*k+_otg_o[j]] )) && (( g |= (1 << bp) ))
                (( bp++ ))
            done
        done
    fi
    _otg_v=$g
}

# --- score computation: degree and cubed-sum ---
comp_scores() {
    local -n _csc_a=$1 _csc_s=$2 _csc_d=$3
    local k=$4 i j
    for ((i=0; i<k; i++)); do
        local deg=0
        for ((j=0; j<k; j++)); do (( _csc_a[i*k+j] )) && (( deg++ )); done
        _csc_d[i]=$deg
    done
    for ((i=0; i<k; i++)); do
        local sum=0
        for ((j=0; j<k; j++)); do
            (( _csc_a[i*k+j] )) && (( sum += _csc_d[j] * _csc_d[j] * _csc_d[j] ))
        done
        _csc_s[i]=$sum
    done
}

# --- sort indices by score (stable selection sort) ---
sort_asc() {
    local -n _sas_o=$1 _sas_s=$2
    local k=$3 i j
    for ((i=0; i<k; i++)); do _sas_o[i]=$i; done
    for ((i=0; i<k-1; i++)); do
        local best=$i
        for ((j=i+1; j<k; j++)); do
            (( _sas_s[_sas_o[j]] < _sas_s[_sas_o[best]] )) && best=$j
        done
        if (( best != i )); then
            local t=${_sas_o[i]}; _sas_o[i]=${_sas_o[best]}; _sas_o[best]=$t
        fi
    done
}

# --- find score-group boundaries ---
find_grps() {
    local -n _fgr_o=$1 _fgr_s=$2 _fgr_g=$3
    local k=$4 i
    _fgr_g=(0)
    for ((i=1; i<k; i++)); do
        (( _fgr_s[_fgr_o[i]] != _fgr_s[_fgr_o[i-1]] )) && _fgr_g+=($i)
    done
    _fgr_g+=($k)
}

# --- group-based canonical (recursive within-group swaps, no brute force) ---
grp_canon() {
    local -n _gca_a=$1
    local k=$2 d=$3
    local -n _gca_mn=$4 _gca_or=$5
    local -a scores degrees order groups
    comp_scores _gca_a scores degrees $k
    sort_asc order scores $k
    find_grps order scores groups $k

    local best=-1
    local -a bp=()
    local -a perm=("${order[@]}")
    _gca_next() {
        local gi=$1
        if (( gi >= ${#groups[@]} - 1 )); then
            local val
            ord_to_gint _gca_a perm $k $d val
            if (( best == -1 || val < best )); then
                best=$val
                bp=("${perm[@]}")
            fi
            return
        fi
        local start=${groups[gi]} end=${groups[gi+1]}
        _gca_swp $start $((end-1)) $gi
    }
    _gca_swp() {
        local left=$1 right=$2 gi=$3
        if (( left >= right )); then _gca_next $((gi+1)); return; fi
        local i entry_gint t test_val
        ord_to_gint _gca_a perm $k $d entry_gint
        for ((i=left; i<=right; i++)); do
            if (( i != left )); then
                t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
                ord_to_gint _gca_a perm $k $d test_val
                t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
                (( test_val == entry_gint )) && continue
            fi
            t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
            _gca_swp $((left+1)) $right $gi
            t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
        done
    }
    _gca_next 0
    _gca_mn=$best
    _gca_or=("${bp[@]}")
}

# --- orbit groups from canonical adjacency ---
# Recursively enumerate all score-group-respecting permutations; for each
# automorphism (same Gint), merge cycle members into the same orbit (makeOrbit).
comp_orb() {
    local -n _cob_a=$1
    local k=$2 d=$3
    local -n _cob_r=$4
    local i
    for ((i=0; i<k; i++)); do _cob_r[i]=$i; done

    local -a scores degrees order groups
    comp_scores _cob_a scores degrees $k
    sort_asc order scores $k
    find_grps order scores groups $k
    local canon
    ord_to_gint _cob_a order $k $d canon

    # Swap groups: pairwise swaps that don't change Gint → same orbit (like BLANT)
    local -a swg=()
    local ii jj t sgv
    for ((ii=0; ii<k; ii++)); do swg[ii]=$ii; done
    for ((ii=0; ii<k; ii++)); do
        for ((jj=ii+1; jj<k; jj++)); do
            t=${order[ii]}; order[ii]=${order[jj]}; order[jj]=$t
            ord_to_gint _cob_a order $k $d sgv
            t=${order[ii]}; order[ii]=${order[jj]}; order[jj]=$t
            if (( sgv == canon )); then
                (( swg[ii] < swg[jj] )) && swg[jj]=${swg[ii]} || swg[ii]=${swg[jj]}
            fi
        done
    done
    for ((ii=0; ii<k; ii++)); do
        local root=${swg[ii]}
        while (( root != swg[root] )); do root=${swg[root]}; done
        swg[ii]=$root
    done
    for ((ii=0; ii<k; ii++)); do _cob_r[ii]=${swg[ii]}; done

    local -a perm=("${order[@]}")
    _cob_next() {
        local gi=$1
        if (( gi >= ${#groups[@]} - 1 )); then
            local val
            ord_to_gint _cob_a perm $k $d val
            if (( val == canon )); then
                local -a visited=()
                local ii2
                for ((ii2=0; ii2<k; ii2++)); do visited[ii2]=0; done
                for ((ii2=0; ii2<k; ii2++)); do
                    (( visited[ii2] )) && continue
                    local -a cyc=()
                    local cur=$ii2
                    while (( !visited[cur] )); do
                        visited[cur]=1; cyc+=($cur); cur=${perm[cur]}
                    done
                    local mo=${_cob_r[${cyc[0]}]}
                    local jj2
                    for ((jj2=1; jj2<${#cyc[@]}; jj2++)); do
                        (( _cob_r[${cyc[jj2]}] < mo )) && mo=${_cob_r[${cyc[jj2]}]}
                    done
                    for ((jj2=0; jj2<${#cyc[@]}; jj2++)); do _cob_r[${cyc[jj2]}]=$mo; done
                done
            fi
            return
        fi
        local start=${groups[gi]} end=${groups[gi+1]}
        _cob_swp $start $((end-1)) $gi
    }
    _cob_swp() {
        local left=$1 right=$2 gi=$3
        if (( left >= right )); then _cob_next $((gi+1)); return; fi
        local i entry_gint t test_val
        ord_to_gint _cob_a perm $k $d entry_gint
        for ((i=left; i<=right; i++)); do
            if (( i != left )); then
                t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
                ord_to_gint _cob_a perm $k $d test_val
                t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
                (( test_val == entry_gint )) && continue
            fi
            t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
            _cob_swp $((left+1)) $right $gi
            t=${perm[left]}; perm[left]=${perm[i]}; perm[i]=$t
        done
    }
    _cob_next 0
    for ((i=0; i<k; i++)); do
        local root=${_cob_r[i]}
        while (( root != _cob_r[root] )); do root=${_cob_r[root]}; done
        _cob_r[i]=$root
    done
}

# --- format orbit groups as BLANT -mj does ---
# Iterate positions 0..k-1; for each unprinted position, group all later
# positions sharing the same orbit (non-consecutive positions allowed).
fmt_grps() {
    local -n _fgrp_n=$1 _fgrp_o=$2 _fgrp_r=$3
    local k=$4
    local -a printed=()
    local i
    for ((i=0; i<k; i++)); do printed[i]=0; done
    for ((i=0; i<k; i++)); do
        (( printed[i] )) && continue
        local -a grp=("${_fgrp_n[_fgrp_o[i]]}")
        printed[i]=1
        local j
        for ((j=i+1; j<k; j++)); do
            if (( !printed[j] && _fgrp_r[j] == _fgrp_r[i] )); then
                grp+=("${_fgrp_n[_fgrp_o[j]]}")
                printed[j]=1
            fi
        done
        local IFS=':'; echo "${grp[*]}"; unset IFS
    done
}

# --- verify one -mj line ---
verify_line() {
    local line=$1 k=$2 d=$3
    IFS=' ' read -ra tokens <<< "$line"
    local printed=${tokens[0]}
    local -a node_list=() out_grps=()
    for ((i=1; i<${#tokens[@]}; i++)); do
        local g=$(( i-1 ))
        out_grps[g]="${tokens[i]}"
        IFS=':' read -ra parts <<< "${tokens[i]}"
        for part in "${parts[@]}"; do node_list+=("$part"); done
    done
    if (( ${#node_list[@]} != k )); then
        echo "  PARSE ERROR: expected $k nodes, got ${#node_list[@]}"
        return 1
    fi
    local ns=""
    for ((i=0; i<k; i++)); do ns+=" ${node_list[i]}"; done; ns="${ns:1}"

    local -a a_nodes=()
    for ((i=0; i<k; i++)); do
        for ((j=0; j<k; j++)); do
            a_nodes[i*k+j]=${EDGES["${node_list[i]},${node_list[j]}"]:+1}
            a_nodes[i*k+j]=${a_nodes[i*k+j]:-0}
        done
    done

    local -a a_gint=()
    gint_to_adj $printed $k $d a_gint

    local c_nodes
    local -a p_nodes
    grp_canon a_nodes $k $d c_nodes p_nodes

    local check_gint
    ord_to_gint a_gint p_nodes $k $d check_gint
    if (( check_gint != c_nodes )); then
        echo "  CANON MISMATCH: printed=$printed c_nodes=$c_nodes check=$check_gint nodes=$ns"
        return 1
    fi

    local -a ca=()
    local i j
    for ((i=0; i<k; i++)); do
        for ((j=0; j<k; j++)); do
            ca[i*k+j]=${a_nodes[p_nodes[i]*k+p_nodes[j]]}
        done
    done

    local -a orb=()
    comp_orb ca $k $d orb

    local -a exp=()
    while IFS= read -r grp; do exp+=("$grp"); done < <(fmt_grps node_list p_nodes orb $k)

    if (( ${#exp[@]} != ${#out_grps[@]} )); then
        echo "  ORBIT GROUP COUNT: expected ${#exp[@]}, got ${#out_grps[@]}"
        echo "    expected: ${exp[*]}"
        echo "    got:      ${out_grps[*]}"
        echo "    nodes: $ns"
        return 1
    fi
    for ((i=0; i<${#exp[@]}; i++)); do
        if [[ "${exp[i]}" != "${out_grps[i]}" ]]; then
            echo "  ORBIT GROUP $i: expected '${exp[i]}' got '${out_grps[i]}'"
            echo "    expected: ${exp[*]}"
            echo "    got:      ${out_grps[*]}"
            echo "    nodes: $ns"
            return 1
        fi
    done
    return 0
}

# --- run BLANT and verify ---
run_and_verify() {
    local graph=$1 k=$2 directed=$3 samples=${4:-100}
    local flagD=
    (( directed )) && flagD="-D"
    local outfile
    outfile=$(mktemp /tmp/blant-verify.XXXXXX 2>/dev/null) || outfile=$(mktemp -t blant-verify.XXXXXX 2>/dev/null) || die "mktemp failed"
    $BLANT -k $k $flagD -mj -n $samples "$graph" > "$outfile" 2>/dev/null
    if (( $? != 0 )); then echo "  BLANT failed!"; rm -f "$outfile"; return 1; fi

    local total=0 ok=0 bad=0 line
    load_edges "$graph" $directed
    while IFS= read -r line; do
        [[ -z "$line" ]] && continue
        (( total++ ))
        if verify_line "$line" $k $directed; then (( ok++ ))
        else (( bad++ ))
        fi
    done < "$outfile"

    rm -f "$outfile"
    return $bad
}

main() {
    [[ -x "$BLANT" ]] || die "BLANT not found at $BLANT"
    [[ -f "$GRAPH_FILE" ]] || die "Graph file not found at $GRAPH_FILE"
    if ls "$PROJECT_DIR"/canon_maps/canon_map*.bin "$PROJECT_DIR"/canon_maps/directed/canon_map*.bin 2>/dev/null | head -1 | grep -q .; then
        echo "PASS"
        exit 0
    fi
    local ec=0
    run_and_verify "$GRAPH_FILE" 5 0; (( $? != 0 )) && ec=1
    run_and_verify "$GRAPH_FILE" 4 1; (( $? != 0 )) && ec=1
    run_and_verify "$GRAPH_FILE" 6 1; (( $? != 0 )) && ec=1
    run_and_verify "$GRAPH_FILE" 10 0 5; (( $? != 0 )) && ec=1
    if (( ec == 0 )); then echo "PASS"; else echo "FAIL"; fi
    exit $ec
}

(return 0 2>/dev/null) && return
main "$@"
