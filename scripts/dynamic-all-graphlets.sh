#!/bin/bash
# verify-dynamic.sh — verify BLANT dynamic-map -mi output against a canon list
# Builds a chained graph containing all connected canonicals, samples BLANT,
# then checks that the number of unique sampled graphlet types matches the
# number of connected canonicals in the list.
# Usage: verify-dynamic.sh K canon_list.txt [directed]

set -euo pipefail

usage() { echo "Usage: $(basename "$0") K canon_list.txt [directed]"; exit 1; }
[[ $# -ge 2 ]] || usage

K=$1
CANON_LIST=$2
DIRECTED=${3:-0}

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/.." || exit 1
PROJECT_DIR=$(pwd)
BLANT="$PROJECT_DIR/blant"

die() { echo "FAIL: $*" >&2; exit 1; }

# ========== main ==========
main() {
    [[ -x "$BLANT" ]] || die "BLANT not found at $BLANT (build it first?)"
    [[ -f "$CANON_LIST" ]] || die "canon list not found: $CANON_LIST"
    echo "# $0: testing BLANT -mi output against $CANON_LIST" >&2

    # Parse canon list — count connected graphlets
    declare -A GINT_TO_ORD
    local -a CONNECTED_EDGES=()
    local ord=0

    while IFS=$'\t' read -r -a fields; do
        local gint="${fields[0]}"
        local rest="${fields[1]}"
        read -r -a parts <<< "$rest"
        local is_grp="${parts[0]}"
        GINT_TO_ORD[$gint]=$ord
        if [[ "$is_grp" == "1" ]]; then
            CONNECTED_EDGES+=("${fields[2]}")  # edges in "u,v u,v ..." format
        fi
        ((++ord))
    done < <(tail -n +2 "$CANON_LIST")

    local NUM_CON=${#CONNECTED_EDGES[@]}
    echo "# $NUM_CON connected canonicals out of $ord total" >&2

    # Shuffle (pure bash, no seq/shuf dependency)
    local -a SHUF
    for ((i=0; i<NUM_CON; i++)); do SHUF[$i]=$i; done
    for ((i=NUM_CON-1; i>0; i--)); do
        local j=$((RANDOM % (i+1)))
        local t=${SHUF[$i]}; SHUF[$i]=${SHUF[$j]}; SHUF[$j]=$t
    done

    # Write chained graph — one redirect for the entire loop
    local EL_FILE="$PROJECT_DIR/chained_k${K}.el"
    local gi idx edges off u v pct_chain=0
    for ((gi=0; gi<NUM_CON; gi++)); do
        idx=${SHUF[$gi]}
        edges="${CONNECTED_EDGES[$idx]}"
        off=$((gi * (K-1)))
        for edge in $edges; do
            u="${edge%,*}"; v="${edge#*,}"
            echo "$((u + off)) $((v + off))"
        done
        local new_pct=$(( (gi + 1) * 100 / NUM_CON / 10 * 10 ))
        if (( new_pct > pct_chain )); then
            pct_chain=$new_pct
            echo "# chaining graphlets: ${pct_chain}%" >&2
        fi
    done > "$EL_FILE"

    local NUM_EDGES
    NUM_EDGES=$(wc -l < "$EL_FILE")
    echo "# chained graph: $NUM_CON graphlets, $NUM_EDGES edges" >&2

    # Sample BLANT in batches
    local BATCH_SIZE=10000
    local MAX_BATCHES=$(( NUM_CON / 10 ))
    (( MAX_BATCHES < 1 )) && MAX_BATCHES=1
    local DIR_FLAG=""
    if (( DIRECTED )); then DIR_FLAG="-D"; fi

    local OUT UNIQUE_FILE
    OUT=$(mktemp) || die "mktemp failed"
    UNIQUE_FILE=$(mktemp) || die "mktemp failed"

    local batch NUM_UNIQUE
    for ((batch=1; batch<=MAX_BATCHES; batch++)); do
        $BLANT -k $K $DIR_FLAG -mi -n $BATCH_SIZE "$EL_FILE" > "$OUT"

        # Extract unique Gints from this batch, merge with accumulated
        cut -d' ' -f1 "$OUT" | sort -u > "${UNIQUE_FILE}.new"
        if [[ -s "$UNIQUE_FILE" ]]; then
            sort -u -m "$UNIQUE_FILE" "${UNIQUE_FILE}.new" > "${UNIQUE_FILE}.tmp"
            mv "${UNIQUE_FILE}.tmp" "$UNIQUE_FILE"
        else
            mv "${UNIQUE_FILE}.new" "$UNIQUE_FILE"
        fi
        rm -f "${UNIQUE_FILE}.new" "${UNIQUE_FILE}.tmp"

        NUM_UNIQUE=$(wc -l < "$UNIQUE_FILE")
        echo "# batch $batch: $NUM_UNIQUE unique graphlets so far" >&2

        if (( NUM_UNIQUE >= NUM_CON )); then
            rm -f "$OUT" "$UNIQUE_FILE" "$EL_FILE"
            echo "PASS ($NUM_UNIQUE/$NUM_CON found after $batch batches, $((batch * BATCH_SIZE)) total samples)"
            exit 0
        fi
    done

    rm -f "$OUT" "$UNIQUE_FILE" "$EL_FILE"
    echo "FAIL ($NUM_UNIQUE/$NUM_CON found after $MAX_BATCHES batches; expected all $NUM_CON)" >&2
    exit 1
}

main "$@"
