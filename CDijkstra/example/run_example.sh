#!/bin/bash
# Run BLANT-Dijkstra on the A. thaliana – H. sapiens example (BioGRID 3.0.064, importance similarity).
# Usage (from the CDijkstra folder): ./example/run_example.sh
set -e
cd "$(dirname "$0")/.."

./build/cdijkstra 0.05 0.1 \
  example/AThaliana-HSapiens-3.0.064.sim \
  example/AThaliana-HSapiens-3.0.064-seed.txt \
  example/AThaliana-3.0.064.el \
  example/HSapiens-3.0.064.el \
  > example/output.txt 2> example/log.txt

echo "Done: $(wc -l < example/output.txt) aligned pairs written to example/output.txt"
