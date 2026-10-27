#!/bin/bash
# Builds and tests the C++ port of scripts/blant-clusters.sh (BlantClusters/).
# It needs only a C++17 compiler, not BLANT itself: the e2e test uses canned
# BLANT output.
cd BlantClusters && make all test
