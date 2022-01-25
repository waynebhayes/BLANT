#ifndef CDIJKSTRA_SIM_HEADER
#define CDIJKSTRA_SIM_HEADER

#include <string>
#include "graph.h"

using std::string;

class SimilarityMatrix {
public:
    SimilarityMatrix(string sim_file, const Graph& graph1, const Graph& graph2);
    ~SimilarityMatrix();

    SimilarityMatrix(const SimilarityMatrix& m) = delete;
    SimilarityMatrix& operator =(SimilarityMatrix& m) = delete;

    double similarity(int n, int m) const;

private:
    unsigned int rows;
    unsigned int cols;
    double** mat;
};

#endif
