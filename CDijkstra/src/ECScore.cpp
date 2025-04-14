#include "ECScore.h"
//induced subgraph
int computeEA(
    const std::vector<int>& alignedG1,
    const std::vector<int>& alignedG2,
    const std::vector<std::vector<bool>>& adj1,
    const std::vector<std::vector<bool>>& adj2)
{
    int EA = 0;
    int n = alignedG1.size();

    for (int i = 0; i < n; ++i) {
        int u1 = alignedG1[i];
        int u2 = alignedG2[i];
        for (int j = i + 1; j < n; ++j) {
            int v1 = alignedG1[j];
            int v2 = alignedG2[j];
            if (adj1[u1][v1] && adj2[u2][v2]) {
                EA = EA + 1;
            }
        }
    }
    return EA;
}

int incrementalEA(
    int newNodeG1,
    int newNodeG2,
    const std::vector<int>& alignedG1,
    const std::vector<int>& alignedG2,
    const std::vector<std::vector<bool>>& adj1,
    const std::vector<std::vector<bool>>& adj2)
{
    int addedEA = 0;
    int n = alignedG1.size();

    for (int i = 0; i < n; ++i) {
        int u1 = alignedG1[i];
        int u2 = alignedG2[i];
        if (adj1[newNodeG1][u1] && adj2[newNodeG2][u2]) {
            ++addedEA;
        }
    }

    return addedEA;
}
//induced graph
int computeE(
    const std::vector<int>& alignedG,
    const std::vector<std::vector<bool>>& adj)
{
    int E = 0;
    int n = alignedG.size();

    for (int i = 0; i < n; ++i) {
        int u1 = alignedG[i];
        for (int j = i + 1; j < n; ++j) {
            int v1 = alignedG[j];
            if (adj[u1][v1] ) {
                ++E;
            }
        }
    }
    return E;
}

int incrementalE(
    int newNodeG,
    const std::vector<int>& alignedG,
    const std::vector<std::vector<bool>>& adj)
{
    int addedE = 0;
    int n = alignedG.size();

    for (int i = 0; i < n; ++i) {
        int u1 = alignedG[i];
        if (adj[newNodeG][u1]) {
            ++addedE;
        }
    }

    return addedE;
}
//EC=EA/E1 or EC = EA/E2
double computeEC(int EA, int E1){
    return (double) EA/E1;
}

double computeS3(int EA, int E1, int E2){
    return (double) EA/(E1+E2-EA);
}