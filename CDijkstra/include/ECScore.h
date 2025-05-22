#ifndef EC_SCORE_H
#define EC_SCORE_H

#include <vector>

int computeEA(
    const std::vector<int>& alignedG1,
    const std::vector<int>& alignedG2,
    const std::vector<std::vector<bool>>& adj1,
    const std::vector<std::vector<bool>>& adj2);

int incrementalEA(
    int newNodeG1,
    int newNodeG2,
    const std::vector<int>& alignedG1,
    const std::vector<int>& alignedG2,
    const std::vector<std::vector<bool>>& adj1,
    const std::vector<std::vector<bool>>& adj2);

int computeE(
    const std::vector<int>& alignedG,
    const std::vector<std::vector<bool>>& adj);

int computeGraphE(
    int numNode,
    const std:: vector<std::vector<bool>>& adj
);
    
int incrementalE(
    int newNodeG,
    const std::vector<int>& alignedG,
    const std::vector<std::vector<bool>>& adj);

double computeEC(
    int EA,
    int E1);

double computeS3(
    int EA, 
    int E1, 
    int E2);

#endif // EC_SCORE_H
