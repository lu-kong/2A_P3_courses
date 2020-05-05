
#include "KnnClassification.hpp"
#include <iostream>
// #include <fstream>
#include <ANN/ANN.h>

 
KnnClassification::KnnClassification(int k, Dataset* dataset, int col_class)
: Classification(dataset, col_class) {
    // TODO Exercise 1
}

KnnClassification::~KnnClassification() {
    // TODO Exercise 1
}

int KnnClassification::Estimate(const Eigen::VectorXd & x, double threshold) {
    // TODO Exercise 2
}

int KnnClassification::getK() {
    return m_k;
}

ANNkd_tree* KnnClassification::getKdTree() {
    return m_kdTree;
}
