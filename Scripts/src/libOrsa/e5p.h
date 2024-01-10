#ifndef E5P_H
#define E5P_H

#include "eigen/Eigen/Eigen"

Eigen::Matrix<double, 10, 20> buildConstraintMatrix(const Eigen::Matrix<double, 1, 4> nullSpace[3][3]);

#endif
