#pragma once
#include <cstdint>

#ifdef _MSC_VER
typedef __int64 int64_t;
#endif

#include <Eigen/Core>


using ScalarArray = Eigen::VectorXd;
using VectorArray = Eigen::MatrixXd;
using IndicesArray = Eigen::MatrixXi;
