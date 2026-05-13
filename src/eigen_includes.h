// DIFFCP SOURCE: diffcp/cpp/include/eigen_includes.h @ v1.1.8 (6a78143)
// Upstream pin: see inst/UPSTREAM.dcf / upstream_info().
//
// We include Eigen via RcppEigen so that the package picks up the
// vendored Eigen shipped with the RcppEigen R package, instead of
// requiring a system-wide Eigen install.

#pragma once

#include <RcppEigen.h>

using Vector       = Eigen::VectorXd;
using Array        = Eigen::Array<double, Eigen::Dynamic, 1>;
using Matrix       = Eigen::MatrixXd;
using MatrixRef    = Eigen::Ref<Eigen::MatrixXd>;
using SparseMatrix = Eigen::SparseMatrix<double>;
