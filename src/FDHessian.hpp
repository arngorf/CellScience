#ifndef FD_HESSIAN_HPP
#define FD_HESSIAN_HPP

#include "types.hpp"

Mat FDHessian(std::function<double(const Vector &)> const& F, Vector const &x, double h = 10e-5);

#endif // FD_HESSIAN_HPP