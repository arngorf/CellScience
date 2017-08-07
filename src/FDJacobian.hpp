#ifndef FDJacobian_HPP
#define FDJacobian_HPP

#include "types.hpp"

Mat FDJacobian(std::function<double(const Vector &)> const& F, Vector const &x, double h = 10e-5);

#endif // FDJacobian_HPP