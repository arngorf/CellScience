#ifndef FD_GRADIENT_HPP
#define FD_GRADIENT_HPP

#include "types.hpp"

Vector FDGradient(std::function<double(const Vector &)> const& F, Vector const &x, double h = 10e-10);

#endif // FD_GRADIENT_HPP