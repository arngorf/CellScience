#include "FDJacobian.hpp"
#include "FDGradient.hpp"

Mat FDJacobian(std::function<double(const Vector &)> const& F, Vector const &x, double h) {
    return FDGradient(F, x, h).transpose();
}
/*
Mat FDJacobian(int dimensions, std::function<Vector(const Vector &)> const& F, Vector const &x, double h) {
    return FDGradient(F, x, h).transpose();
}*/