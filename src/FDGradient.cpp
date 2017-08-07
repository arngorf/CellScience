#include "FDGradient.hpp"

Vector FDGradient(std::function<double(const Vector &)> const& F, Vector const &x, double h) {
    int n = x.size();
    Vector gradient(n);
    for (int i = 0; i < n; ++i) {
        Vector dx = x;
        dx(i) += h;
        gradient(i) = (F(dx) - F(x))/h;
    }
    return gradient;
}

