#include "FDHessian.hpp"
#include <iostream>

Mat FDHessian(std::function<double(const Vector &)> const& F, Vector const &x, double h) {
    int n = x.size();
    Mat Hessian(n,n);
    for (int j = 0; j < n; ++j) {
        for (int i = j; i < n; ++i) {
            if (i == j) {
                Vector a = x; Vector b = x;
                a(i) += h;
                b(i) -= h;
                Hessian(i,i) = (F(a) -2*F(x) + F(b)) / (h*h);
            } else {
                Vector a = x; Vector b = x; Vector c = x; Vector d = x;
                a(i) += h; a(j) += h;
                b(i) += h; b(j) -= h;
                c(i) -= h; c(j) += h;
                d(i) -= h; d(j) -= h;
                Hessian(i,j) = (F(a) - F(b) - F(c) + F(d)) / (4*h*h);
                Hessian(j,i) = Hessian(i,j);
            }
        }
    }
    return Hessian;
}

