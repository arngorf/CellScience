#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include "types.hpp"

class Minimizer {
public:

    Minimizer(double tolerance = 10e-10, int maxIter = 1000, double dxConst = 1.0, double dxMax = 100.0);

    Vector Minimize(std::function<double(const Vector &)> const& F, Vector const &xInit);

    Vector SteepestDecent(std::function<double(const Vector &)> const& F, Vector const &xInit);

    Vector Newtons(std::function<double(const Vector &)> const& F, Vector const &xInit);

    Vector Dogleg(std::function<double(const Vector &)> const& F, Vector const &xInit);

    void setTolerance(double newTolerance) {tolerance = newTolerance;};

    void setMaxIter(double newMaxIter) {maxIter = newMaxIter;};

private:

    double tolerance;
    int maxIter;
    double dxConst;
    double dxMax;

    double backtracking_linesearch(std::function<double(const Vector &)> const& F, Vector const &grad_k, Vector const &xk, Vector const &pk, double c = 0.5);
    inline double dogleg_mk(double fk, Vector gk, Vector pk, Mat Bk);
    inline double dogleg_rk(std::function<double(const Vector &)> const& F, Vector xMin, Vector gk, Vector pk, Mat Bk);
};



#endif // MINIMIZER_HPP