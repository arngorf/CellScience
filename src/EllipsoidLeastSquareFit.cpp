#include "EllipsoidLeastSquareFit.hpp"

#include <exception>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "DistPointHyperellipsoid.hpp"

double EllipsoidSSD(Mat const &X, Vector const &radii, Vector const &center, Mat const &evecs) {

    int N = X.rows();
    Mat rotate = evecs.transpose();

    Vector radiiF(3);
    for (int i = 0; i < 3; ++i) {
        radiiF(i) = std::abs(radii(i));
    }

    Ellipsoid3 ellipsoid(center, rotate, radiiF);

    DCPQuery<Vector, Ellipsoid3> peQuery;

    double SSD = 0.0;
    for (int i = 0; i < N; ++i) {
        double dist = peQuery(X.row(i).transpose(), ellipsoid).distance;
        SSD += std::pow(/*maxValue **/ dist, 2);
        if (dist != dist) {
            throw std::domain_error("Distance was NaN!");
        }
    }

    return SSD;
}

Ellipsoid EllipsoidLeastSquareFit(const Mat &Xin, double &error, bool &errorFlag) {
    Mat X = Xin;
    int N = X.rows();

    /*Mat mu = X.colwise().mean();

    for (int i = 0; i < N; ++i) {
        X.row(i) -= mu;
    }*/

    if (N < 9) {
        std::cout << "Too few points - given " << N << " of 9" << std::endl;
        throw std::invalid_argument("Too few points - given " + std::to_string(N) + " of 9");
    }

    errorFlag = false;

    Array x = X.block(0,0,N,1);
    Array y = X.block(0,1,N,1);
    Array z = X.block(0,2,N,1);

    Mat D = Mat(N,9);

    D.block(0,0,N,1) = x*x + y*y - 2*z*z;
    D.block(0,1,N,1) = x*x + z*z - 2*y*y;
    D.block(0,2,N,1) = 2*x*y;
    D.block(0,3,N,1) = 2*x*z;
    D.block(0,4,N,1) = 2*y*z;
    D.block(0,5,N,1) = 2*x;
    D.block(0,6,N,1) = 2*y;
    D.block(0,7,N,1) = 2*z;
    D.block(0,8,N,1) = Mat::Ones(N,1);

    Vector d = x*x + y*y + z*z;
    Vector u = (D.transpose()*D).colPivHouseholderQr().solve(D.transpose()*d);

    Vector v = Vector::Zero(10);
    v(0) = u(0) +   u(1) - 1;
    v(1) = u(0) - 2*u(1) - 1;
    v(2) = u(1) - 2*u(0) - 1;
    v(3) = u(2); v(4) = u(3); v(5) = u(4); v(6) = u(5); v(7) = u(6);
    v(8) = u(7); v(9) = u(8);

    if (v(0) + v(1) + v(2) < 0) {
        v = -v;
    }

    Mat A(4,4);
    A << v(0), v(3), v(4), v(6),
         v(3), v(1), v(5), v(7),
         v(4), v(5), v(2), v(8),
         v(6), v(7), v(8), v(9);

    Vector center = (-A.block(0,0,3,3)).colPivHouseholderQr().solve(A.block(0,3,3,1));

    Mat T = Mat::Identity(4,4);
    T(3,0) = center(0); T(3,1) = center(1); T(3,2) = center(2);

    Mat R = T*A*T.transpose();
    R /= (-R(3,3));

    Eigen::EigenSolver<Mat> eSolve;
    //eSolve.compute(R.block(0,0,3,3) / (-R(3,3)));
    Mat H = R.block(0,0,3,3);
    eSolve.compute(H);
    Vector radii(3);
    radii(0) = sqrt(1.0 / eSolve.eigenvalues()(0).real());
    radii(1) = sqrt(1.0 / eSolve.eigenvalues()(1).real());
    radii(2) = sqrt(1.0 / eSolve.eigenvalues()(2).real());
    if (eSolve.eigenvalues()(0).real() <= 0 ||
        eSolve.eigenvalues()(1).real() <= 0 ||
        eSolve.eigenvalues()(2).real() <= 0) {
        std::cout << "Error: This is not an ellipsoid we have fittet. Eigenvalues:" << std::endl;
        std::cout << v(0) << "x^2 + " << v(1) << "y^2 + " << v(2) << "z^2 + " << v(3) << "xy + "
                  << v(4) << "xz + " << v(5) << "yz + "  << v(6) << "x + " << v(7) << "y + "
                  << v(8) << "z + " << v(9) << " = 0" << std::endl;
        errorFlag = true;
    }

    /*center(0) += mu(0);
    center(1) += mu(1);
    center(2) += mu(2);

    double cHc = center.transpose() * H * center;

    Vector G = -H * center;

    v(0) = H(0,0);
    v(1) = H(1,1);
    v(2) = H(2,2);
    v(3) = H(0,1);
    v(4) = H(0,2);
    v(5) = H(1,2);
    v(6) = G(0);
    v(7) = G(1);
    v(8) = G(2);
    v(9) = cHc - 1;

    std::cout << H << std::endl;
    std::cout << v << std::endl;*/
    //exit(0);

    Mat evecs(3,3);
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            evecs(i,j) = abs(eSolve.eigenvectors().col(j)(i));
        }
    }

    error = EllipsoidSSD(X, radii, center, evecs);

    return Ellipsoid(center, radii, evecs, v);
}

