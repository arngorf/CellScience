#include "EllipsoidMinimizer.hpp"

#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "DistPointHyperellipsoid.hpp"
#include "ellipsoidHelpers.hpp"
#include "Rotation.hpp"

EllipsoidMinimizer::EllipsoidMinimizer() {

}

Ellipsoid EllipsoidMinimizer::FitEllipsoid(Mat const &XIn, double &error, bool &errorFlag) {

    Mat X = XIn;

    // Store current point-cloud center and translate points to have center in the origin
    Vector centerCorrection = GetCenter(X);
    CenterAtOrigin(X);

    // Approximate the maximum point spread for radii initial value
    double maxd = 0;
    for (int i = 0; i < 3; ++i){
        maxd = std::max(maxd, X.col(i).maxCoeff()-X.col(i).minCoeff());
    }

    Vector centerBest(3);
    Vector radiiBest(3);
    Mat rotateBest = Mat::Identity(3,3);
    double bestAvgResidule = std::numeric_limits<double>::max();

    for (int i = 0; i < 25; ++i) {
        Vector center(3);
        Vector radii(3);
        if (i == 0) {
            radii << maxd/2, maxd/2, maxd/2;
        } else {
            radii << maxd * (Random::randU(Random::generator)-0.5),
                     maxd * (Random::randU(Random::generator)-0.5),
                     maxd * (Random::randU(Random::generator)-0.5);
        }
        center << 0,0,0;
        Mat rotate = RandomRotationMatrix();
        // Calculate minimum
        GetMinimum(X, radii, center, rotate);
        // Setup initial values
        if (latestAvgDist < bestAvgResidule) {
            bestAvgResidule = latestAvgDist;
            centerBest = center;
            radiiBest = radii;
            rotateBest = rotate;
        }
    }

    // Minimize

    error = EllipsoidSSD(X, radiiBest, centerBest, rotateBest);
    errorFlag = false;

    return Ellipsoid(centerBest + centerCorrection, radiiBest, rotateBest);
}

Ellipsoid EllipsoidMinimizer::FitSphere(Mat const &XIn, double &error, bool &errorFlag) {

    Mat X = XIn;

    // Store current point-cloud center and translate points to have center at the origin
    Vector centerCorrection = GetCenter(X);
    CenterAtOrigin(X);

    // Approximate the maximum point spread for radii initial value
    double maxd = 0;
    for (int i = 0; i < 3; ++i){
        maxd = std::max(maxd, X.col(i).maxCoeff()-X.col(i).minCoeff());
    }

    Vector centerBest(3);
    double radiiBest = 0;
    Mat rotateBest = Mat::Identity(3,3);
    double bestAvgResidule = std::numeric_limits<double>::max();

    for (int i = 0; i < 25; ++i) {
        Vector center(3);
        double radii;
        radii = maxd/2;
        center << 0,0,0;

        // Calculate minimum
        GetMinimumSphere(X, radii, center);
        // Setup initial values
        if (latestAvgDist < bestAvgResidule) {
            bestAvgResidule = latestAvgDist;
            centerBest = center;
            radiiBest = radii;
        }
    }

    // Minimize
    Vector radiiAsVector(3);
    radiiAsVector << radiiBest, radiiBest, radiiBest;
    error = EllipsoidSSD(X, radiiAsVector, centerBest, rotateBest);
    errorFlag = false;

    return Ellipsoid(centerBest + centerCorrection, radiiAsVector, rotateBest);
}

void EllipsoidMinimizer::GetMinimum(Mat const &X, Vector &radii, Vector &center, Mat &rotateMat) {

    Vector rotateQuat = MatrixToQuat(rotateMat);

    MinimizeAll(X, radii, center, rotateQuat);

    for (int i = 0; i < 3; ++i) {
        radii(i) = std::abs(radii(i));
    }

    rotateMat = QuatToMatrix(rotateQuat);


    latestAvgDist = AverageResidule(X, radii, center, rotateMat);

    //std::cout << "rotateQuat = " << rotateQuat.transpose() << std::endl;
    //std::cout << "length rotatequat = " << rotateQuat.squaredNorm() << std::endl;
}

void EllipsoidMinimizer::GetMinimumSphere(Mat const &X, double &radii, Vector &center) {

    MinimizeAllSphere(X, radii, center);


    radii = std::abs(radii);

    Vector radiiAsVector(3);
    radiiAsVector << radii, radii, radii;
    latestAvgDist = AverageResidule(X, radiiAsVector, center, Mat::Identity(3,3));

    //std::cout << "rotateQuat = " << rotateQuat.transpose() << std::endl;
    //std::cout << "length rotatequat = " << rotateQuat.squaredNorm() << std::endl;
}

void EllipsoidMinimizer::MinimizeAll(Mat const &X, Vector &radii, Vector &center, Vector &rotate) {
    // input = radii U center U rotate
    std::function<double(Vector const &)> energy =
        [&, X](Vector const &params)
    {
        Vector radiiF = params.segment(0,3);
        Vector centerF = params.segment(3,3);
        Vector rotateF = params.segment(6,4);
        double lambda = 0.0;
        double e1;

        e1 = EllipsoidSSD(X, radiiF, centerF, QuatToMatrix(rotateF));

        double e2 = lambda * std::pow(rotateF.squaredNorm() - 1, 2);

        return e1 + e2;
    };
    Vector params(10);

    params.segment(0,3) = radii;
    params.segment(3,3) = center;
    params.segment(6,4) = rotate;

    params = minimizer.Minimize(energy, params);

    radii = params.segment(0,3);
    center = params.segment(3,3);
    rotate = params.segment(6,4);
}

void EllipsoidMinimizer::MinimizeAllSphere(Mat const &X, double &radii, Vector &center) {
    // input = radii U center U rotate
    std::function<double(Vector const &)> energy =
        [&, X](Vector const &params)
    {
        double radiiF = params(0);
        Vector centerF = params.segment(1,3);
        double e1;
        Vector radiiAsVectorF(3);
        radiiAsVectorF << radiiF, radiiF, radiiF;
        e1 = EllipsoidSSD(X, radiiAsVectorF, centerF, Mat::Identity(3,3));

        return e1;
    };
    Vector params(4);

    params(0) = radii;
    params.segment(1,3) = center;

    params = minimizer.Minimize(energy, params);

    radii = params(0);
    center = params.segment(1,3);
}

Vector EllipsoidMinimizer::MatrixToAngles(Mat const &rotate) {

    AxisAngle<3> aa = Rotation<3>(rotate);
    Vector angle(3);
    if (-1.0 < aa.axis(2)) {
        if (aa.axis(2) < 1.0) {
            angle(0) = std::atan2(aa.axis(1), aa.axis(0));
            angle(1) = std::acos(aa.axis(2));
        } else {
            angle(0) = 0.0;
            angle(1) = 0.0;
        }
    } else {
        angle(0) = 0.0;
        angle(1) = PI;
    }
    return angle;
}

Mat EllipsoidMinimizer::AnglesToMatrix(Vector const &angle) {

    Mat rotate(3,3);
    double cs0 = std::cos(angle(0));
    double sn0 = std::sin(angle(0));
    double cs1 = std::cos(angle(1));
    double sn1 = std::sin(angle(1));
    AxisAngle<3> aa;
    aa.axis(0) = cs0*sn1;
    aa.axis(1) = sn0*sn1;
    aa.axis(2) = cs1;
    aa.angle = angle(2);
    rotate = Rotation<3>(aa);
    return rotate;
}

Vector EllipsoidMinimizer::MatrixToQuat(Mat const &rotate) {
    Vector quaternion(4);
    double trace = rotate(0,0) + rotate(1,1) + rotate(2,2);
    if(trace > 0) {
        double s = 0.5 / std::sqrt(trace + 1.0);
        quaternion(0) = 0.25 / s;
        quaternion(1) = ( rotate(2,1) - rotate(1,2)) * s;
        quaternion(2) = ( rotate(0,2) - rotate(2,0)) * s;
        quaternion(3) = ( rotate(1,0) - rotate(0,1)) * s;
    } else {
        if (rotate(0,0) > rotate(1,1) && rotate(0,0) > rotate(2,2)) {
            double s = 2.0 * std::sqrt(1.0 + rotate(0,0) - rotate(1,1) - rotate(2,2));
            quaternion(0) = (rotate(2,1) - rotate(1,2)) / s;
            quaternion(1) = 0.25 * s;
            quaternion(2) = (rotate(0,1) + rotate(1,0)) / s;
            quaternion(3) = (rotate(0,2) + rotate(2,0)) / s;
        } else if (rotate(1,1) > rotate(2,2)) {
            double s = 2.0 * std::sqrt(1.0 + rotate(1,1) - rotate(0,0) - rotate(2,2));
            quaternion(0) = (rotate(0,2) - rotate(2,0)) / s;
            quaternion(1) = (rotate(0,1) + rotate(1,0)) / s;
            quaternion(2) = 0.25 * s;
            quaternion(3) = (rotate(1,2) + rotate(2,1)) / s;
        } else {
            double s = 2.0 * std::sqrt(1.0 + rotate(2,2) - rotate(0,0) - rotate(1,1));
            quaternion(0) = (rotate(1,0) - rotate(0,1)) / s;
            quaternion(1) = (rotate(0,2) + rotate(2,0)) / s;
            quaternion(2) = (rotate(1,2) + rotate(2,1)) / s;
            quaternion(3) = 0.25 * s;
        }
    }
    return quaternion;
}

Mat EllipsoidMinimizer::QuatToMatrix(Vector const &quaternion) {
    //q = [w,x,y,z] = w + ix + jy + kz
    Mat R(3,3);
    double w = quaternion(0); double x = quaternion(1);
    double y = quaternion(2); double z = quaternion(3);
    double sw = w*w;
    double sx = x*x;
    double sy = y*y;
    double sz = z*z;

    double invs = 1.0 / (sx + sy + sz + sw);
    R(0,0) = ( sx - sy - sz + sw)*invs;
    R(1,1) = (-sx + sy - sz + sw)*invs;
    R(2,2) = (-sx - sy + sz + sw)*invs;

    double tmp1 = x*y;
    double tmp2 = z*w;
    R(1,0) = 2.0 * (tmp1 + tmp2)*invs;
    R(0,1) = 2.0 * (tmp1 - tmp2)*invs;

    tmp1 = x*z;
    tmp2 = y*w;

    R(2,0) = 2.0 * (tmp1 - tmp2)*invs;
    R(0,2) = 2.0 * (tmp1 + tmp2)*invs;
    tmp1 = y*z;
    tmp2 = x*w;
    R(2,1) = 2.0 * (tmp1 + tmp2)*invs;
    R(1,2) = 2.0 * (tmp1 - tmp2)*invs;
    return R;
}

Ellipsoid EllipsoidMinimizer::AlgebraicDistanceEllipsoidFit(const Mat &Xin, double &error, bool &errorFlag) {
    int N = Xin.rows();
    Mat X = Xin;
    Vector centerCorrection = GetCenter(X);
    CenterAtOrigin(X);

    Vector dx(3);
    dx << 30,0,0;
    TranslatePoints(X, dx);

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

    Mat evecs(3,3);
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            evecs(i,j) = abs(eSolve.eigenvectors().col(j)(i));
        }
    }

    error = EllipsoidSSD(X, radii, center, evecs.transpose());
    latestAvgDist = AverageResidule(X, radii, center, evecs.transpose());

    return Ellipsoid(center + centerCorrection, radii, evecs, v);
}

double EllipsoidMinimizer::EllipsoidSSD(Mat const &X, Vector const &radii, Vector const &center, Mat const& rotate) {

    int N = X.rows();

    Vector radiiF(3);
    for (int i = 0; i < 3; ++i) {
        radiiF(i) = std::abs(radii(i));
    }

    //double maxValue = rotate.maxCoeff();
    //double invMax = 1.0 / maxValue;

    //radii *= invMax;
    Ellipsoid3 ellipsoid(center, rotate, radiiF);

    DCPQuery<Vector, Ellipsoid3> peQuery;

    double SSD = 0.0;
    for (int i = 0; i < N; ++i) {
        //Vector diff = X.row(i).transpose();
        //Vector prod = invMax * diff.transpose();
        double dist = peQuery(X.row(i).transpose(), ellipsoid).distance;
        SSD += std::pow(/*maxValue **/ dist, 2);
        if (dist != dist) {
            throw std::domain_error("Distance was NaN!");
        }
    }

    return SSD;
}

double EllipsoidMinimizer::AverageResidule(Mat const &X, Vector const &radii, Vector const &center, Mat const& rotate) {

    int N = X.rows();

    Ellipsoid3 ellipsoid(center, rotate, radii);
    DCPQuery<Vector, Ellipsoid3> peQuery;

    double val = 0;
    for (int i = 0; i < N; ++i) {
        double dist = peQuery(X.row(i).transpose(), ellipsoid).distance;
        val += dist;
    }
    return val / (double)N;
}