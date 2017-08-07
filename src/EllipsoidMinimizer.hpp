#ifndef ELLIPSOID_MINIMIZER_HPP
#define ELLIPSOID_MINIMIZER_HPP

#include "types.hpp"

#include "Ellipsoid.hpp"
#include "Minimizer.hpp"

class EllipsoidMinimizer {
public:

    EllipsoidMinimizer();

    Ellipsoid FitEllipsoid(Mat const &XIn, double &error, bool &errorFlag);

    Ellipsoid FitSphere(Mat const &XIn, double &error, bool &errorFlag);

    Ellipsoid AlgebraicDistanceEllipsoidFit(const Mat &Xin, double &error, bool &errorFlag);

    void GetMinimum(Mat const &X, Vector &radii, Vector &center, Mat &rotate);

    void GetMinimumSphere(Mat const &X, double &radii, Vector &center);

    double getAverageDist() {return latestAvgDist;};

private:

    Minimizer minimizer;

    void MinimizeAll(Mat const &X, Vector &radii, Vector &center, Vector &rotate);

    void MinimizeAllSphere(Mat const &X, double &radii, Vector &center);

    Vector MatrixToAngles(Mat const &rotate);

    Mat AnglesToMatrix(Vector const &angle);

    Vector MatrixToQuat(Mat const &rotate);

    Mat QuatToMatrix(Vector const &quaternion);

    double EllipsoidSSD(Mat const &points, Vector const &radii, Vector const &center, Mat const& rotate);

    double AverageResidule(Mat const &X, Vector const &radii, Vector const &center, Mat const& rotate);

    double latestAvgDist;

};

#endif // ELLIPSOID_MINIMIZER_HPP