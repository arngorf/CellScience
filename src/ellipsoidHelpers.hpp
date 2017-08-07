#ifndef ELLIPSOID_HELPERS_HPP
#define ELLIPSOID_HELPERS_HPP

#include "Ellipsoid.hpp"

#include <random>

struct EllipsoidComparisonResult {
    Ellipsoid flipsoidA;
    Ellipsoid flipsoidB;
    double aMinAngle;
    double bMinAngle;
    double cMinAngle;
    double minAngleAvg;
    double Da;
    double Db;
    double Dc;
    double Aac;
    double Aab;
    double Abc;
    double Bac;
    double Bab;
    double Bbc;
    double Dac;
    double Dab;
    double Dbc;
    double Avol;
    double Bvol;
    double Aatoz;
    double Actoz;
    double Batoz;
    double Bctoz;
    double AAverageError;
    double BAverageError;
    double ASurfaceArea;
    double BSurfaceArea;
};

bool feq(double x, double y, double eps = 10e-10);

namespace Random {
    extern std::default_random_engine generator;
    extern std::uniform_real_distribution<double> randU;
    extern std::normal_distribution<double> randN;
}

bool vecsort(int i,int j);

inline bool fileExists(const std::string& name);

Mat getEllipsoidNormal(Vector v, double x, double y, double z, bool xyflat = false);

std::vector<double> getEllipsoidX(Vector v, double y, double z);

double EllipsoidEnergyFromQC(Vector v, double x, double y, double z);

Mat GetEulerRotationMatrix(Vector const &eulerAngles, bool transposed = false, bool reversed = false);

Mat RandomRotationMatrix();

Mat GetEllipsoidXFromRadii(Vector const &radii, int N = 50);

Mat RandomRotatePoints(Mat &X);

void TranslatePoints(Mat &X, Vector const &dx);

Vector GetCenter(Mat const &X);

void CenterAtOrigin(Mat &X);

void AddNormalDistError(Mat &X, Ellipsoid flipsoid, double variance = 1, bool xyflat = true);

Mat MakeLayeredX(const Mat &X, Ellipsoid flipsoid, double dz, bool noLayering = false);

Ellipsoid GetVerifiedEllipsoid(Mat const &X, Vector const &trueRadii, bool &errorFlag, double eps = 10e-5);

Mat GetEllipsoidPointShell(Ellipsoid flipsoid, int N, double dh);

Mat loadVesicleFromFolder(std::string folderRootPath);

EllipsoidComparisonResult compareEllipsoids(Ellipsoid flipsoidA, Ellipsoid flipsoidB, double eA, double eB);

Vector MatrixToQuat(Mat const &rotate);

Mat QuatToMatrix(Vector const &quaternion);

double EllipsoidEnergy(Mat const &X, Vector const &radii, Vector const &center, Vector const& rotateQuaternion);

double EllipsoidArea(Ellipsoid &flipsoid);

#endif // ELLIPSOID_HELPERS_HPP