#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <random>

#include "DistPointHyperellipsoid.hpp"
#include "ellipsoidHelpers.hpp"
#include "EllipsoidLeastSquareFit.hpp"
#include "Hyperellipsoid.hpp"

namespace Random {
    std::default_random_engine generator = std::default_random_engine(3); //
    std::uniform_real_distribution<double> randU = std::uniform_real_distribution<double>(0, 1); //
    std::normal_distribution<double> randN = std::normal_distribution<double>(0, 1); //
}

bool feq(double x, double y, double eps) {
    return (std::abs(x - y) < eps);
}

bool vecsort(int i,int j) {
    return (i<j);
}

inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

Mat getEllipsoidNormal(Vector v, double x, double y, double z, bool xyflat) {
    double A = v(0); double B = v(1); double C = v(2);
    double D = v(3); double E = v(4); double F = v(5);
    double G = v(6); double H = v(7); double I = v(8);
    Mat n(1,3);
    n(0,0) = 2*A*x + 2*D*y + 2*E*z + 2*G;
    n(0,1) = 2*B*y + 2*D*x + 2*F*z + 2*H;
    if (xyflat) {
        n(0,2) = 0.0;
    } else {
        n(0,2) = 2*C*z + 2*E*x + 2*F*y + 2*I;
    }
    double d = sqrt(n(0,0)*n(0,0) + n(0,1)*n(0,1) + n(0,2)*n(0,2));
    n *= 1.0/d;
    return n;
}

std::vector<double> getEllipsoidX(Vector v, double y, double z) {
    double A = v(0); double B = v(1); double C = v(2);
    double D = v(3); double E = v(4); double F = v(5);
    double G = v(6); double H = v(7); double I = v(8);
    double J = v(9);

    std::vector<double> xes;
    double d = std::pow(2*D*y + 2*E*z + 2*G,2) -4*A*(B*y*y + C*z*z + 2*F*y*z + 2*H*y + 2*I*z + J);
    if (d == 0) {
        double x = - (2*D*y + 2*E*z + 2*G)/(2*A);
        xes.push_back(x);
    } else if (d > 0) {
        double x1 = (  sqrt(d) - (2*D*y + 2*E*z + 2*G))/(2*A);
        double x2 = (- sqrt(d) - (2*D*y + 2*E*z + 2*G))/(2*A);
        xes.push_back(x1);
        xes.push_back(x2);
    }
    return xes;
}

double EllipsoidEnergyFromQC(Vector v, double x, double y, double z) {
    double A = v(0); double B = v(1); double C = v(2);
    double D = v(3); double E = v(4); double F = v(5);
    double G = v(6); double H = v(7); double I = v(8);
    double J = v(9);
    return A*x*x + B*y*y + C*z*z
         + 2*D*x*y + 2*E*x*z + 2*F*y*z
         + 2*G*x + 2*H*y + 2*I*z
         + J;
}

Mat GetEulerRotationMatrix(Vector const &eulerAngles, bool transposed, bool reversed) {
    // Assumes row vectors and right side multiplication;
    Mat Rx(3,3); Mat Ry(3,3); Mat Rz(3,3);

    double alpha = eulerAngles(0);
    double beta = eulerAngles(1);
    double gamma = eulerAngles(2);

    Rx << 1,          0,           0,
          0,          cos(alpha),  -sin(alpha),
          0,          sin(alpha),  cos(alpha);

    Ry << cos(beta),  0,           sin(beta),
          0,          1,           0,
          -sin(beta), 0,           cos(beta);

    Rz << cos(gamma), -sin(gamma), 0,
          sin(gamma), cos(gamma),  0,
          0,          0,           1;

    Mat R;
    if (reversed) {
        R = Rx*Ry*Rz;
    } else {
        R = Rz*Ry*Rx;
    }
    if (transposed) {
        R.transposeInPlace();
    }
    return R;
}

Mat RandomRotationMatrix() {
    /*
     * Rework of http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
     * Author: Jim Arvo, 1991
     */
    Mat R(3,3);
    double theta = Random::randU(Random::generator) * PI * 2;
    double phi   = Random::randU(Random::generator) * PI * 2;
    double z     = Random::randU(Random::generator) * 2.0;

    double r  = std::sqrt(z);
    double Vx = std::sin(phi) * r;
    double Vy = std::cos(phi) * r;
    double Vz = std::sqrt(2.0 - z);

    double st = std::sin(theta);
    double ct = std::cos(theta);
    double Sx = Vx * ct - Vy * st;
    double Sy = Vx * st + Vy * ct;

    /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
    /* is equivalent to V S - R.                                        */

    R(0,0) = Vx * Sx - ct;
    R(0,1) = Vx * Sy - st;
    R(0,2) = Vx * Vz;

    R(1,0) = Vy * Sx + st;
    R(1,1) = Vy * Sy - ct;
    R(1,2) = Vy * Vz;

    R(2,0) = Vz * Sx;
    R(2,1) = Vz * Sy;
    R(2,2) = 1.0 - z;

    return R;
}

Mat GetEllipsoidXFromRadii(Vector const &radii, int N) {
    Mat X = Mat::Zero(N,3);
    for (int i = 0; i < N; ++i) {
        double u = 2 * Random::randU(Random::generator) * PI;
        double v = Random::randU(Random::generator) * PI;
        X(i,0) = radii(0) * cos(u) * sin(v);
        X(i,1) = radii(1) * sin(u) * sin(v);
        X(i,2) = radii(2) * cos(v);
    }
    return X;
}

Mat RandomRotatePoints(Mat &X) {
    Mat R = RandomRotationMatrix();
    X *= R;
    return R;
}

void TranslatePoints(Mat &X, Vector const &dx) {
    int N = X.rows();
    for (int i = 0; i < N; ++i) {
        X.row(i) += dx.transpose();
    }
}

Vector GetCenter(Mat const &X) {
    int N = X.rows();
    Vector center = Vector::Zero(3);
    for (int i = 0; i < N; ++i) {
        center += X.row(i).transpose();
    }
    center /= (double)N;
    return center;
}

void CenterAtOrigin(Mat &X) {
    int N = X.rows();
    Vector center = GetCenter(X);
    for (int i = 0; i < N; ++i) {
        X.row(i) = X.row(i) - center.transpose();
    }
}

void AddNormalDistError(Mat &X, Ellipsoid flipsoid, double variance, bool xyflat) {
    int N = X.rows();
    for (int i = 0; i < N; ++i) {
        double x = X(i,0);
        double y = X(i,1);
        double z = X(i,2);
        Mat Normal = getEllipsoidNormal(flipsoid.getAlgebraicCoefficients(), x, y, z, xyflat);
        X.row(i) += Random::randN(Random::generator) * variance * Normal;
    }
}

Mat MakeLayeredX(const Mat &X, Ellipsoid flipsoid, double dz, bool noLayering) {
    // dh is used here to make sure no points end above or below ellipsoid.
    int N = X.rows();
    Mat LayeredX = X;

    double EnergyTolerance = 10e-10;
    double dt = 0.001;

    double maxZ = X.col(2).maxCoeff();
    double minZ = X.col(2).minCoeff();
    double height = maxZ - minZ;

    if (maxZ < minZ) {
        std::cout << "MakeLayeredX::Error: degenerate ellipsoid thickness." << std::endl;
    }

    int M = height/dz;
    int eIndex = 0;
    double errors[M];
    for (int i = 0; i < M; ++i) {
        errors[i] = Random::randN(Random::generator)*dz/2;
        while (errors[i] >= dz/2 || errors[i] <= -dz/2) {
            errors[i] = Random::randN(Random::generator)*dz/2;
        }
    }
    double dzp = height/(M-1);

    Vector v = flipsoid.getAlgebraicCoefficients();

    for (int i = 0; i < N; ++i) {
        if (!noLayering) {
            LayeredX(i,2) -= minZ;
            LayeredX(i,2) /= dzp;
            LayeredX(i,2) += 0.5;
            LayeredX(i,2) = floor(LayeredX(i,2));
            eIndex = LayeredX(i,2);
            LayeredX(i,2) += 0.5;
            LayeredX(i,2) *= dz;
            LayeredX(i,2) += minZ;
        }
        double prevEnergy = EllipsoidEnergyFromQC(v, LayeredX(i,0), LayeredX(i,1), LayeredX(i,2));
        while (std::abs(prevEnergy) > EnergyTolerance) {

            // Get normal, energy and modify point
            Mat Normal = getEllipsoidNormal(v, LayeredX(i,0), LayeredX(i,1), LayeredX(i,2), true);
            double curEnergy = EllipsoidEnergyFromQC(v, LayeredX(i,0), LayeredX(i,1), LayeredX(i,2));
            LayeredX.block(i,0,1,3) -= dt * curEnergy * Normal;
            // If worse, try other way
            double newEnergy = EllipsoidEnergyFromQC(v, LayeredX(i,0), LayeredX(i,1), LayeredX(i,2));
            if (std::abs(newEnergy) > std::abs(prevEnergy)) {
                LayeredX.block(i,0,1,3) += 2 * dt * curEnergy * Normal;
                // Check if better
                newEnergy = EllipsoidEnergyFromQC(v, LayeredX(i,0), LayeredX(i,1), LayeredX(i,2));
                if (std::abs(newEnergy) > std::abs(prevEnergy)) {
                    LayeredX.block(i,0,1,3) -= dt * curEnergy * Normal;
                    std::cout << "Could not put point on ellipsoid. Out of bounds or too large step." << std::endl;
                }
            }

            prevEnergy = newEnergy;
        }

        if (!noLayering) {
            LayeredX(i,2) += errors[eIndex];
        }

    }
    return LayeredX;

}

Ellipsoid GetVerifiedEllipsoid(Mat const &X, Vector const &trueRadii, bool &errorFlag, double eps) {
    double error;
    Ellipsoid flipsoid = EllipsoidLeastSquareFit(X, error, errorFlag);

    if (errorFlag) {
        std::cout << "Error: Fitted ellipsoid radii did not make it to completion!" << std::endl;
    } else {

        Vector fittedRadii = flipsoid.getRadii();
        std::vector<double> radii1(3);
        std::vector<double> radii2(3);
        for (int i = 0; i < 3; ++i){
            radii1[i] = trueRadii(i);
            radii2[i] = fittedRadii(i);
        }

        std::sort(radii1.begin(), radii1.end(), vecsort);
        std::sort(radii2.begin(), radii2.end(), vecsort);

        for (int i = 0; i < 3; ++i) {
            if (!feq(radii1[i], radii2[i], eps)) {
                std::cout << "Error: Fitted ellipsoid radii not verified!" << std::endl;
                errorFlag = true;
                break;
            }
        }
    }
    return flipsoid;
}

Mat GetEllipsoidPointShell(Ellipsoid flipsoid, int N, double dh) {
    Mat X = Mat::Zero(N,3);

    double maxRadii = flipsoid.getRadii().maxCoeff();

    for (int i = 0; i < N; ++i) {
        Vector y(3);
        for (int i = 0; i < 3; ++i) {
            y(i) = flipsoid.getCenter()(i) + (Random::randU(Random::generator) - 0.5) * 4 * maxRadii;
        }
        if (std::abs(EllipsoidEnergyFromQC(flipsoid.getAlgebraicCoefficients(), y(0), y(1), y(2))) < dh) {
            X.row(i) = y;
            ++i;
        }
    }
    return X;
}

Mat loadVesicleFromFolder(std::string folderRootPath) {

    std::string delim = "/";
    std::string txt = ".txt";
    std::string vesicleLayerFileName = "VesicleSyn1sec";
    std::vector<std::string> vesicleFiles;

    for (int i = 150; i < 300; ++i) {
        std::string layerNumber = std::to_string(i);
        std::string completePath = folderRootPath + delim + vesicleLayerFileName + layerNumber + txt;
        if (fileExists(completePath)) {
            vesicleFiles.push_back(completePath);
        }
    }

    std::vector <std::vector <double> > data;
    double currentZ = 0;

    for (int i = 0; i < vesicleFiles.size(); ++i) {

        std::ifstream infile(vesicleFiles[i]);
        if (infile.good()) {

            while (infile) {
                std::string s;
                if (!getline(infile, s)) break;

                std::istringstream ss(s);
                std::vector<double> record;

                while (ss) {
                    std::string s;
                    if (!getline(ss, s, ',')) break;
                    record.push_back(std::stod(s));
                }

                record.push_back(currentZ);

                data.push_back(record);
            }
            if (!infile.eof()) {
                std::cout << "Fooey!\n";
            }
        } else {
            std::cout << "File does not exist" << std::endl;
        }
        currentZ += 1.0;
    }

    int N = data.size();
    Mat X(N ,3);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < 3; ++i) {
            X(j,i) = data[j][i]*5;
        }
    }
    return X;
}

EllipsoidComparisonResult compareEllipsoids(Ellipsoid flipsoidA, Ellipsoid flipsoidB, double eA, double eB) {

    /*
    double bestValue = 360*3;

    // Find best matching vectors
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            for (int k = 0; k < 3; ++k) {
                if (i == k || j == k) continue;
                double value = std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(0).dot( flipsoidB.getRadiiDirection(i))),2),
                                        std::pow(std::acos(flipsoidA.getRadiiDirection(0).dot(-flipsoidB.getRadiiDirection(i))),2))
                             + std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(1).dot( flipsoidB.getRadiiDirection(j))),2),
                                        std::pow(std::acos(flipsoidA.getRadiiDirection(1).dot(-flipsoidB.getRadiiDirection(j))),2))
                             + std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(2).dot( flipsoidB.getRadiiDirection(k))),2),
                                        std::pow(std::acos(flipsoidA.getRadiiDirection(2).dot(-flipsoidB.getRadiiDirection(k))),2));
                if (value < bestValue) {
                    bestValue = value;
                    indices[0] = i;
                    indices[1] = j;
                    indices[2] = k;
                }
            }
        }
    }
    */

    EllipsoidComparisonResult result;

    result.flipsoidA = flipsoidA;
    result.flipsoidB = flipsoidB;

    result.AAverageError = eA;
    result.BAverageError = eB;

    result.ASurfaceArea = EllipsoidArea(flipsoidA);
    result.BSurfaceArea = EllipsoidArea(flipsoidB);

    Vector zUnit(3);
    zUnit << 0,0,1;

    result.aMinAngle = std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(0).dot(flipsoidB.getRadiiDirection(0))),2),
                                std::pow(std::acos(flipsoidA.getRadiiDirection(0).dot(-flipsoidB.getRadiiDirection(0))),2));
    result.bMinAngle = std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(1).dot(flipsoidB.getRadiiDirection(1))),2),
                                std::pow(std::acos(flipsoidA.getRadiiDirection(1).dot(-flipsoidB.getRadiiDirection(1))),2));
    result.cMinAngle = std::min(std::pow(std::acos(flipsoidA.getRadiiDirection(2).dot(flipsoidB.getRadiiDirection(2))),2),
                                std::pow(std::acos(flipsoidA.getRadiiDirection(2).dot(-flipsoidB.getRadiiDirection(2))),2));
    result.minAngleAvg = (result.aMinAngle + result.bMinAngle + result.cMinAngle)/(6*PI)*360;
    result.Da = -flipsoidA.getRadii(0) + flipsoidB.getRadii(0);
    result.Db = -flipsoidA.getRadii(1) + flipsoidB.getRadii(1);
    result.Dc = -flipsoidA.getRadii(2) + flipsoidB.getRadii(2);
    result.Aac = flipsoidA.getRadii(0) / flipsoidA.getRadii(2);
    result.Aab = flipsoidA.getRadii(0) / flipsoidA.getRadii(1);
    result.Abc = flipsoidA.getRadii(1) / flipsoidA.getRadii(2);
    result.Bac = flipsoidB.getRadii(0) / flipsoidB.getRadii(2);
    result.Bab = flipsoidB.getRadii(0) / flipsoidB.getRadii(1);
    result.Bbc = flipsoidB.getRadii(1) / flipsoidB.getRadii(2);
    result.Dac = -result.Aac + result.Bac;
    result.Dab = -result.Aab + result.Bab;
    result.Dbc = -result.Abc + result.Bbc;
    result.Avol = 4.0*PI*flipsoidA.getRadii(0)*flipsoidA.getRadii(1)*flipsoidA.getRadii(2)/3.0;
    result.Bvol = 4.0*PI*flipsoidB.getRadii(0)*flipsoidB.getRadii(1)*flipsoidB.getRadii(2)/3.0;
    result.Aatoz = std::min(std::acos(flipsoidA.getRadiiDirection(0).dot(zUnit))/(2*PI)*360,
                            std::acos(flipsoidA.getRadiiDirection(0).dot(-zUnit))/(2*PI)*360);
    result.Actoz = std::min(std::acos(flipsoidA.getRadiiDirection(2).dot(zUnit))/(2*PI)*360,
                            std::acos(flipsoidA.getRadiiDirection(2).dot(-zUnit))/(2*PI)*360);
    result.Batoz = std::min(std::acos(flipsoidB.getRadiiDirection(0).dot(zUnit))/(2*PI)*360,
                            std::acos(flipsoidB.getRadiiDirection(0).dot(-zUnit))/(2*PI)*360);
    result.Bctoz = std::min(std::acos(flipsoidB.getRadiiDirection(2).dot(zUnit))/(2*PI)*360,
                            std::acos(flipsoidB.getRadiiDirection(2).dot(-zUnit))/(2*PI)*360);

    return result;
}

Vector MatrixToQuat(Mat const &rotate) {
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

Mat QuatToMatrix(Vector const &quaternion) {
    //q = [w,x,y,z] = w + ix + jy + kz
    Mat R(3,3);
    double w = quaternion(0); double x = quaternion(1);
    double y = quaternion(2); double z = quaternion(3);
    double ww = w*w;
    double xx = x*x;
    double yy = y*y;
    double zz = z*z;
    double wx = w*x;
    double wy = w*y;
    double wz = w*z;
    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double invs = 1.0 / (xx + yy + zz + ww);

    R(0,0) = ( xx - yy - zz + ww)*invs;
    R(1,1) = (-xx + yy - zz + ww)*invs;
    R(2,2) = (-xx - yy + zz + ww)*invs;

    R(1,0) = 2.0 * (xy + wz)*invs;
    R(0,1) = 2.0 * (xy - wz)*invs;

    R(2,0) = 2.0 * (xz - wy)*invs;
    R(0,2) = 2.0 * (xz + wy)*invs;

    R(2,1) = 2.0 * (yz + wx)*invs;
    R(1,2) = 2.0 * (yz - wx)*invs;

    return R;
}

double EllipsoidEnergy(Mat const &X, Vector const &radii, Vector const &center, Vector const& rotateQuaternion) {

    int N = X.rows();
    Mat rotate = QuatToMatrix(rotateQuaternion);

    Vector radiiF(3);
    for (int i = 0; i < 3; ++i) {
        radiiF(i) = std::abs(radii(i)) + 10e-10;
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

double EllipsoidArea(Ellipsoid &flipsoid) {
    Vector radii = flipsoid.getRadii();
    double a = radii(0);
    double b = radii(1);
    double c = radii(2);
    return 4.0*PI*std::pow((std::pow(a*b,1.6)+std::pow(a*c,1.6)+std::pow(b*c,1.6))/3.0, 1.0/1.6);
}
