#include "Minimizer.hpp"

#include <Eigen/LU>
#include "FDGradient.hpp"
#include "FDHessian.hpp"
#include <iostream>


Minimizer::Minimizer(double tolerance, int maxIter, double dxConst, double dxMax) : tolerance(tolerance), maxIter(maxIter), dxConst(dxConst), dxMax(dxMax) {

}

Vector Minimizer::Minimize(std::function<double(const Vector &)> const& F, Vector const &xInit) {
    return SteepestDecent(F, xInit);
    //return Newtons(F, xInit);
    //return Dogleg(F, xInit);
}

Vector Minimizer::SteepestDecent(std::function<double(const Vector &)> const& F, Vector const &xInit) {
    Vector xMin = xInit;

    int i = 0;
    double alpha;
    while (true) {

        Vector pk = -FDGradient(F, xMin);

        //if (pk.norm() < tolerance) {
        //    break;
        //}
        alpha = backtracking_linesearch(F, -pk, xMin, pk, 0.1);
        if (alpha * pk.norm() < tolerance) {
            //std::cout << "alpha * pk.norm() < tolerance" << std::endl;
            break;
        }
        //std::cout << "alpha = " << alpha << std::endl;
        //std::cout << "pk = " << pk.transpose() << std::endl;
        xMin = xMin + alpha*pk;
        ++i;
        if (i >= maxIter) {
            //std::cout << "i >= maxIter" << std::endl;
            break;
        }
    }
    return xMin;
}

Vector Minimizer::Newtons(std::function<double(const Vector &)> const& F, Vector const &xInit) {
    Vector xMin = xInit;

    int i = 0;
    double alpha;
    while (true) {
        Vector grad = FDGradient(F, xMin);
        Mat Bk = FDHessian(F, xMin);
        Vector pk = -(Bk.inverse() * grad);

        /*if (pk.norm() < tolerance) {
            std::cout << "pk.norm() < tolerance" << std::endl;
            break;
        }*/
        alpha = 0.1;// backtracking_linesearch(F, grad, xMin, pk, 0.5);
        if (alpha * pk.norm() < tolerance) {
            std::cout << "alpha * pk.norm() < tolerance" << std::endl;
            break;
        }

        xMin = xMin + alpha*pk;
        ++i;
        if (i >= maxIter) {
            std::cout << "i >= maxIter" << std::endl;
            break;
        }
    }
    return xMin;
}

Vector Minimizer::Dogleg(std::function<double(const Vector &)> const& F, Vector const &xInit) {
    int n = xInit.size();
    Vector xMin = xInit;
    std::cout << "xMin (0): " << xMin.transpose() << std::endl;

    double dx = dxConst;
    double eta = 0.15;

    Mat Bk = Mat::Identity(n,n);
    Mat Hk = Mat::Identity(n,n);
    Vector xMinPrev = xMin;
    int i = 0;
    bool step;
    double rk;

    while (true) {
        Mat J = FDGradient(F, xMin).transpose();
        Vector gk = FDGradient(F, xMin);

        Vector pU;
        Vector pB = -Hk*gk;
        Vector pk;

        step = false;

        if (pB.norm() < dx) {
            std::cout << "Case 1: (pB.norm() < dx)" << std::endl;
            pk = pB;
        } else {
            std::cout << "Case 2: (pB.norm() >= dx)" << std::endl;
            //std::cout << "Here start" << std::endl;
            //std::cout << "gk = " << gk.transpose() << std::endl;
            //std::cout << "Bk = " << std::endl << Bk << std::endl;
            Vector pU = gk;
            pU *= -(gk.transpose() * gk)/(gk.transpose() * Bk * gk);
            //std::cout << "pU = " << pU << std::endl;
            //exit(0);
            double dpU = pU.norm();
            if (dpU >= dx) {
                std::cout << "Case 2.1: (dpU >= dx)" << std::endl;
                pk = pU / dpU * dx;
            } else {
                std::cout << "Case 2.2: (dpU < dx)" << std::endl;
                Vector pD = pB - pU;
                double a = pD.transpose() * pD;
                double b = pU.transpose() * pD;
                double c = pU.transpose() * pU - std::pow(dx, 2);
                double disc = b*b - 4 * a * c;
                double t1 = (-b + std::sqrt(disc)) / (2*a);
                double t2 = (-b - std::sqrt(disc)) / (2*a);
                double t = std::max(t1, t2);
                //std::cout << "pU: " << pU.transpose() << std::endl;
                if (0 <= t && t <= 1) {
                    std::cout << "Case 2.2.1: (0 <= " << t << " <= 1)" << std::endl;
                    pk = t * pU;
                } else if (1 <= t && t <= 2) {
                    std::cout << "Case 2.2.2: (1 <= " << t << " <= 2)" << std::endl;
                    pk = pU + (t - 1) * pD;
                } else {
                    std::cout << "ERROR: tau out of bounds" << std::endl;
                    break;
                }
            }
        }
        double prev_rk = rk;
        if (pk.norm() < 10e-10) {
            std::cout << "STOPPED: (pk.norm() < 10e-10)" << std::endl;
            break;
        }
        rk = dogleg_rk(F, xMin, gk, pk, Bk);

        if (rk < 0.25) {
            std::cout << "dx changed from " << dx << " to " << 0.25*dx << std::endl;
            dx = 0.25*dx;
        } else {
            if (rk > 0.75) {
                std::cout << "dx changed from " << dx << " to " << std::min(2 * dx, dxMax) << std::endl;
                dx = std::min(2 * dx, dxMax);
            }
        }
        //std::cout << "dx = " << dx << " rk = " << rk << " eta = " << eta << std::endl;
        //std::cout << "xMin (1): " << xMin.transpose() << std::endl;
        if (rk > eta) {
            Vector xMinTemp = xMinPrev;
            xMinPrev = xMin;
            xMin = xMin + pk;
            std::cout << "pk = " << pk.transpose() << std::endl;
            if ((xMin - xMinPrev).norm() < 10e-5) {
                //std::cout << "STOPPED: ((xMin - xMinPrev).norm() < 10e-5)" << std::endl;
                //break;
            }
            //std::cout << "pk (2): " << pk.transpose() << std::endl;
            //std::cout << "xMin (2): " << xMin.transpose() << std::endl;
            step = true;
            if (xMinPrev == xMin) {
                xMinPrev = xMinTemp;
            }
        }

        // Update hessian approximation

        if (step) {
            Vector sk = xMin - xMinPrev;

            Vector gkNew   = FDGradient(F, xMin);
            Vector yk = gkNew - gk;

            double r = 1.0/(yk.transpose() * sk + tolerance);

            //Bk = (Mat::Identity(n,n) - r*yk*sk.transpose())*Hk*(Mat::Identity(n,n) - r*sk*yk.transpose()) + r*yk*yk.transpose();
            Bk = (Mat::Identity(n,n) - r*yk*sk.transpose())*Bk*(Mat::Identity(n,n) - r*sk*yk.transpose()) + r*yk*yk.transpose();
            Hk = (Mat::Identity(n,n) - r*sk*yk.transpose())*Hk*(Mat::Identity(n,n) - r*yk*sk.transpose()) + r*sk*sk.transpose();
        }

        i += 1;
        std::cout << "i = " << i << " XMin = " << xMin.transpose() << " rk = " << rk << std::endl;
        if (i >= maxIter) {
            std::cout << "STOPPED: (i >= maxIter)" << std::endl;
            break;
        }

        if (i != 0 && prev_rk == rk && (xMinPrev - xMin).norm() < 10e-5) {
            std::cout << "STOPPED: (i != 0 && prev_rk == rk && (xMinPrev - xMin).norm() < 10e-5)" << std::endl;
            break;
        }
    }
    return xMin;
}

double Minimizer::backtracking_linesearch(std::function<double(const Vector &)> const& F, Vector const &grad_k, Vector const &xk, Vector const &pk, double c) {
    double alpha = 1;
    double p = 0.1;
    bool run = true;
    while (run) {
        double Fstep;
        try {
            Fstep = F(xk + alpha * pk);
        } catch (const std::invalid_argument &e) {
            alpha = alpha*p;
            continue;
        }
        if (Fstep > F(xk) + c*alpha*grad_k.transpose()*pk) {
            alpha = alpha*p;
        } else {
            run = false;
        }
    }
    return alpha;
}

inline double Minimizer::dogleg_mk(double fk, Vector gk, Vector pk, Mat Bk) {
    return fk + gk.transpose() * pk + 0.5 * pk.transpose() * Bk * pk;
}

inline double Minimizer::dogleg_rk(std::function<double(const Vector &)> const& F, Vector xMin, Vector gk, Vector pk, Mat Bk) {
    double fk = F(xMin);
    double fkp = F(xMin + pk);
    double mk = dogleg_mk(fk, gk, pk, Bk);
    return (fk - fkp)/(fk - mk);
}