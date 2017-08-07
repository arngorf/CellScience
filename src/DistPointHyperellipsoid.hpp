#ifndef DIST_POINT_HYPERELLIPSOID_HPP
#define DIST_POINT_HYPERELLIPSOID_HPP

// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#include <algorithm>
#include <exception>
#include <iostream>
#include "types.hpp"
#include "Hyperellipsoid.hpp"


// Compute the distance from a point to a hyperellipsoid.  In 2D, this is a
// point-ellipse distance query.  In 3D, this is a point-ellipsoid distance
// query.  The following document describes the algorithm.
//   http://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
// The hyperellipsoid can have arbitrary center and orientation; that is, it
// does not have to be axis-aligned with center at the origin.
//
// For the 2D query,
//   Vector2<Real> point;  // initialized to something
//   Ellipse2<Real> ellipse;  // initialized to something
//   DCPPoint2Ellipse2<Real> query;
//   auto result = query(point, ellipse);
//   Real distance = result.distance;
//   Vector2<Real> closestEllipsePoint = result.closest;
//
// For the 3D query,
//   Vector3<Real> point;  // initialized to something
//   Ellipsoid3<Real> ellipsoid;  // initialized to something
//   DCPPoint3Ellipsoid3<Real> query;
//   auto result = query(point, ellipsoid);
//   Real distance = result.distance;
//   Vector3<Real> closestEllipsoidPoint = result.closest;

// Distance and closest-point queries.
template <typename Type0, typename Type1>
class DCPQuery {

public:
    struct Result {
        // A DCPQuery-base class B must define a B::Result struct with member
        // 'Real distance'.  A DCPQuery-derived class D must also derive a
        // D::Result from B:Result but may have no members.  The idea is to
        // allow Result to store closest-point information in addition to the
        // distance.  The operator() is non-const to allow DCPQuery to store
        // and modify private state that supports the query.
    };
    Result operator()(Type0 const& primitive0, Type1 const& primitive1);
};


template <int N>
class DCPQuery<Vector, Hyperellipsoid<N>> {

public:
    struct Result
    {
        double distance, sqrDistance;
        Vector closest;
    };

    // The query for any hyperellipsoid.
    Result operator()(Vector const &point, Hyperellipsoid<N> const &hyperellipsoid);

    // The 'hyperellipsoid' is assumed to be axis-aligned and centered at the
    // origin , so only the extent[] values are used.
    Result operator()(Vector const &point, Vector const &extent);

private:
    // The hyperellipsoid is sum_{d=0}^{N-1} (x[d]/e[d])^2 = 1 with no
    // constraints on the orderind of the e[d].  The query point is
    // (y[0],...,y[N-1]) with no constraints on the signs of the components.
    // The function returns the squared distance from the query point to the
    // hyperellipsoid.   It also computes the hyperellipsoid point
    // (x[0],...,x[N-1]) that is closest to (y[0],...,y[N-1]).
    double SqrDistance(Vector const &e, Vector const &y, Vector &x);

    // The hyperellipsoid is sum_{d=0}^{N-1} (x[d]/e[d])^2 = 1 with the e[d]
    // positive and nonincreasing:  e[d] >= e[d + 1] for all d.  The query
    // point is (y[0],...,y[N-1]) with y[d] >= 0 for all d.  The function
    // returns the squared distance from the query point to the
    // hyperellipsoid.  It also computes the hyperellipsoid point
    // (x[0],...,x[N-1]) that is closest to (y[0],...,y[N-1]), where
    // x[d] >= 0 for all d.
    double SqrDistanceSpecial(Vector const &e, Vector const &y, Vector &x);

    // The bisection algorithm to find the unique root of F(t).
    double Bisector(int numComponents, Vector const &e, Vector const &y, Vector &x);

    double LengthRobust(Vector const &v);
};

// Template aliases for convenience.
template <int N>
using DCPPointHyperellipsoid =
DCPQuery<Vector, Hyperellipsoid<N>>;

using DCPPoint2Ellipse2 = DCPPointHyperellipsoid<2>;

using DCPPoint3Ellipsoid3 = DCPPointHyperellipsoid<3>;

template <int N>
typename DCPQuery<Vector, Hyperellipsoid<N>>::Result
DCPQuery<Vector, Hyperellipsoid<N>>::operator()(Vector const &point, Hyperellipsoid<N> const &hyperellipsoid) {
    for (int i = 0; i < N; ++i) {
        if (hyperellipsoid.extent(i) <= 0) {
            std::cout << "extents = " << hyperellipsoid.extent.transpose() << std::endl;
            throw std::invalid_argument("Hyperellipsoid extents must be positive.");
        }
    }

    Result result;

    // Compute the coordinates of Y in the hyperellipsoid coordinate system.
    Vector diff = point - hyperellipsoid.center;
    //std::cout << "diff = " << diff.transpose() << std::endl;
    //std::cout << "point = " << point.transpose() << std::endl;
    //std::cout << "h.center = " << hyperellipsoid.center.transpose() << std::endl;
    Vector y(N);
    for (int i = 0; i < N; ++i) {
        y(i) = diff.transpose() * hyperellipsoid.axis.row(i).transpose();
    }
    //std::cout << "h.extent = " << hyperellipsoid.extent.transpose() << std::endl;

    // Compute the closest hyperellipsoid point in the axis-aligned
    // coordinate system.
    Vector x = Vector::Zero(N);
    result.sqrDistance = SqrDistance(hyperellipsoid.extent, y, x);
    //std::cout << "result.sqrDistance = " << hyperellipsoid.extent << ", " << y << ", " << x << std::endl;
    result.distance = std::sqrt(result.sqrDistance);

    // Convert back to the original coordinate system.
    result.closest = hyperellipsoid.center;
    for (int i = 0; i < N; ++i) {
        result.closest += x(i) * hyperellipsoid.axis.row(i);
    }

    return result;
}

template <int N>
typename DCPQuery<Vector, Hyperellipsoid<N>>::Result
DCPQuery<Vector, Hyperellipsoid<N>>::operator()(Vector const &point, Vector const &extent) {
    for (int i = 0; i < N; ++i) {
        if (extent(i) <= 0) {
            throw std::invalid_argument("Hyperellipsoid extents must be positive.");
        }
    }
    Result result;
    result.sqrDistance = SqrDistance(extent, point, result.closest);
    result.distance = std::sqrt(result.sqrDistance);
    return result;
}

template <int N>
double DCPQuery<Vector, Hyperellipsoid<N>>::SqrDistance(Vector const &e, Vector const &y, Vector &x) {
    //std::cout << "pre pre e = " << e.transpose() << std::endl;
    // Determine negations for y to the first octant.
    //std::cout << "e = " << e.transpose() << std::endl;
    //std::cout << "y = " << y.transpose() << std::endl;
    std::array<bool, N> negate;
    for (int i = 0; i < N; ++i) {
        negate[i] = (y(i) < 0.0);
    }

    // Determine the axis order for decreasing extents.
    std::array<std::pair<double, int>, N> permute;
    for (int i = 0; i < N; ++i) {
        permute[i].first = -e(i);
        permute[i].second = i;
    }
    std::sort(permute.begin(), permute.end());

    std::array<int, N> invPermute;
    for (int i = 0; i < N; ++i) {
        invPermute[permute[i].second] = i;
    }

    Vector locE = Vector::Zero(N);
    Vector locY = Vector::Zero(N);
    for (int i = 0; i < N; ++i) {
        int j = permute[i].second;
        locE(i) = e(j);
        locY(i) = std::abs(y(j));
    }

    Vector locX = Vector::Zero(N);
    double sqrDistance = SqrDistanceSpecial(locE, locY, locX);

    // Restore the axis order and reflections.
    for (int i = 0; i < N; ++i) {
        int j = invPermute[i];
        if (negate[i]) {
            locX(j) = -locX(j);
        }
        x(i) = locX(j);
    }

    return sqrDistance;
}

template <int N>
double DCPQuery<Vector, Hyperellipsoid<N>>::
SqrDistanceSpecial(Vector const &e, Vector const &y, Vector& x)
{
    //std::cout << "e = " << e.transpose() << std::endl;
    //std::cout << "y = " << y.transpose() << std::endl;
    double sqrDistance = 0.0;
    Vector ePos = Vector::Zero(N);
    Vector yPos = Vector::Zero(N);
    Vector xPos = Vector::Zero(N);
    int numPos = 0;
    for (int i = 0; i < N; ++i) {
        if (y(i) > 0.0) {
            ePos(numPos) = e(i);
            yPos(numPos) = y(i);
            ++numPos;
        }
        else {
            x(i) = 0.0;
        }
    }
    if (y(N - 1) > 0.0) {
        //std::cout << "aaa bisector" << std::endl;
        sqrDistance = Bisector(numPos, ePos, yPos, xPos);
    } else {
        // y[N-1] = 0

        Vector numer = Vector::Zero(N - 1);
        Vector denom = Vector::Zero(N - 1);
        double eNm1Sqr = e(N - 1) * e(N - 1);
        for (int i = 0; i < numPos; ++i) {
            numer(i) = ePos(i) * yPos(i);
            denom(i) = ePos(i) * ePos(i) - eNm1Sqr;
        }

        bool inSubHyperbox = true;
        for (int i = 0; i < numPos; ++i) {
            if (numer(i) >= denom(i)) {
                inSubHyperbox = false;
                break;
            }
        }

        bool inSubHyperellipsoid = false;
        if (inSubHyperbox) {
            // yPos[] is inside the axis-aligned bounding box of the
            // subhyperellipsoid.  This intermediate test is designed to guard
            // against the division by zero when ePos[i] == e[N-1] for some i.
            Vector xde = Vector::Zero(N - 1);
            double discr = 1.0;
            for (int i = 0; i < numPos; ++i) {
                if (denom(i) == 0) {
                    throw std::invalid_argument("E: DistPointHyperellipsoid::SqrDistanceSpecial: division by zero error");
                }
                xde(i) = numer(i) / denom(i);
                discr -= xde(i) * xde(i);
            }

            if (discr > 0.0) {
                // yPos[] is inside the subhyperellipsoid.  The closest
                // hyperellipsoid point has x[N-1] > 0.
                sqrDistance = 0.0;
                for (int i = 0; i < numPos; ++i) {
                    xPos(i) = ePos(i) * xde(i);
                    double diff = xPos(i) - yPos(i);
                    sqrDistance += diff * diff;
                }
                x(N - 1) = e(N - 1) * std::sqrt(discr);
                sqrDistance += x(N - 1) * x(N - 1);
                inSubHyperellipsoid = true;
            }
        }
        if (!inSubHyperellipsoid) {
            // yPos[] is outside the subhyperellipsoid.  The closest
            // hyperellipsoid point has x[N-1] == 0 and is on the
            // domain-boundary hyperellipsoid.
            x(N - 1) = 0.0;
            //std::cout << "bbb bisector" << std::endl;
            sqrDistance = Bisector(numPos, ePos, yPos, xPos);

        }
    }
    // Fill in those x[] values that were not zeroed out initially.
    for (int i = 0, numPos = 0; i < N; ++i) {
        if (y(i) > 0.0) {
            x(i) = xPos(numPos);
            ++numPos;
        }
    }

    return sqrDistance;
}

template <int N>
double DCPQuery<Vector, Hyperellipsoid<N>>::Bisector(int numComponents, Vector const &e, Vector const &y, Vector &x) {
    //std::cout << "begin" << std::endl;
    //std::cout << "numComponents = " << numComponents << std::endl;
    //std::cout << "e = " << e.transpose() << std::endl;
    //std::cout << "y = " << y.transpose() << std::endl;
    //std::cout << "x = " << x.transpose() << std::endl;
    Vector z = Vector::Zero(N);
    double sumZSqr = 0.0;
    int i;
    for (i = 0; i < numComponents; ++i) {
        if (e(i) == 0) {
            throw std::invalid_argument("E: DistPointHyperellipsoid::Bisector: division by zero error (e(i))");
        }
        z(i) = y(i) / e(i);
        //std::cout << "z = y / e :: " << z(i) << " = " << y(i) << " / " << e(i) << std::endl;
        sumZSqr += z(i) * z(i);
    }
    //std::cout << "z = " << z.transpose() << std::endl;
    if (sumZSqr == 1.0) {
        // The point is on the hyperellipsoid.
        for (i = 0; i < numComponents; ++i) {
            x(i) = y(i);
        }
        return 0.0;
    }

    double emin = e(numComponents - 1);
    Vector pSqr = Vector::Zero(N);
    Vector numerator = Vector::Zero(N);
    for (i = 0; i < numComponents; ++i) {
        if (emin == 0) {
            throw std::invalid_argument("E: DistPointHyperellipsoid::Bisector: division by zero error (emin)");
        }
        double p = e(i) / emin;
        pSqr(i) = p * p;
        numerator(i) = pSqr(i) * z(i);
        //std::cout << "num = pSqr * z :: " << numerator(i) << " = " << pSqr(i) << " * " << z(i) << std::endl;
    }

    double s = 0.0;
    double smin = z(numComponents - 1) - 1.0;
    double smax;
    if (sumZSqr < 1.0) {
        // The point is strictly inside the hyperellipsoid.
        smax = 0.0;
    } else {
        // The point is strictly outside the hyperellipsoid.
        smax = LengthRobust(numerator) - 1.0;
    }
    //std::cout << "sumZSqr = " << sumZSqr << std::endl;
    //std::cout << "smin = " << smin << std::endl;
    //std::cout << "smax = " << smax << std::endl;

    // The use of 'double' is intentional in case Real is a BSNumber or
    // BSRational type.  We want the bisections to terminate in a reasonable
    // amount of time.
    unsigned int const jmax = (unsigned int)(3 + std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent);
    for (unsigned int j = 0; j < jmax; ++j) {
        s = (smin + smax) * 0.5;
        if (s == smin || s == smax) {
            break;
        }

        double g = -1.0;
        for (i = 0; i < numComponents; ++i) {
            if ((s + pSqr(i)) == 0) {
                throw std::invalid_argument("E: DistPointHyperellipsoid::Bisector: division by zero error ((s + pSqr(i))) (1)");
            }
            double ratio = numerator(i) / (s + pSqr(i));
            g += ratio * ratio;
        }

        if (g > 0.0) {
            smin = s;
        } else if (g < 0.0) {
            smax = s;
        }
        else {
            break;
        }
    }

    double sqrDistance = 0.0;
    //std::cout << "1) x = " << x.transpose() << std::endl;
    //std::cout << "1) y = " << y.transpose() << std::endl;
    //std::cout << "1) pSqr = " << pSqr.transpose() << std::endl;
    //std::cout << "1) smin = " << smin << std::endl;
    //std::cout << "1) smax = " << smax << std::endl;
    //std::cout << "1) s = " << s << std::endl;
    for (i = 0; i < numComponents; ++i) {
        if ((s + pSqr(i)) == 0) {
            throw std::invalid_argument("E: Bisector: division by zero error ((s + pSqr(i))) (2)");
        }
        x(i) = pSqr(i) * y(i) / (s + pSqr(i));
        double diff = x(i) - y(i);
        sqrDistance += diff * diff;
    }
    //std::cout << "2) x = " << x.transpose() << std::endl;
    return sqrDistance;
}

template <int N>
double DCPQuery<Vector, Hyperellipsoid<N>>::LengthRobust(Vector const &v) {
    double maxAbsComp = std::abs(v(0));
    for (int i = 1; i < v.size(); ++i) {
        double absComp = std::abs(v(i));
        if (absComp > maxAbsComp) {
            maxAbsComp = absComp;
        }
    }

    double length;
    if (maxAbsComp > 0.0) {
        if (maxAbsComp == 0) {
            throw std::invalid_argument("E: DistPointHyperellipsoid::Bisector: division by zero error (maxAbsComp)");
        }
        Vector scaled = v / maxAbsComp;
        length = maxAbsComp*scaled.norm();
    } else {
        length = 0.0;
    }
    return length;
};

#endif // DIST_POINT_HYPERELLIPSOID_HPP