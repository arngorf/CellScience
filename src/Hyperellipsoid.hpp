// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#ifndef HYPEREELIPSOID_HPP
#define HYPEREELIPSOID_HPP

// A hyperellipsoid has center K; axis directions U[0] through U[N-1], all
// unit-length vectors; and extents e[0] through e[N-1], all positive numbers.
// A point X = K + sum_{d=0}^{N-1} y[d]*U[d] is on the hyperellipsoid whenever
// sum_{d=0}^{N-1} (y[d]/e[d])^2 = 1.  An algebraic representation for the
// hyperellipsoid is (X-K)^T * M * (X-K) = 1, where M is the NxN symmetric
// matrix M = sum_{d=0}^{N-1} U[d]*U[d]^T/e[d]^2, where the superscript T
// denotes transpose.  Observe that U[i]*U[i]^T is a matrix, not a scalar dot
// product.  The hyperellipsoid is also represented by a quadratic equation
// 0 = C + B^T*X + X^T*A*X, where C is a scalar, B is an Nx1 vector, and A is
// an NxN symmetric matrix with positive eigenvalues.  The coefficients can be
// stored from lowest degree to highest degree,
//   C = k[0]
//   B = k[1], ..., k[N]
//   A = k[N+1], ..., k[(N+1)(N+2)/2 - 1]
// where the A-coefficients are the upper-triangular elements of A listed in
// row-major order.  For N = 2, X = (x[0],x[1]) and
//   0 = k[0] +
//       k[1]*x[0] + k[2]*x[1] +
//       k[3]*x[0]*x[0] + k[4]*x[0]*x[1]
//                      + k[5]*x[1]*x[1]
// For N = 3, X = (x[0],x[1],x[2]) and
//   0 = k[0] +
//       k[1]*x[0] + k[2]*x[1] + k[3]*x[2] +
//       k[4]*x[0]*x[0] + k[5]*x[0]*x[1] + k[6]*x[0]*x[2] +
//                      + k[7]*x[1]*x[1] + k[8]*x[1]*x[2] +
//                                       + k[9]*x[2]*x[2]
// This equation can be factored to the form (X-K)^T * M * (X-K) = 1, where
// K = -A^{-1}*B/2, M = A/(B^T*A^{-1}*B/4-C).

#include "types.hpp"

template <int N>
class Hyperellipsoid {
public:

    // Construction and destruction.  The default constructor sets the center
    // to Vector<N,double>::Zero(), the axes to Vector<N,double>::Unit(d), and all
    // extents to 1.
    Hyperellipsoid();
    Hyperellipsoid(Vector const &inCenter, Mat const &inAxis, Vector const &inExtent);

    // Compute M = sum_{d=0}^{N-1} U[d]*U[d]^T/e[d]^2.
    void GetM(Mat &M) const;

    // Compute M^{-1} = sum_{d=0}^{N-1} U[d]*U[d]^T*e[d]^2.
    void GetMInverse(Mat &MInverse) const;

    // Public member access.
    Vector center;
    Mat axis;
    Vector extent;

    // Comparisons to support sorted containers.
    bool operator==(Hyperellipsoid const& hyperellipsoid) const;
    bool operator!=(Hyperellipsoid const& hyperellipsoid) const;
    /*bool operator< (Hyperellipsoid const& hyperellipsoid) const;
    bool operator<=(Hyperellipsoid const& hyperellipsoid) const;
    bool operator> (Hyperellipsoid const& hyperellipsoid) const;
    bool operator>=(Hyperellipsoid const& hyperellipsoid) const;*/
};

// Template aliases for convenience.
using Ellipse2 = Hyperellipsoid<2>;

using Ellipsoid3 = Hyperellipsoid<3>;


template <int N>
Hyperellipsoid<N>::Hyperellipsoid() {
    center = Mat::Zero(N, N);
    for (int d = 0; d < N; ++d) {
        Vector newVector = Vector::Zero(N);
        newVector(d) = 1;
        axis.row(d) = newVector;
        extent(d) = 1.0;
    }
}

template <int N>
Hyperellipsoid<N>::Hyperellipsoid(Vector const &inCenter, Mat const &inAxis, Vector const &inExtent)
    : center(inCenter), axis(inAxis), extent(inExtent) {

}

template <int N>
void Hyperellipsoid<N>::GetM(Mat &M) const {
    M = Mat::Zero(N, N);
    for (int d = 0; d < N; ++d) {
        Vector ratio = axis.row(d) / extent(d);
        M += ratio * ratio.transpose();
    }
}

template <int N>
void Hyperellipsoid<N>::GetMInverse(Mat &MInverse) const {
    MInverse = Mat::Zero(N, N);
    for (int d = 0; d < N; ++d) {
        Vector product = axis.row(d) * extent[d];
        MInverse += product * product.transpose();
    }
}

template <int N>
bool Hyperellipsoid<N>::operator==(Hyperellipsoid const& hyperellipsoid) const {
    return center == hyperellipsoid.center && axis == hyperellipsoid.axis
        && extent == hyperellipsoid.extent;
}

template <int N>
bool Hyperellipsoid<N>::operator!=(Hyperellipsoid const& hyperellipsoid) const {
    return !operator==(hyperellipsoid);
}

/*template <int N>
bool Hyperellipsoid<N>::operator<(Hyperellipsoid const& hyperellipsoid) const {
    if (LexLess(center, hyperellipsoid.center)) {
        return true;
    }

    if (LexGreat(center, hyperellipsoid.center)) {
        return false;
    }

    if (LexLess(axis, hyperellipsoid.axis)) {
        return true;
    }

    if (LexGreat(axis, hyperellipsoid.axis)) {
        return false;
    }

    return LexLess(extent, hyperellipsoid.extent);
}

template <int N>
bool Hyperellipsoid<N>::operator<=(Hyperellipsoid const& hyperellipsoid) const {
    return operator<(hyperellipsoid) || operator==(hyperellipsoid);
}

template <int N>
bool Hyperellipsoid<N>::operator>(Hyperellipsoid const& hyperellipsoid) const {
    return !operator<=(hyperellipsoid);
}

template <int N>
bool Hyperellipsoid<N>::operator>=(Hyperellipsoid const& hyperellipsoid) const {
    return !operator<(hyperellipsoid);
}

bool LexLess(Vector v1, Vector v2) {
    for (int i = 0; i < v1.size(); ++i) {
        if (v1(i) >= v2(i)) {
            return false;
        }
    }
    return true;
}

bool LexGreat(Vector v1, Vector v2) {
    return LexLess(v2, v1);
}*/

#endif // HYPEREELIPSOID_HPP