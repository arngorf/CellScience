#ifndef ROTATION_HPP
#define ROTATION_HPP

// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

// Conversions among various representations of rotations.  The value of
// N must be 3 or 4.  The latter case supports affine algebra when you use
// 4-tuple vectors (w-component is 1 for points and 0 for vector) and 4x4
// matrices for affine transformations.  Rotation axes must be unit length.
// The angles are in radians.  The Euler angles are in world coordinates;
// we have not yet added support for body coordinates.

#include "types.hpp"
#include <iostream> // debug

template <int N>
class AxisAngle
{
public:
    AxisAngle();
    AxisAngle(Vector const &inAxis, double inAngle);

    Vector axis;
    double angle;
};


template <int N>
AxisAngle<N>::AxisAngle()
{
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");
    axis = Vector::Zero(3);
}

template <int N>
AxisAngle<N>::AxisAngle(Vector const &inAxis, double inAngle) : axis(inAxis), angle(inAngle) {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");
}


template <int N>
class Rotation {
public:
    // Create rotations from various representations.
    Rotation(Mat const &matrix);
    Rotation(AxisAngle<N> const& axisAngle);

    // Convert one representation to another.
    operator Mat() const;
    operator AxisAngle<N>() const;

private:
    enum RepresentationType {
        IS_MATRIX,
        IS_AXIS_ANGLE,
    };

    RepresentationType mType;
    mutable Mat mMatrix;
    mutable AxisAngle<N> mAxisAngle;

    // Convert a rotation matrix to an axis-angle pair.  Let (x0,x1,x2) be the
    // axis let t be an angle of rotation.  The rotation matrix is
    // [GTE_USE_MAT_VEC]
    //   R = I + sin(t)*S + (1-cos(t))*S^2
    // or
    // [GTE_USE_VEC_MAT]
    //   R = I - sin(t)*S + (1-cos(t))*S^2
    // where I is the identity and S = {{0,-x2,x1},{x2,0,-x0},{-x1,x0,0}}
    // where the inner-brace triples are the rows of the matrix.  If t > 0,
    // R represents a counterclockwise rotation; see the comments for the
    // constructor Matrix3x3(axis,angle).  It may be shown that cos(t) =
    // (trace(R)-1)/2 and R - Transpose(R) = 2*sin(t)*S.  As long as sin(t) is
    // not zero, we may solve for S in the second equation, which produces the
    // axis direction U = (S21,S02,S10).  When t = 0, the rotation is the
    // identity, in which case any axis direction is valid; we choose (1,0,0).
    // When t = pi, it must be that R - Transpose(R) = 0, which prevents us
    // from extracting the axis.  Instead, note that (R+I)/2 = I+S^2 = U*U^T,
    // where U is a unit-length axis direction.
    static void Convert(Mat const &r, AxisAngle<N> &a);

    // Convert an axis-angle pair to a rotation matrix.  Assuming (x0,x1,x2)
    // is for a right-handed world (x0 to right, x1 up, x2 out of plane of
    // page), a positive angle corresponds to a counterclockwise rotation from
    // the perspective of an observer looking at the origin of the plane of
    // rotation and having view direction the negative of the axis direction.
    // The coordinate-axis rotations are the following, where
    // unit(0) = (1,0,0), unit(1) = (0,1,0), unit(2) = (0,0,1),
    // [GTE_USE_MAT_VEC]
    //   R(unit(0),t) = {{ 1, 0, 0}, { 0, c,-s}, { 0, s, c}}
    //   R(unit(1),t) = {{ c, 0, s}, { 0, 1, 0}, {-s, 0, c}}
    //   R(unit(2),t) = {{ c,-s, 0}, { s, c, 0}, { 0, 0, 1}}
    // or
    // [GTE_USE_VEC_MAT]
    //   R(unit(0),t) = {{ 1, 0, 0}, { 0, c, s}, { 0,-s, c}}
    //   R(unit(1),t) = {{ c, 0,-s}, { 0, 1, 0}, { s, 0, c}}
    //   R(unit(2),t) = {{ c, s, 0}, {-s, c, 0}, { 0, 0, 1}}
    // where c = cos(t), s = sin(t), and the inner-brace triples are rows of
    // the matrix.  The general matrix is
    // [GTE_USE_MAT_VEC]
    //       +-                                                          -+
    //   R = | (1-c)*x0^2  + c     (1-c)*x0*x1 - s*x2  (1-c)*x0*x2 + s*x1 |
    //       | (1-c)*x0*x1 + s*x2  (1-c)*x1^2  + c     (1-c)*x1*x2 - s*x0 |
    //       | (1-c)*x0*x2 - s*x1  (1-c)*x1*x2 + s*x0  (1-c)*x2^2  + c    |
    //       +-                                                          -+
    // [GTE_USE_VEC_MAT]
    //       +-                                                          -+
    //   R = | (1-c)*x0^2  + c     (1-c)*x0*x1 + s*x2  (1-c)*x0*x2 - s*x1 |
    //       | (1-c)*x0*x1 - s*x2  (1-c)*x1^2  + c     (1-c)*x1*x2 + s*x0 |
    //       | (1-c)*x0*x2 + s*x1  (1-c)*x1*x2 - s*x0  (1-c)*x2^2  + c    |
    //       +-                                                          -+
    static void Convert(AxisAngle<N> const &a, Mat &r);
};


template <int N>
Rotation<N>::Rotation(Mat const &matrix) : mType(IS_MATRIX), mMatrix(matrix) {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");
}

template <int N>
Rotation<N>::Rotation(AxisAngle<N> const &axisAngle) : mType(IS_AXIS_ANGLE), mAxisAngle(axisAngle) {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");
}

template <int N>
Rotation<N>::operator Mat() const {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");

    switch (mType)
    {
    case IS_MATRIX:
        break;
    case IS_AXIS_ANGLE:
        Convert(mAxisAngle, mMatrix);
        break;
    }

    return mMatrix;
}

template <int N>
Rotation<N>::operator AxisAngle<N>() const {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");

    switch (mType) {
    case IS_MATRIX:
        Convert(mMatrix, mAxisAngle);
        break;
    case IS_AXIS_ANGLE:
        break;
    }

    return mAxisAngle;
}

template <int N>
void Rotation<N>::Convert(Mat const &r, AxisAngle<N> &a) {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");

    double trace = r(0, 0) + r(1, 1) + r(2, 2);
    double cs = 0.5 * (trace - 1.0);
    cs = std::max(std::min(cs, 1.0), -1.0);
    a.angle = std::acos(cs);  // The angle is in [0,pi].
    a.axis = Vector::Zero(N);

    if (a.angle > 0.0) {
        if (a.angle < PI) {
            // The angle is in (0,pi).
            a.axis(0) = r(1, 2) - r(2, 1);
            a.axis(1) = r(2, 0) - r(0, 2);
            a.axis(2) = r(0, 1) - r(1, 0);
            a.axis /= a.axis.norm();
        } else {
            // The angle is pi, in which case R is symmetric and
            // R+I = 2*(I+S^2) = 2*U*U^T, where U = (u0,u1,u2) is the
            // unit-length direction of the rotation axis.  Determine the
            // largest diagonal entry of R+I and normalize the
            // corresponding row to produce U.  It does not matter the
            // sign on u[d] for chosen diagonal d, because R(U,pi) = R(-U,pi).
            if (r(0, 0) >= r(1, 1)) {
                if (r(0, 0) >= r(2, 2)) {
                    // r00 is maximum diagonal term
                    a.axis(0) = r(0, 0) + 1.0;
                    a.axis(1) = 0.5 * (r(0, 1) + r(1, 0));
                    a.axis(2) = 0.5 * (r(0, 2) + r(2, 0));
                } else {
                    // r22 is maximum diagonal term
                    a.axis(0) = 0.5 * (r(2, 0) + r(0, 2));
                    a.axis(1) = 0.5 * (r(2, 1) + r(1, 2));
                    a.axis(2) = r(2, 2) + 1.0;
                }
            } else {
                if (r(1, 1) >= r(2, 2)) {
                    // r11 is maximum diagonal term
                    a.axis(0) = 0.5*(r(1, 0) + r(0, 1));
                    a.axis(1) = r(1, 1) + 1.0;
                    a.axis(2) = 0.5*(r(1, 2) + r(2, 1));
                } else {
                    // r22 is maximum diagonal term
                    a.axis(0) = 0.5*(r(2, 0) + r(0, 2));
                    a.axis(1) = 0.5*(r(2, 1) + r(1, 2));
                    a.axis(2) = r(2, 2) + 1.0;
                }
            }
            a.axis /= a.axis.norm();
        }
    } else {
        // The angle is 0 and the matrix is the identity.  Any axis will
        // work, so choose the Unit(0) axis.
        a.axis(0) = 1.0;
    }
}

template <int N>
void Rotation<N>::Convert(AxisAngle<N> const &a, Mat &r) {
    static_assert(N == 3 || N == 4, "Dimension must be 3 or 4.");

    r = Mat::Identity(N,N);

    double cs = std::cos(a.angle);
    double sn = std::sin(a.angle);
    double oneMinusCos = 1.0 - cs;
    double x0sqr = a.axis(0) * a.axis(0);
    double x1sqr = a.axis(1) * a.axis(1);
    double x2sqr = a.axis(2) * a.axis(2);
    double x0x1m = a.axis(0) * a.axis(1) * oneMinusCos;
    double x0x2m = a.axis(0) * a.axis(2) * oneMinusCos;
    double x1x2m = a.axis(1) * a.axis(2) * oneMinusCos;
    double x0Sin = a.axis(0) * sn;
    double x1Sin = a.axis(1) * sn;
    double x2Sin = a.axis(2) * sn;

    r(0, 0) = x0sqr * oneMinusCos + cs;
    r(1, 0) = x0x1m - x2Sin;
    r(2, 0) = x0x2m + x1Sin;
    r(0, 1) = x0x1m + x2Sin;
    r(1, 1) = x1sqr * oneMinusCos + cs;
    r(2, 1) = x1x2m - x0Sin;
    r(0, 2) = x0x2m - x1Sin;
    r(1, 2) = x1x2m + x0Sin;
    r(2, 2) = x2sqr * oneMinusCos + cs;
}

#endif // ROTATION_HPP