// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#include "MinimizerN.hpp"
#include <cstring>
#include <iostream>

// The Cartesian-product domain provided to GetMinimum(*) has minimum values
// stored in t0[0..d-1] and maximum values stored in t1[0..d-1], where d is
// 'dimensions'.  The domain is searched along lines through the current
// estimate of the minimum location.  Each such line is searched for a minimum
// using a Minimize1<double> object.  This is called "Powell's Direction Set
// Method".  The parameters 'maxLevel' and 'maxBracket' are used by
// Minimize1<double>, so read the documentation for that class (in its header
// file) to understand what these mean.  The input 'maxIterations' is the
// number of iterations for the direction-set method.

MinimizeN::MinimizeN(int dimensions, std::function<double(const Vector &)> const& F, int maxLevel, int maxBracket, int maxIterations)
    : mDimensions(dimensions), mFunction(F), mMaxIterations(maxIterations), mDirections(dimensions + 1),
    mDConjIndex(dimensions), mDCurrIndex(0), mTCurr(dimensions), mTSave(dimensions), mMinimizer(
        [this](double t) {
            return mFunction(mTCurr + t*mDirections[mDCurrIndex]);
        },
        maxLevel, maxBracket) {
            for (auto& direction : mDirections) {
                direction.resize(dimensions);
            }
        }

void MinimizeN::GetMinimum(Vector const &t0, Vector const &t1, Vector const &tInitial, Vector &tMin, double& fMin) {
    // The initial guess.
    mFCurr = mFunction(tInitial);
    mTCurr = tInitial;
    mTSave = tInitial;

    std::cout << "mTCurr: " << mTCurr.transpose() << std::endl;

    // Initialize the direction set to the standard Euclidean basis.
    for (int i = 0; i < mDimensions; ++i) {
        Vector newDirection = Vector::Zero(mDimensions);
        newDirection(i) = 1;
        mDirections[i] = newDirection;
    }

    double ell0, ell1, ellMin;
    for (int iter = 0; iter < mMaxIterations; ++iter) {
        // Find minimum in each direction and update current location.
        for (int i = 0; i < mDimensions; ++i) {
            mDCurrIndex = i;
            ComputeDomain(t0, t1, ell0, ell1);
            mMinimizer.GetMinimum(ell0, ell1, 0.0, ellMin, mFCurr);
            mTCurr += ellMin*mDirections[i];
            std::cout << "mTCurr (inner): " << mTCurr.transpose() << std::endl;
        }
        std::cout << "mTCurr (outer1): " << mTCurr.transpose() << std::endl;

        // Estimate a unit-length conjugate direction.  TODO: Expose
        // epsilon to the caller.
        mDirections[mDConjIndex] = mTCurr - mTSave;
        std::cout << "mTCurr (outer2): " << mTCurr.transpose() << std::endl;
        double length = mDirections[mDConjIndex].norm();
        double const epsilon = 1e-06;
        if (length < epsilon) {
            // New position did not change significantly from old one.
            // Should there be a better convergence criterion here?
            break;
        }

        mDirections[mDConjIndex] /= length;

        // Minimize in conjugate direction.
        mDCurrIndex = mDConjIndex;
        ComputeDomain(t0, t1, ell0, ell1);
        mMinimizer.GetMinimum(ell0, ell1, 0.0, ellMin, mFCurr);
        mTCurr += ellMin*mDirections[mDCurrIndex];
        std::cout << "mTCurr (outer3): " << mTCurr.transpose() << std::endl;

        // Cycle the directions and add conjugate direction to set.
        mDConjIndex = 0;
        for (int i = 0; i < mDimensions; ++i) {
            mDirections[i] = mDirections[i + 1];
        }

        // Set parameters for next pass.
        mTSave = mTCurr;
    }

    tMin = mTCurr;
    fMin = mFCurr;
    std::cout << "tMin (inside): " << tMin.transpose() << std::endl;

}

void MinimizeN::ComputeDomain(Vector const &t0, Vector const &t1, double &ell0, double &ell1) {
    ell0 = -std::numeric_limits<double>::max();
    ell1 = +std::numeric_limits<double>::max();

    for (int i = 0; i < mDimensions; ++i) {
        double value = mDirections[mDCurrIndex][i];
        if (value != 0.0) {
            double b0 = t0(i) - mTCurr(i);
            double b1 = t1(i) - mTCurr(i);
            double inv = 1.0 / value;
            if (value > 0.0) {
                // The valid t-interval is [b0,b1].
                b0 *= inv;
                if (b0 > ell0) {
                    ell0 = b0;
                }
                b1 *= inv;
                if (b1 < ell1) {
                    ell1 = b1;
                }
            } else {
                // The valid t-interval is [b1,b0].
                b0 *= inv;
                if (b0 < ell1) {
                    ell1 = b0;
                }
                b1 *= inv;
                if (b1 > ell0) {
                    ell0 = b1;
                }
            }
        }
    }

    // Correction if numerical errors lead to values nearly zero.
    if (ell0 > 0.0) {
        ell0 = 0.0;
    }

    if (ell1 < 0.0) {
        ell1 = 0.0;
    }
}
