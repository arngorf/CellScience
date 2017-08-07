// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#ifndef MINIMIZER_N_HPP
#define MINIMIZER_N_HPP

#include "types.hpp"
#include "Minimizer1.hpp"

// The Cartesian-product domain provided to GetMinimum(*) has minimum values
// stored in t0[0..d-1] and maximum values stored in t1[0..d-1], where d is
// 'dimensions'.  The domain is searched along lines through the current
// estimate of the minimum location.  Each such line is searched for a minimum
// using a Minimize1<double> object.  This is called "Powell's Direction Set
// Method".  The parameters 'maxLevel' and 'maxBracket' are used by
// Minimize1<double>, so read the documentation for that class (in its header
// file) to understand what these mean.  The input 'maxIterations' is the
// number of iterations for the direction-set method.


class MinimizeN {
public:
    // Construction.
    MinimizeN(int dimensions, std::function<double(const Vector &)> const& F, int maxLevel, int maxBracket, int maxIterations);

    // Find the minimum on the Cartesian-product domain whose minimum values
    // are stored in t0[0..d-1] and whose maximum values are stored in
    // t1[0..d-1], where d is 'dimensions'.  An initial guess is specified in
    // tInitial[0..d-1].  The location of the minimum is tMin[0..d-1] and
    // the value of the minimum is 'fMin'.
    void GetMinimum(Vector const &t0, Vector const &t1, Vector const &tInitial, Vector &tMin, double& fMin);

private:
    // The current estimate of the minimum location is mTCurr[0..d-1].  The
    // direction of the current line to search is mDCurr[0..d-1].  This line
    // must be clipped against the Cartesian-product domain, a process
    // implemented in this function.  If the line is mTCurr+s*mDCurr, the
    // clip result is the s-interval [ell0,ell1].
    void ComputeDomain(Vector const &t0, Vector const &t1, double& ell0, double& ell1);

    int mDimensions;
    std::function<double(Vector)> mFunction;
    int mMaxIterations;
    std::vector<Vector> mDirections;
    int mDConjIndex;
    int mDCurrIndex;
    Vector mTCurr;
    Vector mTSave;
    double mFCurr;
    Minimize1 mMinimizer;
};

#endif // MINIMIZER_N_HPP