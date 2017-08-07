// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#include "Minimizer1.hpp"
#include <cmath>
#include <limits>

// The interval [t0,t1] provided to GetMinimum(double,double,double,double&,double&)
// is processed by examining subintervals.  On each subinteral [a,b], the
// values f0 = F(a), f1 = F((a+b)/2), and f2 = F(b) are examined.  If
// {f0,f1,f2} is monotonic, then [a,b] is subdivided and processed.  The
// maximum depth of the recursion is limited by 'maxLevel'.  If {f0,f1,f2}
// is not monotonic, then two cases arise.  First, if f1 = min{f0,f1,f2},
// then {f0,f1,f2} is said to "bracket a minimum" and GetBracketedMinimum(*)
// is called to locate the function minimum.  The process uses a form of
// bisection called "parabolic interpolation" and the maximum number of
// bisection steps is 'maxBracket'.  Second, if f1 = max{f0,f1,f2}, then
// {f0,f1,f2} brackets a maximum.  The minimum search continues recursively
// as before on [a,(a+b)/2] and [(a+b)/2,b].

Minimize1::Minimize1(std::function<double(double)> const& F, int maxLevel, int maxBracket)
    : mFunction(F), mMaxLevel(maxLevel), mMaxBracket(maxBracket) {
}

void Minimize1::GetMinimum(double t0, double t1, double tInitial, double& tMin, double& fMin) {

    mTMin = std::numeric_limits<double>::max();
    mFMin = std::numeric_limits<double>::max();

    double f0 = mFunction(t0);
    if (f0 < mFMin)
    {
        mTMin = t0;
        mFMin = f0;
    }

    double fInitial = mFunction(tInitial);
    if (fInitial < mFMin)
    {
        mTMin = tInitial;
        mFMin = fInitial;
    }

    double f1 = mFunction(t1);
    if (f1 < mFMin)
    {
        mTMin = t1;
        mFMin = f1;
    }

    GetMinimum(t0, f0, tInitial, fInitial, t1, f1, mMaxLevel);

    tMin = mTMin;
    fMin = mFMin;
}

void Minimize1::GetMinimum(double t0, double f0, double tm, double fm, double t1, double f1, int level)
{
    if (level-- == 0)
    {
        return;
    }

    if ((t1 - tm)*(f0 - fm) > (tm - t0)*(fm - f1))
    {
        // The quadratic fit has positive second derivative at the midpoint.
        if (f1 > f0)
        {
            if (fm >= f0)
            {
                // Increasing, repeat on [t0,tm].
                GetMinimum(t0, f0, tm, fm, level);
            }
            else
            {
                // Not monotonic, have a bracket.
                GetBracketedMinimum(t0, f0, tm, fm, t1, f1, level);
            }
        }
        else if (f1 < f0)
        {
            if (fm >= f1)
            {
                // Decreasing, repeat on [tm,t1].
                GetMinimum(tm, fm, t1, f1, level);
            }
            else
            {
                // Not monotonic, have a bracket.
                GetBracketedMinimum(t0, f0, tm, fm, t1, f1, level);
            }
        }
        else
        {
            // Constant, repeat on [t0,tm] and [tm,t1].
            GetMinimum(t0, f0, tm, fm, level);
            GetMinimum(tm, fm, t1, f1, level);
        }
    }
    else
    {
        // The quadratic fit has a nonpositive second derivative at the
        // midpoint.
        if (f1 > f0)
        {
            // Repeat on [t0,tm].
            GetMinimum(t0, f0, tm, fm, level);
        }
        else if (f1 < f0)
        {
            // Repeat on [tm,t1].
            GetMinimum(tm, fm, t1, f1, level);
        }
        else
        {
            // Repeat on [t0,tm] and [tm,t1].
            GetMinimum(t0, f0, tm, fm, level);
            GetMinimum(tm, fm, t1, f1, level);
        }
    }
}

void Minimize1::GetMinimum(double t0, double f0, double t1, double f1, int level) {
    if (level-- == 0) {
        return;
    }

    double tm = 0.5*(t0 + t1);
    double fm = mFunction(tm);
    if (fm < mFMin) {
        mTMin = tm;
        mFMin = fm;
    }

    if (f0 - 2.0*fm + f1 > 0.0) {
        // The quadratic fit has positive second derivative at the midpoint.
        if (f1 > f0) {
            if (fm >= f0) {
                // Increasing, repeat on [t0,tm].
                GetMinimum(t0, f0, tm, fm, level);
            } else {
                // Not monotonic, have a bracket.
                GetBracketedMinimum(t0, f0, tm, fm, t1, f1, level);
            }
        } else if (f1 < f0) {
            if (fm >= f1) {
                // Decreasing, repeat on [tm,t1].
                GetMinimum(tm, fm, t1, f1, level);
            } else {
                // Not monotonic, have a bracket.
                GetBracketedMinimum(t0, f0, tm, fm, t1, f1, level);
            }
        } else {
            // Constant, repeat on [t0,tm] and [tm,t1].
            GetMinimum(t0, f0, tm, fm, level);
            GetMinimum(tm, fm, t1, f1, level);
        }
    } else {
        // The quadratic fit has nonpositive second derivative at the
        // midpoint.
        if (f1 > f0) {
            // Repeat on [t0,tm].
            GetMinimum(t0, f0, tm, fm, level);
        } else if (f1 < f0) {
            // Repeat on [tm,t1].
            GetMinimum(tm, fm, t1, f1, level);
        } else {
            // Repeat on [t0,tm] and [tm,t1].
            GetMinimum(t0, f0, tm, fm, level);
            GetMinimum(tm, fm, t1, f1, level);
        }
    }
}

void Minimize1::GetBracketedMinimum(double t0, double f0, double tm, double fm, double t1, double f1, int level) {
    for (int i = 0; i < mMaxBracket; ++i) {
        // Update minimum value.
        if (fm < mFMin) {
            mTMin = tm;
            mFMin = fm;
        }

        // Test for convergence.  TODO: Expose the epsilon and tolerance
        // parameters to the caller.
        double const epsilon = 1e-08;
        double const tolerance = 1e-04;
        if (std::abs(t1 - t0) <= 2*tolerance*std::abs(tm) + epsilon) {
            break;
        }

        // Compute vertex of interpolating parabola.
        double dt0 = t0 - tm;
        double dt1 = t1 - tm;
        double df0 = f0 - fm;
        double df1 = f1 - fm;
        double tmp0 = dt0*df1;
        double tmp1 = dt1*df0;
        double denom = tmp1 - tmp0;
        if (std::abs(denom) < epsilon) {
            return;
        }

        double tv = tm + 0.5*(dt1*tmp1 - dt0*tmp0) / denom;
        double fv = mFunction(tv);
        if (fv < mFMin) {
            mTMin = tv;
            mFMin = fv;
        }

        if (tv < tm) {
            if (fv < fm) {
                t1 = tm;
                f1 = fm;
                tm = tv;
                fm = fv;
            } else {
                t0 = tv;
                f0 = fv;
            }
        } else if (tv > tm) {
            if (fv < fm) {
                t0 = tm;
                f0 = fm;
                tm = tv;
                fm = fv;
            } else {
                t1 = tv;
                f1 = fv;
            }
        } else {
            // The vertex of parabola is already at middle sample point.
            GetMinimum(t0, f0, tm, fm, level);
            GetMinimum(tm, fm, t1, f1, level);
        }
    }
}
