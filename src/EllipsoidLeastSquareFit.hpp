#ifndef ELLIPSOID_LEAST_SQUARE_FIT_HPP
#define ELLIPSOID_LEAST_SQUARE_FIT_HPP

#include "types.hpp"

#include "Ellipsoid.hpp"

Ellipsoid EllipsoidLeastSquareFit(const Mat &Xin, double &error, bool &errorFlag);

#endif // ELLIPSOID_LEAST_SQUARE_FIT_HPP