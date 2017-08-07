#ifndef GAUSS_SMOOTH_HPP
#define GAUSS_SMOOTH_HPP

#include <cmath>

void GaussSmooth(float *image,
                 const int width,
                 const int height,
                 const int depth,
                 float sigma,
                 int kernelWidth);

#endif // GAUSS_SMOOTH_HPP