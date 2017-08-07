#include "SignedDistanceMapResult.hpp"

#include "Debug.hpp"

SignedDistanceMapResult::SignedDistanceMapResult(int width,
                                                 int height,
                                                 int depth)
                                                 :
                                                 width(width),
                                                 height(height),
                                                 depth(depth)
{
    int N = width * height * depth;

    sdm = (double *) malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        sdm[i] = 0;
    }
}