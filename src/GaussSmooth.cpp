
#include "GaussSmooth.hpp"

#include <cmath>

#include "Debug.hpp"
#include "UtilityFunctions.hpp"

void GaussSmooth(float *image,
                 const int width,
                 const int height,
                 const int depth,
                 float sigma,
                 int kernelWidth)
{
    float kernel[kernelWidth];

    int c = kernelWidth / 2;
    float sum = 0;

    for (int i = 0; i < kernelWidth; ++i)
    {
        float val = 1.0 / (std::sqrt(2*M_PI) * sigma) * std::exp(-(std::pow(i-c, 2)/(2*sigma*sigma)));
        kernel[i] = val;
        sum += val;
    }

    for (int i = 0; i < kernelWidth; ++i)
    {
        kernel[i] /= sum;
    }

    int N = width*height*depth;

    float *tmp = (float*) malloc(N * sizeof(float));

    for (int i = 0; i < N; ++i) tmp[i] = image[i];

    // Do 1D convolution in x-direction
    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                float kernelSum = 0;

                for (int ki = 0; ki < kernelWidth; ++ki)
                {
                    int ii = std::min(std::max(i - c + ki, 0), width - 1);
                    kernelSum += tmp[gIndex(ii,j,k,height,width,0)] * kernel[ki];

                }

                image[gIndex(i,j,k,height,width,0)] = kernelSum;
            }
        }
    }

    // Do 1D convolution in y-direction
    for (int k = 0; k < depth; ++k)
    {
        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                float kernelSum = 0;

                for (int ki = 0; ki < kernelWidth; ++ki)
                {
                    int jj = std::min(std::max(j - c + ki, 0), height - 1);
                    kernelSum += image[gIndex(i,jj,k,height,width,0)] * kernel[ki];

                }

                tmp[gIndex(i,j,k,height,width,0)] = kernelSum;
            }
        }
    }

    // Do 1D convolution in z-direction
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            for (int k = 0; k < depth; ++k)
            {
                float kernelSum = 0;

                for (int ki = 0; ki < kernelWidth; ++ki)
                {
                    int kk = std::min(std::max(k - c + ki, 0), depth - 1);
                    kernelSum += tmp[gIndex(i,j,kk,height,width,0)] * kernel[ki];

                }

                image[gIndex(i,j,k,height,width,0)] = kernelSum;
            }
        }
    }
}