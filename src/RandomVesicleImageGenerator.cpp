#include "Debug.hpp"

#include <cmath>

#include "ellipsoidHelpers.hpp"
#include "GaussSmooth.hpp"
#include "RandomVesicleImageGenerator.hpp"

RandomVesicleImageGenerator::RandomVesicleImageGenerator() :
                                                         minRadiiA(3),
                                                         maxRadiiA(3),
                                                         minRadiiB(3),
                                                         maxRadiiB(3),
                                                         minRadiiC(3),
                                                         maxRadiiC(3),
                                                         randomRotation(false)
{

}
RandomVesicleImageGenerator::RandomVesicleImageGenerator(float minRadiiA,
                                                         float maxRadiiA,
                                                         float minRadiiB,
                                                         float maxRadiiB,
                                                         float minRadiiC,
                                                         float maxRadiiC,
                                                         bool randomRotation):
                                                         minRadiiA(minRadiiA),
                                                         maxRadiiA(maxRadiiA),
                                                         minRadiiB(minRadiiB),
                                                         maxRadiiB(maxRadiiB),
                                                         minRadiiC(minRadiiC),
                                                         maxRadiiC(maxRadiiC),
                                                         randomRotation(randomRotation)
{

}

void RandomVesicleImageGenerator::generateRandomVesicleImage(float *image,
                                                             unsigned int &label,
                                                             int width,
                                                             int height,
                                                             int depth,
                                                             bool dense,
                                                             float noise)
{
    float ci = (float) (width  - 1) / 2.0;
    float cj = (float) (height - 1) / 2.0;
    float ck = (float) (depth  - 1) / 2.0;

    float *R = (float *) malloc(3*3 * sizeof(float));
    float *Rt = (float *) malloc(3*3 * sizeof(float));
    float *D = (float *) malloc(3*3 * sizeof(float));
    float *RD = (float *) malloc(3*3 * sizeof(float));
    float *RDRt = (float *) malloc(3*3 * sizeof(float));
    float *RDRt_inside = (float *) malloc(3*3 * sizeof(float));

    // Generate random rotation matrix
    RandomRotationMatrix(R);

    // Make a transpose of the rotation matrix
    for (int i = 0; i < 9; ++i) Rt[i] = R[i];
    Transpose(Rt);

    // Generate radii

    float a = Random::randU(Random::generator) * (maxRadiiA - minRadiiA) + minRadiiA;
    float b = Random::randU(Random::generator) * (maxRadiiB - minRadiiB) + minRadiiB;
    float c = Random::randU(Random::generator) * (maxRadiiC - minRadiiC) + minRadiiC;

    // Set the diagonal radii squared matrix
    for (int i = 0; i < 9; ++i) D[i] = 0;

    D[0] = 1.0/(a * a);
    D[4] = 1.0/(b * b);
    D[8] = 1.0/(c * c);

    MatrixMult(R, D, RD);
    MatrixMult(RD, Rt, RDRt);

    float a_inside = a - 2.5;
    float b_inside = b - 2.5;
    float c_inside = c - 2.5;

    bool makeInside = (not dense) and (a_inside > 0)
                      and (b_inside > 0) and (c_inside > 0);

    if (makeInside) {
        D[0] = 1.0/(a_inside * a_inside);
        D[4] = 1.0/(b_inside * b_inside);
        D[8] = 1.0/(c_inside * c_inside);

        MatrixMult(R, D, RD);
        MatrixMult(RD, Rt, RDRt_inside);
    }

    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                float val = GetAlgebraicFormValue(i, j, k, ci, cj, ck, RDRt);

                if (val < 0) image[gIndex(i, j, k, height, width)] = 110;
                else image[gIndex(i, j, k, height, width)] = 173;

                if (makeInside)
                {
                    if (GetAlgebraicFormValue(i, j, k, ci, cj, ck, RDRt_inside) < 0)
                    {
                        image[gIndex(i, j, k, height, width)] = 150;
                    }
                }
            }
        }
    }

    GaussSmooth(image, width, height, depth, 0.5, 5);

    if (noise)
    {
        addNoise(image, width, height, depth, 10);
    }

    // Set label
    if (a < b) std::swap(a, b);
    if (a < c) std::swap(a, c);
    if (b < c) std::swap(b, c);

    if (c / a > 0.9)
    {
        // Spheroid
        if (dense) label = 3;
        else label = 0;
    }
    else
    {
        if (a - b < b - c)
        {
            // Oblate
            if (dense) label = 4;
            else label = 1;
        }
        else
        {
            //prolate
            if (dense) label = 5;
            else label = 2;
        }
    }

    free(R);
    free(Rt);
    free(D);
    free(RD);
    free(RDRt);
    free(RDRt_inside);
}

void RandomVesicleImageGenerator::generateResponses(float *responseKernels, int width, int height, int depth, int featureSize)
{
    float ci = (float) (width  - 1) / 2.0;
    float cj = (float) (height - 1) / 2.0;
    float ck = (float) (depth  - 1) / 2.0;

    float *R = (float *) malloc(3*3 * sizeof(float));
    float *Rt = (float *) malloc(3*3 * sizeof(float));
    float *D = (float *) malloc(3*3 * sizeof(float));
    float *RD = (float *) malloc(3*3 * sizeof(float));
    float *RDRt = (float *) malloc(3*3 * sizeof(float));
    float *RDRt_inside = (float *) malloc(3*3 * sizeof(float));
    float *image = (float *) malloc(width*height*depth * sizeof(float));

    for (int fIdx = 0; fIdx < featureSize; ++fIdx)
    {
        // Generate random rotation matrix
        RandomRotationMatrix(R);

        // Make a transpose of the rotation matrix
        for (int i = 0; i < 9; ++i) Rt[i] = R[i];
        Transpose(Rt);

        // Generate radii

        float a = Random::randU(Random::generator) * (maxRadiiA - minRadiiA) + minRadiiA;
        float b = Random::randU(Random::generator) * (maxRadiiB - minRadiiB) + minRadiiB;
        float c = Random::randU(Random::generator) * (maxRadiiC - minRadiiC) + minRadiiC;

        // Set the diagonal radii squared matrix
        for (int i = 0; i < 9; ++i) D[i] = 0;

        D[0] = 1.0/(a * a);
        D[4] = 1.0/(b * b);
        D[8] = 1.0/(c * c);

        MatrixMult(R, D, RD);
        MatrixMult(RD, Rt, RDRt);

        float a_inside = a - 2.5;
        float b_inside = b - 2.5;
        float c_inside = c - 2.5;

        bool dense = fIdx % 2 == 0;

        bool makeInside = (not dense) and (a_inside > 0)
                          and (b_inside > 0) and (c_inside > 0);

        if (makeInside) {
            D[0] = 1.0/(a_inside * a_inside);
            D[4] = 1.0/(b_inside * b_inside);
            D[8] = 1.0/(c_inside * c_inside);

            MatrixMult(R, D, RD);
            MatrixMult(RD, Rt, RDRt_inside);
        }

        for (int k = 0; k < depth; ++k)
        {
            for (int j = 0; j < height; ++j)
            {
                for (int i = 0; i < width; ++i)
                {
                    float val = GetAlgebraicFormValue(i, j, k, ci, cj, ck, RDRt);

                    if (val < 0) image[gIndex(i, j, k, height, width)] = 110;
                    else image[gIndex(i, j, k, height, width)] = 173;

                    if (makeInside)
                    {
                        if (GetAlgebraicFormValue(i, j, k, ci, cj, ck, RDRt_inside) < 0)
                        {
                            image[gIndex(i, j, k, height, width)] = 150;
                        }
                    }
                }
            }
        }

        GaussSmooth(image, width, height, depth, 0.5, 5);

        Normalize(image, width, height, depth);

        for (int k = 0; k < depth; ++k)
        {
            for (int j = 0; j < height; ++j)
            {
                for (int i = 0; i < width; ++i)
                {
                    int rIdx = gIndex(i, j, k, fIdx, 17, 17, 17);
                    int iIdx = gIndex(i, j, k, 17, 17);

                    responseKernels[rIdx] = image[iIdx];
                }
            }
        }
    }

    free(R);
    free(Rt);
    free(D);
    free(RD);
    free(RDRt);
    free(RDRt_inside);
    free(image);
}

void RandomVesicleImageGenerator::MatrixMult(float *A, float *B, float *C)
{
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            C[gIndex(i,j,3)] = 0;

            for (int k = 0; k < 3; ++k)
            {
                C[gIndex(i,j,3)] += A[gIndex(k,j,3)] * B[gIndex(i,k,3)];
            }
        }
    }
}

void RandomVesicleImageGenerator::RandomRotationMatrix(float *R) {
    /*
     * Rework of http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
     * Author: Jim Arvo, 1991
     */
    float theta = Random::randU(Random::generator) * M_PI * 2;
    float phi   = Random::randU(Random::generator) * M_PI * 2;
    float z     = Random::randU(Random::generator) * 2.0;

    float r  = std::sqrt(z);
    float Vx = std::sin(phi) * r;
    float Vy = std::cos(phi) * r;
    float Vz = std::sqrt(2.0 - z);

    float st = std::sin(theta);
    float ct = std::cos(theta);
    float Sx = Vx * ct - Vy * st;
    float Sy = Vx * st + Vy * ct;

    /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
    /* is equivalent to V S - R.                                        */

    R[gIndex(0,0,3)] = Vx * Sx - ct;
    R[gIndex(0,1,3)] = Vx * Sy - st;
    R[gIndex(0,2,3)] = Vx * Vz;

    R[gIndex(1,0,3)] = Vy * Sx + st;
    R[gIndex(1,1,3)] = Vy * Sy - ct;
    R[gIndex(1,2,3)] = Vy * Vz;

    R[gIndex(2,0,3)] = Vz * Sx;
    R[gIndex(2,1,3)] = Vz * Sy;
    R[gIndex(2,2,3)] = 1.0 - z;
}

void RandomVesicleImageGenerator::Transpose(float *A)
{
    std::swap(A[gIndex(1,0,3)], A[gIndex(0,1,3)]);
    std::swap(A[gIndex(2,0,3)], A[gIndex(0,2,3)]);
    std::swap(A[gIndex(2,1,3)], A[gIndex(1,2,3)]);
}

float RandomVesicleImageGenerator::GetAlgebraicFormValue(float i,
                                                         float j,
                                                         float k,
                                                         float ci,
                                                         float cj,
                                                         float ck,
                                                         float *RDRt)
{
    float px = i - ci;
    float py = j - cj;
    float pz = k - ck;

    float tmpA = px * RDRt[gIndex(0,0,3)]
               + py * RDRt[gIndex(0,1,3)]
               + pz * RDRt[gIndex(0,2,3)];

    float tmpB = px * RDRt[gIndex(1,0,3)]
               + py * RDRt[gIndex(1,1,3)]
               + pz * RDRt[gIndex(1,2,3)];

    float tmpC = px * RDRt[gIndex(2,0,3)]
               + py * RDRt[gIndex(2,1,3)]
               + pz * RDRt[gIndex(2,2,3)];

    return px * tmpA + py * tmpB + pz * tmpC - 1;
}

void RandomVesicleImageGenerator::addNoise(float *image, int width, int height, int depth, float stddev)
{
    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                int gIdx = gIndex(i, j, k, height, width);

                float error = Random::randU(Random::generator) * stddev;

                float newVal = image[gIdx] + error;

                image[gIdx] = std::max((float) 0, std::min(newVal, (float) 255));
            }
        }
    }
}

void RandomVesicleImageGenerator::Normalize(float *image, int width, int height, int depth)
{
    int N = width * height * depth;

    float min = 255.0;
    float max = 0;

    for (int i = 0; i < N; ++i)
    {
        min = std::min(min, image[i]);
        max = std::max(max, image[i]);
    }

    float d = max - min;

    for (int i = 0; i < N; ++i)
    {
        image[i] = (image[i] - min) / d;
    }
}