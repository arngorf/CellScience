#ifndef SIGNEDDISTANCEMAP_HPP
#define SIGNEDDISTANCEMAP_HPP

#include "types.hpp"
#include "Debug.hpp"
#include <cmath>
#include <algorithm>

inline float pos(float x)
{
    if (x > 0) return x;
    else return 0;
}

inline float neg(float x)
{
    if (x < 0) return x;
    else return 0;
}

inline float max(float x, float y)
{
    if (x > y) return x;
    else return y;
}

void SignedDistanceMap(float * const Phi,
                       const int nx,
                       const int ny,
                       const int nz,
                       const int numIter)
{
    int N = nz*ny*nx;

    float dx = 1.0;

    float *Phi0 = (float*) malloc(N * sizeof(float));

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                Phi0[gIndex(x, y, z, ny, nx, 0)] =
                    Phi[gIndex(x, y, z, ny, nx, 1)];
            }
        }
    }

    for (int n = 0; n < numIter; ++n) {

        float *V = (float*) malloc(N * sizeof(float));

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {

                    float Phixyz = Phi[gIndex(x, y, z, ny, nx, 1)];
                    float Phi0xyz = Phi0[gIndex(x, y, z, ny, nx, 0)];

                    float a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
                    if (x > 0)      a = ( Phixyz
                                        - Phi[gIndex(x-1, y, z, ny, nx, 1)])
                                        / dx;
                    if (x < nx - 1) b = ( Phi[gIndex(x+1, y, z, ny, nx, 1)]
                                        - Phixyz)
                                        / dx;
                    if (y > 0)      c = ( Phixyz
                                        - Phi[gIndex(x, y-1, z, ny, nx, 1)])
                                        / dx;
                    if (y < ny - 1) d = ( Phi[gIndex(x, y+1, z, ny, nx, 1)]
                                        - Phixyz)
                                        / dx;
                    if (z > 0)      e = ( Phixyz
                                        - Phi[gIndex(x, y, z-1, ny, nx, 1)])
                                        / dx;
                    if (z < nz - 1) f = ( Phi[gIndex(x, y, z+1, ny, nx, 1)]
                                        - Phixyz)
                                        / dx;

                    float S = Phi0xyz / sqrt(Phi0xyz*Phi0xyz + 1);

                    float G = 0.0;

                    if (Phi0xyz > 0) {
                        G = 1. - std::sqrt(
                                  std::max(std::pow(pos(a),2),std::pow(neg(b),2))
                                + std::max(std::pow(pos(c),2),std::pow(neg(d),2))
                                + std::max(std::pow(pos(e),2),std::pow(neg(f),2))
                                    );
                    } else if (Phi0xyz < 0) {
                        G = 1. - std::sqrt(
                                  std::max(std::pow(neg(a),2),std::pow(pos(b),2))
                                + std::max(std::pow(neg(c),2),std::pow(pos(d),2))
                                + std::max(std::pow(neg(e),2),std::pow(pos(f),2))
                                    );
                    }

                    V[gIndex(x, y, z, ny, nx, 0)] = S * G;

                }
            }
        }

        float vMax = 0;
        for (int i = 0; i < nx*ny*nz; ++i) {
            vMax = std::max(vMax, (float) fabs(V[i]));
        }

        if (vMax < 1) {
            Debug::Info("SignedDistanceMap: Iterations stopped at iteration: "
                       + STR(n+1) + " with vMax: " + STR(vMax));
            free(V);
            break;
        }

        float dt = 0.5*dx/vMax;

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    Phi[gIndex(x,y,z,ny,nx,1)] = Phi[gIndex(x,y,z,ny,nx,1)] + dt * V[gIndex(x,y,z,ny,nx,0)];
                    if (Phi[gIndex(x,y,z,ny,nx,1)] * Phi0[gIndex(x,y,z,ny,nx,0)] < 0) Phi[gIndex(x,y,z,ny,nx,1)] = -Phi[gIndex(x,y,z,ny,nx,1)];
                }
            }
        }

        free(V);
    }

    free(Phi0);
}

#endif // SIGNEDDISTANCEMAP_HPP