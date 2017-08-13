#include "ChanVese3D.hpp"

#include "Debug.hpp"
#include "types.hpp"
#include "UtilityFunctions.hpp"

#include <cmath>
#include <algorithm>

ChanVese3D::ChanVese3D()
{
    dx = 1.0;
    numericEps = 0.000000000001;
    maxCurvature = 1.0/dx;
}

void ChanVese3D::segmentImage(const float * const image,
                              float * const Phi,
                              const float mu,
                              const float nu,
                              int nx,
                              int ny,
                              int nz,
                              bool initPhi,
                              int maxIter,
                              float &c1,
                              float &c2,
                              bool isConstc1,
                              bool isConstc2,
                              MainWindow *mw)
{
    int N = nz * ny * nx;

    float eps = 1.0;
    float lam1 = 1.0;
    float lam2 = 1.0;

    if (initPhi) initPhiBubble(Phi, nx, ny, nz, 0, 0, 0, 10);

    float *C = (float*) malloc(N * sizeof(float));
    float *V = (float*) malloc(N * sizeof(float));

    for (int n = 0; n < maxIter; ++n)
    {

        if (n % 50 == 0) mw->setMessage(("Segmenting: " +STR(((n+1) * 100) / maxIter) + " %").c_str());

        if (not isConstc1)
        {
            c1 = areaInside(image, Phi, nx, ny, nz, eps);
        }

        if (not isConstc2)
        {
            c2 = areaOutside(image, Phi, nx, ny, nz, eps);
        }

        if (mu > 0) curvature(Phi, C, nz, ny, nx);

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z)
        {
            for (int y = 0; y < ny; ++y)
            {
                for (int x = 0; x < nx; ++x)
                {
                    if (mu > 0)
                    {
                        V[gIndex(x,y,z,ny,nx,0)] = delta(Phi[gIndex(x,y,z,ny,nx,1)], eps)
                            * ( mu * C[gIndex(x,y,z,ny,nx,0)] - nu
                            - lam1 * pow(image[gIndex(x,y,z,ny,nx,0)] - c1, 2)
                            + lam2 * pow(image[gIndex(x,y,z,ny,nx,0)] - c2, 2)
                            );
                    }
                    else
                    {
                        V[gIndex(x,y,z,ny,nx,0)] = delta(Phi[gIndex(x,y,z,ny,nx,1)], eps)
                            * ( - nu
                            - lam1 * pow(image[gIndex(x,y,z,ny,nx,0)] - c1, 2)
                            + lam2 * pow(image[gIndex(x,y,z,ny,nx,0)] - c2, 2)
                            );
                    }
                }
            }
        }

        float Vmax = 0;
        for (int i = 0; i < nx*ny*nz; ++i) {
            Vmax = std::max(Vmax, (float) fabs(V[i]));
        }

        float dt = 0.5*dx/Vmax;

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    int bIndex = gIndex(x,y,z,ny,nx,1);
                    int sIndex = gIndex(x,y,z,ny,nx,0);
                    Phi[bIndex] = Phi[bIndex] + dt * V[sIndex];
                }
            }
        }

        applyBoundaryCondition(Phi, nz, ny, nx);

        if (n % 50 == 0) Debug::Info("Segmenting - iteration: " + STR(n)
                                    + " dt: " + STR(dt)
                                    + " vMax: " + STR(Vmax)
                                    + " averages: " + STR(c1)
                                    + ", " + STR(c2));
    }

    free(C);
    free(V);
}

void ChanVese3D::setConstant_c1(bool isConstc1, float val)
{
    this->isConstc1 = isConstc1;
    this->constc1 = val;
}

void ChanVese3D::setConstant_c2(bool isConstc2, float val)
{
    this->isConstc2 = isConstc2;
    this->constc2 = val;
}

void ChanVese3D::initPhiBubble(float * const out, const int nx, const int ny, const int nz, const float sz, const float sy, const float sx, const float bubblePeriod)
{
    #pragma omp parallel for collapse(3)
    for (int z = -1; z < nz + 1; ++z) {
        for (int y = -1; y < ny + 1; ++y) {
            for (int x = -1; x < nx + 1; ++x) {
                if (z == -1 or y == -1 or x == -1 or z == nz or y == ny or x == nx) {
                    out[gIndex(x,y,z,ny,nx,1)] = 0;
                } else {
                    out[gIndex(x,y,z,ny,nx,1)] = 5
                        * cos(PI*float(x - sx)/bubblePeriod)
                        * cos(PI*float(y - sy)/bubblePeriod)
                        * cos(PI*float(z - sz)/bubblePeriod);
                }
            }
        }
    }
}

void ChanVese3D::initPhiSphere(float * const out, const int nz, const int ny, const int nx, const float sz, const float sy, const float sx, const float radius)
{
    #pragma omp parallel for collapse(3)
    for (int z = -1; z < nz + 1; ++z) {
        for (int y = -1; y < ny + 1; ++y) {
            for (int x = -1; x < nx + 1; ++x) {
                out[gIndex(x,y,z,ny,nx,1)] = pow(float(x - sx), 2.0)
                                           + pow(float(y - sy), 2.0)
                                           + pow(float(z - sz), 2.0)
                                           - radius * radius;
            }
        }
    }
}

void ChanVese3D::curvature(const float * const Phi, float * const C, const int nz, const int ny, const int nx)
{
    #pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {

                float Dx  = ( Phi[gIndex(x+1,y  ,z  ,ny,nx,1)]
                           -   Phi[gIndex(x-1,y  ,z  ,ny,nx,1)])/(2*dx);
                float Dy  = ( Phi[gIndex(x  ,y+1,z  ,ny,nx,1)]
                           -   Phi[gIndex(x  ,y-1,z  ,ny,nx,1)])/(2*dx);
                float Dz  = ( Phi[gIndex(x  ,y  ,z+1,ny,nx,1)]
                           -   Phi[gIndex(x  ,y  ,z-1,ny,nx,1)])/(2*dx);

                float Dxx = ( Phi[gIndex(x+1,y  ,z  ,ny,nx,1)]
                           - 2*Phi[gIndex(x  ,y  ,z  ,ny,nx,1)]
                           +   Phi[gIndex(x-1,y  ,z  ,ny,nx,1)])/(dx*dx);
                float Dyy = ( Phi[gIndex(x  ,y+1,z  ,ny,nx,1)]
                           - 2*Phi[gIndex(x  ,y  ,z  ,ny,nx,1)]
                           +   Phi[gIndex(x  ,y-1,z  ,ny,nx,1)])/(dx*dx);
                float Dzz = ( Phi[gIndex(x  ,y  ,z+1,ny,nx,1)]
                           - 2*Phi[gIndex(x  ,y  ,z  ,ny,nx,1)]
                           +   Phi[gIndex(x  ,y  ,z-1,ny,nx,1)])/(dx*dx);

                float Dxy = (Phi[gIndex(x+1,y+1,z  ,ny,nx,1)]
                           -  Phi[gIndex(x-1,y+1,z  ,ny,nx,1)]
                           -  Phi[gIndex(x+1,y-1,z  ,ny,nx,1)]
                           +  Phi[gIndex(x-1,y-1,z  ,ny,nx,1)])/(4*dx*dx);
                float Dxz = (Phi[gIndex(x+1,y  ,z+1,ny,nx,1)]
                           -  Phi[gIndex(x-1,y  ,z+1,ny,nx,1)]
                           -  Phi[gIndex(x+1,y  ,z-1,ny,nx,1)]
                           +  Phi[gIndex(x-1,y  ,z-1,ny,nx,1)])/(4*dx*dx);
                float Dyz = (Phi[gIndex(x  ,y+1,z+1,ny,nx,1)]
                           -  Phi[gIndex(x  ,y-1,z+1,ny,nx,1)]
                           -  Phi[gIndex(x  ,y+1,z-1,ny,nx,1)]
                           +  Phi[gIndex(x  ,y-1,z-1,ny,nx,1)])/(4*dx*dx);

                float G = pow(sqrt(Dx*Dx + Dy*Dy + Dz*Dz), 3);

                C[gIndex(x,y,z,ny,nx,0)] = ( Dx*Dx*(Dyy + Dzz)
                                           + Dy*Dy*(Dxx + Dzz)
                                           + Dz*Dz*(Dxx + Dyy)
                                           - 2*Dx*Dy*Dxy
                                           - 2*Dx*Dz*Dxz
                                           - 2*Dy*Dz*Dyz
                                           ) / (G + numericEps);

                // Clamp curvature to prevent excessive curvature of Phi
                C[gIndex(x,y,z,ny,nx,0)] =
                    std::max(-maxCurvature,
                        std::min(C[gIndex(x,y,z,ny,nx,0)],
                            maxCurvature
                            )
                        );
            }
        }
    }
}

inline float ChanVese3D::delta(float x, float eps)
{
    return eps / (PI*(x*x + eps*eps));
}

inline float ChanVese3D::heaviside(float x, float eps)
{
    return 0.5 * (1.0 + 2.0/PI * atan( x / eps));
}

float ChanVese3D::areaInside(const float * const image, const float * const Phi, int nx, int ny, int nz, float eps)
{
    float pixelWeight = 0.0;
    float totalWeight = 0.0;

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                pixelWeight += image[gIndex(x,y,z,ny,nx,0)]
                            * heaviside(Phi[gIndex(x,y,z,ny,nx,1)], eps);
                totalWeight += heaviside(Phi[gIndex(x,y,z,ny,nx,1)], eps);
            }
        }
    }

    return pixelWeight / totalWeight;
}

float ChanVese3D::areaOutside(const float * const image, const float * const Phi,
                  int nx, int ny, int nz, float eps)
{

    float pixelWeight = 0.0;
    float totalWeight = 0.0;

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                pixelWeight += image[gIndex(x,y,z,ny,nx,0)]
                            * (1 - heaviside(Phi[gIndex(x,y,z,ny,nx,1)], eps));
                totalWeight += 1 - heaviside(Phi[gIndex(x,y,z,ny,nx,1)], eps);
            }
        }
    }

    return pixelWeight / totalWeight;
}

void ChanVese3D::applyBoundaryCondition(float * const Phi, int nz, int ny, int nx) {

    #pragma omp parallel for collapse(3)
    for (int z = -1; z < nz+1; ++z) {
        for (int y = -1; y < ny+1; ++y) {
            for (int x = -1; x < nx+1; ++x) {

                int cx = std::min( std::max(x,0), nx-1 );
                int cy = std::min( std::max(y,0), ny-1 );
                int cz = std::min( std::max(z,0), nz-1 );

                Phi[gIndex(x,y,z,ny,nx,1)] = Phi[gIndex(cx,cy,cz,ny,nx,1)];
            }
        }
    }
}
