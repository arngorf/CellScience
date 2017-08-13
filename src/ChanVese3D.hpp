#ifndef CHANVESE3D_HPP
#define CHANVESE3D_HPP

class ChanVese3D {
public:

    ChanVese3D();

    void segmentImage(const float * const image,
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
                      bool isConstc2);

    void setConstant_c1(bool isConstc1, float val = -1);

    void setConstant_c2(bool isConstc2, float val = -1);

    void initPhiBubble(float * const out, const int nx, const int ny, const int nz, const float sz, const float sy, const float sx, const float bubblePeriod);

    void initPhiSphere(float * const out, const int nz, const int ny, const int nx, const float sz, const float sy, const float sx, const float radius);

    void curvature(const float * const Phi, float * const C, const int nz, const int ny, const int nx);

private:

    float dx;
    float numericEps;
    float maxCurvature;
    bool isConstc1;
    bool isConstc2;
    float constc1;
    float constc2;

    inline float delta(float x, float eps);

    inline float heaviside(float x, float eps);

    float areaInside(const float * const image, const float * const Phi, int nx, int ny, int nz, float eps);

    float areaOutside(const float * const image, const float * const Phi,
                      int nx, int ny, int nz, float eps);

    void applyBoundaryCondition(float * const Phi, int nz, int ny, int nx);
};

#endif // CHANVESE3D_HPP