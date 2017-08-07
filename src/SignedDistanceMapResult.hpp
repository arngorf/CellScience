#ifndef SIGNEDDISTANCEMAPRESULT_HPP
#define SIGNEDDISTANCEMAPRESULT_HPP




class SignedDistanceMapResult
{
public:

    SignedDistanceMapResult(int width, int height, int depth);

    void addSource();

    void addBoundary();

private:

    int width;
    int height;
    int depth;

    double *sdm;
};

#endif // SIGNEDDISTANCEMAPRESULT_HPP