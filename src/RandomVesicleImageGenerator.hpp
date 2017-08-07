#ifndef RANDOM_VESICLE_IMAGE_GENERATOR_HPP
#define RANDOM_VESICLE_IMAGE_GENERATOR_HPP

class RandomVesicleImageGenerator
{
public:
    RandomVesicleImageGenerator();

    RandomVesicleImageGenerator(float minRadiiA, float maxRadiiA,
                                float minRadiiB, float maxRadiiB,
                                float minRadiiC, float maxRadiiC,
                                bool randomRotation);

    void generateRandomVesicleImage(float *image,
                                    unsigned int &label,
                                    int width,
                                    int height,
                                    int depth,
                                    bool dense,
                                    float noise);

    void generateResponses(float *responseKernels, int width, int height, int depth, int featureSize);

private:
    void RandomRotationMatrix(float *R);

    void Transpose(float *A);

    void MatrixMult(float *A, float *B, float *C);

    float GetAlgebraicFormValue(float i, float j, float k,
                                float ci, float cj, float ck,
                                float *RdRt);

    void addNoise(float *image, int width, int height, int depth, float stddev);

    void Normalize(float *image, int width, int height, int depth);

    float minRadiiA;
    float maxRadiiA;
    float minRadiiB;
    float maxRadiiB;
    float minRadiiC;
    float maxRadiiC;
    bool randomRotation;
};

#endif // RANDOM_VESICLE_IMAGE_GENERATOR_HPP