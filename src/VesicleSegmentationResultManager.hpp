#ifndef VESICLE_SEGMENTATION_RESULT_MANAGER_HPP
#define VESICLE_SEGMENTATION_RESULT_MANAGER_HPP

#include "TiffReader.hpp"

#include "CellObject.hpp"
#include "VesicleSegmentationResult.hpp"

#include <Eigen/Dense>

#include <string>

class VesicleSegmentationResultManager
{
public:
    VesicleSegmentationResultManager(char const path[], int width, int height);

    ~VesicleSegmentationResultManager();

    /*
     * Grabs next vesicle + rotation combination */
    void NextSlice();

    float GetCompletionValue();

    float GetMinCompletionValue();

    float *GetImageSlice();

    void SetSegmentationResult(std::vector<float> x, std::vector<float> y);

    void AddVesicle(CellObject *object);

private:
    std::string imagePath;

    TiffImage *tiffImage;

    std::vector<VesicleSegmentationResult*> vesicleResultArray;

    int completedSlices;

    int numRotations;

    int curVesicleIndex;
    int curDirectionIndex;

    Mat currentR;
    Mat currentRInverse;

    int width;
    int height;

    int cxImg;
    int cyImg;
    int czImg;

    int x_center;
    int y_center;
    int z_center;

    int xBegin;
    int yBegin;
    int zBegin;
    void loadSegmentationResultsFromFile();

    void saveSegmentationResultsToFile();

    void generateRotation(int direction);

    void RotatePoint(float &x, float &y, float &z, const Mat &R);

    void RotatePointAroundPoint(float &x, float &y, float &z, const float cx, const float cy, const float cz, const Mat &R);

    void TranslatePointsToGlobal(float &x, float &y, float &z);
};

#endif // VESICLE_SEGMENTATION_RESULT_MANAGER_HPP