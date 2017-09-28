#ifndef VESICLE_SEGMENTATION_RESULT_HPP
#define VESICLE_SEGMENTATION_RESULT_HPP

#include "UtilityFunctions.hpp"
#include "ellipsoidHelpers.hpp"
#include "CellObject.hpp"

#include <vector>
#include <string>

#include <QDataStream>

class VesicleSegmentationResult
{
public:
    VesicleSegmentationResult();

    VesicleSegmentationResult(CellObject *object);

    void addSegmentingPoint(float x, float y, float z);

    bool isEqual(CellObject *object);

    std::string getNameID()
    {
        return name_ID;
    }

    void getCenter(int &x, int &y, int &z);

    void markSegmentationDirectionComplete(int direction);

    void OutputSegmentation();

    float GetCompletionValue()
    {
        return (float) completedSlices / 9.0;
    };

    Mat getBaseRotationMatrix()
    {
        return R;
    };

    // ostream, << overloading
    friend QDataStream &operator<<(QDataStream &out, const VesicleSegmentationResult &object)
    {
        out << object.x_center << object.y_center << object.z_center;

        out << object.completedSlices;

        for (int j = 0; j < 3; ++j)
        {
            for (int i = 0; i < 3; ++i)
            {
                out << object.R(j,i);
            }
        }

        out << object.n_segmentations;

        for (int i = 0; i < object.n_segmentations; ++i)
        {
            out << object.X[i];
        }
        for (int i = 0; i < object.n_segmentations; ++i)
        {
            out << object.Y[i];
        }
        for (int i = 0; i < object.n_segmentations; ++i)
        {
            out << object.Z[i];
        }

        return out;
    }

    // istream, >> overloading
    friend QDataStream &operator>>(QDataStream &in, VesicleSegmentationResult &object)
    {
        in >> object.x_center >> object.y_center >> object.z_center;


        object.name_ID = "V_" + STR(object.x_center)
                       + "_" + STR(object.y_center)
                       + "_" + STR(object.z_center);

        in >> object.completedSlices;

        Mat Rin(3,3);

        for (int j = 0; j < 3; ++j)
        {
            for (int i = 0; i < 3; ++i)
            {
                in >> Rin(j,i);
            }
        }

        object.R = Rin;

        in >> object.n_segmentations;

        float val;

        for (int i = 0; i < object.n_segmentations; ++i)
        {
            in >> val;
            object.X.push_back(val);
        }

        for (int i = 0; i < object.n_segmentations; ++i)
        {
            in >> val;
            object.Y.push_back(val);
        }

        for (int i = 0; i < object.n_segmentations; ++i)
        {
            in >> val;
            object.Z.push_back(val);
        }

        return in;
    }

private:
    std::string name_ID;

    int x_center;
    int y_center;
    int z_center;

    int n_segmentations;

    std::vector<float> X;
    std::vector<float> Y;
    std::vector<float> Z;

    int completedSlices;

    Mat R;
};

#endif // VESICLE_SEGMENTATION_RESULT_HPP