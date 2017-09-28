#include "VesicleSegmentationResult.hpp"

#include "Debug.hpp"

#include <stdlib.h>
#include <iostream>

VesicleSegmentationResult::VesicleSegmentationResult()
                                            : n_segmentations(0)
                                            , X(0)
                                            , Y(0)
                                            , completedSlices(0)
{

}

VesicleSegmentationResult::VesicleSegmentationResult(CellObject *object)
                                            : n_segmentations(0)
                                            , X(0)
                                            , Y(0)
                                            , completedSlices(0)
{
    object->getCenter(x_center, y_center, z_center);

    R = RandomRotationMatrix();

    name_ID = "V_" + STR(x_center) + "_" + STR(y_center) + "_" + STR(z_center);
}

void VesicleSegmentationResult::addSegmentingPoint(float x, float y, float z)
{
    ++n_segmentations;
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);

    Debug::Info("Adding segmenting point: " + STR(x) + ", " + STR(y) + ", " + STR(z));
}

bool VesicleSegmentationResult::isEqual(CellObject *object)
{
    int other_center_x, other_center_y, other_center_z;

    object->getCenter(other_center_x, other_center_y, other_center_z);

    return other_center_x == x_center and
           other_center_y == y_center and
           other_center_z == z_center;
}

void VesicleSegmentationResult::getCenter(int &x, int &y, int &z)
{
    x = x_center;
    y = y_center;
    z = z_center;
}

void VesicleSegmentationResult::markSegmentationDirectionComplete(int direction)
{

}

void VesicleSegmentationResult::OutputSegmentation()
{
    std::cout << name_ID << " = [[";
    for (int i = 0; i < n_segmentations; ++i)
    {
        std::cout << X[i];
        if (i == n_segmentations - 1) std::cout << "]";
        else std::cout << ",";
    }
    std::cout << ",[";
    for (int i = 0; i < n_segmentations; ++i)
    {
        std::cout << Y[i];
        if (i == n_segmentations - 1) std::cout << "]";
        else std::cout << ",";
    }
    std::cout << ",[";
    for (int i = 0; i < n_segmentations; ++i)
    {
        std::cout << Z[i];
        if (i == n_segmentations - 1) std::cout << "]";
        else std::cout << ",";
    }
    std::cout << "]" << std::endl;
}