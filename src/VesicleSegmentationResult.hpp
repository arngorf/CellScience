#ifndef VESICLE_SEGMENTATION_RESULT_HPP
#define VESICLE_SEGMENTATION_RESULT_HPP

#include <vector>

class VesicleSegmentationResult
{
public:
    VesicleSegmentationResult();

    // ostream, << overloading
    /*friend QDataStream &operator<<(QDataStream &out, const CellObject &object)
    {
        int typeN = typeToInt(object.type);

        out << typeN << object.n;
        out << object.xMin << object.xMax
            << object.yMin << object.yMax
            << object.zMin << object.zMax;

        for (int i = 0; i < object.n; ++i)
        {
            out << object.x[i] << object.y[i] << object.z[i];
        }

        out << object.phiN;

        for (int i = 0; i < object.phiN; ++i)
        {
            out << (double) object.phi[i];
        }

        return out;
    }*/

    // istream, >> overloading
    /*friend QDataStream &operator>>(QDataStream &in, CellObject &object)
    {
        object.isRootItem = false;

        int typeN;

        in >> typeN >> object.n;
        in >> object.xMin >> object.xMax
           >> object.yMin >> object.yMax
           >> object.zMin >> object.zMax;

        object.x = (unsigned int*) malloc(object.n * sizeof(unsigned int));
        object.y = (unsigned int*) malloc(object.n * sizeof(unsigned int));
        object.z = (unsigned int*) malloc(object.n * sizeof(unsigned int));

        for (int i = 0; i < object.n; ++i)
        {
            in >> object.x[i] >> object.y[i] >> object.z[i];
        }

        in >> object.phiN;

        object.phi = (float *) malloc(object.phiN * sizeof(float));

        if (object.phi == NULL)
        {
            Debug::Error("CellObject::operator>>: Could not malloc phi");
        }

        for (int i = 0; i < object.phiN; ++i)
        {
            double tmp;
            in >> tmp;
            object.phi[i] = (float) tmp;
        }

        object.type = intToType(typeN);

        return in;
    }*/

private:
    int x_center;
    int y_center;

    int n_directions;
    std::vector<float*> X;
    std::vector<int> n_segmentations;
};

#endif // VESICLE_SEGMENTATION_RESULT_HPP