#ifndef CELLOBJECT_HPP
#define CELLOBJECT_HPP

#include <QDataStream>
#include <QFile>
#include <QTreeWidget>

#include "Debug.hpp"

enum CellObjectType
{
    COT_UNCATEGORIZED = 0,
    COT_MEMBRANE      = 1,
    COT_VESICLE       = 2,
    COT_FILAMENT      = 3,
    COT_SYNAPSE       = 4,
    COT_ENDOSOME      = 5,
    COT_ER            = 6
};

int typeToInt(CellObjectType type);

CellObjectType intToType(int i);


class CellObject : public QTreeWidgetItem
{
    Q_GADGET

public:
    CellObject(CellObjectType type = COT_UNCATEGORIZED,
               bool isRootItem = false);

    CellObject(QTreeWidget *parent,
               CellObjectType type = COT_UNCATEGORIZED,
               bool isRootItem = false);

    CellObject(CellObjectType type,
               bool isRootItem,
               unsigned int *x,
               unsigned int *y,
               unsigned int *z,
               unsigned int n,
               float *phasefield,
               int totalXDiff,
               int totalYDiff,
               int totalZDiff,
               int activeRegionDepth,
               int activeRegionHeight,
               int activeRegionWidth);

    CellObject(QTreeWidget *parent,
               CellObjectType type,
               bool isRootItem,
               unsigned int *x,
               unsigned int *y,
               unsigned int *z,
               unsigned int n,
               float *phasefield,
               int totalXDiff,
               int totalYDiff,
               int totalZDiff,
               int activeRegionDepth,
               int activeRegionHeight,
               int activeRegionWidth);

    void getCenter(float &center_x, float &center_y, float &center_z);

    void getCenter(int &center_x, int &center_y, int &center_z);

    void getBounds(int &beginX, int &endX,
                   int &beginY, int &endY,
                   int &beginZ, int &endZ);

    unsigned int *getX()
    {
        return x;
    }

    unsigned int *getY()
    {
        return y;
    }

    unsigned int *getZ()
    {
        return z;
    }

    unsigned int getN()
    {
        return n;
    }

    CellObjectType getType()
    {
        return type;
    }

    float *getPhiPtr()
    {
        return phi;
    }

    bool isRootObject()
    {
        return isRootItem;
    }

    // ostream, << overloading
    friend QDataStream &operator<<(QDataStream &out, const CellObject &object)
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
    }

    // istream, >> overloading
    friend QDataStream &operator>>(QDataStream &in, CellObject &object)
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
    }

private:
    CellObjectType type;

    bool isRootItem;

    unsigned int n;

    unsigned int *x;
    unsigned int *y;
    unsigned int *z;

    unsigned int xMin;
    unsigned int xMax;
    unsigned int yMin;
    unsigned int yMax;
    unsigned int zMin;
    unsigned int zMax;

    int phiN;
    float *phi;

    void setBounds();

    void setPhi(float *phasefield,
                int totalXDiff,
                int totalYDiff,
                int totalZDiff,
                int activeRegionDepth,
                int activeRegionHeight,
                int activeRegionWidth);

public:
    void getPhi(float *phasefield,
                int totalXDiff,
                int totalYDiff,
                int totalZDiff,
                int activeRegionDepth,
                int activeRegionHeight,
                int activeRegionWidth);
};

#endif // CELLOBJECT_HPP