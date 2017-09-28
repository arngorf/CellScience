#include "CellObject.hpp"

#include "UtilityFunctions.hpp"

int typeToInt(CellObjectType type)
{
    int i = 0;

    switch (type)
    {
    case COT_UNCATEGORIZED:
        i = 0;
        break;
    case COT_MEMBRANE:
        i = 1;
        break;
    case COT_VESICLE:
        i = 2;
        break;
    case COT_FILAMENT:
        i = 3;
        break;
    case COT_SYNAPSE:
        i = 4;
        break;
    case COT_ENDOSOME:
        i = 5;
        break;
    case COT_ER:
        i = 6;
        break;
    }

    return i;
}

CellObjectType intToType(int i)
{
    CellObjectType type = COT_UNCATEGORIZED;

    switch (i)
    {
    case 0:
        type = COT_UNCATEGORIZED;
        break;
    case 1:
        type = COT_MEMBRANE;
        break;
    case 2:
        type = COT_VESICLE;
        break;
    case 3:
        type = COT_FILAMENT;
        break;
    case 4:
        type = COT_SYNAPSE;
        break;
    case 5:
        type = COT_ENDOSOME;
        break;
    case 6:
        type = COT_ER;
        break;
    }

    return type;
}

CellObject::CellObject(CellObjectType type,
                       bool isRootItem)
                       : QTreeWidgetItem()
                       , type(type)
                       , isRootItem(isRootItem)
{

}

CellObject::CellObject(QTreeWidget *parent,
                       CellObjectType type,
                       bool isRootItem)
                       : QTreeWidgetItem(parent)
                       , type(type)
                       , isRootItem(isRootItem)
{

}

CellObject::CellObject(CellObjectType type,
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
                       int activeRegionWidth)
                       : QTreeWidgetItem()
                       , type(type)
                       , isRootItem(isRootItem)
                       , x(x)
                       , y(y)
                       , z(z)
                       , n(n)
                       , phiN(0)
                       , phi(NULL)
{
    setBounds();

    isRootItem = false;

    if (phasefield != NULL)
    {
    setPhi(phasefield,
           totalXDiff,
           totalYDiff,
           totalZDiff,
           activeRegionDepth,
           activeRegionHeight,
           activeRegionWidth);
    }
}

CellObject::CellObject(QTreeWidget *parent,
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
                       int activeRegionWidth)
                       : QTreeWidgetItem(parent)
                       , type(type)
                       , isRootItem(isRootItem)
                       , x(x)
                       , y(y)
                       , z(z)
                       , n(n)
                       , phi(NULL)
{
    setBounds();

    isRootItem = false;

    if (phasefield != NULL)
    {
        setPhi(phasefield,
               totalXDiff,
               totalYDiff,
               totalZDiff,
               activeRegionDepth,
               activeRegionHeight,
               activeRegionWidth);
    }

}

void CellObject::getCenter(float &center_x, float &center_y, float &center_z)
{
    if (n <= 0) return;

    center_x = 0;
    center_y = 0;
    center_z = 0;

    for (int i = 0; i < n; ++i)
    {
        center_x += x[i];
        center_y += y[i];
        center_z += z[i];
    }

    center_x /= n;
    center_y /= n;
    center_z /= n;
}

void CellObject::getCenter(int &center_x, int &center_y, int &center_z)
{
    if (n <= 0) return;

    float float_center_x, float_center_y, float_center_z;

    getCenter(float_center_x, float_center_y, float_center_z);

    center_x = (int) float_center_x;
    center_y = (int) float_center_y;
    center_z = (int) float_center_z;
}

void CellObject::getBounds(int &beginX, int &endX,
                           int &beginY, int &endY,
                           int &beginZ, int &endZ)
{
    //Debug::Info("CellObject::getBounds: Entering");

    beginX = xMin;
    endX = xMax;
    beginY = yMin;
    endY = yMax;
    beginZ = zMin;
    endZ = zMax;

    //Debug::Info("CellObject::getBounds: Leaving");
}

void CellObject::setBounds()
{
    Debug::Info("CellObject::setBounds: Entering");

    xMin = x[0];
    xMax = x[0];
    yMin = y[0];
    yMax = y[0];
    zMin = z[0];
    zMax = z[0];

    for (int i = 0; i < n; ++i)
    {
        xMin = std::min(xMin, x[i]);
        xMax = std::max(xMax, x[i]);
        yMin = std::min(yMin, y[i]);
        yMax = std::max(yMax, y[i]);
        zMin = std::min(zMin, z[i]);
        zMax = std::max(zMax, z[i]);
    }

    Debug::Info("CellObject::setBounds: Leaving");
}

void CellObject::setPhi(float *phasefield,
                        int totalXDiff,
                        int totalYDiff,
                        int totalZDiff,
                        int activeRegionDepth,
                        int activeRegionHeight,
                        int activeRegionWidth)
{
    Debug::Info("CellObject::setPhi: Entering");

    int width = xMax - xMin + 3;
    int height = yMax - yMin + 3;
    int depth = zMax - zMin + 3;

    phiN = width * height * depth;

    if (phi != NULL)
    {
        Debug::Warning("CellObject::setPhi: phi is being overwritten");

        free(phi);
    }

    if (phasefield == NULL)
    {
        Debug::Error("CellObject::setPhi: phasefield is NULL");

        return;
    }

    phi = (float *) malloc(phiN * sizeof(float));

    for (int z = 0; z < depth; ++z)
    {
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                int gx = x + xMin - totalXDiff - 1;
                int gy = y + yMin - totalYDiff - 1;
                int gz = z + zMin - totalZDiff - 1;

                int localIndex = gIndex(x, y, z, height, width, 0);
                int globalIndex = gIndex(gx,
                                         gy,
                                         gz,
                                         activeRegionHeight,
                                         activeRegionWidth,
                                         1);

                if (gx < 0 or gy < 0 or gz < 0 or
                    gx >= activeRegionWidth or
                    gy >= activeRegionHeight or
                    gz >= activeRegionDepth)
                {
                    phi[localIndex] = -1;
                }
                else
                {
                    phi[localIndex] = phasefield[globalIndex];
                }
            }
        }
    }

    Debug::Info("CellObject::setPhi: Leaving");
}

void CellObject::getPhi(float *phasefield,
                        int totalXDiff,
                        int totalYDiff,
                        int totalZDiff,
                        int activeRegionDepth,
                        int activeRegionHeight,
                        int activeRegionWidth)
{
    Debug::Info("CellObject::getPhi: Entering");

    int width = xMax - xMin + 3;
    int height = yMax - yMin + 3;
    int depth = zMax - zMin + 3;

    if (phi == NULL)
    {
        Debug::Error("CellObject::getPhi: phi is NULL");

        return;
    }
    if (phasefield == NULL)
    {
        Debug::Error("CellObject::getPhi: phasefield is NULL");

        return;
    }

    for (int z = 0; z < depth; ++z)
    {
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                int gx = x + xMin - totalXDiff - 1;
                int gy = y + yMin - totalYDiff - 1;
                int gz = z + zMin - totalZDiff - 1;

                int localIndex = gIndex(x, y, z, height, width, 0);
                int globalIndex = gIndex(gx,
                                         gy,
                                         gz,
                                         activeRegionHeight,
                                         activeRegionWidth,
                                         1);

                phasefield[globalIndex] = phi[localIndex];
            }
        }
    }

    Debug::Info("CellObject::getPhi: Leaving");
}