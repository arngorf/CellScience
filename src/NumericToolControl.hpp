#ifndef NUMERICTOOLCONTROL_HPP
#define NUMERICTOOLCONTROL_HPP

#include <QDockWidget>
#include <QObject>
#include <QTreeWidget>

#include "CellObject.hpp"
#include "ChanVese3D.hpp"
#include "MainWindow.h"
#include "TiffReader.hpp"
#include "UtilityFunctions.hpp"
#include "VesicleSegmentationWindow.h"
#include "VesicleSegmentationResult.hpp"
#include "VesicleSegmentationResultManager.hpp"

class QAction;

void sort_indexes(const std::vector<float> &v, std::vector<int> &idx);

class NumericToolControl : public QObject
{
    Q_OBJECT
public:
    NumericToolControl(char const path[]);

    ~NumericToolControl();

    int run(MainWindow *mainWin, int argc, char *argv[]);

private:
    void updateViewImage(bool activeRegionOnly = false);

    void updateViewSegmented();

    void updateViewPhi();

    void updateViewSignedMap();

    void createCellTree();

    void loadObjectsFromFile();

    void saveObjectsToFile();

    void updateSegmented();

    QTreeWidgetItem *addTreeRoot(QString name, CellObjectType type);

    TiffImageRef *tiffImage;
    ChanVese3D cv3d;
    MainWindow *mainWin;
    VesicleSegmentationWindow *vesicleSegmentationWindow;

    QTreeWidget *cellTree;

    int imageBeginX;
    int imageEndX;
    int imageBeginY;
    int imageEndY;
    int imageBeginZ;
    int imageEndZ;

    // Cell Tree root items
    QTreeWidgetItem *cellMembranes;
    QTreeWidgetItem *cellVesicles;
    QTreeWidgetItem *cellFilaments;
    QTreeWidgetItem *cellUncategorized;
    QTreeWidgetItem *cellSynapses;
    QTreeWidgetItem *cellEndosomes;
    QTreeWidgetItem *cellERs;

    int cellMembraneCount;
    int cellVesicleCount;
    int cellFilamentCount;
    int cellUncategorizedCount;
    int cellSynapseCount;
    int cellEndosomeCount;
    int cellERCount;

    std::vector<CellObject*> sdBoundaryObjects;
    std::vector<CellObject*> sdSourceObjects;
    std::vector<CellObject*> sdPassThroughObjects;

    VesicleSegmentationResultManager *vesicleSegmentationResultManager;

    const int signedDistanceIterationFactor;

    float cytosolI;

    float travelWeight;

    int currentPage;

    bool activeRegion;

    int activeRegionOriginX;
    int activeRegionOriginY;
    int activeRegionOriginZ;
    int activeRegionWidth;
    int activeRegionHeight;
    int activeRegionDepth;

    int sdOriginX;
    int sdOriginY;
    int sdOriginZ;
    int sdWidth;
    int sdHeight;
    int sdDepth;

    float *activeImage;
    float *activeImageOriginal;
    bool *activeFlagArray;
    bool *segmented;
    float *phi;
    float *signedDistanceMap;
    int *filamentContour;

public slots:

    /*************************************************************************
    *****  MainWindow Slots  *************************************************
    *************************************************************************/

    void nextImage();

    void prevImage();

    void gotoPageSlot(int newZ);

    void newActiveRegion(int beginX, int endX, int beginY, int endY, int beginZ, int endZ, float sigma);

    void setActiveRegion(int beginX, int endX, int beginY, int endY, int beginZ, int endZ, float sigma);

    void runChanVeseSlot(int numIterations,
                         float mu, float nu,
                         float c1, float c2,
                         bool constc1, bool constc2);

    void markVesicleLocationSlot(int x, int y);

    void markNeurofilamentSlot(int x, int y);

    void finalizeNeurofilamentSlot();

    void setPosSlot(int beginX, int endX, int beginY, int endY, int beginZ, int endZ);

    void setNegSlot(int beginX, int endX, int beginY, int endY, int beginZ, int endZ);

    void setDisallowedSlot(int beginX, int endX, int beginY, int endY, int beginZ, int endZ);

    void swapPhiSignSlot();

    void reinitPhiSlot();

    void renderIn3DSlot();

    void storeObjectSlot();

    void loadObjectSlot();

    void deleteObjectSlot();

    void startVesicleSegmentationSlot();

    void addVesiclesToSegmentationPoolSlot();

    void blightSlot(int imageX, int imageY, float sigma);

    void restoreSlot(int imageX, int imageY, float sigma);

    void SDSetSourcesSlot();

    void SDSetBoundariesSlot();

    void SDSetPassThroughSlot(float travelCost);

    void SDRunBoundedSignedDistanceMapSlot();

    void SDReadDataSlot();

    void SDReadFullDataSlot();

    void ripleysKFunctionSlot();

    void runCustomFunctionSlot();

    /*************************************************************************
    *****  VesicleSegmentationWindow Slots  **********************************
    *************************************************************************/

    void storeSegmentationSlot(std::vector<float> x, std::vector<float> y);

    void nextSlot();

private:

    void cleanSignedDistanceMap(int width,
                                int height,
                                int depth,
                                float *signedDistanceMap);

    void normalizeImageToRange(float *image,
                               float &phiZero,
                               const int width,
                               const int height,
                               const float imageMin,
                               const float imageMax,
                               const float boundMin,
                               const float boundMax);

    void cosineModulation(float *image, int width, int height, float period);

    void restoreImageLine(float beginX,
                          float endX,
                          float beginY,
                          float endY,
                          float beginZ,
                          float endZ);

    void loadMatlabMembraneSegmentation();

    void trainAndPredictForVesicleCalssification();

    void trainAndPredictForVesicleCalssification2();

    void ffmTest();
};

#endif // NUMERICTOOLCONTROL_HPP