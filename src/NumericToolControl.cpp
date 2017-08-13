#include "NumericToolControl.hpp"

#include <QApplication>
#include <QDataStream>
#include <QFile>
#include <QtWidgets>

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "CellPlot.hpp"
#include "Debug.hpp"
#include "ellipsoidHelpers.hpp"
#include "FastMarching.hpp"
#include "GaussSmooth.hpp"
#include "RandomVesicleImageGenerator.hpp"
#include "SignedDistanceMap.hpp"

bool cellObjectLexiographicCompare (QTreeWidgetItem *a,QTreeWidgetItem *b)
{
    return a->text(0) < b->text(0);
}

void sort_indexes(const std::vector<float> &v, std::vector<int> &idx)
{
    bool flag = false;
    iota(idx.begin(), idx.end(), 0);
    for (int i = 0; i < idx.size(); ++i) {
        if (idx[i] < 0 or idx[i] >= v.size()) {
            flag = true;
        }
    }
    if (flag) for (int i = 0; i < idx.size(); ++i) Debug::Info(STR(idx[i]));
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    if (flag) for (int i = 0; i < idx.size(); ++i) Debug::Info(STR(idx[i]));
    if (flag) exit(0);
}

NumericToolControl::NumericToolControl(char const path[]) :
                                       mainWin(NULL)
                                     , cellTree(new QTreeWidget)
                                     , imageBeginX(0)
                                     , imageEndX(2174) // 1000
                                     , imageBeginY(0) // 100
                                     , imageEndY(1854) // 600
                                     , imageBeginZ(0)
                                     , imageEndZ(1065) // 800
                                     , cellMembraneCount(0)
                                     , cellVesicleCount(0)
                                     , cellFilamentCount(0)
                                     , cellUncategorizedCount(0)
                                     , cellSynapseCount(0)
                                     , cellEndosomeCount(0)
                                     , cellERCount(0)
                                     , sdBoundaryObjects(0)
                                     , sdSourceObjects(0)
                                     , sdPassThroughObjects(0)
                                     , signedDistanceIterationFactor(4)
                                     , cytosolI(172)
                                     , travelWeight(1)
                                     , currentPage(0)
                                     , activeRegion(false)
                                     , activeImage(NULL)
                                     , activeImageOriginal(NULL)
                                     , activeFlagArray(NULL)
                                     , segmented(NULL)
                                     , phi(NULL)
                                     , signedDistanceMap(NULL)
                                     , filamentContour(NULL)
{
    Debug::Info("NumericToolControl::NumericToolControl: Entering");

    tiffImage = new TiffImageRef(path, imageEndX, imageEndY, imageEndZ);

    cellTree->setSelectionMode(QAbstractItemView::MultiSelection);

    Debug::Info("NumericToolControl::NumericToolControl: Leaving");
}

NumericToolControl::~NumericToolControl()
{
    Debug::Info("NumericToolControl::~NumericToolControl: Entering");

    saveObjectsToFile();

    Debug::Info("NumericToolControl::~NumericToolControl: Leaving");
}

int NumericToolControl::run(MainWindow *mainWin, int argc, char *argv[])
{
    Debug::Info("NumericToolControl::run: Entering");

    for (int i = 1; i < argc; ++i)
    {
        Debug::Warning("Warning: running program with argument '"
                       + std::string(argv[i]) + "' is not supported");
    }

    this->mainWin = mainWin;

    unsigned int width = tiffImage->getWidth();
    unsigned int height = tiffImage->getHeight();
    unsigned int depth = tiffImage->getDepth();

    int N = width * height * depth;

    float *imageSlice = tiffImage->getSlice(0);
    //segmented = (bool *) malloc(N * sizeof(bool));
    segmented = (bool *) malloc(width*height * sizeof(bool));

    if (segmented == NULL)
    {
        Debug::Error("NumericToolControl::run: Malloc failed (segmented)");
    }

    for (int i = 0; i < N; ++i)
    {
        segmented[i] = false;
    }

    mainWin->InitViewArea(width, height, depth);
    mainWin->setFillImage(width, height, 0, imageSlice);

    createCellTree();

    QObject::connect(mainWin, SIGNAL(nextImageSig()),
                     this,    SLOT(nextImage()));

    QObject::connect(mainWin, SIGNAL(prevImageSig()),
                     this,    SLOT(prevImage()));

    QObject::connect(mainWin, SIGNAL(newActiveRegionSig(int, int, int, int, int, int, float)),
                     this,    SLOT(newActiveRegion(int, int, int, int, int, int, float)));

    QObject::connect(mainWin, SIGNAL(runChanVeseSig(int, float, float, float, float, bool, bool)),
                     this,    SLOT(runChanVeseSlot(int, float, float, float, float, bool, bool)));

    QObject::connect(mainWin, SIGNAL(markVesicleLocationSig(int, int)),
                     this,    SLOT(markVesicleLocationSlot(int, int)));

    QObject::connect(mainWin, SIGNAL(markNeurofilamentSig(int, int)),
                     this,    SLOT(markNeurofilamentSlot(int, int)));

    QObject::connect(mainWin, SIGNAL(setPosSig(int, int, int, int, int, int)),
                     this,    SLOT(setPosSlot(int, int, int, int, int, int)));

    QObject::connect(mainWin, SIGNAL(setNegSig(int, int, int, int, int, int)),
                     this,    SLOT(setNegSlot(int, int, int, int, int, int)));

    QObject::connect(mainWin, SIGNAL(setDisallowedSig(int, int, int, int, int, int)),
                     this,    SLOT(setDisallowedSlot(int, int, int, int, int, int)));

    QObject::connect(mainWin, SIGNAL(swapPhiSignSig()),
                     this,    SLOT(swapPhiSignSlot()));

    QObject::connect(mainWin, SIGNAL(reinitPhiSig()),
                     this,    SLOT(reinitPhiSlot()));

    QObject::connect(mainWin, SIGNAL(renderIn3DSig()),
                     this,    SLOT(renderIn3DSlot()));

    QObject::connect(mainWin, SIGNAL(storeObjectSig()),
                     this,    SLOT(storeObjectSlot()));

    QObject::connect(mainWin, SIGNAL(loadObjectSig()),
                     this,    SLOT(loadObjectSlot()));

    QObject::connect(mainWin, SIGNAL(deleteObjectSig()),
                     this,    SLOT(deleteObjectSlot()));

    QObject::connect(mainWin, SIGNAL(blightSig(int, int, float)),
                     this,    SLOT(blightSlot(int, int, float)));

    QObject::connect(mainWin, SIGNAL(restoreSig(int, int, float)),
                     this,    SLOT(restoreSlot(int, int, float)));

    QObject::connect(mainWin, SIGNAL(SDSetSourcesSig()),
                     this,    SLOT(SDSetSourcesSlot()));

    QObject::connect(mainWin, SIGNAL(SDSetBoundariesSig()),
                     this,    SLOT(SDSetBoundariesSlot()));

    QObject::connect(mainWin, SIGNAL(SDSetPassThroughSig(float)),
                     this,    SLOT(SDSetPassThroughSlot(float)));

    QObject::connect(mainWin, SIGNAL(SDRunBoundedSignedDistanceMapSig()),
                     this,    SLOT(SDRunBoundedSignedDistanceMapSlot()));

    QObject::connect(mainWin, SIGNAL(SDReadDataSig()),
                     this,    SLOT(SDReadDataSlot()));

    QObject::connect(mainWin, SIGNAL(SDReadFullDataSig()),
                     this,    SLOT(SDReadFullDataSlot()));

    QObject::connect(mainWin, SIGNAL(ripleysKFunctionSig()),
                     this,    SLOT(ripleysKFunctionSlot()));

    QObject::connect(mainWin, SIGNAL(runCustomFunctionSig()),
                     this,    SLOT(runCustomFunctionSlot()));

    QObject::connect(mainWin, SIGNAL(finalizeNeurofilamentSig()),
                     this,    SLOT(finalizeNeurofilamentSlot()));

    loadObjectsFromFile();

    Debug::Info("NumericToolControl::run: Leaving");

    return 0;
}

void NumericToolControl::updateViewImage(bool activeRegionOnly)
{
    Debug::Info("NumericToolControl::updateViewImage: Entering");

    unsigned int width = tiffImage->getWidth(); // imageEndX - imageBeginX - 1;
    unsigned int height = tiffImage->getHeight(); // imageEndY - imageBeginY - 1;

    float *imageSlice = NULL;
    if (not activeRegionOnly)
    {
        imageSlice = tiffImage->getSlice(currentPage);
    }

    if (activeRegion)
    {
        int k = currentPage - activeRegionOriginZ;
        float *activeImageSlice = &activeImage[gIndex(0, 0, k,
                                                       activeRegionHeight,
                                                       activeRegionWidth,
                                                       0)];

        mainWin->setFillImage(width, height, currentPage, imageSlice,
                              activeImageSlice, activeRegionOnly);
    }
    else
    {
        if (not activeRegionOnly)
        {
            mainWin->setFillImage(width, height, currentPage, imageSlice);
        }
    }

    Debug::Info("NumericToolControl::updateViewImage: Leaving");
}

void NumericToolControl::updateViewSegmented()
{
    Debug::Info("NumericToolControl::updateViewSegmented: Entering");

    if (segmented != NULL)
    {
        bool *segmentedSlice = &segmented[0];

        mainWin->updateSegmented(segmentedSlice);
    }
    else
    {
        Debug::Warning("NumericToolControl::updateViewSegmented: segmented flag array is NULL");
    }

    Debug::Info("NumericToolControl::updateViewSegmented: Leaving");
}

void NumericToolControl::updateViewPhi()
{
    Debug::Info("NumericToolControl::updateViewPhi: Entering");

    if (phi != NULL) {
        if (currentPage >= activeRegionOriginZ and currentPage < activeRegionOriginZ + activeRegionDepth) {

            float *phiSlice = (float*) malloc(activeRegionWidth * activeRegionHeight * sizeof(float));

            if (phiSlice == NULL)
            {
                Debug::Error("NumericToolControl::updateViewPhi: Malloc failed (phiSlice)");
            }

            #pragma omp parallel for collapse(2)
            for (int j = 0; j < activeRegionHeight; ++j)
            {
                for (int i = 0; i < activeRegionWidth; ++i)
                {
                    int k = currentPage - activeRegionOriginZ;

                    phiSlice[j*activeRegionWidth + i] = phi[gIndex(i, j, k, activeRegionHeight, activeRegionWidth, 1)];
                }
            }

            float phiMax = phi[0];
            float phiMin = phi[0];

            int phiN = (activeRegionWidth + 2)
                     * (activeRegionHeight + 2)
                     * (activeRegionDepth + 2);

            for (int i = 1; i < phiN; ++i)
            {
                phiMin = std::min(phiMin, phi[i]);
                phiMax = std::max(phiMax, phi[i]);
            }

            float phiZero;
            normalizeImageToRange(phiSlice, phiZero, activeRegionWidth, activeRegionHeight, phiMin, phiMax, 0, 255);

            mainWin->updatePhi(phiSlice, phiZero);

            free(phiSlice);
        } else {
            mainWin->updatePhi(NULL);
        }
    }

    Debug::Info("NumericToolControl::updateViewPhi: Leaving");
}

void NumericToolControl::updateViewSignedMap()
{
    Debug::Info("NumericToolControl::updateViewSignedMap: Entering");

    if (signedDistanceMap != NULL) {

        if (currentPage >= sdOriginZ and currentPage < sdOriginZ + sdDepth) {

            float *signedSlice = (float*) malloc(sdWidth * sdHeight * sizeof(float));

            if (signedSlice == NULL) {
                Debug::Error("NumericToolControl::updateViewSignedMap: Malloc failed (signedSlice)");
            }

            #pragma omp parallel for collapse(2)
            for (int j = 0; j < sdHeight; ++j)
            {
                for (int i = 0; i < sdWidth; ++i)
                {
                    int k = currentPage - sdOriginZ;
                    signedSlice[j*sdWidth + i] = signedDistanceMap[gIndex(i, j, k, sdHeight, sdWidth, 0)];
                }
            }

            float signedMax = signedSlice[0];
            float signedMin = signedSlice[0];

            int signedN = (sdWidth)
                        * (sdHeight)
                        * (sdDepth);

            for (int i = 1; i < signedN; ++i)
            {
                signedMin = std::min(signedMin, signedDistanceMap[i]);
                signedMax = std::max(signedMax, signedDistanceMap[i]);
            }

            float signedZero;

            normalizeImageToRange(signedSlice, signedZero, sdWidth, sdHeight, signedMin, signedMax, 0, 255);

            mainWin->updateSignedMap(signedSlice, sdOriginX, sdOriginY, sdWidth, sdHeight);

            free(signedSlice);
        } else {
            mainWin->updateSignedMap(NULL, sdOriginX, sdOriginY, sdWidth, sdHeight);
        }
    }

    Debug::Info("NumericToolControl::updateViewSignedMap: Leaving");
}

void NumericToolControl::createCellTree()
{
    Debug::Info("NumericToolControl::createCellTree: Entering");

    cellTree->setColumnCount(2);
    cellTree->setHeaderLabels(QStringList() << "Name" << "Voxels");

    // Add root nodes
    cellUncategorized = addTreeRoot("Uncategorized", COT_UNCATEGORIZED);
    cellMembranes = addTreeRoot("Membrane", COT_MEMBRANE);
    cellVesicles = addTreeRoot("Vesicles", COT_VESICLE);
    cellFilaments = addTreeRoot("Filaments", COT_FILAMENT);
    cellSynapses = addTreeRoot("Synapses", COT_SYNAPSE);
    cellEndosomes = addTreeRoot("Endosomes", COT_ENDOSOME);
    cellERs = addTreeRoot("ERs", COT_ER);

    if (mainWin != NULL) {
        mainWin->createTree(cellTree);
    } else {
        Debug::Error("NumericToolControl::createCellTree: mainWin not "
                     "initialized, could not set cell tree structure.");
    }

    Debug::Info("NumericToolControl::createCellTree: Leaving");
}

void NumericToolControl::loadObjectsFromFile()
{
    Debug::Info("NumericToolControl::loadObjectsFromFile: Entering");

    QString filename = "cellObjects.data";
    QFile file(filename);

    if(!file.open(QIODevice::ReadOnly))
    {
        Debug::Warning("NumericToolControl::loadObjectsFromFile: Could not "
                       "open cell object file");
    }

    QDataStream in(&file);

    int n;

    in >> n;

    Debug::Info("NumericToolControl::loadObjectsFromFile: Loading " + STR(n) +
                " objects");

    for (int i = 0; i < n; ++i)
    {
        CellObject *object = new CellObject();

        in >> (*object);

        object->setText(1, QString::number(object->getN()));

        switch (object->getType())
        {
        case COT_UNCATEGORIZED:

            object->setText(0, QString("U")
                               + QString::number(cellUncategorizedCount));
            cellUncategorized->addChild(object);
            ++cellUncategorizedCount;
            break;

        case COT_MEMBRANE:

            object->setText(0, QString("M")
                               + QString::number(cellMembraneCount));
            cellMembranes->addChild(object);
            ++cellMembraneCount;
            break;

        case COT_VESICLE:

            object->setText(0, QString("V")
                               + QString::number(cellVesicleCount));
            cellVesicles->addChild(object);
            ++cellVesicleCount;
            break;

        case COT_FILAMENT:

            object->setText(0, QString("F")
                               + QString::number(cellFilamentCount));
            cellFilaments->addChild(object);
            ++cellFilamentCount;
            break;

        case COT_SYNAPSE:

            object->setText(0, QString("S")
                               + QString::number(cellSynapseCount));
            cellSynapses->addChild(object);
            ++cellSynapseCount;
            break;

        case COT_ENDOSOME:

            object->setText(0, QString("En")
                               + QString::number(cellEndosomeCount));
            cellEndosomes->addChild(object);
            ++cellEndosomeCount;
            break;

        case COT_ER:

            object->setText(0, QString("ER")
                               + QString::number(cellERCount));
            cellERs->addChild(object);
            ++cellERCount;
            break;
        }
    }


    file.close();

    updateSegmented();

    Debug::Info("NumericToolControl::loadObjectsFromFile: Leaving");
}

void NumericToolControl::saveObjectsToFile()
{
    Debug::Info("NumericToolControl::saveObjectsToFile: Entering");

    QString filename = "cellObjects.data";
    QFile file(filename);

    if(!file.open(QIODevice::WriteOnly))
    {
        Debug::Warning("NumericToolControl::saveObjectsToFile:Could not "
                       "open cell object file");
    }

    int n = 0;

    for (int j = 0; j < cellTree->topLevelItemCount(); ++j)
    {
        n += cellTree->topLevelItem(j)->childCount();
    }

    QDataStream out(&file);

    out << n;

    for (int j = 0; j < cellTree->topLevelItemCount(); ++j)
    {
        CellObject *topLevelObject = (CellObject *) cellTree->topLevelItem(j);

        for (int i = 0; i < topLevelObject->childCount(); ++i)
        {
            CellObject *object = (CellObject *) topLevelObject->child(i);
            out << (*object);
        }
    }

    file.flush();
    file.close();

    Debug::Info("NumericToolControl::saveObjectsToFile: Leaving");
}

void NumericToolControl::updateSegmented()
{
    Debug::Info("NumericToolControl::updateSegmented: Entering");

    unsigned int width = tiffImage->getWidth();
    unsigned int height = tiffImage->getHeight();

    int N = width * height;

    for (int i = 0; i < N; ++i)
    {
        segmented[i] = false;
    }

    for (int j = 0; j < cellTree->topLevelItemCount(); ++j)
    {
        CellObject *topLevelObject = (CellObject *) cellTree->topLevelItem(j);

        for (int i = 0; i < topLevelObject->childCount(); ++i)
        {
            CellObject *object = (CellObject *) topLevelObject->child(i);

            int beginX, endX, beginY, endY, beginZ, endZ;

            object->getBounds(beginX, endX,
                              beginY, endY,
                              beginZ, endZ);

            if (currentPage < beginZ or currentPage > endZ) {
                continue;
            }

            unsigned int n = object->getN();
            unsigned int *x = object->getX();
            unsigned int *y = object->getY();
            unsigned int *z = object->getZ();

            for (int k = 0; k < n; ++k)
            {
                unsigned int xx = x[k] - imageBeginX;
                unsigned int yy = y[k] - imageBeginY;
                unsigned int zz = z[k] - imageBeginZ;

                //segmented[zz*width*height + yy*width + xx] = true;
                if (zz == currentPage)
                {
                    segmented[yy*width + xx] = true;
                }
                //segmented[zz*width*height + yy*width + xx] = true;
            }
        }
    }

    updateViewSegmented();

    Debug::Info("NumericToolControl::updateSegmented: Leaving");
}

QTreeWidgetItem *NumericToolControl::addTreeRoot(QString name, CellObjectType type)
{
    Debug::Info("NumericToolControl::addTreeRoot: Entering");

    QTreeWidgetItem *treeItem = new CellObject(cellTree, type, true);

    treeItem->setText(0, name);

    Debug::Info("NumericToolControl::addTreeRoot: Leaving");

    return treeItem;
}

void NumericToolControl::nextImage()
{
    Debug::Info("NumericToolControl::nextImage: Entering");

    unsigned int depth = tiffImage->getDepth();

    currentPage = std::min(currentPage + 1, (int) depth - 1);

    updateViewImage();
    updateViewPhi();
    updateViewSignedMap();
    updateSegmented();
    //updateViewSegmented();

    Debug::Info("NumericToolControl::nextImage: Leaving");
}

void NumericToolControl::prevImage()
{
    Debug::Info("NumericToolControl::prevImage: Entering");

    currentPage = std::max(currentPage - 1, 0);

    updateViewImage();
    updateViewPhi();
    updateViewSignedMap();
    updateSegmented();

    Debug::Info("NumericToolControl::prevImage: Leaving");
}

void NumericToolControl::newActiveRegion(int beginX, int endX, int beginY, int endY, int beginZ, int endZ, float sigma)
{
    Debug::Info("NumericToolControl::newActiveRegion: Entering");

    setActiveRegion(beginX, endX, beginY, endY, beginZ, endZ, sigma);

    float c1 = 0;
    float c2 = 0;

    cv3d.segmentImage(activeImage,
                      phi,
                      1.0,
                      12.5,
                      activeRegionWidth,
                      activeRegionHeight,
                      activeRegionDepth,
                      true,
                      0,
                      c1,
                      c2,
                      false,
                      false);

    updateViewImage();
    updateViewPhi();
    updateViewSegmented();

    Debug::Info("NumericToolControl::newActiveRegion: Leaving");
}

void NumericToolControl::setActiveRegion(int beginX, int endX, int beginY, int endY, int beginZ, int endZ, float sigma)
{
    Debug::Info("NumericToolControl::setActiveRegion: Entering");

    activeRegionOriginX = beginX;
    activeRegionOriginY = beginY;
    activeRegionOriginZ = beginZ;
    activeRegionWidth   = endX - beginX + 1;
    activeRegionHeight  = endY - beginY + 1;
    activeRegionDepth   = endZ - beginZ + 1;

    if (activeRegionWidth <= 0)
    {
        Debug::Error("NumericToolControl::setActiveRegion: activeRegionWidth "
                     "is " + STR(activeRegionWidth) + " from " + STR(endX) +
                     " - " + STR(beginX));
    }
    if (activeRegionHeight <= 0)
    {
        Debug::Error("NumericToolControl::setActiveRegion: activeRegionHeight "
                     "is " + STR(activeRegionHeight) + " from " + STR(endY) +
                     " - " + STR(beginY));
    }
    if (activeRegionDepth <= 0)
    {
        Debug::Error("NumericToolControl::setActiveRegion: activeRegionDepth "
                     "is " + STR(activeRegionDepth) + " from " + STR(endZ) +
                     " - " + STR(beginZ));
    }

    int N = activeRegionWidth * activeRegionHeight * activeRegionDepth;
    int paddedN = (activeRegionWidth + 2) * (activeRegionHeight + 2) * (activeRegionDepth + 2);

    if (activeImage != NULL) free(activeImage);

    if (activeImageOriginal != NULL) free(activeImageOriginal);

    if (activeFlagArray != NULL) free(activeFlagArray);

    if (phi != NULL) free(phi);

    if (filamentContour != NULL) free(filamentContour);

    activeImage = (float*) malloc(N * sizeof(float));
    if (activeImage == NULL) {
        Debug::Error("NumericToolControl::newActiveRegion: Malloc failed (activeImage)");
    }

    activeImageOriginal = (float*) malloc(N * sizeof(float));
    if (activeImageOriginal == NULL) {
        Debug::Error("NumericToolControl::newActiveRegion: Malloc failed (activeImageOriginal)");
    }

    activeFlagArray = (bool*) malloc(N * sizeof(bool));
    if (activeFlagArray == NULL) {
        Debug::Error("NumericToolControl::newActiveRegion: Malloc failed (activeFlagArray)");
    }

    phi = (float*) malloc(paddedN * sizeof(float));
    if (phi == NULL) {
        Debug::Error("NumericToolControl::newActiveRegion: Malloc failed (phi)");
    }

    filamentContour = NULL;

    activeRegion = true;

    for (int k = 0; k < activeRegionDepth; ++k)
    {
        for (int j = 0; j < activeRegionHeight; ++j)
        {
            for (int i = 0; i < activeRegionWidth; ++i)
            {
                activeImage[gIndex(i, j, k, activeRegionHeight, activeRegionWidth)]
                    = tiffImage->getValue(i + activeRegionOriginX, j + activeRegionOriginY, k + activeRegionOriginZ);
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        activeFlagArray[i] = true;
    }

    if (sigma > 0.0001)
    {
        GaussSmooth(activeImage, activeRegionWidth, activeRegionHeight, activeRegionDepth, sigma, 5);
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        activeImageOriginal[i] = activeImage[i];
    }

    Debug::Info("NumericToolControl::setActiveRegion: Leaving");
}

void NumericToolControl::runChanVeseSlot(int numIterations,
                                         float mu, float nu,
                                         float c1, float c2,
                                         bool constc1, bool constc2)
{
    Debug::Info("NumericToolControl::runChanVeseSlot: Entering");

    if (activeRegion)
    {
        float newc1 = c1;
        float newc2 = c2;
        cv3d.segmentImage(activeImage,
                          phi,
                          mu,
                          nu,
                          activeRegionWidth,
                          activeRegionHeight,
                          activeRegionDepth,
                          false,
                          numIterations,
                          newc1,
                          newc2,
                          constc1,
                          constc2);

        updateViewPhi();
    }

    Debug::Info("NumericToolControl::runChanVeseSlot: Leaving");
}

void NumericToolControl::markVesicleLocationSlot(int vesicleCenterX,
                                                 int vesicleCenterY)
{
    Debug::Info("NumericToolControl::markVesicleLocationSlot: Entering");

    int beginX = vesicleCenterX - 10;
    int endX = vesicleCenterX + 10;
    int beginY = vesicleCenterY - 10;
    int endY = vesicleCenterY + 10;
    int beginZ = currentPage - 10;
    int endZ = currentPage + 10;

    unsigned int width = tiffImage->getWidth();
    unsigned int height = tiffImage->getHeight();
    unsigned int depth = tiffImage->getDepth();

    beginX = std::max(beginX, 0         );
    endX   = std::min(endX, (int) width  - 1);
    beginY = std::max(beginY, 0         );
    endY   = std::min(endY, (int) height - 1);
    beginZ = std::max(beginZ, 0         );
    endZ   = std::min(endZ, (int) depth  - 1);

    setActiveRegion(beginX, endX, beginY, endY, beginZ, endZ, 0.1);

    for (int z = 0; z < activeRegionDepth+2; ++z)
    {
        for (int y = 0; y < activeRegionHeight+2; ++y)
        {
            for (int x = 0; x < activeRegionWidth+2; ++x)
            {
                int index = gIndex(x,
                                   y,
                                   z,
                                   activeRegionHeight+2,
                                   activeRegionWidth+2);

                phi[index] = 5.5 - std::sqrt(std::pow(x-11, 2)
                                           + std::pow(y-11, 2)
                                           + std::pow(z-11, 2));
            }
        }
    }

    mainWin->forceActiveRegionChange(beginX, endX,
                                     beginY, endY,
                                     beginZ, endZ);

    updateViewImage();
    updateViewPhi();

    Debug::Info("NumericToolControl::markVesicleLocationSlot: Leaving");
}

void NumericToolControl::markNeurofilamentSlot(int neurofilamentGlobalCenterX,
                                               int neurofilamentGlobalCenterY)
{
    Debug::Info("NumericToolControl::markNeurofilamentSlot: Entering");

    int centerX = neurofilamentGlobalCenterX - activeRegionOriginX;
    int centerY = neurofilamentGlobalCenterY - activeRegionOriginY;
    int centerZ = currentPage - activeRegionOriginZ;

    if (filamentContour == NULL)
    {
        filamentContour = (int *) malloc(2*(activeRegionDepth+2) * sizeof(int));

        for (int i = 0; i < activeRegionDepth+2; ++i)
        {
            filamentContour[gIndex(0, i, 2)] = -1;
            filamentContour[gIndex(1, i, 2)] = -1;
        }
    }

    filamentContour[gIndex(0, centerZ, 2)] = centerX;
    filamentContour[gIndex(1, centerZ, 2)] = centerY;

    for (int y = 0; y < activeRegionHeight + 2; ++y)
    {
        for (int x = 0; x < activeRegionWidth + 2; ++x)
        {
            int index = gIndex(x, y, centerZ + 1, activeRegionHeight + 2, activeRegionWidth + 2);

            phi[index] = std::max(2.5 - std::sqrt(std::pow(x - (centerX + 1), 2)
                                                + std::pow(y - (centerY + 1), 2))
                                 ,-1.0);
        }
    }

    updateViewImage(true);
    updateViewPhi();

    Debug::Info("NumericToolControl::markNeurofilamentSlot: Leaving");
}

void NumericToolControl::finalizeNeurofilamentSlot()
{
    Debug::Info("NumericToolControl::finalizeNeurofilamentSlot: Entering");

    if (filamentContour == NULL)
    {
        return;
    }

    int N = activeRegionWidth * activeRegionHeight * activeRegionDepth;

    for (int i = 0; i < N; ++i)
    {
        activeImage[i] = cytosolI;
    }

    int prevX = -1;
    int prevY = -1;
    int prevZ = -1;

    for (int k = 1; k < activeRegionDepth+1; ++k)
    {
        if (filamentContour[gIndex(0,k,2)] == -1)
        {
            continue;
        }

        int curX = filamentContour[gIndex(0,k,2)];
        int curY = filamentContour[gIndex(1,k,2)];
        int curZ = k;

        if (prevZ != -1)
        {
            restoreImageLine(prevX, curX, prevY, curY, prevZ, curZ);
        }

        prevX = curX;
        prevY = curY;
        prevZ = curZ;
    }

    Debug::Info("NumericToolControl::finalizeNeurofilamentSlot: Leaving");
}

void NumericToolControl::setPosSlot(int beginX, int endX,
                                    int beginY, int endY,
                                    int beginZ, int endZ)
{
    Debug::Info("NumericToolControl::setPosSlot: Entering");

    int width = endX - beginX + 1;
    int height = endY - beginY + 1;
    int depth = endZ - beginZ + 1;

    int x0 = beginX - activeRegionOriginX;
    int y0 = beginY - activeRegionOriginY;
    int z0 = beginZ - activeRegionOriginZ;

    for (int k = z0; k < z0 + depth; ++k)
    {
        for (int j = y0; j < y0 + height; ++j)
        {
            for (int i = x0; i < x0 + width; ++i)
            {
                int phiIndex = gIndex(i,
                                      j,
                                      k,
                                      activeRegionHeight,
                                      activeRegionWidth,
                                      1);
                if (phi[phiIndex] < 0) phi[phiIndex] = -phi[phiIndex];
            }
        }
    }

    updateViewPhi();

    Debug::Info("NumericToolControl::setPosSlot: Leaving");
}

void NumericToolControl::setNegSlot(int beginX, int endX, int beginY, int endY, int beginZ, int endZ)
{
    Debug::Info("NumericToolControl::setNegSlot: Entering");

    int width = endX - beginX + 1;
    int height = endY - beginY + 1;
    int depth = endZ - beginZ + 1;

    int x0 = beginX - activeRegionOriginX;
    int y0 = beginY - activeRegionOriginY;
    int z0 = beginZ - activeRegionOriginZ;

    for (int k = z0; k < z0 + depth; ++k)
    {
        for (int j = y0; j < y0 + height; ++j)
        {
            for (int i = x0; i < x0 + width; ++i)
            {
                if (i < 0 or j < 0 or k < 0 or i >= activeRegionWidth
                    or j >= activeRegionHeight or k >= activeRegionDepth)
                {
                    continue;
                }
                int phiIndex = gIndex(i,
                                      j,
                                      k,
                                      activeRegionHeight,
                                      activeRegionWidth,
                                      1);
                if (phi[phiIndex] > 0) phi[phiIndex] = -phi[phiIndex];
                int imageIndex = gIndex(i,
                                        j,
                                        k,
                                        activeRegionHeight,
                                        activeRegionWidth);

                //distance from XY border
                float dFXYB = std::min(
                                        std::min(i-x0,
                                                 j-y0),
                                        std::min(x0 + width - 1 - i,
                                                 y0 + height - 1 - j)
                                        );

                float alpha = std::min(dFXYB / 3.0, 3.0) / 6.0;
                activeImage[imageIndex] =
                    clamp(( alpha       * cytosolI
                          + (1 - alpha) * activeImage[imageIndex]),
                          0, 255);
            }
        }
    }

    updateViewImage();
    updateViewPhi();
    updateViewSegmented();

    Debug::Info("NumericToolControl::setNegSlot: Leaving");
}

void NumericToolControl::setDisallowedSlot(int beginX, int endX, int beginY, int endY, int beginZ, int endZ)
{
    Debug::Info("NumericToolControl::setDisallowedSlot: Entering");

    int width = endX - beginX + 1;
    int height = endY - beginY + 1;
    int depth = endZ - beginZ + 1;

    int x0 = beginX - activeRegionOriginX;
    int y0 = beginY - activeRegionOriginY;
    int z0 = beginZ - activeRegionOriginZ;

    for (int k = z0; k < z0 + depth; ++k)
    {
        for (int j = y0; j < y0 + height; ++j)
        {
            for (int i = x0; i < x0 + width; ++i)
            {
                int index = gIndex(i, j, k, activeRegionHeight, activeRegionWidth);
                activeFlagArray[index] = false;
            }
        }
    }

    Debug::Info("NumericToolControl::setDisallowedSlot: Leaving");
}

void NumericToolControl::swapPhiSignSlot()
{
    Debug::Info("NumericToolControl::swapPhiSignSlot: Entering");

    for (int k = 0; k < activeRegionDepth; ++k)
    {
        for (int j = 0; j < activeRegionHeight; ++j)
        {
            for (int i = 0; i < activeRegionWidth; ++i)
            {
                int index = gIndex(i, j, k, activeRegionHeight, activeRegionWidth, 1);
                phi[index] = -phi[index];
            }
        }
    }

    updateViewPhi();

    Debug::Info("NumericToolControl::swapPhiSignSlot: Leaving");
}

void NumericToolControl::reinitPhiSlot()
{
    Debug::Info("NumericToolControl::reinitPhiSlot: Entering");

    if (not activeRegion)
    {
        Debug::Info("NumericToolControl::reinitPhiSlot: Leaving");

        return;
    }

    //int diagonalDistance = std::sqrt(activeRegionWidth * activeRegionWidth
    //                               + activeRegionHeight * activeRegionHeight
    //                              + activeRegionDepth * activeRegionDepth);

    SignedDistanceMap(phi,
                      activeRegionWidth,
                      activeRegionHeight,
                      activeRegionDepth,
                      /*diagonalDistance * signedDistanceIterationFactor*/ 50);

    updateViewPhi();

    Debug::Info("NumericToolControl::reinitPhiSlot: Leaving");
}

void NumericToolControl::renderIn3DSlot()
{
    Debug::Info("NumericToolControl::renderIn3DSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    std::vector<CellObject*> objects;

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            objects.push_back(object);
        }
    }

    if (objects.size() > 0)
    {
        CellPlot cp;
        cp.MultiObjectInteractiveImplicitSurface(objects);
    }

    Debug::Info("NumericToolControl::renderIn3DSlot: Leaving");
}

void NumericToolControl::storeObjectSlot()
{
    Debug::Info("NumericToolControl::storeObjectSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    if (selectedObjects.length() != 1) return;

    CellObject *selectedObject = (CellObject *) selectedObjects[0];

    if (not selectedObject->isRootObject()) return;

    CellObjectType type = selectedObject->getType();

    int n = 0;

    for (int z = 0; z < activeRegionDepth; ++z)
    {
        for (int y = 0; y < activeRegionHeight; ++y)
        {
            for (int x = 0; x < activeRegionWidth; ++x)
            {
                int index = gIndex(x, y, z, activeRegionHeight,
                                   activeRegionWidth, 1);
                if (phi[index] > 0) n += 1;
            }
        }
    }

    int totalXDiff = activeRegionOriginX + imageBeginX;
    int totalYDiff = activeRegionOriginY + imageBeginY;
    int totalZDiff = activeRegionOriginZ + imageBeginZ;

    unsigned int *X = (unsigned int*) malloc(n * sizeof(unsigned int));
    unsigned int *Y = (unsigned int*) malloc(n * sizeof(unsigned int));
    unsigned int *Z = (unsigned int*) malloc(n * sizeof(unsigned int));

    int curGlobalIndex = 0;

    for (int z = 0; z < activeRegionDepth; ++z)
    {
        for (int y = 0; y < activeRegionHeight; ++y)
        {
            for (int x = 0; x < activeRegionWidth; ++x)
            {
                int index = gIndex(x, y, z, activeRegionHeight,
                                   activeRegionWidth, 1);
                if (phi[index] > 0)
                {
                    X[curGlobalIndex] = x + totalXDiff;
                    Y[curGlobalIndex] = y + totalYDiff;
                    Z[curGlobalIndex] = z + totalZDiff;

                    ++curGlobalIndex;
                }
            }
        }
    }

    CellObject *treeItem = new CellObject(type,
                                          false,
                                          X,
                                          Y,
                                          Z,
                                          n,
                                          phi,
                                          totalXDiff,
                                          totalYDiff,
                                          totalZDiff,
                                          activeRegionDepth,
                                          activeRegionHeight,
                                          activeRegionWidth);

    treeItem->setText(1, QString::number(n));
    QString name;

    switch (type)
    {
    case COT_UNCATEGORIZED:

        treeItem->setText(0, QString("U")
                           + QString::number(cellUncategorizedCount));
        cellUncategorized->addChild(treeItem);
        ++cellUncategorizedCount;
        break;

    case COT_MEMBRANE:

        treeItem->setText(0, QString("M")
                           + QString::number(cellMembraneCount));
        cellMembranes->addChild(treeItem);
        ++cellMembraneCount;
        break;

    case COT_VESICLE:

        treeItem->setText(0, QString("V")
                           + QString::number(cellVesicleCount));
        cellVesicles->addChild(treeItem);
        ++cellVesicleCount;
        break;

    case COT_FILAMENT:

        treeItem->setText(0, QString("F")
                           + QString::number(cellFilamentCount));
        cellFilaments->addChild(treeItem);
        ++cellFilamentCount;
        break;

    case COT_SYNAPSE:

        treeItem->setText(0, QString("S")
                           + QString::number(cellSynapseCount));
        cellSynapses->addChild(treeItem);
        ++cellSynapseCount;
        break;

    case COT_ENDOSOME:

        treeItem->setText(0, QString("En")
                           + QString::number(cellEndosomeCount));
        cellEndosomes->addChild(treeItem);
        ++cellEndosomeCount;
        break;

    case COT_ER:

        treeItem->setText(0, QString("ER")
                           + QString::number(cellERCount));
        cellERs->addChild(treeItem);
        ++cellERCount;
        break;
    }

    updateSegmented();

    Debug::Info("NumericToolControl::storeObjectSlot: Leaving");
}

void NumericToolControl::loadObjectSlot()
{
    Debug::Info("NumericToolControl::loadObjectSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    if (selectedObjects.length() <= 0) return;

    CellObject *object = (CellObject *) selectedObjects[0];

    if (object->isRootObject()) return;

    int margin = 3;

    int beginX, endX, beginY, endY, beginZ, endZ;

    object->getBounds(beginX, endX, beginY, endY, beginZ, endZ);

    unsigned int width = tiffImage->getWidth();
    unsigned int height = tiffImage->getHeight();
    unsigned int depth = tiffImage->getDepth();

    beginX -= imageBeginX;
    endX -= imageBeginX;
    beginY -= imageBeginY;
    endY -= imageBeginY;
    beginZ -= imageBeginZ;
    endZ -= imageBeginZ;

    beginX = std::max(beginX - margin, 0         );
    endX   = std::min(endX   + margin, (int) width  - 1);
    beginY = std::max(beginY - margin, 0         );
    endY   = std::min(endY   + margin, (int) height - 1);
    beginZ = std::max(beginZ - margin, 0         );
    endZ   = std::min(endZ   + margin, (int) depth  - 1);

    currentPage = (beginZ + endZ) / 2;

    setActiveRegion(beginX, endX, beginY, endY, beginZ, endZ, 0.0);

    float c1 = 0;
    float c2 = 0;

    cv3d.segmentImage(activeImage,
                      phi,
                      1.0,
                      12.5,
                      activeRegionWidth,
                      activeRegionHeight,
                      activeRegionDepth,
                      true,
                      0,
                      c1,
                      c2,
                      false,
                      false);

    for (int z = 0; z < activeRegionDepth; ++z)
    {
        for (int y = 0; y < activeRegionHeight; ++y)
        {
            for (int x = 0; x < activeRegionWidth; ++x)
            {
                int index = gIndex(x, y, z, activeRegionHeight,
                                   activeRegionWidth, 1);
                phi[index] = -1;
            }
        }
    }

    unsigned int n = object->getN();
    unsigned int *X = object->getX();
    unsigned int *Y = object->getY();
    unsigned int *Z = object->getZ();

    int totalXDiff = activeRegionOriginX + imageBeginX;
    int totalYDiff = activeRegionOriginY + imageBeginY;
    int totalZDiff = activeRegionOriginZ + imageBeginZ;

    for (int i = 0; i < n; ++i)
    {

        int index = gIndex(X[i] - totalXDiff,
                           Y[i] - totalYDiff,
                           Z[i] - totalZDiff,
                           activeRegionHeight,
                           activeRegionWidth,
                           1);

        phi[index] = 1;
    }

    int diagonalDistance = (int) std::sqrt(
                                  activeRegionWidth * activeRegionWidth
                                + activeRegionHeight * activeRegionHeight
                                + activeRegionDepth * activeRegionDepth
                                );

    Debug::Info("Running Signed distance map numeric function using: "
                + STR(diagonalDistance * signedDistanceIterationFactor)
                + " iterations");

    SignedDistanceMap(phi,
                      activeRegionWidth,
                      activeRegionHeight,
                      activeRegionDepth,
                      diagonalDistance * signedDistanceIterationFactor);

    mainWin->forceActiveRegionChange(beginX, endX,
                                     beginY, endY,
                                     beginZ, endZ);

    object->getPhi(phi,
                   totalXDiff,
                   totalYDiff,
                   totalZDiff,
                   activeRegionDepth,
                   activeRegionHeight,
                   activeRegionWidth);

    updateViewImage();
    updateViewPhi();
    updateViewSegmented();

    Debug::Info("NumericToolControl::loadObjectSlot: Leaving");
}

void NumericToolControl::deleteObjectSlot()
{
    Debug::Info("NumericToolControl::deleteObjectSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];

        if (not object->isRootObject())
        {
            delete selectedObjects[i];
        }
    }

    Debug::Info("NumericToolControl::deleteObjectSlot: Leaving");
}

void NumericToolControl::blightSlot(int imageX, int imageY, float sigma)
{
    Debug::Info("NumericToolControl::blightSlot: Entering");

    int kernelWidth = ((int) sigma + 1) * 4 + 1;

    if (kernelWidth % 2 == 0) kernelWidth += 1;

    int c = kernelWidth / 2;

    int activeImageX = imageX - activeRegionOriginX;
    int activeImageY = imageY - activeRegionOriginY;

    for (int y = 0; y < kernelWidth; ++y)
    {
        for (int x = 0; x < kernelWidth; ++x)
        {
            float kernelxy = 1.0 / 2.0
                            * std::exp(-((std::pow(x-c, 2)
                            + std::pow(y-c, 2)) / (2*sigma*sigma)));

            int ix = activeImageX + x - c;
            int iy = activeImageY + y - c;

            int iz = currentPage - activeRegionOriginZ;
            if (ix >= 0 and ix < activeRegionWidth and
                iy >= 0 and iy < activeRegionHeight and
                iz >= 0 and iz < activeRegionDepth)
            {
                int actImgIdx = gIndex(ix, iy, iz,
                                       activeRegionHeight,
                                       activeRegionWidth);

                activeImage[actImgIdx] =      kernelxy  * cytosolI
                                       + (1 - kernelxy) * activeImage[actImgIdx];
            }
        }
    }

    updateViewImage(true);

    Debug::Info("NumericToolControl::blightSlot: Leaving");
}

void NumericToolControl::restoreSlot(int imageX, int imageY, float sigma)
{
    Debug::Info("NumericToolControl::blightSlot: Entering");

    int kernelWidth = ((int) sigma + 1) * 4 + 1;

    if (kernelWidth % 2 == 0) kernelWidth += 1;

    int c = kernelWidth / 2;

    int activeImageX = imageX - activeRegionOriginX;
    int activeImageY = imageY - activeRegionOriginY;

    for (int y = 0; y < kernelWidth; ++y)
    {
        for (int x = 0; x < kernelWidth; ++x)
        {
            float kernelxy = 1.0 / 2.0
                            * std::exp(-((std::pow(x-c, 2)
                            + std::pow(y-c, 2)) / (2*sigma*sigma)));

            int ix = activeImageX + x - c;
            int iy = activeImageY + y - c;

            int iz = currentPage - activeRegionOriginZ;
            if (ix >= 0 and ix < activeRegionWidth and
                iy >= 0 and iy < activeRegionHeight and
                iz >= 0 and iz < activeRegionDepth)
            {
                int actImgIdx = gIndex(ix, iy, iz,
                                       activeRegionHeight,
                                       activeRegionWidth);

                activeImage[actImgIdx] =      kernelxy  * activeImageOriginal[actImgIdx]
                                       + (1 - kernelxy) * activeImage[actImgIdx];
            }
        }
    }

    updateViewImage(true);

    Debug::Info("NumericToolControl::blightSlot: Leaving");
}

void NumericToolControl::SDSetSourcesSlot()
{
    Debug::Info("NumericToolControl::SDSetSourcesSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    sdSourceObjects.clear();

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            sdSourceObjects.push_back(object);
        }
    }

    Debug::Info("NumericToolControl::SDSetSourcesSlot: Length: " + STR(sdSourceObjects.size()));

    Debug::Info("NumericToolControl::SDSetSourcesSlot: Leaving");
}

void NumericToolControl::SDSetBoundariesSlot()
{
    Debug::Info("NumericToolControl::SDSetBoundariesSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    sdBoundaryObjects.clear();

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            sdBoundaryObjects.push_back(object);
        }
    }

    Debug::Info("NumericToolControl::SDSetBoundariesSlot: Length: " + STR(sdBoundaryObjects.size()));

    Debug::Info("NumericToolControl::SDSetBoundariesSlot: Leaving");
}

void NumericToolControl::SDSetPassThroughSlot(float travelCost)
{
    Debug::Info("NumericToolControl::SDSetPassThroughSlot: Entering");

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();

    sdPassThroughObjects.clear();
    travelWeight = travelCost;

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            sdPassThroughObjects.push_back(object);
        }
    }

    Debug::Info("NumericToolControl::SDSetPassThroughSlot: Length: " + STR(sdPassThroughObjects.size()));

    Debug::Info("NumericToolControl::SDSetPassThroughSlot: Leaving");
}

void NumericToolControl::SDRunBoundedSignedDistanceMapSlot()
{
    Debug::Info("NumericToolControl::SDRunBoundedSignedDistanceMapSlot:"
                " Entering");

    // Run through all sources and boundaries and make a 3d array float
    int xMin = INT_MAX, xMax = INT_MIN;
    int yMin = INT_MAX, yMax = INT_MIN;
    int zMin = INT_MAX, zMax = INT_MIN;

    int beginX, endX, beginY, endY, beginZ, endZ;

    if (sdSourceObjects.size() == 0)
    {
        Debug::Info("NumericToolControl::SDRunBoundedSignedDistanceMapSlot:"
            " Leaving (Early)");
        return;
    }

    // Get bounds from all objects involved

    for (int i = 0; i < sdSourceObjects.size(); ++i)
    {
        CellObject *object = sdSourceObjects[i];

        object->getBounds(beginX, endX, beginY, endY, beginZ, endZ);

        xMin = std::min(xMin, beginX); xMax = std::max(xMax, endX);
        yMin = std::min(yMin, beginY); yMax = std::max(yMax, endY);
        zMin = std::min(zMin, beginZ); zMax = std::max(zMax, endZ);
    }

    for (int i = 0; i < sdBoundaryObjects.size(); ++i)
    {
        CellObject *object = sdBoundaryObjects[i];

        object->getBounds(beginX, endX, beginY, endY, beginZ, endZ);

        xMin = std::min(xMin, beginX); xMax = std::max(xMax, endX);
        yMin = std::min(yMin, beginY); yMax = std::max(yMax, endY);
        zMin = std::min(zMin, beginZ); zMax = std::max(zMax, endZ);
    }

    for (int i = 0; i < sdPassThroughObjects.size(); ++i)
    {
        CellObject *object = sdPassThroughObjects[i];

        object->getBounds(beginX, endX, beginY, endY, beginZ, endZ);

        xMin = std::min(xMin, beginX); xMax = std::max(xMax, endX);
        yMin = std::min(yMin, beginY); yMax = std::max(yMax, endY);
        zMin = std::min(zMin, beginZ); zMax = std::max(zMax, endZ);
    }

    // Calculate final bounds

    sdWidth = xMax - xMin + 1;
    sdHeight = yMax - yMin + 1;
    sdDepth = zMax - zMin + 1;

    sdOriginX = xMin - imageBeginX;
    sdOriginY = yMin - imageBeginY;
    sdOriginZ = zMin - imageBeginZ;

    // Reinit signed distance map

    int N = sdWidth * sdHeight * sdDepth;

    if (signedDistanceMap != NULL) free(signedDistanceMap);

    signedDistanceMap = (float *) malloc(N * sizeof(float));

    if (signedDistanceMap == NULL) {
        Debug::Error("NumericToolControl::SDRunBoundedSignedDistanceMapSlot: "
                     "Malloc failed (signedDistanceMap)");
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) signedDistanceMap[i] = 1;

    for (int i = 0; i < sdBoundaryObjects.size(); ++i)
    {
        CellObject *object = sdBoundaryObjects[i];

        unsigned int n = object->getN();

        unsigned int* X = object->getX();
        unsigned int* Y = object->getY();
        unsigned int* Z = object->getZ();

        #pragma omp parallel for
        for (int j = 0; j < n; ++j)
        {
            int x = X[j] - xMin;
            int y = Y[j] - yMin;
            int z = Z[j] - zMin;

            int idx = gIndex(x, y, z, sdHeight, sdWidth);

            signedDistanceMap[idx] = -1;
        }
    }

    for (int i = 0; i < sdSourceObjects.size(); ++i)
    {
        CellObject *object = sdSourceObjects[i];

        unsigned int n = object->getN();

        unsigned int* X = object->getX();
        unsigned int* Y = object->getY();
        unsigned int* Z = object->getZ();

        #pragma omp parallel for
        for (int j = 0; j < n; ++j)
        {
            int x = X[j] - xMin;
            int y = Y[j] - yMin;
            int z = Z[j] - zMin;

            int idx = gIndex(x, y, z, sdHeight, sdWidth);

            signedDistanceMap[idx] = -2;
        }
    }

    for (int i = 0; i < sdPassThroughObjects.size(); ++i)
    {
        CellObject *object = sdPassThroughObjects[i];

        unsigned int n = object->getN();

        unsigned int* X = object->getX();
        unsigned int* Y = object->getY();
        unsigned int* Z = object->getZ();

        #pragma omp parallel for
        for (int j = 0; j < n; ++j)
        {
            int x = X[j] - xMin;
            int y = Y[j] - yMin;
            int z = Z[j] - zMin;

            int idx = gIndex(x, y, z, sdHeight, sdWidth);

            signedDistanceMap[idx] = travelWeight;
        }
    }

    FastMarching(sdWidth, sdHeight, sdDepth, signedDistanceMap);

    cleanSignedDistanceMap(sdWidth, sdHeight, sdDepth, signedDistanceMap);

    updateViewSignedMap();

    Debug::Info("NumericToolControl::SDRunBoundedSignedDistanceMapSlot:"
                " Leaving");
}

void NumericToolControl::SDReadDataSlot()
{
    Debug::Info("NumericToolControl::SDReadDataSlot: Entering");

    if (signedDistanceMap == NULL)
    {
        return;
    }

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();
    std::sort(selectedObjects.begin(),
              selectedObjects.end(),
              cellObjectLexiographicCompare);

    std::string distanceMassString = "";
    std::string dataString = "";
    std::string dataString2 = "";

    float maxDist = std::sqrt(std::pow(sdWidth, 2) + std::pow(sdHeight, 2) + std::pow(sdDepth, 2));

    float globalMaxDist = 0;

    float stepSize = 2.5; // change this accordingly
    int N = (int) (maxDist / stepSize) + 1;

    int *bins = (int *) malloc(N * sizeof(int));

    for (int i = 0; i < N; ++i) bins[i] = 0;

    for (int i = 0; i < sdWidth * sdHeight * sdDepth; ++i)
    {
        float val = signedDistanceMap[i];
        globalMaxDist = std::max(globalMaxDist, signedDistanceMap[i]);
        if (val == 0) continue;
        bins[(int) (val / stepSize)] += 1;
    }

    for (int i = 0; i < N; ++i) distanceMassString += STR(bins[i]) + ", ";

    free(bins);

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            unsigned int n = object->getN();
            unsigned int *X = object->getX();
            unsigned int *Y = object->getY();
            unsigned int *Z = object->getZ();

            int xCenterOfMass = 0;
            int yCenterOfMass = 0;
            int zCenterOfMass = 0;

            for (int j = 0; j < n; ++j)
            {
                xCenterOfMass += X[j];
                yCenterOfMass += Y[j];
                zCenterOfMass += Z[j];
            }

            xCenterOfMass /= n;
            yCenterOfMass /= n;
            zCenterOfMass /= n;

            xCenterOfMass = xCenterOfMass - sdOriginX - imageBeginX;
            yCenterOfMass = yCenterOfMass - sdOriginY - imageBeginY;
            zCenterOfMass = zCenterOfMass - sdOriginZ - imageBeginZ;

            if (xCenterOfMass >= 0 and xCenterOfMass < sdWidth and
                yCenterOfMass >= 0 and yCenterOfMass < sdHeight and
                zCenterOfMass >= 0 and zCenterOfMass < sdDepth)
            {
                dataString += STR(signedDistanceMap[gIndex(xCenterOfMass,
                                                           yCenterOfMass,
                                                           zCenterOfMass,
                                                           sdHeight,
                                                           sdWidth,
                                                           0)]) + ", ";
                Debug::Info(object->text(0).toStdString() + ", " + STR(signedDistanceMap[gIndex(xCenterOfMass,
                                                           yCenterOfMass,
                                                           zCenterOfMass,
                                                           sdHeight,
                                                           sdWidth,
                                                           0)]));
            }

            float minDist = std::sqrt(std::pow(sdWidth, 2) + std::pow(sdHeight, 2) + std::pow(sdDepth, 2));

            for (int j = 0; j < n; ++j)
            {
                int x = X[j] - sdOriginX - imageBeginX;
                int y = Y[j] - sdOriginY - imageBeginY;
                int z = Z[j] - sdOriginZ - imageBeginZ;

                if (x >= 0 and x < sdWidth and
                    y >= 0 and y < sdHeight and
                    z >= 0 and z < sdDepth)
                {
                    float newDist = signedDistanceMap[gIndex(x,
                                                             y,
                                                             z,
                                                             sdHeight,
                                                             sdWidth)];

                    if (newDist > 0 and newDist < minDist)
                    {
                        minDist = newDist;
                    }
                }
            }

            dataString2 += STR(minDist) + ", ";
        }
    }

    Debug::Info("Data: distance mass (step size: " + STR(stepSize) + "):\n" + distanceMassString.substr(0, distanceMassString.size()-2));

    Debug::Info("Data: object average distances:\n" + dataString.substr(0, dataString.size()-2));

    Debug::Info("Data: object minimum distances:\n" + dataString2.substr(0, dataString2.size()-2));

    Debug::Info("Data: Max Distance on distance map: " + STR(globalMaxDist));

    Debug::Info("NumericToolControl::SDReadDataSlot: Leaving");
}

void NumericToolControl::SDReadFullDataSlot()
{
    Debug::Info("NumericToolControl::SDReadFullDataSlot: Entering");

    if (signedDistanceMap == NULL)
    {
        return;
    }

    QList<QTreeWidgetItem *> selectedObjects = cellTree->selectedItems();
    std::sort(selectedObjects.begin(),
              selectedObjects.end(),
              cellObjectLexiographicCompare);

    std::string distanceMassString = "";
    std::string dataString = "";
    std::string dataString2 = "";

    float maxDist = std::sqrt(std::pow(sdWidth, 2) + std::pow(sdHeight, 2) + std::pow(sdDepth, 2));

    float globalMaxDist = 0;

    float stepSize = 2.5; // change this accordingly
    int N = (int) (maxDist / stepSize) + 1;

    int *bins = (int *) malloc(N * sizeof(int));

    for (int i = 0; i < N; ++i) bins[i] = 0;

    for (int i = 0; i < sdWidth * sdHeight * sdDepth; ++i)
    {
        float val = signedDistanceMap[i];
        globalMaxDist = std::max(globalMaxDist, signedDistanceMap[i]);
        if (val == 0) continue;
        bins[(int) (val / stepSize)] += 1;
    }

    for (int i = 0; i < N; ++i) distanceMassString += STR(bins[i]) + ", ";

    free(bins);

    for (int i = 0; i < selectedObjects.length(); ++i)
    {
        CellObject *object = (CellObject *) selectedObjects[i];
        if (not object->isRootObject())
        {
            unsigned int n = object->getN();
            unsigned int *X = object->getX();
            unsigned int *Y = object->getY();
            unsigned int *Z = object->getZ();

            for (int j = 0; j < n; ++j)
            {
                int x = X[j];
                int y = Y[j];
                int z = Z[j];


                x = x - sdOriginX - imageBeginX;
                y = y - sdOriginY - imageBeginY;
                z = z - sdOriginZ - imageBeginZ;

                if (x >= 0 and x < sdWidth and
                    y >= 0 and y < sdHeight and
                    z >= 0 and z < sdDepth)
                {
                    dataString += STR(signedDistanceMap[gIndex(x,
                                                               y,
                                                               z,
                                                               sdHeight,
                                                               sdWidth,
                                                               0)]) + ", ";
                }
            }
        }
    }

    Debug::Info("Data: distance mass (step size: " + STR(stepSize) + "):\n" + distanceMassString.substr(0, distanceMassString.size()-2));

    Debug::Info("Data: object full set of point distances:\n" + dataString.substr(0, dataString.size()-2));

    Debug::Info("Data: Max Distance on distance map: " + STR(globalMaxDist));

    Debug::Info("NumericToolControl::SDReadFullDataSlot: Leaving");
}

void NumericToolControl::runCustomFunctionSlot()
{
    //loadMatlabMembraneSegmentation();
    //trainAndPredictForVesicleCalssification();
    ffmTest();
}

void NumericToolControl::cleanSignedDistanceMap(int width,
                                                int height,
                                                int depth,
                                                float *sdm)
{
    int N = width * height * depth;

    bool *clearFlagArray = (bool *) malloc(N * sizeof(bool));

    for (int i = 0; i < N; ++i)
    {
        clearFlagArray[i] = false;
    }

    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                int gIdx = gIndex(i, j, k, height, width, 0);
                if (sdm[gIdx] > 0) {
                    clearFlagArray[gIdx] = true;
                } else {
                    break;
                }
            }
            for (int i = width - 1; i >= 0; --i)
            {
                int gIdx = gIndex(i, j, k, height, width, 0);
                if (sdm[gIdx] > 0) {
                    clearFlagArray[gIdx] = true;
                } else {
                    break;
                }
            }
        }

        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                int gIdx = gIndex(i, j, k, height, width, 0);
                if (sdm[gIdx] > 0) {
                    clearFlagArray[gIdx] = true;
                } else {
                    break;
                }
            }
            for (int j = height - 1; j >= 0; --j)
            {
                int gIdx = gIndex(i, j, k, height, width, 0);
                if (sdm[gIdx] > 0) {
                    clearFlagArray[gIdx] = true;
                } else {
                    break;
                }
            }
        }
    }

    for (int i = 0; i < N; ++i)
    {
        if (clearFlagArray[i])
        {
            sdm[i] = 0;
        }
    }

    free(clearFlagArray);
}

void NumericToolControl::normalizeImageToRange(float *image,
                                               float &phiZero,
                                               const int width,
                                               const int height,
                                               const float imageMin,
                                               const float imageMax,
                                               const float boundMin,
                                               const float boundMax)
{
    Debug::Info("NumericToolControl::normalizeImageToRange: Entering");

    const int N = width * height;

    float boundRange = boundMax - boundMin;
    float imageRange = imageMax - imageMin;

    for (int i = 0; i < N; ++i)
    {
        float val = (image[i] - imageMin) / imageRange;
        image[i] = std::max(std::min(val * boundRange + boundMin, boundMax),boundMin);
    }

    phiZero = -imageMin / imageRange * boundRange + boundMin;

    Debug::Info("NumericToolControl::normalizeImageToRange: Leaving");
}

void NumericToolControl::cosineModulation(float *image, int width, int height, float period)
{
    Debug::Info("NumericToolControl::cosineModulation: Entering");

    const int N = width * height;

    for (int i = 0; i < N; ++i)
    {
        image[i] = std::pow(2.0,std::cos(image[i]*2*PI / period)) * (255.0 / 2.0);
    }

    Debug::Info("NumericToolControl::cosineModulation: Leaving");
}

void NumericToolControl::restoreImageLine(float beginX,
                                          float endX,
                                          float beginY,
                                          float endY,
                                          float beginZ,
                                          float endZ)
{
    float maxStepSize = 1.0;
    float lw = 1.5;

    float d = std::sqrt(std::pow(endX - beginX, (float) 2) +
                        std::pow(endY - beginY, (float) 2) +
                        std::pow(endZ - beginZ, (float) 2));

    int numSteps = (int) (d / maxStepSize) + 1;
    float stepSize = d / numSteps;

    float dx = (endX - beginX) / d;
    float dy = (endY - beginY) / d;
    float dz = (endZ - beginZ) / d;

    auto valueFunction = [lw](float d)
    {
        float val;
        if (d < lw)
        {
            val = (float) 1.0;
        }
        else if (d > 2*lw)
        {
            val = (float) 0.0;
        }
        else
        {
            val =  (2*lw - d) / lw;
        }
        return val;
    };

    float x = (float) beginX;
    float y = (float) beginY;
    float z = (float) beginZ;

    for (int step = 0; step <= numSteps; ++step)
    {
        int i = std::round(x);
        int j = std::round(y);
        int k = std::round(z);

        #pragma omp parallel for collapse(3)
        for (int ii = i-11; ii <= i+11; ++ii)
        {
            for (int jj = j-11; jj <= j+11; ++jj)
            {
                for (int kk = k-11; kk <= k+11; ++kk)
                {
                    if (ii < 0 or ii >= activeRegionWidth) continue;
                    if (jj < 0 or jj >= activeRegionHeight) continue;
                    if (kk < 0 or kk >= activeRegionDepth) continue;

                    // indexing boundaries here
                    float pd = std::sqrt(std::pow(ii-x, 2.0) +
                                         std::pow(jj-y, 2.0) +
                                         std::pow(kk-z, 2.0));

                    float val = valueFunction(pd);

                    int gIdx = gIndex(ii,
                                      jj,
                                      kk,
                                      activeRegionHeight,
                                      activeRegionWidth,
                                      0);

                    val =    val  * activeImageOriginal[gIdx]
                        + (1-val) * cytosolI;

                    activeImage[gIdx] = std::min(activeImage[gIdx], val);
                }
            }
        }

        x += stepSize * dx;
        y += stepSize * dy;
        z += stepSize * dz;
    }
}

void NumericToolControl::loadMatlabMembraneSegmentation()
{
    // Load text file
    std::string path = "../data/volumedata_registered_4.txt";

    std::ifstream infile(path);

    int width = tiffImage->getWidth();
    int height = tiffImage->getHeight();
    int depth = tiffImage->getDepth();

    int k = 0;
    int lineIndex = 0;

    std::vector<float> x[depth];
    std::vector<float> y[depth];

    while (infile)
    {
        std::string s;

        if (!getline(infile, s)) break;

        std::stringstream ss(s);

        float val;

        if (lineIndex % 2 == 0)
        {
            while (ss >> val)
            {
                y[k].push_back(val - imageBeginX);
                if (ss.peek() == ',') ss.ignore();
            }
        }
        else
        {
            while (ss >> val)
            {
                x[k].push_back(val - imageBeginY);
                if (ss.peek() == ',') ss.ignore();
            }
        }

        if (lineIndex % 2 == 1)
        {
            if (x[k].size() != y[k].size())
            {
                Debug::Error("Misaligned x, y vectors read from file");
            }

            for (int i = 0; i < x[k].size(); ++i)
            {
                Debug::Info(STR(x[k][i]) + ", " + STR(y[k][i]) + ", " + STR(k));
            }
        }

        k += lineIndex % 2;
        ++lineIndex;

    }

    int margin = 10;

    if (!infile.eof())
    {
        Debug::Warning("NumericToolControl::loadMatlabMembraneSegmentationSlot: "
                       "File read error");
    }

    int beginX = width-1;
    int endX = 0;
    int beginY = height-1;
    int endY = 0;
    int beginZ = depth-1;
    int endZ = 0;

    for (int k = 0; k < depth; ++k)
    {
        for (int i = 0; i < x[k].size(); ++i)
        {
            if (i == 0)
            {
                beginZ = std::min(beginZ, k);
                endZ = std::max(endZ, k);
            }

            beginX = std::min(beginX, (int) x[k][i]);
            endX = std::max(endX, (int) x[k][i] + 1);
            beginY = std::min(beginY, (int) y[k][i]);
            endY = std::max(endY, (int) y[k][i] + 1);
        }
    }

    // Add extra margin
    beginX = std::max(beginX - margin, 0);
    endX = std::min(endX + margin, width - 1);
    beginY = std::max(beginY - margin, 0);
    endY = std::min(endY + margin, height - 1);
    beginZ = std::max(beginZ - margin, 0);
    endZ = std::min(endZ + margin, depth - 1);

    // Copy/paste start

    setActiveRegion(beginX, endX, beginY, endY, beginZ, endZ, 0.0);

    float c1 = 0;
    float c2 = 0;

    cv3d.segmentImage(activeImage,
                      phi,
                      1.0,
                      12.5,
                      activeRegionWidth,
                      activeRegionHeight,
                      activeRegionDepth,
                      true,
                      0,
                      c1,
                      c2,
                      false,
                      false);

    int N = activeRegionWidth * activeRegionHeight * activeRegionDepth;

    for (int i = 0; i < N; ++i)
    {
        activeImage[i] = cytosolI;
    }

    // Transform points to activeArea
    for (int k = 0; k < depth; ++k)
    {
        int n = x[k].size();

        if (n == 0) continue;

        for (int i = 0; i < n; ++i)
        {
            x[k][i] -= activeRegionOriginX;
            y[k][i] -= activeRegionOriginY;
        }
    }

    for (int k = 0; k < depth; ++k)
    {
        int n = x[k].size();

        Debug::Info("Drawing contour of z = " + STR(k));

        if (n == 0) continue;

        float xPrev = x[k][n-1];
        float yPrev = y[k][n-1];

        for (int i = 0; i < n; ++i)
        {
            float xCur = x[k][i];
            float yCur = y[k][i];

            restoreImageLine(xPrev,
                             xCur,
                             yPrev,
                             yCur,
                             k,
                             k);

            if (k > 0 and x[k-1].size() > 0)
            {
                int m = x[k-1].size();

                int closestIdx = -1;
                int nextClosestIdx = -1;

                float closestDist = std::sqrt(activeRegionWidth*activeRegionWidth +
                                              activeRegionHeight*activeRegionHeight +
                                              activeRegionDepth*activeRegionDepth);

                float nextClosestDist = closestDist;

                for (int j = 0; j < m; ++j)
                {
                    float xAbove = x[k-1][j];
                    float yAbove = y[k-1][j];

                    float newDist = std::sqrt(std::pow(xAbove - xCur, 2.0) +
                                              std::pow(yAbove - yCur, 2.0) +
                                              1);

                    if (newDist < nextClosestDist)
                    {
                        nextClosestDist = newDist;
                        nextClosestIdx = j;

                        if (nextClosestDist < closestDist)
                        {
                            std::swap(closestDist, nextClosestDist);
                            std::swap(closestIdx, nextClosestIdx);
                        }
                    }
                }

                restoreImageLine(x[k-1][closestIdx],
                                 xCur,
                                 y[k-1][closestIdx],
                                 yCur,
                                 k-1,
                                 k);

                restoreImageLine(x[k-1][nextClosestIdx],
                                 xCur,
                                 y[k-1][nextClosestIdx],
                                 yCur,
                                 k-1,
                                 k);

            }


            xPrev = xCur;
            yPrev = yCur;
        }
    }

    mainWin->forceActiveRegionChange(beginX, endX,
                                     beginY, endY,
                                     beginZ, endZ);

    currentPage = (beginZ + endZ) / 2;

    updateViewImage();
    updateViewPhi();
    updateViewSegmented();
}

void NumericToolControl::ripleysKFunctionSlot()
{
    Debug::Info("NumericToolControl::ripleysKFunctionSlot: Entering");

    float maxRipleySdVal = 200.0;
    float intensityEstimateRadius = 30.0;

    QList<QTreeWidgetItem *> unfilteredSelectedObjects = cellTree->selectedItems();

    int N = unfilteredSelectedObjects.length();

    bool GeneratedPoints = (N == 0);

    QList<QTreeWidgetItem *> selectedObjects;

    // Filter objects start

    for (int i = 0; i < N; ++i)
    {
        CellObject *object = (CellObject *) unfilteredSelectedObjects[i];

        if (object->isRootObject())
        {
            Debug::Warning("NumericToolControl::computeRipleysKFunctionSDRead"
                           "DataSlot: You have selected a root object - "
                           "terminating function call");

            return;
        }
        else
        {
            unsigned int n = object->getN();
            unsigned int *X = object->getX();
            unsigned int *Y = object->getY();
            unsigned int *Z = object->getZ();

            float xCenterOfMass = 0;
            float yCenterOfMass = 0;
            float zCenterOfMass = 0;

            for (int j = 0; j < n; ++j)
            {
                xCenterOfMass += X[j];
                yCenterOfMass += Y[j];
                zCenterOfMass += Z[j];
            }

            float x = std::round(xCenterOfMass / (float) n) - sdOriginX - imageBeginX;
            float y = std::round(yCenterOfMass / (float) n) - sdOriginY - imageBeginY;
            float z = std::round(zCenterOfMass / (float) n) - sdOriginZ - imageBeginZ;

            float sdVal = signedDistanceMap[gIndex((int) x,
                                                   (int) y,
                                                   (int) z,
                                                   sdHeight,
                                                   sdWidth)];

            if (sdVal != 0 and sdVal < maxRipleySdVal)
            {
                selectedObjects.append(object);
            }
        }
    }

    N = selectedObjects.length();

    if (GeneratedPoints)
    {
        N = 350*sdSourceObjects.size();
    }

    // Filter objects done

    float *distanceMatrix = (float *) malloc(N*N * sizeof(float));

    float *p = (float *) malloc(N*3 * sizeof(float));

    float A = 0.0;

    for (int k = 0; k < sdDepth; ++k)
    {
        for (int j = 0; j < sdHeight; ++j)
        {
            for (int i = 0; i < sdWidth; ++i)
            {
                int gIdx = gIndex(i,j,k,sdHeight,sdWidth);

                float sdVal = signedDistanceMap[gIdx];

                if (sdVal != 0 and sdVal < maxRipleySdVal)
                {
                    A += 1.0;
                }
            }
        }
    }



    float maxDist = 0;

    if (GeneratedPoints)
    {
        int n = 0;
        while (n < N)
        {
            p[gIndex(0,n,3)] = Random::randU(Random::generator) * sdWidth;
            p[gIndex(1,n,3)] = Random::randU(Random::generator) * sdHeight;
            p[gIndex(2,n,3)] = Random::randU(Random::generator) * sdDepth;

            float sdVal = signedDistanceMap[gIndex((int) p[gIndex(0,n,3)],
                                                   (int) p[gIndex(1,n,3)],
                                                   (int) p[gIndex(2,n,3)],
                                                   sdHeight,
                                                   sdWidth)];

            if (sdVal != 0 and sdVal < maxRipleySdVal) ++n;
        }
    }
    else
    {
        for (int i = 0; i < N; ++i)
        {
            CellObject *object = (CellObject *) selectedObjects[i];

            if (object->isRootObject())
            {
                Debug::Warning("NumericToolControl::computeRipleysKFunctionSDRead"
                               "DataSlot: You have selected a root object - "
                               "terminating function call");

                free(distanceMatrix);
                free(p);

                return;
            }
            else
            {
                unsigned int n = object->getN();
                unsigned int *X = object->getX();
                unsigned int *Y = object->getY();
                unsigned int *Z = object->getZ();

                float xCenterOfMass = 0;
                float yCenterOfMass = 0;
                float zCenterOfMass = 0;

                for (int j = 0; j < n; ++j)
                {
                    xCenterOfMass += X[j];
                    yCenterOfMass += Y[j];
                    zCenterOfMass += Z[j];
                }

                p[gIndex(0,i,3)] = std::round(xCenterOfMass / (float) n);
                p[gIndex(1,i,3)] = std::round(yCenterOfMass / (float) n);
                p[gIndex(2,i,3)] = std::round(zCenterOfMass / (float) n);
            }
        }
    }

    auto D = [&p](const int i, const int j)
    {
        return std::sqrt(std::pow(p[gIndex(0,i,3)] - p[gIndex(0,j,3)], 2.0) +
                         std::pow(p[gIndex(1,i,3)] - p[gIndex(1,j,3)], 2.0) +
                         std::pow(p[gIndex(2,i,3)] - p[gIndex(2,j,3)], 2.0));
    };

    for (int j = 0; j < N; ++j)
    {
        for (int i = j; i < N; ++i)
        {
            if (i == j)
            {
                int gIdx = gIndex(i,j,N);

                distanceMatrix[gIdx] = 0;
            }
            else
            {
                int gIdx = gIndex(i,j,N);
                int gIdxFlip = gIndex(j,i,N);

                float newDist = D(i,j);

                maxDist = std::max(maxDist, newDist);

                distanceMatrix[gIdx] = newDist;
                distanceMatrix[gIdxFlip] = newDist;
            }
        }
    }

    float stepSize = 0.1;
    int numSteps = (int) (maxDist / stepSize) + 1;

    float *allRipleysFunctions = (float *) malloc(N*numSteps * sizeof(float));
    float *intensity = (float *) malloc(N * sizeof(float));

    std::vector<std::vector<int>> idxMatrix(N, std::vector<int>(N));
    std::vector<std::vector<float>> distances(N, std::vector<float>(N));

    for (int j = 0; j < N; ++j)
    {

        for (int i = 0; i < N; ++i)
        {
            distances[j][i] = distanceMatrix[gIndex(i,j,N)];
        }

        sort_indexes(distances[j], idxMatrix[j]);

        int pIdx = 0;

        while (pIdx < N and distances[j][idxMatrix[j][pIdx]] < intensityEstimateRadius)
        {
            ++pIdx;
        }

        //intensity[j] = ((float) pIdx) / (4.0 / 3.0 * M_PI * std::pow(intensityEstimateRadius, 3.0)-4.0 / 3.0 * M_PI * std::pow(10.0, 3.0));
        intensity[j] = ((float) N) / A;

        //std::sort(distances.begin(), distances.end());
        //intensityEstimateRadius
    }

    for (int j = 0; j < N; ++j)
    {

        int pIdx = 1;
        float val = 0.0;

        for (int i = 0; i < numSteps; ++i)
        {
            float t = i * stepSize;

            while (pIdx < N and distances[j][idxMatrix[j][pIdx]] < t)
            {
                val += 1.0 / intensity[pIdx];//A / ((float) N); // * 1.0
                //val += 1.0 / (intensity[pIdx] * intensity[j]);
                ++pIdx;
            }

            allRipleysFunctions[gIndex(i,j,numSteps)] = val;
        }
    }

    std::string dataString = "";

    for (int i = 0; i < numSteps; ++i)
    {
        float val = 0.0;

        for (int j = 0; j < N; ++j)
        {
            val += allRipleysFunctions[gIndex(i,j,numSteps)];
        }

        val = val / ((float) (N));

        dataString += STR(val);
        if (i < numSteps - 1) dataString += ", ";
    }

    free(allRipleysFunctions);
    free(intensity);
    free(distanceMatrix);
    free(p);

    Debug::Info("Data:\n" + dataString);

    Debug::Info("NumericToolControl::ripleysKFunctionSlot: Leaving");
}

void NumericToolControl::ffmTest()
{
    // -2 source
    // -1 bound
    // 1 normal speed
    // 0 'infinite speed'

    int width = 100;
    int height = 100;
    float radius = 5.000001;
    std::string sdm_data_string = "sdm data:\n";
    std::string error_data_string = "error data:\n";

    float *sdm = (float *) malloc(width*height * sizeof(float));
    float *groundTruth = (float *) malloc(width*height * sizeof(float));

    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            //if (i == 20 and j == 20)
            if (std::sqrt( std::pow((float) (i-20), 2.0) + std::pow((float) (j-20), 2.0) ) <= radius*radius)
            {
                sdm[gIndex(i,j,width)] = -2;
            }
            else
            {
                sdm[gIndex(i,j,width)] = 1;
            }
            groundTruth[gIndex(i,j,width)] = std::max( std::sqrt( std::pow((float) (i-20), 2.0) + std::pow((float) (j-20), 2.0) ) - radius*radius, 0.0 );
        }
    }

    FastMarching(width, height, 1, sdm);

    sdm_data_string += "[";

    for (int j = 0; j < height; ++j)
    {
        sdm_data_string += "[";
        for (int i = 0; i < width; ++i)
        {
            sdm_data_string += STR(sdm[gIndex(i,j,width)]);

            if (i < width - 1)
            {
                sdm_data_string += ", ";
            }
        }

        sdm_data_string += "]";

        if (j < height - 1)
        {
            sdm_data_string += ",\n";
        }
    }

    sdm_data_string += "]";

    error_data_string += "[";

    for (int j = 0; j < height; ++j)
    {
        error_data_string += "[";
        for (int i = 0; i < width; ++i)
        {
            error_data_string += STR(sdm[gIndex(i,j,width)] - groundTruth[gIndex(i,j,width)]);

            if (i < width - 1)
            {
                error_data_string += ", ";
            }
        }

        error_data_string += "]";

        if (j < height - 1)
        {
            error_data_string += ",\n";
        }
    }

    error_data_string += "]";

    Debug::Info(sdm_data_string);
    Debug::Info(error_data_string);
}