#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QActionGroup>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QTreeWidget>

#include "CellObject.hpp"
#include "MouseoverScrollArea.hpp"

class QAction;
class QMenu;
class QPlainTextEdit;
class QSessionManager;
class QScrollArea;
class QScrollBar;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

    void loadFile(const QString &fileName);

    void setFillImage(const int width, const int height, const int z,
                      const float *arrayImage,
                      const float *activeRegionImage = NULL,
                      bool activeRegionOnly = false);

    void updateSegmented(const int width, const bool *segmented);

    void updatePhi(const float *arrayImage, float newPhiZero = 0);

    void updateSignedMap(const float *arrayImage, int sdOriginX,
                         int sdOriginY, int sdWidth, int sdHeight,
                         float signedZero = 0);

    void InitViewArea(const int width, const int height, const int depth);

    void createTree(QTreeWidget *treeWidget);

    void forceActiveRegionChange(int beginX, int endX,
                                 int beginY, int endY,
                                 int beginZ, int endZ);

protected:
    void closeEvent(QCloseEvent *event) override;

public slots:
    void setMessage(char const message[]);

    void scrollAreaMouseClick(float x, float y);

private slots:
    void newFile();
    void open();
    bool save();
    bool saveAs();
    void about();
    void documentWasModified();
    void zoomIn();
    void zoomOut();
    void nextImage();
    void prevImage();
    void SDsetPassThrough();
    void runCV();
    void setPos();
    void setNeg();
    void runBoundedSignedDistanceMap();
    void setDisallowed();
    void swapPhiSign();
    void reinitPhi();
    void showImage(bool activeRegionOnly = false);
    void showPhi();
    void showSigned();
    void renderIn3D();
    void toggleBoundary();
    void toggleSegmented();
    void setActiveRegion();
    void storeObject();
    void loadObject();
    void deleteObject();


#ifndef QT_NO_SESSIONMANAGER
    void commitData(QSessionManager &);
#endif

private:
    void redrawBackground(bool activeRegionOnly = false);
    void createActions();
    void createStatusBar();
    void readSettings();
    void writeSettings();
    bool maybeSave();
    bool saveFile(const QString &fileName);
    void setCurrentFile(const QString &fileName);
    void scaleImage(float factor);
    void adjustScrollBar(QScrollBar *scrollBar, float factor);
    void generatePhiBoundary();
    void drawSelectionBorder();
    void drawActiveRegionBorder();
    void updateCheckedActions();
    void forceSelectionIntoActiveRegion();

    QString strippedName(const QString &fullFileName);
    QImage ArrayImageToQImage(const int width, const int height,
                              const float *image);
    QPlainTextEdit *textEdit;
    QString curFile;

    QPixmap pixmap;

    QImage image;
    QImage segmented;
    QImage phi;
    QImage signedMap;
    QImage boundary;
    QImage selection;
    QImage active;

    float phiZero;

    int imageWidth;
    int imageHeight;
    int imageDepth;

    int currentZ;

    bool activeRegionPresent;
    bool imageShown;
    bool phiShown;
    bool signedShown;
    bool selectionShown;

    bool signedMapCalculated;

    QLabel *imageLabel;

    QLineEdit *blightLineEdit;
    QLineEdit *iterLineEdit;
    QLineEdit *muLineEdit;
    QLineEdit *nuLineEdit;
    QLineEdit *c1LineEdit;
    QLineEdit *c2LineEdit;

    MouseoverScrollArea *scrollArea;
    float scaleFactor;
    QDockWidget *cellDockWidget;

    bool selectClickedOnce;
    int selectionBeginX;
    int selectionEndX;
    int selectionBeginY;
    int selectionEndY;
    int selectionBeginZ;
    int selectionEndZ;

    int activeRegionBeginX;
    int activeRegionEndX;
    int activeRegionBeginY;
    int activeRegionEndY;
    int activeRegionBeginZ;
    int activeRegionEndZ;

    QActionGroup *selectGroup;
    QActionGroup *viewGroup;

    QAction *selectAct;
    QAction *vesicleCenterAct;
    QAction *neurofilamentCenterAct;
    QAction *blightSelectAct;
    QAction *restoreSelectAct;
    QAction *boundaryAct;
    QAction *segmentedAct;

    QAction *c1lockAct;
    QAction *c2lockAct;

    QAction *imageAct;
    QAction *phiAct;
    QAction *signedAct;

    QAction *zoomInAct;
    QAction *zoomOutAct;
    QAction *sd_runAct;
    QAction *sd_setSourcesAct;
    QAction *sd_setBoundariesAct;
    QAction *sd_setPassThroughAct;
    QAction *sd_readDataAct;
    QAction *sd_readFullDataAct;
    QAction *ripleyKAct;

    QAction *nextImageAct;
    QAction *prevImageAct;

signals:
    void nextImageSig();

    void prevImageSig();

    void newActiveRegionSig(int beginX, int endX,
                            int beginY, int endY,
                            int beginZ, int endZ,
                            float sigma);

    void runChanVeseSig(int numIterations,
                        float mu, float nu,
                        float c1, float c2,
                        bool constc1, bool constc2);

    void markVesicleLocationSig(int x, int y);

    void markNeurofilamentSig(int x, int y);

    void finalizeNeurofilamentSig();

    void setPosSig(int beginX, int endX,
                   int beginY, int endY,
                   int beginZ, int endZ);

    void setNegSig(int beginX, int endX,
                   int beginY, int endY,
                   int beginZ, int endZ);

    void setDisallowedSig(int beginX, int endX,
                          int beginY, int endY,
                          int beginZ, int endZ);

    void swapPhiSignSig();

    void reinitPhiSig();

    void renderIn3DSig();

    void storeObjectSig();

    void loadObjectSig();

    void deleteObjectSig();

    void blightSig(int imageX, int imageY, float sigma);

    void restoreSig(int imageX, int imageY, float sigma);

    void SDSetSourcesSig();

    void SDSetBoundariesSig();

    void SDSetPassThroughSig(float travelCost);

    void SDRunBoundedSignedDistanceMapSig();

    void SDReadDataSig();

    void SDReadFullDataSig();

    void ripleysKFunctionSig();

    void runCustomFunctionSig();
};

#endif