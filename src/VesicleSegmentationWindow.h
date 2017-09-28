#ifndef VESICLE_SEGMENTATION_WINDOW_H
#define VESICLE_SEGMENTATION_WINDOW_H

#include <QActionGroup>
#include <QImage>
#include <QLabel>
#include <QMainWindow>

#include "MouseoverScrollArea.hpp"

class QAction;
class QMenu;
class QPlainTextEdit;
class QSessionManager;
class QScrollArea;
class QScrollBar;

class VesicleSegmentationWindow : public QMainWindow
{
    Q_OBJECT

public:
    VesicleSegmentationWindow();

    void setFillImage(const float *arrayImage);

    void updateSegmented();

    void InitViewArea(const int width, const int height);

protected:
    void closeEvent(QCloseEvent *event) override;

public slots:
    void setMessage(char const message[]);

    void setCompletion(float completionValue, float minCompletionValue);

    void scrollAreaMouseClick(float x, float y);

private slots:
    void zoomIn();
    void zoomOut();
    void showImage();
    void undo();
    void nextSlice();

private:
    void redrawBackground();
    void createActions();
    void createStatusBar();
    void scaleImage(float factor);
    void adjustScrollBar(QScrollBar *scrollBar, float factor);

    QImage ArrayImageToQImage(const int width, const int height,
                              const float *image);

    QPixmap pixmap;

    QImage image;
    QImage segmented;

    int imageWidth;
    int imageHeight;

    QLabel *imageLabel;
    QLabel *completionLabel;

    MouseoverScrollArea *scrollArea;
    float scaleFactor;

    std::vector<float> segmentingPointX;
    std::vector<float> segmentingPointY;

    QActionGroup *segmentationGroup;

    QAction *addAct;
    QAction *moveAct;
    QAction *removeAct;

    QAction *zoomInAct;
    QAction *zoomOutAct;

signals:

    void storeSegmentationSig(std::vector<float> x, std::vector<float> y);

    void nextSig();
};

#endif // VESICLE_SEGMENTATION_WINDOW_H