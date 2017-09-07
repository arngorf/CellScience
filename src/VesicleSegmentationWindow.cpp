#include "VesicleSegmentationWindow.h"

#include <QtWidgets>

#include "Debug.hpp"

VesicleSegmentationWindow::VesicleSegmentationWindow() :
     segmentedFlag(NULL)
   , imageLabel(new QLabel)
   , scrollArea(new MouseoverScrollArea)
   , scaleFactor(1)
{

    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setVisible(false);
    setCentralWidget(scrollArea);

    createActions();
    createStatusBar();

    connect(scrollArea, SIGNAL( scrollAreaMouseClick(float, float) ),
            this,       SLOT  ( scrollAreaMouseClick(float, float) )
            );

}

void VesicleSegmentationWindow::closeEvent(QCloseEvent *event)
{
    event->accept();
}

void VesicleSegmentationWindow::createActions()
{
    Debug::Info("VesicleSegmentationWindow::createActions: Entering");

    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    QToolBar *topToolBar = addToolBar(tr("Top"));
    QMenu *viewMenu = menuBar()->addMenu(tr("&View"));
    QMenu *toolsMenu = menuBar()->addMenu(tr("&Tools"));

    segmentationGroup = new QActionGroup(this);

    // Add icon
    QIcon addIcon = QIcon::fromTheme("document-new", QIcon(":/images/select.png"));
    addAct = new QAction(addIcon, tr("&Select"), this);
    addAct->setStatusTip(tr("Add boundary point"));
    addAct->setCheckable(true);
    addAct->setChecked(true);
    segmentationGroup->addAction(addAct);
    topToolBar->addAction(addAct);

    // Vesicle center icon
    QIcon moveIcon = QIcon::fromTheme("document-new", QIcon(":/images/selectVesicle.png"));
    moveAct = new QAction(moveIcon, tr("&Select vesicle"), this);
    moveAct->setStatusTip(tr("Annotate the center of a vesicle"));
    moveAct->setCheckable(true);
    segmentationGroup->addAction(moveAct);
    topToolBar->addAction(moveAct);

    // Vesicle center icon
    QIcon removeIcon = QIcon::fromTheme("document-new", QIcon(":/images/markNeurofilament.png"));
    removeAct = new QAction(removeIcon, tr("&Select neurofilament"), this);
    removeAct->setStatusTip(tr("Annotate the center of a neurofilament"));
    removeAct->setCheckable(true);
    segmentationGroup->addAction(removeAct);
    topToolBar->addAction(removeAct);

    Debug::Info("VesicleSegmentationWindow::createActions: Leaving");
}

void VesicleSegmentationWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void VesicleSegmentationWindow::redrawBackground()
{
    showImage();
}

void VesicleSegmentationWindow::showImage()
{
    QPainter p(&pixmap);
    p.setOpacity(1.0);
    p.setCompositionMode(QPainter::CompositionMode_SourceOver);

    p.drawImage(QPoint(0, 0), image);

    p.drawImage(QPoint(0, 0), segmented);

    imageLabel->setPixmap(pixmap);
}

void VesicleSegmentationWindow::setFillImage(const float *arrayImage)
{
    Debug::Info("VesicleSegmentationWindow::setFillImage: Entering");

    image = ArrayImageToQImage(imageWidth, imageHeight, arrayImage);

    QPixmap newPixmap(image.size());
    pixmap = newPixmap;

    scrollArea->setVisible(true);

    redrawBackground();

    scaleImage(1.0);

    Debug::Info("VesicleSegmentationWindow::setFillImage: Leaving");
}

void VesicleSegmentationWindow::updateSegmented()
{
    Debug::Info("VesicleSegmentationWindow::updateSegmented: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    newQImage.fill(QColor(0,0,0,0));

    if (segmentedFlag != NULL)
    {
        for (int j = 0; j < imageHeight; ++j)
        {
            for (int i = 0; i < imageWidth; ++i)
            {
                if (segmentedFlag[j*imageWidth + i])
                {
                    newQImage.setPixelColor(i, j, QColor(75, 50, 25, 100));
                }
            }
        }
    }

    segmented = newQImage;

    redrawBackground();

    Debug::Info("VesicleSegmentationWindow::updateSegmented: Leaving");
}

void VesicleSegmentationWindow::InitViewArea(const int width, const int height)
{
    Debug::Info("VesicleSegmentationWindow::InitViewArea: Entering");

    scaleFactor = 1.0;

    imageWidth = width;
    imageHeight = height;

    bool *segmentedFlag = (bool *) malloc(width*height * sizeof(bool));

    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            segmentedFlag[gIndex(i,j,height)] = false;
        }
    }

    Debug::Info("VesicleSegmentationWindow::InitViewArea: Leaving");
}

QImage VesicleSegmentationWindow::ArrayImageToQImage(const int width, const int height, const float *arrayImage)
{
    Debug::Info("VesicleSegmentationWindow::ArrayImageToQImage: Entering");

    QImage newQImage(width, height, QImage::Format_ARGB32);

    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            int pixelVal = (int) arrayImage[j*width + i];
            newQImage.setPixelColor(i, j, QColor(pixelVal, pixelVal, pixelVal, 255));
        }
    }

    Debug::Info("VesicleSegmentationWindow::ArrayImageToQImage: Leaving");
    return newQImage;
}

void VesicleSegmentationWindow::zoomIn()
{
    scaleImage(1.25);
}

void VesicleSegmentationWindow::zoomOut()
{
    scaleImage(0.8);
}

void VesicleSegmentationWindow::setMessage(char const message[])
{
    statusBar()->showMessage(tr(message));
}

void VesicleSegmentationWindow::scrollAreaMouseClick(float x, float y)
{
    float dx = (float) scrollArea->horizontalScrollBar()->value();
    float dy = (float) scrollArea->verticalScrollBar()->value();

    int imageX = int((x+dx)/imageLabel->width() * imageWidth);
    int imageY = int((y+dy)/imageLabel->height() * imageHeight);

    if (imageX >= 0 and imageY >= 0 and imageX < imageWidth and imageY < imageHeight)
    {
        if (addAct->isChecked())
        {

            //redrawBackground();
        }
        else if (moveAct->isChecked())
        {

        }
        else if (removeAct->isChecked())
        {

        }
    }

}

void VesicleSegmentationWindow::adjustScrollBar(QScrollBar *scrollBar, float factor)
{
    scrollBar->setValue(int(factor * scrollBar->value()
                            + ((factor - 1) * scrollBar->pageStep()/2)));
}

void VesicleSegmentationWindow::scaleImage(float factor)
{
    Q_ASSERT(imageLabel->pixmap());
    scaleFactor *= factor;
    imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);

    scrollArea->setScale(scaleFactor);
}

