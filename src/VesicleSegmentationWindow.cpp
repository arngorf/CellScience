#include "VesicleSegmentationWindow.h"

#include <QtWidgets>

#include "Debug.hpp"

VesicleSegmentationWindow::VesicleSegmentationWindow() :
     imageLabel(new QLabel)
   , completionLabel(new QLabel)
   , scrollArea(new MouseoverScrollArea)
   , scaleFactor(10.0)
   , segmentingPointX(0)
   , segmentingPointY(0)
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

    QMenu *vesicleMenu = menuBar()->addMenu(tr("&VesSeg"));
    //QMenu *viewMenu = menuBar()->addMenu(tr("&View"));
    //QMenu *toolsMenu = menuBar()->addMenu(tr("&Tools"));
    QToolBar *topToolBar = addToolBar(tr("Top"));

    segmentationGroup = new QActionGroup(this);

    QAction *saveAct = new QAction(tr("&Save (and next)"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save and get next vesicle"));
    connect(saveAct, &QAction::triggered, this, &VesicleSegmentationWindow::nextSlice);
    vesicleMenu->addAction(saveAct);

    QAction *undoAct = new QAction(tr("Undo&z..."), this);
    undoAct->setShortcuts(QKeySequence::Undo);
    undoAct->setStatusTip(tr("Undo last add action"));
    connect(undoAct, &QAction::triggered, this, &VesicleSegmentationWindow::undo);
    vesicleMenu->addAction(undoAct);

    const QIcon exitIcon = QIcon::fromTheme("application-exit");
    QAction *exitAct = vesicleMenu->addAction(exitIcon, tr("&QuitSeg"), this, &QWidget::close);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("End segmentation"));

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

    topToolBar->addSeparator();

    completionLabel = new QLabel(tr("0.0"));
    topToolBar->addWidget(completionLabel);

    topToolBar->addSeparator();

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

void VesicleSegmentationWindow::undo()
{
    int n = segmentingPointX.size();

    if (n > 0)
    {
        segmentingPointX.resize(n-1);
        segmentingPointY.resize(n-1);

        updateSegmented();
    }
}

void VesicleSegmentationWindow::nextSlice()
{
    emit storeSegmentationSig(segmentingPointX, segmentingPointY);

    segmentingPointX.resize(0);
    segmentingPointY.resize(0);

    updateSegmented();

    emit nextSig();
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

    QString shownName = "Vesicle Segmentation Window";
    setWindowFilePath(shownName);

    Debug::Info("VesicleSegmentationWindow::setFillImage: Leaving");
}

void VesicleSegmentationWindow::updateSegmented()
{
    Debug::Info("VesicleSegmentationWindow::updateSegmented: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    newQImage.fill(QColor(0,0,0,0));

    for (int n = 0; n < segmentingPointX.size(); ++n)
    {
        int i = (int) segmentingPointX[n];
        int j = (int) segmentingPointY[n];
        float di = segmentingPointX[n] - i;
        float dj = segmentingPointY[n] - j;

        int v00 = (int) ((1.0-di)*(1.0-dj)*255.0);
        int v10 = (int) ((    di)*(1.0-dj)*255.0);
        int v01 = (int) ((1.0-di)*(    dj)*255.0);
        int v11 = (int) ((    di)*(    dj)*255.0);

        newQImage.setPixelColor(i  , j  , QColor(255, 0, 0, v00));
        newQImage.setPixelColor(i+1, j  , QColor(255, 0, 0, v10));
        newQImage.setPixelColor(i  , j+1, QColor(255, 0, 0, v01));
        newQImage.setPixelColor(i+1, j+1, QColor(255, 0, 0, v11));
    }

    segmented = newQImage;

    redrawBackground();

    Debug::Info("VesicleSegmentationWindow::updateSegmented: Leaving");
}

void VesicleSegmentationWindow::InitViewArea(const int width, const int height)
{
    Debug::Info("VesicleSegmentationWindow::InitViewArea: Entering");

    imageWidth = width;
    imageHeight = height;

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

void VesicleSegmentationWindow::setCompletion(float completionValue, float minCompletionValue)
{
    QString text = QString(STR(completion))
                 + QString(" (")
                 + QString(STR(minCompletionValue))
                 + QString(")");

    completionLabel->setText(text);
}

void VesicleSegmentationWindow::scrollAreaMouseClick(float x, float y)
{
    float dx = (float) scrollArea->horizontalScrollBar()->value();
    float dy = (float) scrollArea->verticalScrollBar()->value();

    float imageX = (x+dx)/imageLabel->width() * imageWidth;
    float imageY = (y+dy)/imageLabel->height() * imageHeight;

    if (imageX >= 0 and imageY >= 0 and imageX < imageWidth and imageY < imageHeight)
    {
        if (addAct->isChecked())
        {
            Debug::Info("Add: " + STR(imageX) + ", " + STR(imageY));

            segmentingPointX.push_back(imageX);
            segmentingPointY.push_back(imageY);

            updateSegmented();
            //redrawBackground();
        }
        else if (moveAct->isChecked())
        {
            Debug::Info("Move: " + STR(imageX) + ", " + STR(imageY));
        }
        else if (removeAct->isChecked())
        {
            Debug::Info("Remove: " + STR(imageX) + ", " + STR(imageY));
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
    Debug::Info("VesicleSegmentationWindow::scaleImage: Entering");

    Q_ASSERT(imageLabel->pixmap());
    scaleFactor *= factor;
    imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);

    scrollArea->setScale(scaleFactor);

    Debug::Info("VesicleSegmentationWindow::scaleImage: Leaving");
}

