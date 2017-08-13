#include "MainWindow.h"

#include <QtWidgets>

#include "Debug.hpp"

MainWindow::MainWindow() :
     textEdit(new QPlainTextEdit)
   , phiZero(0)
   , activeRegionPresent(false)
   , imageShown(true)
   , phiShown(false)
   , signedShown(false)
   , selectionShown(false)
   , signedMapCalculated(false)
   , imageLabel(new QLabel)
   , blightLineEdit(new QLineEdit("2.5"))
   , iterLineEdit(new QLineEdit("2500"))
   , muLineEdit(new QLineEdit("1.0"))
   , nuLineEdit(new QLineEdit("15.0"))
   , c1LineEdit(new QLineEdit("130.0"))
   , c2LineEdit(new QLineEdit("173.0"))
   , scrollArea(new MouseoverScrollArea)
   , scaleFactor(1)
   , selectClickedOnce(false)
{
    //setCentralWidget(textEdit);

    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setVisible(false);
    setCentralWidget(scrollArea);

    createActions();
    createStatusBar();

    readSettings();

    connect(scrollArea, SIGNAL( scrollAreaMouseClick(float, float) ),
            this,       SLOT  ( scrollAreaMouseClick(float, float) )
            );

    connect(textEdit->document(), &QTextDocument::contentsChanged,
            this, &MainWindow::documentWasModified);

#ifndef QT_NO_SESSIONMANAGER
    QGuiApplication::setFallbackSessionManagementEnabled(false);
    connect(qApp, &QGuiApplication::commitDataRequest,
            this, &MainWindow::commitData);
#endif

    setCurrentFile(QString());
    setUnifiedTitleAndToolBarOnMac(true);
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    if (maybeSave()) {
        writeSettings();
        event->accept();
    } else {
        event->ignore();
    }
}

void MainWindow::newFile()
{
    if (maybeSave()) {
        textEdit->clear();
        setCurrentFile(QString());
    }
}

void MainWindow::open()
{
    if (maybeSave()) {
        QString fileName = QFileDialog::getOpenFileName(this);
        if (!fileName.isEmpty())
            loadFile(fileName);
    }
}

bool MainWindow::save()
{
    if (curFile.isEmpty()) {
        return saveAs();
    } else {
        return saveFile(curFile);
    }
}

bool MainWindow::saveAs()
{
    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if (dialog.exec() != QDialog::Accepted)
        return false;
    return saveFile(dialog.selectedFiles().first());
}

void MainWindow::about()
{
   QMessageBox::about(this, tr("About Application"),
            tr("The <b>Application</b> example demonstrates how to "
               "write modern GUI applications using Qt, with a menu bar, "
               "toolbars, and a status bar."));
}

void MainWindow::documentWasModified()
{
    setWindowModified(textEdit->document()->isModified());
}

void MainWindow::createActions()
{
    Debug::Info("MainWindow::createActions: Entering");

    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    QToolBar *topToolBar = addToolBar(tr("Top"));
    QMenu *viewMenu = menuBar()->addMenu(tr("&View"));
    QMenu *toolsMenu = menuBar()->addMenu(tr("&Tools"));


    selectGroup = new QActionGroup(this);

    // Select icon
    QIcon selectIcon = QIcon::fromTheme("document-new", QIcon(":/images/select.png"));
    selectAct = new QAction(selectIcon, tr("&Select"), this);
    selectAct->setStatusTip(tr("Select region of image"));
    selectAct->setCheckable(true);
    selectAct->setChecked(true);
    selectGroup->addAction(selectAct);
    topToolBar->addAction(selectAct);

    // Vesicle center icon
    QIcon vesicleCenterIcon = QIcon::fromTheme("document-new", QIcon(":/images/selectVesicle.png"));
    vesicleCenterAct = new QAction(vesicleCenterIcon, tr("&Select vesicle"), this);
    vesicleCenterAct->setStatusTip(tr("Annotate the center of a vesicle"));
    vesicleCenterAct->setCheckable(true);
    selectGroup->addAction(vesicleCenterAct);
    topToolBar->addAction(vesicleCenterAct);

    // Vesicle center icon
    QIcon neurofilamentCenterIcon = QIcon::fromTheme("document-new", QIcon(":/images/markNeurofilament.png"));
    neurofilamentCenterAct = new QAction(neurofilamentCenterIcon, tr("&Select neurofilament"), this);
    neurofilamentCenterAct->setStatusTip(tr("Annotate the center of a neurofilament"));
    neurofilamentCenterAct->setCheckable(true);
    selectGroup->addAction(neurofilamentCenterAct);
    topToolBar->addAction(neurofilamentCenterAct);

    // blight icon
    QIcon blightSelectIcon = QIcon::fromTheme("document-new", QIcon(":/images/blight.png"));
    blightSelectAct = new QAction(blightSelectIcon, tr("&remove"), this);
    blightSelectAct->setStatusTip(tr("Blight region of image"));
    blightSelectAct->setCheckable(true);
    selectGroup->addAction(blightSelectAct);
    topToolBar->addAction(blightSelectAct);

    // restore icon
    QIcon restoreSelectIcon = QIcon::fromTheme("document-new", QIcon(":/images/restore.png"));
    restoreSelectAct = new QAction(restoreSelectIcon, tr("&restore"), this);
    restoreSelectAct->setStatusTip(tr("restore region of image"));
    restoreSelectAct->setCheckable(true);
    selectGroup->addAction(restoreSelectAct);
    topToolBar->addAction(restoreSelectAct);

    QLabel *blightLabel = new QLabel(tr("Tool Sigma"));
    blightLineEdit->setMaximumWidth(54);
    blightLineEdit->setValidator( new QDoubleValidator(0, 255, 1, this) );
    topToolBar->addWidget(blightLabel);
    topToolBar->addWidget(blightLineEdit);

    topToolBar->addSeparator();

    // Set Chan Vese region icon
    QIcon setcvIcon = QIcon::fromTheme("document-new", QIcon(":/images/setcv.png"));
    QAction *setcvAct = new QAction(setcvIcon, tr("&Set active region"), this);
    setcvAct->setStatusTip(tr("Set active region for segementation"));
    connect(setcvAct, &QAction::triggered, this, &MainWindow::setActiveRegion);
    topToolBar->addAction(setcvAct);

    QLabel *iterLabel = new QLabel(tr("Iter"));
    iterLineEdit->setMaximumWidth(54);
    iterLineEdit->setValidator( new QIntValidator(0, 100000, this) );
    topToolBar->addWidget(iterLabel);
    topToolBar->addWidget(iterLineEdit);

    QLabel *muLabel = new QLabel(tr("mu"));
    muLineEdit->setMaximumWidth(54);
    muLineEdit->setValidator( new QDoubleValidator(0, 10000, 2, this) );
    topToolBar->addWidget(muLabel);
    topToolBar->addWidget(muLineEdit);

    QLabel *nuLabel = new QLabel(tr("nu"));
    nuLineEdit->setMaximumWidth(54);
    nuLineEdit->setValidator( new QDoubleValidator(0, 10000, 2, this) );
    topToolBar->addWidget(nuLabel);
    topToolBar->addWidget(nuLineEdit);

    QLabel *c1Label = new QLabel(tr("c1"));
    c1LineEdit->setMaximumWidth(54);
    c1LineEdit->setValidator( new QDoubleValidator(0, 255, 1, this) );
    topToolBar->addWidget(c1Label);
    topToolBar->addWidget(c1LineEdit);

    // Lock c1 icon
    QIcon c1lockIcon = QIcon::fromTheme("document-new", QIcon(":/images/lock.png"));
    c1lockAct = new QAction(c1lockIcon, tr("&Lock c1"), this);
    c1lockAct->setStatusTip(tr("Lock c1"));
    c1lockAct->setCheckable(true);
    topToolBar->addAction(c1lockAct);

    QLabel *c2Label = new QLabel(tr("c2"));
    c2LineEdit->setMaximumWidth(54);
    c2LineEdit->setValidator( new QDoubleValidator(0, 255, 1, this) );
    topToolBar->addWidget(c2Label);
    topToolBar->addWidget(c2LineEdit);

    // Lock c2 icon
    QIcon c2lockIcon = QIcon::fromTheme("document-new", QIcon(":/images/lock.png"));
    c2lockAct = new QAction(c2lockIcon, tr("&Lock c2"), this);
    c2lockAct->setStatusTip(tr("Lock c2"));
    c2lockAct->setCheckable(true);
    topToolBar->addAction(c2lockAct);

    // Run Chan Vese icon
    QIcon runcvIcon = QIcon::fromTheme("document-new", QIcon(":/images/runcv.png"));
    QAction *runcvAct = new QAction(runcvIcon, tr("&Run Chan-Vese"), this);
    runcvAct->setStatusTip(tr("Run Chan Vese"));
    connect(runcvAct, &QAction::triggered, this, &MainWindow::runCV);
    topToolBar->addAction(runcvAct);

    // Store Object icon
    QIcon storeIcon = QIcon::fromTheme("document-new", QIcon(":/images/store.png"));
    QAction *storeAct = new QAction(storeIcon, tr("&Store"), this);
    storeAct->setStatusTip(tr("Store cell object"));
    connect(storeAct, &QAction::triggered, this, &MainWindow::storeObject);
    topToolBar->addAction(storeAct);

    // Load Object icon
    QIcon loadIcon = QIcon::fromTheme("document-new", QIcon(":/images/load.png"));
    QAction *loadAct = new QAction(loadIcon, tr("&Load"), this);
    loadAct->setStatusTip(tr("Load cell object"));
    connect(loadAct, &QAction::triggered, this, &MainWindow::loadObject);
    topToolBar->addAction(loadAct);

    // Delete Object icon
    QIcon deleteIcon = QIcon::fromTheme("document-new", QIcon(":/images/delete.png"));
    QAction *deleteAct = new QAction(deleteIcon, tr("&Delete"), this);
    deleteAct->setStatusTip(tr("Delete cell object"));
    connect(deleteAct, &QAction::triggered, this, &MainWindow::deleteObject);
    topToolBar->addAction(deleteAct);

    // Set region positive icon
    QIcon setposIcon = QIcon::fromTheme("document-new", QIcon(":/images/setpos.png"));
    QAction *setposAct = new QAction(setposIcon, tr("&Set positive"), this);
    setposAct->setStatusTip(tr("Set region positive"));
    connect(setposAct, &QAction::triggered, this, &MainWindow::setPos);
    topToolBar->addAction(setposAct);

    // Set region positive icon
    QIcon setnegIcon = QIcon::fromTheme("document-new", QIcon(":/images/setneg.png"));
    QAction *setnegAct = new QAction(setnegIcon, tr("&Set negative"), this);
    setnegAct->setStatusTip(tr("Set region negative"));
    connect(setnegAct, &QAction::triggered, this, &MainWindow::setNeg);
    topToolBar->addAction(setnegAct);

    // Set region disallowed icon
    QIcon setdisallowedIcon = QIcon::fromTheme("document-new", QIcon(":/images/setdisallowed.png"));
    QAction *setdisallowedAct = new QAction(setdisallowedIcon, tr("&Set disallowed"), this);
    setdisallowedAct->setStatusTip(tr("Set region disallowed"));
    connect(setdisallowedAct, &QAction::triggered, this, &MainWindow::setDisallowed);
    topToolBar->addAction(setdisallowedAct);

    // Swap sign icon
    QIcon swapSignIcon = QIcon::fromTheme("document-new", QIcon(":/images/swapsign.png"));
    QAction *swapSignAct = new QAction(swapSignIcon, tr("&Set disallowed"), this);
    swapSignAct->setStatusTip(tr("Swap sign of phasefield"));
    connect(swapSignAct, &QAction::triggered, this, &MainWindow::swapPhiSign);
    topToolBar->addAction(swapSignAct);

    // re-init icon
    QIcon reinitIcon = QIcon::fromTheme("document-new", QIcon(":/images/reinit.png"));
    QAction *reinitAct = new QAction(reinitIcon, tr("&Re-init phi"), this);
    reinitAct->setStatusTip(tr("Reinitialize the phasefield to signed distance map phasefield"));
    connect(reinitAct, &QAction::triggered, this, &MainWindow::reinitPhi);
    topToolBar->addAction(reinitAct);

    topToolBar->addSeparator();

    // Zoom in icon
    const QIcon zoomInIcon = QIcon::fromTheme("document-open", QIcon(":/images/zoomin.png"));
    zoomInAct = viewMenu->addAction(zoomInIcon ,tr("Zoom &In (25%)"), this, &MainWindow::zoomIn);
    zoomInAct->setShortcut(Qt::Key_Plus);
    zoomInAct->setStatusTip(tr("Zoom in"));
    topToolBar->addAction(zoomInAct);

    // Zoom out icon
    const QIcon zoomOutIcon = QIcon::fromTheme("document-open", QIcon(":/images/zoomout.png"));
    zoomOutAct = viewMenu->addAction(zoomOutIcon ,tr("Zoom &Out (25%)"), this, &MainWindow::zoomOut);
    zoomOutAct->setShortcut(Qt::Key_Minus);
    zoomOutAct->setStatusTip(tr("Zoom out"));
    topToolBar->addAction(zoomOutAct);

    // Signed Distance: Run
    sd_runAct = toolsMenu->addAction(tr("Run sdm"), this, &MainWindow::runBoundedSignedDistanceMap);

    // Signed Distance: Set sources
    sd_setSourcesAct = toolsMenu->addAction(tr("Set sources"), this, &MainWindow::SDSetSourcesSig);

    // Signed Distance: Set boundaries
    sd_setBoundariesAct = toolsMenu->addAction(tr("Set boundaries"), this, &MainWindow::SDSetBoundariesSig);

    // Signed Distance: Set pass through
    sd_setPassThroughAct = toolsMenu->addAction(tr("Set pass through"), this, &MainWindow::SDsetPassThrough);

    // Signed Distance: Read off distances to objects
    sd_readDataAct = toolsMenu->addAction(tr("Get distance to objects"), this, &MainWindow::SDReadDataSig);

    // Signed Distance: Read off full distances to objects
    sd_readFullDataAct = toolsMenu->addAction(tr("Get full distance to objects"), this, &MainWindow::SDReadFullDataSig);

    // Signed Distance: Read off distances to objects
    ripleyKAct = toolsMenu->addAction(tr("Compute Ripleys K-function"), this, &MainWindow::ripleysKFunctionSig);

    toolsMenu->addSeparator();

    customFunctionAct = toolsMenu->addAction(tr("Run custom function"), this, &MainWindow::runCustomFunctionSig);

    toolsMenu->addSeparator();

    finalizeNeurofilamentAct = toolsMenu->addAction(tr("Finalize Neurofilament annotation"), this, &MainWindow::finalizeNeurofilamentSig);

    topToolBar->addSeparator();

    // Image icon
    QIcon imageIcon = QIcon::fromTheme("document-new", QIcon(":/images/image.png"));
    imageAct = new QAction(imageIcon, tr("&Show image"), this);
    imageAct->setStatusTip(tr("Show image as background"));
    imageAct->setCheckable(true);
    imageAct->setChecked(true);
    connect(imageAct, &QAction::triggered, this, &MainWindow::showImage);
    topToolBar->addAction(imageAct);

    // Phi icon
    QIcon phiIcon = QIcon::fromTheme("document-new", QIcon(":/images/phi.png"));
    phiAct = new QAction(phiIcon, tr("&Show phasefield"), this);
    phiAct->setStatusTip(tr("Show phasefield as background"));
    phiAct->setCheckable(true);
    connect(phiAct, &QAction::triggered, this, &MainWindow::showPhi);
    topToolBar->addAction(phiAct);

    // signed distance function icon
    QIcon signedIcon = QIcon::fromTheme("document-new", QIcon(":/images/signed.png"));
    signedAct = new QAction(signedIcon, tr("&Show signed distance function"), this);
    signedAct->setStatusTip(tr("Show signed distance function as background"));
    signedAct->setCheckable(true);
    connect(signedAct, &QAction::triggered, this, &MainWindow::showSigned);
    topToolBar->addAction(signedAct);

    topToolBar->addSeparator();

    // boundary icon
    QIcon boundaryIcon = QIcon::fromTheme("document-new", QIcon(":/images/boundary.png"));
    boundaryAct = new QAction(boundaryIcon, tr("&Show level-set boundary"), this);
    boundaryAct->setStatusTip(tr("Show or hide level-set boundary"));
    boundaryAct->setCheckable(true);
    connect(boundaryAct, &QAction::triggered, this, &MainWindow::toggleBoundary);
    topToolBar->addAction(boundaryAct);

    // segmented icon
    QIcon segmentedIcon = QIcon::fromTheme("document-new", QIcon(":/images/segmented.png"));
    segmentedAct = new QAction(segmentedIcon, tr("&Show segmented objects"), this);
    segmentedAct->setStatusTip(tr("Show or hide segmented objects"));
    segmentedAct->setCheckable(true);
    connect(segmentedAct, &QAction::triggered, this, &MainWindow::toggleSegmented);
    topToolBar->addAction(segmentedAct);

    // 3d icon
    QIcon renderIn3DIcon = QIcon::fromTheme("document-new", QIcon(":/images/3d.png"));
    QAction *renderIn3DAct = new QAction(renderIn3DIcon, tr("&Show level-set in 3D"), this);
    renderIn3DAct->setStatusTip(tr("Show level-set in 3D"));
    connect(renderIn3DAct, &QAction::triggered, this, &MainWindow::renderIn3D);
    topToolBar->addAction(renderIn3DAct);

    topToolBar->addSeparator();

    QAction *openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::open);
    fileMenu->addAction(openAct);

    const QIcon saveAsIcon = QIcon::fromTheme("document-save-as");
    QAction *saveAsAct = fileMenu->addAction(saveAsIcon, tr("Save &As..."), this, &MainWindow::saveAs);
    saveAsAct->setShortcuts(QKeySequence::SaveAs);
    saveAsAct->setStatusTip(tr("Save the document under a new name"));

    fileMenu->addSeparator();

    const QIcon exitIcon = QIcon::fromTheme("application-exit");
    QAction *exitAct = fileMenu->addAction(exitIcon, tr("E&xit"), this, &QWidget::close);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit the application"));

    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    QAction *aboutAct = helpMenu->addAction(tr("&About"), this, &MainWindow::about);
    aboutAct->setStatusTip(tr("Show the application's About box"));

    QAction *aboutQtAct = helpMenu->addAction(tr("About &Qt"), qApp, &QApplication::aboutQt);
    aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));

    nextImageAct = viewMenu->addAction(tr("Next image"), this, &MainWindow::nextImage);
    nextImageAct->setShortcut(QKeySequence::ZoomIn);

    prevImageAct = viewMenu->addAction(tr("Previous image"), this, &MainWindow::prevImage);
    prevImageAct->setShortcut(QKeySequence::ZoomOut);

    Debug::Info("MainWindow::createActions: Leaving");
}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings()
{
    QSettings settings(QCoreApplication::organizationName(), QCoreApplication::applicationName());
    const QByteArray geometry = settings.value("geometry", QByteArray()).toByteArray();
    if (geometry.isEmpty()) {
        const QRect availableGeometry = QApplication::desktop()->availableGeometry(this);
        resize(availableGeometry.width() / 3, availableGeometry.height() / 2);
        move((availableGeometry.width() - width()) / 2,
             (availableGeometry.height() - height()) / 2);
    } else {
        restoreGeometry(geometry);
    }
}

void MainWindow::writeSettings()
{
    QSettings settings(QCoreApplication::organizationName(), QCoreApplication::applicationName());
    settings.setValue("geometry", saveGeometry());
}

bool MainWindow::maybeSave()
{
    if (!textEdit->document()->isModified())
        return true;
    const QMessageBox::StandardButton ret
        = QMessageBox::warning(this, tr("Application"),
                               tr("The document has been modified.\n"
                                  "Do you want to save your changes?"),
                               QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    switch (ret) {
    case QMessageBox::Save:
        return save();
    case QMessageBox::Cancel:
        return false;
    default:
        break;
    }
    return true;
}

void MainWindow::runCV()
{
    setMessage("Running Chan Vese segmentation");

    int numIterations = iterLineEdit->text().toInt();
    float mu = muLineEdit->text().toFloat();
    float nu = nuLineEdit->text().toFloat();

    float c1 = c1LineEdit->text().toFloat();
    float c2 = c2LineEdit->text().toFloat();

    bool constc1 = c1lockAct->isChecked();
    bool constc2 = c2lockAct->isChecked();

    emit runChanVeseSig(numIterations, mu, nu, c1, c2, constc1, constc2);

    setMessage("Done");
}

void MainWindow::setPos()
{
    if (selectClickedOnce == false)
    {
        forceSelectionIntoActiveRegion();

        emit setPosSig(selectionBeginX, selectionEndX,
                       selectionBeginY, selectionEndY,
                       selectionBeginZ, selectionEndZ);
    }
}

void MainWindow::setNeg()
{
    if (selectClickedOnce == false)
    {
        //forceSelectionIntoActiveRegion();

        emit setNegSig(selectionBeginX, selectionEndX,
                       selectionBeginY, selectionEndY,
                       selectionBeginZ, selectionEndZ);
    }
}

void MainWindow::runBoundedSignedDistanceMap()
{
    emit SDRunBoundedSignedDistanceMapSig();
}

void MainWindow::setDisallowed()
{
    if (selectClickedOnce == false)
    {
        forceSelectionIntoActiveRegion();

        emit setDisallowedSig(selectionBeginX, selectionEndX,
                              selectionBeginY, selectionEndY,
                              selectionBeginZ, selectionEndZ);
    }
}

void MainWindow::swapPhiSign()
{
    emit swapPhiSignSig();
}

void MainWindow::reinitPhi()
{
    emit reinitPhiSig();
}

void MainWindow::redrawBackground(bool activeRegionOnly)
{
    if (imageShown) {
        showImage(activeRegionOnly);
    } else if (phiShown) {
        showPhi();
    } else if (signedShown) {
        showSigned();
    }
}

void MainWindow::showImage(bool activeRegionOnly)
{
    imageShown = true;
    phiShown = false;
    signedShown = false;
    updateCheckedActions();
    QPainter p(&pixmap);
    p.setOpacity(1.0);
    p.setCompositionMode(QPainter::CompositionMode_SourceOver);
    p.drawImage(QPoint(0, 0), image);
    if (segmentedAct->isChecked() and not activeRegionOnly) p.drawImage(QPoint(0, 0), segmented);
    p.setOpacity(0.75);
    if (boundaryAct->isChecked() and not activeRegionOnly) p.drawImage(QPoint(0, 0), boundary);
    if (currentZ >= activeRegionBeginZ and currentZ <= activeRegionEndZ)
    {
        p.drawImage(QPoint(activeRegionBeginX - 1, activeRegionBeginY - 1), active);
    }
    if (selectionShown and currentZ >= selectionBeginZ and currentZ <= selectionEndZ and not activeRegionOnly) p.drawImage(QPoint(0, 0), selection);
    imageLabel->setPixmap(pixmap);
}

void MainWindow::showPhi()
{
    imageShown = false;
    phiShown = true;
    signedShown = false;
    updateCheckedActions();

    QPainter p(&pixmap);
    p.setOpacity(1.0);
    p.setCompositionMode(QPainter::CompositionMode_SourceOver);
    p.drawImage(QPoint(0, 0), phi);
    if (selectionShown and currentZ >= selectionBeginZ and currentZ <= selectionEndZ) p.drawImage(QPoint(0, 0), selection);
    imageLabel->setPixmap(pixmap);
}

void MainWindow::showSigned()
{
    imageShown = false;
    phiShown = false;
    signedShown = true;
    updateCheckedActions();

    QPainter p(&pixmap);
    p.setOpacity(1.0);
    p.setCompositionMode(QPainter::CompositionMode_SourceOver);
    p.drawImage(QPoint(0, 0), signedMap);
    if (selectionShown and currentZ >= selectionBeginZ and currentZ <= selectionEndZ) p.drawImage(QPoint(0, 0), selection);
    imageLabel->setPixmap(pixmap);
}

void MainWindow::renderIn3D()
{
    Debug::Info("MainWindow::renderIn3D: Entering");

    emit renderIn3DSig();

    Debug::Info("MainWindow::renderIn3D: Leaving");
}

void MainWindow::toggleBoundary()
{
    redrawBackground();
}

void MainWindow::toggleSegmented()
{
    redrawBackground();
}

void MainWindow::loadFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName), file.errorString()));
        return;
    }

    QTextStream in(&file);
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    textEdit->setPlainText(in.readAll());
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif

    setCurrentFile(fileName);
}

bool MainWindow::saveFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName),
                                  file.errorString()));
        return false;
    }

    QTextStream out(&file);
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    out << textEdit->toPlainText();
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif

    setCurrentFile(fileName);
    return true;
}

void MainWindow::setCurrentFile(const QString &fileName)
{
    curFile = fileName;
    textEdit->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty())
        shownName = "untitled.txt";
    setWindowFilePath(shownName);
}

QString MainWindow::strippedName(const QString &fullFileName)
{
    return QFileInfo(fullFileName).fileName();
}
#ifndef QT_NO_SESSIONMANAGER
void MainWindow::commitData(QSessionManager &manager)
{
    if (manager.allowsInteraction()) {
        if (!maybeSave())
            manager.cancel();
    } else {
        // Non-interactive: save without asking
        if (textEdit->document()->isModified())
            save();
    }
}
#endif

void MainWindow::setFillImage(const int width,
                              const int height,
                              const int z,
                              const float *arrayImage,
                              const float *activeRegionImage,
                              bool activeRegionOnly)
{
    Debug::Info("MainWindow::setFillImage: Entering");

    currentZ = z;

    QImage newQImage;
    if (activeRegionOnly)
    {
        newQImage = pixmap.toImage();
    }
    else
    {
        newQImage = ArrayImageToQImage(width, height, arrayImage);
    }

    if (activeRegionImage != NULL and currentZ >= activeRegionBeginZ and currentZ <= activeRegionEndZ)
    {

        int activeRegionWidth = activeRegionEndX - activeRegionBeginX + 1;

        for (int j = activeRegionBeginY; j < activeRegionEndY + 1; ++j)
        {
            for (int i = activeRegionBeginX; i < activeRegionEndX + 1; ++i)
            {
                int ii = i - activeRegionBeginX;
                int jj = j - activeRegionBeginY;

                int pixelVal = (int) activeRegionImage[jj*activeRegionWidth + ii];

                newQImage.setPixelColor(i, j, QColor(pixelVal, pixelVal, pixelVal, 255));
            }
        }
    }

    image = newQImage;

    QPixmap newPixmap(image.size());
    pixmap = newPixmap;

    scrollArea->setVisible(true);

    redrawBackground(activeRegionOnly);

    Debug::Info("MainWindow::setFillImage: Leaving");
}

void MainWindow::updateSegmented(const bool *segmentedFlag)
{
    Debug::Info("MainWindow::updateSegmented: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    newQImage.fill(QColor(0,0,0,0));

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

    segmented = newQImage;

    redrawBackground();

    Debug::Info("MainWindow::updateSegmented: Leaving");
}

void MainWindow::updatePhi(const float *arrayImage, float newPhiZero)
{
    Debug::Info("MainWindow::updatePhi: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);
    // (1) init phi to be black
    newQImage.fill(QColor(0,0,0,255));

    if (arrayImage != NULL) {

        phiZero = newPhiZero;

        // (2) calculate width, height and depth of sub-image
        int activeRegionWidth = activeRegionEndX - activeRegionBeginX + 1;
        // (3) copy in relevant parts of the incoming array
        for (int j = activeRegionBeginY; j < activeRegionEndY + 1; ++j)
        {
            for (int i = activeRegionBeginX; i < activeRegionEndX + 1; ++i)
            {
                int ii = i - activeRegionBeginX;
                int jj = j - activeRegionBeginY;

                int pixelVal = (int) arrayImage[jj*activeRegionWidth + ii];

                newQImage.setPixelColor(i, j, QColor(pixelVal, pixelVal, pixelVal, 255));
            }
        }
    }

    phi = newQImage;
    signedMapCalculated = false;

    generatePhiBoundary();
    redrawBackground();

    Debug::Info("MainWindow::updatePhi: Leaving");
}

void MainWindow::updateSignedMap(const float *arrayImage, int sdOriginX,
                                 int sdOriginY, int sdWidth, int sdHeight)
{
    Debug::Info("MainWindow::updateSignedMap: Entering");


    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    // (1) init phi to be black
    newQImage.fill(QColor(0,0,0,255));

    if (arrayImage != NULL)
    {
        // (3) copy in relevant parts of the incoming array
        for (int j = sdOriginY; j < sdOriginY + sdHeight; ++j)
        {
            for (int i = sdOriginX; i < sdOriginX + sdWidth; ++i)
            {
                int ii = i - sdOriginX;
                int jj = j - sdOriginY;

                float pixelVal = arrayImage[jj*sdWidth + ii];

                QColor color(0, (int) pixelVal, (int) pixelVal, 255);

                if (color.isValid())
                {
                    newQImage.setPixelColor(i, j, color);
                }
                else
                {
                    newQImage.setPixelColor(i, j, QColor(255, 0, 0, 255));
                    Debug::Error("Color is not valid at (" + STR(i) + ", "
                                 + STR(j) + ") with value " + STR(pixelVal)
                                 + " coming from gidx " + STR(jj*sdWidth + ii)
                                 + " of " + STR(sdWidth*sdHeight));
                }
            }
        }
    }

    signedMap = newQImage;
    signedMapCalculated = true;

    redrawBackground();

    Debug::Info("MainWindow::updateSignedMap: Leaving");
}

void MainWindow::InitViewArea(const int width, const int height, const int depth)
{
    Debug::Info("MainWindow::InitViewArea: Entering");

    scaleFactor = 1.0;

    imageWidth = width;
    imageHeight = height;
    imageDepth = depth;

    activeRegionBeginX = 0;
    activeRegionEndX = width - 1;
    activeRegionBeginY = 0;
    activeRegionEndY = height - 1;
    activeRegionBeginZ = 0;
    activeRegionEndZ = depth - 1;

    Debug::Info("MainWindow::InitViewArea: Leaving");
}

void MainWindow::createTree(QTreeWidget *treeWidget)
{
    cellDockWidget = new QDockWidget(tr("Shapes"));
    cellDockWidget->setObjectName("cellDockWidget");
    cellDockWidget->setAllowedAreas(Qt::LeftDockWidgetArea
                                   | Qt::RightDockWidgetArea);
    cellDockWidget->setWidget(treeWidget);
    addDockWidget(Qt::RightDockWidgetArea, cellDockWidget);

}

void MainWindow::forceActiveRegionChange(int beginX, int endX,
                                         int beginY, int endY,
                                         int beginZ, int endZ)
{
    activeRegionPresent = true;

    activeRegionBeginX = beginX;
    activeRegionEndX = endX;
    activeRegionBeginY = beginY;
    activeRegionEndY = endY;
    activeRegionBeginZ = beginZ;
    activeRegionEndZ = endZ;

    drawActiveRegionBorder();
    redrawBackground();
}

void MainWindow::generatePhiBoundary()
{
    Debug::Info("MainWindow::generatePhiBoundary: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);
    // (1) init phi to be black
    newQImage.fill(QColor(0,0,0,0));

    // (2) calculate width, height and depth of sub-image
    // (3) copy in relevant parts of the incoming array
    for (int j = activeRegionBeginY; j < activeRegionEndY + 1; ++j)
    {
        for (int i = activeRegionBeginX; i < activeRegionEndX + 1; ++i)
        {
            float ori = phi.pixelColor(i,j).red() - phiZero;
            float a = ori, b = ori, c = ori, d = ori;
            if (i < activeRegionEndX) a = phi.pixelColor(i+1,j).red() - phiZero;
            if (i > activeRegionBeginX) b = phi.pixelColor(i-1,j).red() - phiZero;
            if (j < activeRegionEndY) c = phi.pixelColor(i,j+1).red() - phiZero;
            if (j > activeRegionBeginY) d = phi.pixelColor(i,j-1).red() - phiZero;

            if (ori < 0 and (ori * a < 0 or ori * b < 0 or ori * c < 0 or ori * d < 0))
            {
                newQImage.setPixelColor(i, j, QColor(255, 0, 0, 255));
            }
        }
    }


    boundary = newQImage;

    Debug::Info("MainWindow::generatePhiBoundary: Leaving");
}

QImage MainWindow::ArrayImageToQImage(const int width, const int height, const float *arrayImage)
{
    Debug::Info("MainWindow::ArrayImageToQImage: Entering");

    QImage newQImage(width, height, QImage::Format_ARGB32);

    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            int pixelVal = (int) arrayImage[j*width + i];
            newQImage.setPixelColor(i, j, QColor(pixelVal, pixelVal, pixelVal, 255));
        }
    }

    Debug::Info("MainWindow::ArrayImageToQImage: Leaving");
    return newQImage;
}

void MainWindow::zoomIn()
{
    scaleImage(1.25);
}

void MainWindow::zoomOut()
{
    scaleImage(0.8);
}

void MainWindow::nextImage()
{
    emit nextImageSig();
}

void MainWindow::prevImage()
{
    emit prevImageSig();
}

void MainWindow::SDsetPassThrough()
{
    float sigma = blightLineEdit->text().toFloat();

    emit SDSetPassThroughSig(sigma);
}

void MainWindow::setMessage(char const message[])
{
    statusBar()->showMessage(tr(message));
}

void MainWindow::scrollAreaMouseClick(float x, float y)
{
    float dx = (float) scrollArea->horizontalScrollBar()->value();
    float dy = (float) scrollArea->verticalScrollBar()->value();

    int imageX = int((x+dx)/imageLabel->width() * imageWidth);
    int imageY = int((y+dy)/imageLabel->height() * imageHeight);
    if (imageX >= 0 and imageY >= 0 and imageX < imageWidth and imageY < imageHeight)
    {
        if (selectAct->isChecked()) {
            if (selectClickedOnce) {
                selectionEndX = imageX;
                selectionEndY = imageY;
                selectionEndZ = currentZ;

                if (selectionEndX < selectionBeginX) std::swap(selectionBeginX, selectionEndX);
                if (selectionEndY < selectionBeginY) std::swap(selectionBeginY, selectionEndY);
                if (selectionEndZ < selectionBeginZ) std::swap(selectionBeginZ, selectionEndZ);

                selectClickedOnce = false;

                drawSelectionBorder();
                selectionShown = true;
            } else {
                selectionBeginX = imageX;
                selectionBeginY = imageY;
                selectionBeginZ = currentZ;

                selectClickedOnce = true;
                selectionShown = false;
            }
            redrawBackground();
        }
        else if (vesicleCenterAct->isChecked())
        {
            Debug::Info("Vesicle: " + STR(imageX) + ", " + STR(imageY));
            emit markVesicleLocationSig(imageX, imageY);
        }
        else if (neurofilamentCenterAct->isChecked())
        {
            if (    imageX >= activeRegionBeginX
                and imageX <= activeRegionEndX
                and imageY >= activeRegionBeginY
                and imageY <= activeRegionEndY)
            {
                emit markNeurofilamentSig(imageX, imageY);
            }
        }
        else if (blightSelectAct->isChecked())
        {
            if (    imageX >= activeRegionBeginX
                and imageX <= activeRegionEndX
                and imageY >= activeRegionBeginY
                and imageY <= activeRegionEndY)
            {
                float sigma = blightLineEdit->text().toFloat();
                emit blightSig(imageX, imageY, sigma);
            }
        }
        else if (restoreSelectAct->isChecked())
        {
            if (    imageX >= activeRegionBeginX
                and imageX <= activeRegionEndX
                and imageY >= activeRegionBeginY
                and imageY <= activeRegionEndY)
            {
                float sigma = blightLineEdit->text().toFloat();
                emit restoreSig(imageX, imageY, sigma);
            }
        }
    }

}

void MainWindow::adjustScrollBar(QScrollBar *scrollBar, float factor)
{
    scrollBar->setValue(int(factor * scrollBar->value()
                            + ((factor - 1) * scrollBar->pageStep()/2)));
}

void MainWindow::scaleImage(float factor)
{
    Q_ASSERT(imageLabel->pixmap());
    scaleFactor *= factor;
    imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);

    scrollArea->setScale(scaleFactor);
}

void MainWindow::drawSelectionBorder()
{
    Debug::Info("MainWindow::drawSelectionBorder: Entering");

    QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    for (int j = 0; j < imageHeight; ++j)
    {
        for (int i = 0; i < imageWidth; ++i)
        {
            if (((i >= selectionBeginX - 1 and (j == selectionBeginY - 1 or j == selectionEndY + 1)) and
                (i <= selectionEndX   + 1 and (j == selectionBeginY - 1 or j == selectionEndY + 1))) or
                ((j >= selectionBeginY - 1 and (i == selectionBeginX - 1 or i == selectionEndX + 1)) and
                (j <= selectionEndY   + 1 and (i == selectionBeginX - 1 or i == selectionEndX + 1))))
            {
                newQImage.setPixelColor(i, j, QColor(0, 255, 0, 255));
            }
            else
            {
                newQImage.setPixelColor(i, j, QColor(0, 0, 0, 0));
            }
        }
    }
    selection = newQImage;

    Debug::Info("MainWindow::drawSelectionBorder: Leaving");
}

void MainWindow::drawActiveRegionBorder()
{
    Debug::Info("MainWindow::drawActiveRegionBorder: Entering");

    //QImage newQImage(imageWidth, imageHeight, QImage::Format_ARGB32);

    int activeBorderRegionWidth = activeRegionEndX - activeRegionBeginX + 3;
    int activeBorderRegionHeight = activeRegionEndY - activeRegionBeginY + 3;

    QImage newQImage(activeBorderRegionWidth, activeBorderRegionHeight, QImage::Format_ARGB32);


    /*for (int j = 0; j < imageHeight; ++j)
    {
        for (int i = 0; i < imageWidth; ++i)
        {
            if (((i >= activeRegionBeginX - 1 and (j == activeRegionBeginY - 1 or j == activeRegionEndY + 1)) and
                (i <= activeRegionEndX   + 1 and (j == activeRegionBeginY - 1 or j == activeRegionEndY + 1))) or
                ((j >= activeRegionBeginY - 1 and (i == activeRegionBeginX - 1 or i == activeRegionEndX + 1)) and
                (j <= activeRegionEndY   + 1 and (i == activeRegionBeginX - 1 or i == activeRegionEndX + 1))))
            {
                newQImage.setPixelColor(i, j, QColor(0, 0, 255, 255));
            }
            else
            {
                newQImage.setPixelColor(i, j, QColor(0, 0, 0, 0));
            }
        }
    }*/

    newQImage.fill(QColor(0,0,0,0));

    for (int j = 0; j < activeBorderRegionHeight; ++j)
    {
        for (int i = 0; i < activeBorderRegionWidth; ++i)
        {
            if (j == 0 or j == activeBorderRegionHeight - 1 or
                i == 0 or i == activeBorderRegionWidth - 1)
            {
                newQImage.setPixelColor(i, j, QColor(0, 0, 255, 255));
            }
        }
    }
    active = newQImage;

    Debug::Info("MainWindow::drawActiveRegionBorder: Leaving");
}

void MainWindow::setActiveRegion()
{
    Debug::Info("MainWindow::setActiveRegion: Entering");

    setMessage("Setting new active region");

    if (selectClickedOnce == false) {

        bool commitToMakingNewActiveRegion = not activeRegionPresent;

        if (activeRegionPresent)
        {
            QMessageBox::StandardButton reply;

            reply = QMessageBox::question(this,
                                          "Test", "Active region present! \n"
                                          "Are you sure you want to make a "
                                          "new one?",
                                          QMessageBox::Yes|QMessageBox::No);

            commitToMakingNewActiveRegion = (reply == QMessageBox::Yes);
        }

        if (commitToMakingNewActiveRegion)
        {
            activeRegionPresent = true;

            activeRegionBeginX = selectionBeginX;
            activeRegionEndX = selectionEndX;
            activeRegionBeginY = selectionBeginY;
            activeRegionEndY = selectionEndY;
            activeRegionBeginZ = selectionBeginZ;
            activeRegionEndZ = selectionEndZ;

            selectionShown = false;

            drawActiveRegionBorder();
            redrawBackground();

            float sigma = blightLineEdit->text().toFloat();

            emit newActiveRegionSig(activeRegionBeginX, activeRegionEndX,
                                    activeRegionBeginY, activeRegionEndY,
                                    activeRegionBeginZ, activeRegionEndZ,
                                    sigma);
        }
    }

    setMessage("Done");

    Debug::Info("MainWindow::setActiveRegion: Leaving");
}

void MainWindow::storeObject()
{
    Debug::Info("MainWindow::storeObject: Entering");

    emit storeObjectSig();

    Debug::Info("MainWindow::storeObject: Leaving");
}

void MainWindow::loadObject()
{
    Debug::Info("MainWindow::loadObject: Entering");

    emit loadObjectSig();

    Debug::Info("MainWindow::loadObject: Leaving");
}

void MainWindow::deleteObject()
{
    Debug::Info("MainWindow::deleteObject: Entering");

    QMessageBox::StandardButton reply;

    reply = QMessageBox::question(this,
                                  "Test", "Are you sure you want to delete "
                                  "this object",
                                  QMessageBox::Yes|QMessageBox::No);

    if (reply == QMessageBox::Yes)
    {
        emit deleteObjectSig();
    }

    Debug::Info("MainWindow::deleteObject: Leaving");
}

void MainWindow::updateCheckedActions()
{
    imageAct->setChecked(imageShown);
    phiAct->setChecked(phiShown);
    signedAct->setChecked(signedShown);
}

void MainWindow::forceSelectionIntoActiveRegion()
{
    selectionBeginX = std::max(std::min(selectionBeginX, activeRegionEndX),
                               activeRegionBeginX);
    selectionEndX   = std::max(std::min(selectionEndX, activeRegionEndX),
                               activeRegionBeginX);
    selectionBeginY = std::max(std::min(selectionBeginY, activeRegionEndY),
                               activeRegionBeginY);
    selectionEndY   = std::max(std::min(selectionEndY, activeRegionEndY),
                               activeRegionBeginY);
    selectionBeginZ = std::max(std::min(selectionBeginZ, activeRegionEndZ),
                               activeRegionBeginZ);
    selectionEndZ   = std::max(std::min(selectionEndZ, activeRegionEndZ),
                               activeRegionBeginZ);
}