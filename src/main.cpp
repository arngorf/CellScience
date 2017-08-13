#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "ellipsoidHelpers.hpp"
#include "EllipsoidLeastSquareFit.hpp"
#include "EllipsoidMinimizer.hpp"
#include "Minimizer.hpp"
#include "MinimizerN.hpp"
#include "CellPlot.hpp"

#include "ChanVese3D.hpp"
#include "TiffReader.hpp"
#include "types.hpp"
#include <stdlib.h>

#include "NumericToolControl.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QCommandLineOption>

#include "MainWindow.h"
#include "Debug.hpp"

/* mu parameters
 *
 * Vesicle, Endosome: 1000 (lock 130/160)
 * Membrane, Synapse: 2000
 * Filament: 1500 (lock 160/178)
 * FIB-SEM image size (450x450)
 */

int main(int argc, char *argv[])
{
    static Debug debug;

    Debug::Info("Main: Entering");

    Q_INIT_RESOURCE(application);

    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("QtProject");
    QCoreApplication::setApplicationName("Application Example");
    QCoreApplication::setApplicationVersion(QT_VERSION_STR);

    NumericToolControl ntCtrl("/"
                              "home/"
                              "dith/"
                              "MahdiehVesicle/"
                              "3dvesicleSnakes/"
                              "LaussaneDataNAnnotations/"
                              "volumedata_registered.tif");
    MainWindow *mainWin = new MainWindow();
    ntCtrl.run(mainWin, argc, argv);
    mainWin->show();

    return app.exec();

    Debug::Info("Main: Leaving");
}