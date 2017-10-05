#include "VesicleSegmentationResultManager.hpp"

#include "UtilityFunctions.hpp"
#include "Debug.hpp"
#include "ellipsoidHelpers.hpp"

#include <random>

#include <QDirIterator>

bool vesicleResultCompletionCmp (VesicleSegmentationResult *a,VesicleSegmentationResult *b)
{
    return a->GetCompletionValue() < b->GetCompletionValue();
}

VesicleSegmentationResultManager::VesicleSegmentationResultManager(
                                                char const path[],
                                                int width,
                                                int height
                                                )
                                                : imagePath(path)
                                                , tiffImage(NULL)
                                                , vesicleResultArray(0)
                                                , width(width)
                                                , height(height)
                                                , cxImg(25)
                                                , cyImg(25)
                                                , czImg(25)
{
    Debug::Info("VesicleSegmentationResultManager::VesicleSegmentationResultManager: Entering");

    loadSegmentationResultsFromFile();

    Debug::Info("VesicleSegmentationResultManager::VesicleSegmentationResultManager: Leaving");
}

VesicleSegmentationResultManager::~VesicleSegmentationResultManager()
{
    Debug::Info("VesicleSegmentationResultManager::~VesicleSegmentationResultManager: Entering");

    saveSegmentationResultsToFile();

    Debug::Info("VesicleSegmentationResultManager::~VesicleSegmentationResultManager: Leaving");
}

void VesicleSegmentationResultManager::NextSlice()
{
    //vesicleSegmentationWindow->setBigMessage("Loading vesicle slice"); vesicleSegmentationWindow->showBigMessage();

    // If no vesicles vesicleSegmentationWindow->setBigMessage("No Vesicles found"); vesicleSegmentationWindow->showBigMessage(); return;

    std::random_device rd;
    std::mt19937 gen(rd());

    int cellVesicleCount = vesicleResultArray.size();

    std::uniform_int_distribution<> randI(0, cellVesicleCount*9 - 1);
    int randInt = randI(gen);

    curVesicleIndex = randInt / 9;
    curDirectionIndex = randInt % 9;

    vesicleResultArray[curVesicleIndex]->getCenter(x_center, y_center, z_center);

    TiffReader reader;

    xBegin = x_center - cxImg;
    yBegin = y_center - cyImg;
    zBegin = z_center - czImg;

    tiffImage = reader.readImageToTiffClass(imagePath.c_str(),
                                            xBegin, xBegin + cxImg*2+1,
                                            yBegin, yBegin + cyImg*2+1,
                                            zBegin, zBegin + czImg*2+1);
}

float VesicleSegmentationResultManager::GetCompletionValue()
{
    return ((float) completedSlices)
           / ((float) (vesicleResultArray.size() * numRotations));
}

float VesicleSegmentationResultManager::GetMinCompletionValue()
{
    if (vesicleResultArray.size() <= 0)
    {
        return 0;
    }

    float minValue = vesicleResultArray[0]->GetCompletionValue();

    for (int i = 1; i < vesicleResultArray.size(); ++i)
    {
        minValue = std::min(vesicleResultArray[i]->GetCompletionValue(),
                            minValue);
    }

    return minValue;
}

float *VesicleSegmentationResultManager::GetImageSlice()
{
    generateRotation(curDirectionIndex);

    float* imageSlice = (float *) malloc(width*height * sizeof(float));

    int imgWidth = tiffImage->getWidth();
    int imgHeight = tiffImage->getHeight();
    int imgDepth = tiffImage->getDepth();

    float *rotatedImage = (float *) malloc(imgWidth*imgHeight*imgDepth * sizeof(float));

    for (int k = 0; k < imgDepth; ++k)
    {
        for (int j = 0; j < imgHeight; ++j)
        {
            for (int i = 0; i < imgWidth; ++i)
            {
                float x = (float) i;
                float y = (float) j;
                float z = (float) k;

                RotatePointAroundPoint(x, y, z, cxImg, cyImg, czImg, currentRInverse);

                int ii = (int) x;
                int jj = (int) y;
                int kk = (int) z;

                int idx = gIndex(i,j,k,imgHeight,imgWidth);

                if (ii >= 0 and jj >= 0 and kk >= 0 and
                    ii+1 < imgWidth and jj+1 < imgHeight and kk+1 < imgDepth)
                {
                    float dx = x - ii;
                    float dy = y - jj;
                    float dz = z - kk;

                    float v000 = tiffImage->getValue(ii  , jj  , kk  );
                    float v100 = tiffImage->getValue(ii+1, jj  , kk  );
                    float v010 = tiffImage->getValue(ii  , jj+1, kk  );
                    float v110 = tiffImage->getValue(ii+1, jj+1, kk  );
                    float v001 = tiffImage->getValue(ii  , jj  , kk+1);
                    float v101 = tiffImage->getValue(ii+1, jj  , kk+1);
                    float v011 = tiffImage->getValue(ii  , jj+1, kk+1);
                    float v111 = tiffImage->getValue(ii+1, jj+1, kk+1);

                    float value = (1-dx)*(1-dy)*(1-dz)*v000
                                + (  dx)*(1-dy)*(1-dz)*v100
                                + (1-dx)*(  dy)*(1-dz)*v010
                                + (  dx)*(  dy)*(1-dz)*v110
                                + (1-dx)*(1-dy)*(  dz)*v001
                                + (  dx)*(1-dy)*(  dz)*v101
                                + (1-dx)*(  dy)*(  dz)*v011
                                + (  dx)*(  dy)*(  dz)*v111;

                    rotatedImage[idx] = value;
                }
                else
                {
                    rotatedImage[idx] = 0;
                }
            }
        }
    }

    int di = (int) ((imgWidth - width) / 2.0);
    int dj = (int) ((imgHeight - height) / 2.0);

    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            int idx = gIndex(i + di, j + dj, czImg, imgHeight, imgWidth);
            imageSlice[gIndex(i, j, height)] = rotatedImage[idx];
        }
    }

    free(rotatedImage);

    return imageSlice;
}

void VesicleSegmentationResultManager::SetSegmentationResult(std::vector<float> x, std::vector<float> y)
{
    int n = x.size();

    int imgWidth = tiffImage->getWidth();
    int imgHeight = tiffImage->getHeight();
    int imgDepth = tiffImage->getDepth();

    int di = (int) ((imgWidth - width) / 2.0);
    int dj = (int) ((imgHeight - height) / 2.0);

    for (int i = 0; i < n; ++i)
    {
        float rx = x[i] + (float) di;
        float ry = y[i] + (float) dj;
        float rz = czImg;

        RotatePointAroundPoint(rx, ry, rz, cxImg, cyImg, czImg, currentR);

        TranslatePointsToGlobal(rx, ry, rz);

        vesicleResultArray[curVesicleIndex]->addSegmentingPoint(rx, ry, rz);
    }
}

void VesicleSegmentationResultManager::AddVesicle(CellObject *object)
{
    bool found = false;

    for (int i = 0; i < vesicleResultArray.size(); ++i)
    {
        if (vesicleResultArray[i]->isEqual(object))
        {
            found = true;
            break;
        }
    }

    if (!found)
    {
        VesicleSegmentationResult *newVesRes = new VesicleSegmentationResult(object);
        vesicleResultArray.push_back(newVesRes);
    }
}

void VesicleSegmentationResultManager::loadSegmentationResultsFromFile()
{
    Debug::Info("VesicleSegmentationResultManager::loadSegmentationResultsFromFile: Entering");

    std::vector<QString> filenames;

    QDirIterator it(".", QDirIterator::Subdirectories);
    while (it.hasNext()) {
        QString filename = it.next();
        std::cout << "(1) " << filename.toStdString() << std::endl;
        if (QString("./V") == filename.left(3) and
            QString(".data") == filename.right(5))
        {
            std::cout << "(2) " << filename.toStdString() << std::endl;
            filenames.push_back(filename);
        }
    }

    for (int i = 0; i < filenames.size(); ++i)
    {

        QString filename = filenames[i];
        QFile file(filename);

        if(!file.open(QIODevice::ReadOnly))
        {
            Debug::Warning("VesicleSegmentationResultManager::"
                           "loadSegmentationResultsFromFile: Could not open "
                           "cell object file");
        }

        QDataStream in(&file);

        VesicleSegmentationResult *result = new VesicleSegmentationResult();
        in >> (*result);

        vesicleResultArray.push_back(result);

        vesicleResultArray[vesicleResultArray.size() - 1]->OutputSegmentation();

        file.close();
    }

    Debug::Info("VesicleSegmentationResultManager::loadSegmentationResultsFromFile: Leaving");
}

void VesicleSegmentationResultManager::saveSegmentationResultsToFile()
{
    Debug::Info("VesicleSegmentationResultManager::saveSegmentationResultsToFile: Entering");

    Debug::Info("VesicleSegmentationResultManager::"
                "saveSegmentationResultsToFile: Storing "
                + STR(vesicleResultArray.size())
                + " vesicle segmentation results");

    for (int i = 0; i < vesicleResultArray.size(); ++i)
    {
        std::string nameID = vesicleResultArray[i]->getNameID();
        QString filename = QString(nameID.c_str()) + QString(".data");
        QFile file(filename);

        if(!file.open(QIODevice::WriteOnly))
        {
            Debug::Warning("VesicleSegmentationResultManager::"
                           "saveSegmentationResultsToFile: Could not open "
                           "vesicle segmentation result file");
        }

        QDataStream out(&file);

        VesicleSegmentationResult *result = vesicleResultArray[i];
        out << (*result);

        file.flush();
        file.close();
    }

    Debug::Info("VesicleSegmentationResultManager::saveSegmentationResultsToFile: Leaving");
}

void VesicleSegmentationResultManager::generateRotation(int direction)
{
    currentR = Mat::Identity(3,3);

    Mat baseR = vesicleResultArray[curVesicleIndex]->getBaseRotationMatrix();

    switch (direction)
    {
        case 0:
        {
            break;
        }
        case 1:
        {
            currentR = GetAxisAngleRotationMatrix(0,1,0,PI/4.0);
            break;
        }
        case 2:
        {
            currentR = GetAxisAngleRotationMatrix(0,1,0,PI/2.0);
            break;
        }
        case 3:
        {
            currentR = GetAxisAngleRotationMatrix(0,1,0,3.0*PI/4.0);
            break;
        }
        case 4:
        {
            currentR = GetAxisAngleRotationMatrix(1,0,0,PI/4.0);
            break;
        }
        case 5:
        {
            currentR = GetAxisAngleRotationMatrix(1,0,0,3.0*PI/2.0);
            break;
        }
        case 6:
        {
            currentR = GetAxisAngleRotationMatrix(1,0,0,3.0*PI/4.0);
            break;
        }
        case 7:
        {
            Mat R1 = GetAxisAngleRotationMatrix(0,1,0,PI/2.0);
            Mat R2 = GetAxisAngleRotationMatrix(0,0,1,PI/4.0);
            currentR = R2 * R1;
            break;
        }
        case 8:
        {
            Mat R1 = GetAxisAngleRotationMatrix(0,1,0,PI/2.0);
            Mat R2 = GetAxisAngleRotationMatrix(0,0,1,3.0*PI/4.0);
            currentR = R2 * R1;
            break;
        }
        default:
        {
            Debug::Error("VesicleSegmentationResultManager::generateRotation:"
                         " Rotation direction not matching valid rotation");
        }
    }

    currentR = currentR * baseR;

    float vx = 0, vy = 0, vz = 1;
    float theta = Random::randU(Random::generator) * PI * 2;
    RotatePoint(vx, vy, vz, currentR);
    Mat finalR = GetAxisAngleRotationMatrix(0,1,0, theta);

    currentR = finalR * currentR;

    currentRInverse = currentR.inverse();
}

void VesicleSegmentationResultManager::RotatePoint(float &x, float &y, float &z, const Mat &R)
{
    Vector v(3);
    v(0) = x; v(1) = y; v(2) = z;
    Vector u = R*v;
    x = u(0); y = u(1); z = u(2);
}

void VesicleSegmentationResultManager::RotatePointAroundPoint(float &x, float &y, float &z, const float cx, const float cy, const float cz, const Mat &R)
{
    x -= cx;
    y -= cy;
    z -= cz;

    RotatePoint(x, y, z, R);

    x += cx;
    y += cy;
    z += cz;
}

void VesicleSegmentationResultManager::TranslatePointsToGlobal(float &x, float &y, float &z)
{
    x += xBegin;
    y += yBegin;
    z += zBegin;
}