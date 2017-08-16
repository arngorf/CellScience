#include "TiffReader.hpp"

#include <sstream>
#include <string>

#include "UtilityFunctions.hpp"

void DummyHandler(const char* module, const char* fmt, va_list ap) {
    // ignore errors and warnings (or handle them your own way)
    module = module; // I just want
    fmt = fmt; // the warnings to go away
    ap = ap; // until I come up with proper handling
}

void grabImageSlice(const std::string path,
                    const int width,
                    const int height,
                    const int beginZ,
                    const int endZ,
                    bool &grabbing,
                    TiffImage *&imagePointer,
                    std::mutex &tiffImageMutex)
{
    TiffReader tiffReader;

    TiffImage *newImagePointer = tiffReader.readImageToTiffClass(path.c_str(),
                                                                 0,
                                                                 width,
                                                                 0,
                                                                 height,
                                                                 beginZ,
                                                                 endZ);

    Debug::Warning("Grabber thread: getting lock");
    //std::lock_guard<std::mutex> guard(tiffImageMutex);
    tiffImageMutex.lock();

    if (imagePointer != NULL)
    {
        Debug::Warning("Deleting imagePointer");
        delete imagePointer;
    }

    imagePointer = newImagePointer;

    grabbing = false;
    tiffImageMutex.unlock();
    Debug::Warning("Grabber thread: done");
}

TiffImage::TiffImage() :_width(0), _height(0), _depth(0), _zOffset(0)
{
    Debug::Info("TiffImage::TiffImage (default): Entering");

    _image = NULL;
    _histogram = NULL;

    Debug::Info("TiffImage::TiffImage (default): Leaving");
}

TiffImage::TiffImage(const int width,
                     const int height,
                     const int depth,
                     const int zOffset) :
                     _width(width), _height(height), _depth(depth), _zOffset(zOffset)
{
    //Debug::Info("TiffImage::TiffImage: Entering");

    _image = NULL;
    _histogram = NULL;

    _image = (float*) malloc(width*height*depth * sizeof(float));

    if (_image == NULL)
    {
        Debug::Error("TiffImage::TiffImage: malloc of _image did not "
                     "succeed!");
    }

    _histogram = (int*) malloc(256 * sizeof(int));

    if (_histogram == NULL)
    {
        Debug::Error("TiffImage::TiffImage: malloc of _histogram did not "
                     "succeed!");
    }

    for (int i = 0; i < 256; ++i)
    {
        _histogram[i] = 0;
    }

    //Debug::Info("TiffImage::TiffImage: Leaving");
}

TiffImage::~TiffImage()
{
    if (_image != NULL) free(_image);
    if (_histogram != NULL) free(_histogram);
}

int TiffReader::read(char const path[], float *image,
                     int _xOffset, int _yOffset, int _zOffset,
                     int _width, int _height, int _depth)
{
    Debug::Info("TiffImage::read: Entering");

    // disable warnings
    TIFFSetWarningHandler(DummyHandler);


    TIFF *tif = TIFFOpen(path, "r");
    if (tif) {

        unsigned int width, height;
        // get the size of the tiff

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        if (_xOffset + _width > width)
        {
            Debug::Error("TiffImage::read: width of image exceeded");
        }

        if (_yOffset + _height > height)
        {
            Debug::Error("TiffImage::read: height of image exceeded");
        }

        size_t pageCounter = 0;

        do {

            uint32* raster;
            pageCounter++;
            //std::cout << pageCounter << std::endl;

            uint npixels = width*height; // get the total number of pixels

             // allocate temp memory (must use the tiff library malloc)
            raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));

            if (raster == NULL)
            {
                TIFFClose(tif);

                Debug::Error("TiffImage::read: Could not allocate memory for "
                             "raster of TIFF image");
            }

            // Check the tif read to the raster correctly
            if (!TIFFReadRGBAImage(tif, width, height, raster, 0)) {

                TIFFClose(tif);

                Debug::Error("TiffImage::read: Could not read raster of TIFF "
                             "image");

                 _TIFFfree(raster);
                free(image);
                return -1;
            }

            uint subImageSize = _height * _width;

            // iterate through all the pixels of the tif

            for (uint x = 0; x < width; x++) {
                for (uint y = 0; y < height; y++) {

                    int subImageZ = pageCounter - _zOffset;
                    int subImageY = (height - y) - _yOffset;
                    int subImageX = x - _xOffset;

                    if (subImageZ >= 0 and subImageZ < _depth  and
                        subImageY >= 0 and subImageY < _height and
                        subImageX >= 0 and subImageX < _width) {

                        // read the current pixel of the TIF
                        uint32& TiffPixel = raster[y*width + x];

                        // Take last 8 bits
                        unsigned low8bits = TiffPixel & 0xFF;

                        int index = subImageSize*subImageZ
                                  + subImageY*_width
                                  + subImageX;

                        image[index] = (float) low8bits;
                    }
                }
            }
            _TIFFfree(raster); // release temp memory

        } while (TIFFReadDirectory(tif)); // get the next tif

    if (_zOffset + _depth > pageCounter)
    {
        Debug::Error("TiffImage::read: depth of image exceeded");

        free(image);
    }

    TIFFClose(tif);

    _depth = pageCounter;
    _height = height;
    _width = width;


    } else {
        Debug::Error("TiffImage::read: could not read file "
                     + std::string(path) + ")");
    }

    Debug::Info("TiffImage::read: Leaving");

    return 0;
}

int TiffReader::readToImageHistogram(char const path[],
                                     std::vector<int> &IFreq)
{
    Debug::Info("TiffImage::readToImageHistogram: Entering");

    // disable warnings
    TIFFSetWarningHandler(DummyHandler);


    TIFF *tif = TIFFOpen(path, "r");
    if (tif) {

        unsigned int width, height;
        // get the size of the tiff

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        do {

            uint32* raster;

            uint npixels = width*height; // get the total number of pixels

            // allocate temp memory (must use the tiff library malloc)
            raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));

            if (raster == NULL)
            {
                TIFFClose(tif);

                Debug::Error("TiffImage::readToImageHistogram: Could not"
                             "allocate memory for raster of TIFF image");
            }

            // Check the tif read to the raster correctly
            if (!TIFFReadRGBAImage(tif, width, height, raster, 0)) {

                TIFFClose(tif);

                Debug::Error("TiffImage::readToImageHistogram: Could not read"
                             "raster of TIFF image");

                 _TIFFfree(raster);
                return -1;
            }

            // iterate through all the pixels of the tif

            for (uint x = 0; x < width; x++) {
                for (uint y = 0; y < height; y++) {

                    // read the current pixel of the TIF
                    uint32& TiffPixel = raster[y*width + x];

                    // Take last 8 bits
                    unsigned low8bits = TiffPixel & 0xFF;

                    IFreq[low8bits] += 1;

                }
            }
            _TIFFfree(raster); // release temp memory

        } while (TIFFReadDirectory(tif)); // get the next tif

    TIFFClose(tif);

    } else {
        Debug::Error("TiffImage::readToImageHistogram: could not read file "
                     + std::string(path) + ")");
    }

    Debug::Info("TiffImage::readToImageHistogram: Leaving");

    return 0;
}

TiffImage *TiffReader::readImageToTiffClass(char const path[],
                                           const int x_start, const int x_end,
                                           const int y_start, const int y_end,
                                           const int z_start, const int z_end)
{
    //Debug::Info("TiffImage::readImageToTiffClass: Entering");

    // disable warnings
    TIFFSetWarningHandler(DummyHandler);

    // Page count cannot be read directly. Its value is hacked here.
    size_t pageCounter = 0;
    unsigned int width, height;
        // get the size of the tiff

    TIFF *tif = TIFFOpen(path, "r");
    if (tif)
    {
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        do {
            pageCounter++;
        } while (TIFFReadDirectory(tif)); // get the next tif
    }
    TIFFClose(tif);

    const int windowWidth = x_end - x_start;
    const int windowHeight = y_end - y_start;
    const int windowDepth = z_end - z_start;

    TiffImage *tiffImage = new TiffImage(windowWidth,
                                         windowHeight,
                                         windowDepth,
                                         z_start);

    // Open and read image
    tif = TIFFOpen(path, "r");
    int z = -1;
    if (tif)
    {
        do {
            ++z;
            int window_Z = z - z_start;
            if (not (window_Z >= 0 and window_Z < windowDepth)) continue;
            uint32* raster;

            //std::cout << pageCounter << std::endl;

            uint npixels = width*height; // get the total number of pixels

            // allocate temp memory (must use the tiff library malloc)
            raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));

            if (raster == NULL)
            {
                TIFFClose(tif);

                Debug::Error("TiffImage::readImageToTiffClass: Could not"
                             "allocate memory for raster of TIFF image");
            }

            // Check the tif read to the raster correctly
            if (!TIFFReadRGBAImage(tif, width, height, raster, 0)) {

                TIFFClose(tif);

                Debug::Error("TiffImage::readImageToTiffClass: Could not read"
                             "raster of TIFF image");

                _TIFFfree(raster);
                break;
            }

            // iterate through all the pixels of the tif

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {

                    //int window_Z = z - z_start;
                    // y is always flipped
                    int window_Y = (height - y - 1) - y_start;
                    int window_X = x - x_start;

                    if (/*window_Z >= 0 and window_Z < windowDepth  and*/
                        window_Y >= 0 and window_Y < windowHeight and
                        window_X >= 0 and window_X < windowWidth) {
                        // read the current pixel of the TIF
                        uint32& TiffPixel = raster[y*width + x];
                        // Take last 8 bits
                        unsigned low8bits = TiffPixel & 0xFF;
                        tiffImage->setValue(window_X,
                                            window_Y,
                                            window_Z,
                                            low8bits);
                    }
                }
            }
            _TIFFfree(raster); // release temp memory

        } while (TIFFReadDirectory(tif)); // get the next tif
    TIFFClose(tif);

    } else {
        Debug::Error("TiffImage::readImageToTiffClass: could not read file "
        + std::string(path) + ")");
    }

    //Debug::Info("TiffImage::readImageToTiffClass: Leaving");

    return tiffImage;
}

TiffImageRef::TiffImageRef()
                           : path("")
                           , tiffImage(NULL)
                           , currentZ(0)
                           , _width(0)
                           , _height(0)
                           , _depth(0)
{
    Debug::Info("TiffImage::TiffImage (default): Entering");

    _histogram = NULL;

    Debug::Info("TiffImage::TiffImage (default): Leaving");
}

TiffImageRef::TiffImageRef(char const pathIn[],
                           const int width,
                           const int height,
                           const int depth)
                           : path(std::string(pathIn))
                           , tiffImage(NULL)
                           , imageSlicePointer(NULL)
                           , currentZ(0)
                           , _width(width)
                           , _height(height)
                           , _depth(depth)
                           , tiffImageRetrievalInProgress(false)
    {
    Debug::Info("TiffImage::TiffImage: Entering");

    TiffReader tiffReader;
    tiffImage = tiffReader.readImageToTiffClass(pathIn,
                                                0,
                                                width,
                                                0,
                                                height,
                                                currentZ,
                                                currentZ+40);

    _histogram = (int*) malloc(256 * sizeof(int));

    if (_histogram == NULL)
    {
        Debug::Error("TiffImage::TiffImage: malloc of _histogram did not "
                     "succeed!");
    }

    for (int i = 0; i < 256; ++i)
    {
        _histogram[i] = 0;
    }

    Debug::Info("TiffImage::TiffImage: Leaving");
}

TiffImageRef::~TiffImageRef()
{
    if (tiffImage != NULL) delete tiffImage;
    if (_histogram != NULL) free(_histogram);
    if (imageSlicePointer != NULL) free(imageSlicePointer);
}

float TiffImageRef::getValue(const int x,
                             const int y,
                             const int z,
                             const int subWindowWidth,
                             const int subWindowHeight)
{
    //Debug::Info("TiffImageRef::getValue: Entering");
    updateImageSlice(z);

    float value = 0;

    while (true)
    {
        tiffImageMutex.lock();

        int subImageZ = z - tiffImage->getZOffset();

        bool availableSlice = tiffImage->availableSlice(subImageZ);

        if (availableSlice)
        {
            value = tiffImage->getValue(x,y,subImageZ,subWindowWidth, subWindowHeight);
            tiffImageMutex.unlock();
            break;
        }
        tiffImageMutex.unlock();
    }
    //Debug::Info("TiffImageRef::getValue: Leaving");
    return value;
};

float *TiffImageRef::getSlice(const int z)
{
    //Debug::Info("TiffImageRef::getSlice: Entering");
    updateImageSlice(z);

    if (imageSlicePointer == NULL)
    {
        imageSlicePointer = (float *) malloc(_width*_height * sizeof(float));
    }

    float *tmp_imageSlicePointer;

    while (true)
    {
        tiffImageMutex.lock();

        int subImageZ = z - tiffImage->getZOffset();

        bool availableSlice = tiffImage->availableSlice(subImageZ);

        if (availableSlice)
        {
            tmp_imageSlicePointer = tiffImage->getSlice(subImageZ);
            // Don't unlock yet, now we have the slice we need to keep the
            // lock lest it gets invalid in the meantime
            break;
        }
        tiffImageMutex.unlock();
    }

    for (int i = 0; i < _width*_height; ++i)
    {
        imageSlicePointer[i] = tmp_imageSlicePointer[i];
    }
    tiffImageMutex.unlock();

    //Debug::Info("TiffImageRef::getSlice: Leaving");
    return imageSlicePointer;
};

void TiffImageRef::logHistogram()
{
    Debug::Info("## Pixel Histogram ##");
    Debug::Info("Value - Count:");
    for (int i = 0; i < 256; ++i)
    {
        Debug::Info(STR(i) + " - " + STR(_histogram[i]));
    }
}

void TiffImageRef::updateImageSlice(const int z)
{
    //Debug::Info("TiffImageRef::updateImageSlice: Entering");

    if (z < 0 or z >= _depth)
    {
        return;
    }

    tiffImageMutex.lock();

    if (tiffImageRetrievalInProgress)
    {
        tiffImageMutex.unlock();
        return;
    }

    tiffImageRetrievalInProgress = true;
    tiffImageMutex.unlock();

    int curBeginZ;
    int curEndZ;

    int curDepth;

    tiffImageMutex.lock();

    curBeginZ = tiffImage->getZOffset();
    curDepth = tiffImage->getDepth();
    curEndZ = curBeginZ + curDepth;

    tiffImageMutex.unlock();

    int curCenterZ = (curBeginZ + curEndZ) / 2;

    // We are still in the safe range, don't do anything.
    if (z >= curCenterZ - curDepth / 4 and z < curCenterZ + curDepth / 4)
    {
        tiffImageMutex.lock();
        tiffImageRetrievalInProgress = false;
        tiffImageMutex.unlock();
        return;
    }

    // We are at the top of the image, don't do anything.
    if (z < curCenterZ - curDepth / 4 and curBeginZ == 0)
    {
        tiffImageMutex.lock();
        tiffImageRetrievalInProgress = false;
        tiffImageMutex.unlock();
        return;
    }

    // We are at the bottom of the image, don't do anything.
    if (z >= curCenterZ + curDepth / 4 and curEndZ == _depth)
    {
        tiffImageMutex.lock();
        tiffImageRetrievalInProgress = false;
        tiffImageMutex.unlock();
        return;
    }

    int newBeginZ = z - 20;
    int newEndZ = z + 20;

    if (newBeginZ < 0)
    {
        newBeginZ = 0;
        newEndZ = std::min(40, curDepth);
    }

    if (newEndZ > _depth)
    {
        newEndZ = _depth;
        newBeginZ = newEndZ - std::min(40, curDepth);
    }

    std::thread t(grabImageSlice,
                  path,
                  _width,
                  _height,
                  newBeginZ,
                  newEndZ,
                  std::ref(tiffImageRetrievalInProgress),
                  std::ref(tiffImage),
                  std::ref(tiffImageMutex));

    t.detach();

    Debug::Info("TiffImageRef::updateImageSlice: Leaving");
}