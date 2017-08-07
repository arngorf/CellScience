#include "TiffReader.hpp"

#include <sstream>
#include <string>

void DummyHandler(const char* module, const char* fmt, va_list ap) {
    // ignore errors and warnings (or handle them your own way)
}

TiffImage::TiffImage() :_width(0), _height(0), _depth(0)
{
    Debug::Info("TiffImage::TiffImage (default): Entering");

    _image = NULL;
    _histogram = NULL;

    Debug::Info("TiffImage::TiffImage (default): Leaving");
}

TiffImage::TiffImage(const unsigned int width,
                     const unsigned int height,
                     const unsigned int depth) :
                     _width(width), _height(height), _depth(depth)
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
                           const unsigned int width,
                           const unsigned int height,
                           const unsigned int depth)
                           : path(std::string(pathIn))
                           , tiffImage(NULL)
                           , currentZ(0)
                           , _width(width)
                           , _height(height)
                           , _depth(depth)
    {
    Debug::Info("TiffImage::TiffImage: Entering");

    TiffReader tiffReader;
    tiffImage = tiffReader.readImageToTiffClass(pathIn,
                                                0,
                                                width,
                                                0,
                                                height,
                                                currentZ,
                                                currentZ+1);

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
                                         windowDepth);

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