#ifndef TIFFREADER_HPP
#define TIFFREADER_HPP

#include <tiffio.h>
#include <vector>
#include "types.hpp"
#include <iostream>

#include "Debug.hpp"

void DummyHandler(const char* module,
                  const char* fmt,
                  va_list ap);

class TiffImage {
public:

    TiffImage();

    TiffImage(const unsigned int width,
              const unsigned int height,
              const unsigned int depth);

    ~TiffImage();

    unsigned int getWidth()
    {
        return _width;
    };

    unsigned int getHeight()
    {
        return _height;
    };

    unsigned int getDepth()
    {
        return _depth;
    };

    void setValue(const unsigned int i,
                  const unsigned int value)
    {
        _image[i] = (float) value;

        if (value < 255)
        {
            _histogram[value] += 1;
        }

        else
        {
            Debug::Warning("TiffImage::setValue: pixel value " + STR(value) +
                           " out of range");
        }
    };

    void setValue(const unsigned int x,
                  const unsigned int y,
                  const unsigned int z,
                  const unsigned int value)
    {
        _image[globalIndex(x, y, z, _width, _height)] = (float) value;

        if (value < 255)
        {
            _histogram[value] += 1;
        }
        else
        {
            Debug::Warning("TiffImage::setValue: pixel value " + STR(value) +
                           " out of range");
        }
    };

    float getValue(const unsigned int i)
    {
        return _image[i];
    };

    float getValue(const unsigned int x,
                    const unsigned int y,
                    const unsigned int z,
                    const int subWindowWidth = -1,
                    const int subWindowHeight = -1)
    {
        if (subWindowWidth == -1 or subWindowHeight == -1) {
            return _image[globalIndex(x, y, z, _width, _height)];
        }
        return _image[globalIndex(x, y, z, subWindowWidth, subWindowHeight)];
    };

    float *getSlice(const unsigned int z)
    {
        return &_image[globalIndex(0, 0, z, _width, _height)];
    };

    void logHistogram()
    {
        Debug::Info("## Pixel Histogram ##");
        Debug::Info("Value - Count:");
        for (int i = 0; i < 256; ++i)
        {
            Debug::Info(STR(i) + " - " + STR(_histogram[i]));
        }
    }

private:
    float *_image;
    int *_histogram;

    unsigned int _width;
    unsigned int _height;
    unsigned int _depth;

    unsigned int globalIndex(const unsigned int x,
                             const unsigned int y,
                             const unsigned int z,
                             const int width,
                             const int height)
    {
        return z * width*height + y * width + x;
    }
};

class TiffReader {

public:

    int read(char const path[], float *image,
             int _xOffset, int _yOffset, int _zOffset,
             int _width, int _height, int _depth);

    int readToImageHistogram(char const path[], std::vector<int> &IFreq);

    TiffImage *readImageToTiffClass(char const path[],
                                    const int x_start, const int x_end,
                                    const int y_start, const int y_end,
                                    const int z_start, const int z_end);


};

class TiffImageRef {
public:

    TiffImageRef();

    TiffImageRef(char const path[],
                 const unsigned int width,
                 const unsigned int height,
                 const unsigned int depth);

    ~TiffImageRef();

    unsigned int getWidth()
    {
        return _width;
    };

    unsigned int getHeight()
    {
        return _height;
    };

    unsigned int getDepth()
    {
        return _depth;
    };

    float getValue(const unsigned int x,
                    const unsigned int y,
                    const unsigned int z,
                    const int subWindowWidth = -1,
                    const int subWindowHeight = -1)
    {
        updateImageSlice(z);

        return tiffImage->getValue(x,y,0,subWindowWidth, subWindowHeight);
    };

    float *getSlice(const unsigned int z)
    {
        updateImageSlice(z);

        return tiffImage->getSlice(0);
    };

    void logHistogram()
    {
        Debug::Info("## Pixel Histogram ##");
        Debug::Info("Value - Count:");
        for (int i = 0; i < 256; ++i)
        {
            Debug::Info(STR(i) + " - " + STR(_histogram[i]));
        }
    }

private:
    std::string path;
    TiffImage *tiffImage;
    int *_histogram;

    unsigned int currentZ;

    unsigned int _width;
    unsigned int _height;
    unsigned int _depth;

    unsigned int globalIndex(const unsigned int x,
                             const unsigned int y,
                             const unsigned int z,
                             const int width,
                             const int height)
    {
        return z * width*height + y * width + x;
    }

    void updateImageSlice(const unsigned int z)
    {
        if (z != currentZ and z < _depth)
        {
            if (tiffImage != NULL) delete tiffImage;
            currentZ = z;
            TiffReader tiffReader;
            tiffImage = tiffReader.readImageToTiffClass(path.c_str(),
                                                        0,
                                                        _width,
                                                        0,
                                                        _height,
                                                        currentZ,
                                                        currentZ+1);
        }
    }
};

#endif // TIFFREADER_HPP