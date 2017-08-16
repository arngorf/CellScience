#ifndef TIFFREADER_HPP
#define TIFFREADER_HPP

#include <tiffio.h>

#include <vector>
#include <iostream>
#include <thread>
#include <mutex>

#include "types.hpp"

#include "Debug.hpp"

void DummyHandler(const char* module,
                  const char* fmt,
                  va_list ap);

class TiffImage
{
public:

    TiffImage();

    TiffImage(const int width,
              const int height,
              const int depth,
              const int zOffset);

    ~TiffImage();

    inline int getWidth()
    {
        return _width;
    };

    inline int getHeight()
    {
        return _height;
    };

    inline int getDepth()
    {
        return _depth;
    };

    inline int getZOffset()
    {
        return _zOffset;
    };

    inline bool availableSlice(int z)
    {
        return z >= 0 and z < _depth;
    }

    void setValue(const int i,
                  const int value)
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

    void setValue(const int x,
                  const int y,
                  const int z,
                  const int value)
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

    float getValue(const int i)
    {
        return _image[i];
    };

    float getValue(const int x,
                   const int y,
                   const int z,
                   const int subWindowWidth = -1,
                   const int subWindowHeight = -1)
    {
        if (subWindowWidth == -1 or subWindowHeight == -1) {
            return _image[globalIndex(x, y, z, _width, _height)];
        }
        return _image[globalIndex(x, y, z, subWindowWidth, subWindowHeight)];
    };

    float *getSlice(const int z)
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
    };

private:
    float *_image;
    int *_histogram;

    int _width;
    int _height;
    int _depth;

    int _zOffset;

    int globalIndex(const int x,
                    const int y,
                    const int z,
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

void grabImageSlice(const std::string path,
                    const int width,
                    const int height,
                    const int beginZ,
                    const int endZ,
                    bool &grabbing,
                    TiffImage *&imagePointer,
                    std::mutex &tiffImageMutex);

class TiffImageRef {
public:

    TiffImageRef();

    TiffImageRef(char const path[],
                 const int width,
                 const int height,
                 const int depth);

    ~TiffImageRef();

    inline int getWidth()
    {
        return _width;
    }

    inline int getHeight()
    {
        return _height;
    }

    inline int getDepth()
    {
        return _depth;
    }

    float getValue(const int x,
                   const int y,
                   const int z,
                   const int subWindowWidth = -1,
                   const int subWindowHeight = -1);

    float *getSlice(const int z);

    void logHistogram();

private:
    std::string path;
    TiffImage *tiffImage;
    float *imageSlicePointer;
    int *_histogram;

    int currentZ;

    int _width;
    int _height;
    int _depth;

    bool tiffImageRetrievalInProgress;

    std::mutex tiffImageMutex;

    int globalIndex(const int x,
                             const int y,
                             const int z,
                             const int width,
                             const int height)
    {
        return z * width*height + y * width + x;
    }

    void updateImageSlice(const int z);
};

#endif // TIFFREADER_HPP