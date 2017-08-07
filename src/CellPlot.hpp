#ifndef CELL_PLOT_HPP
#define CELL_PLOT_HPP

#include <vector>

#include "CellObject.hpp"

class CellPlot {
public:
    CellPlot();
    ~CellPlot();
    void InterativeCubeField(float *phaseField, int nz, int ny, int nx, float levelSet = 0, bool invertPhasefield = false);
    void InteractiveGuessMesh(float *phaseField, int nz, int ny, int nx);
    void InteractiveImplicitSurface(float *phaseField, int nz, int ny, int nx);
    void MultiObjectInteractiveImplicitSurface(std::vector<CellObject*> objects);
};

#endif // CELL_PLOT_HPP