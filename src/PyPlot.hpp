//////////////////////////////////////////////////////////////////////////////
// Simple C++ API for the python matplotlib library                         //
//                                                                          //
// Author: Hans Jacob Teglbjaerg Stephensen 2016                            //
//////////////////////////////////////////////////////////////////////////////

#ifndef PYPLOT_HPP
#define PYPLOT_HPP

#include <vector>
#include <string>
#include <Python.h>
#include <Eigen/Core>

#include "types.hpp"
#include "UtilityFunctions.hpp"

/*
 * Singleton 2D Graph Class
 * This class initializes a python program managing variables and plots.
 * This class is not safe to use along other python embedding. The immediate
 * problem is variable name overlaps.
 *
 * Use: 1) get pointer to Graph2D instance by calling pyGraph().
 *      2) use plot functions like plot, hist or scatter.
 *      3) call show() to display graph and clear stored plots.
 */
class Graph2D {
public:
    /*
     * Instance function. Returns pointer to current singleton graph instance.
     */
    static Graph2D* Instance();
private:
    /*
     * Singleton class constructor. Class holds own pointer.
     * Use the "new" Instance() to get pointer to graph instance.
     */
    Graph2D();

    /*
     * Destructor, this is to never be called manually.
     */

public:
    ~Graph2D();

    /**
     * given start point, end point and number of points, generate a vector
     * of values in between with equal spacing.
     * @Param: s: double starting value.
     * @Param: e: double ending value.
     * @Param: n: int number of values
     */
    std::vector<double> linspace(double s, double e, int n);

    /*
     * Simple plot function for plotting a graph of a 2d function.
     * @Param: x_vals: Vector of x-values.
     * @Param: y_vals: Vector of y-values.
     */
    void plot(std::vector<double> x_vals, std::vector<double> y_vals, std::string color = "'blue'", double lw = 1.0);

    /*
     * Simple plot function for plotting a graph of a 2d function.
     * @Param: x_vals: Eigen Matrix.
     * @Param: y_vals: Eigen Matrix.
     */
    void plot(Mat x_vals, Mat y_vals, std::string color = "'blue'", double lw = 1.0);

    /*
     * Simple plot function for plotting a graph of a 2d step function.
     * @Param: x_vals: Vector of x-values.
     * @Param: y_vals: Vector of y-values.
     */
    void step(std::vector<double> x_vals, std::vector<double> y_vals);

    /*
     * Histogram plotting function.
     * @Param: x_vals: Vector of x-values.
     * @Param: bins: Number of bins.
     * @Param: normed: Set if histogram is to be normalized.
     */
    void hist(std::vector<double> values, int bins = 10, bool normed = false, double alpha = 1.0);

    /*
     * Histogram plotting function.
     * @Param: x_vals: Vector of x-values.
     * @Param: bins: Alternative bin parameter
     * @Param: normed: Set if histogram is to be normalized.
     */
    void hist(std::vector<double> values, std::string bins, bool normed = false, double alpha = 1.0);

    /*
     * Scatter plot function.
     * @Param: x_vals: Vector of x-values.
     * @Param: y_vals: Vector of y-values.
     * @Param: color: color of dots.
     * @Param: size: Size^2 of dots.
     */
    void scatter(std::vector<double> x_vals, std::vector<double> y_vals, std::string color = "'blue'", double size = 20.0);

    /*
     * Scatter plot function.
     * @Param: x_vals: Eigen Matrix of x-values.
     * @Param: y_vals: Eigen Matrix of y-values.
     * @Param: color: color of dots.
     * @Param: size: Size^2 of dots.
     */
    void scatter(Mat x_vals, Mat y_vals, std::string color = "'blue'", double size = 20.0);

    void imshow(Mat image, float vmin = -1, float vmax = -1);

    void imshow(Mat image, Mat contour, float vmin = -1, float vmax = -1);

    /*
     * @Param: Nx3 matrix of N points in 3D.
     */
    void scatter3d(Mat X, std::string color = "blue");

    void plotSphere3d();

    void plotTriSurfVesicle(Mat X);

    /*
     * Set current subplot parameters
     * @Param: numRows: The row plot spacing parameter
     * @Param: numCols: The column plot spacing parameter
     * @Param: plotID: The placing of next plot. (Row major order)
     */
    void subplot(int numRows, int numCols, int plotID);

    /*
     * Show plot (clearing anything plotted).
     */
    void show();

    /*
     * Set x limits.
     */
    void x_lim(double x_min, double x_max);

    /*
     * Set y limits.
     */
    void y_lim(double y_min, double y_max);

    /*
     * Set x label.
     */
    void xlabel(std::string axisName);

    /*
     * Set y label.
     */
    void ylabel(std::string axisName);

    /*
     * Set plot title
     */
    void title(std::string t);

    /*
     * Set plot legend
     */
    void legend(std::vector<std::string> labels, int location);

    /*
     * Set if a grid on tick points should be shown.
     */
    void grid(bool on);

    /*
     * Stores an input vector v as a new python list with name 'varName' in
     * the running python environment.
     * Don't do something stupid with this. It was never meant for the public
     * use.
     */
    void setValues(std::vector<double> &v, std::string varName);

    /*
     * Stores an input Matrix v as a new python list with name 'varName' in
     * the running python environment.
     * Don't do something stupid with this. It was never meant for the public
     * use.
     */
    void setValues(Mat &v, std::string varName, bool conserveSize = false);

private:

    void loadAxes3D();

    static Graph2D* _instance;
    void initPython();
    std::string pyBool(bool b);
    bool axes3dloaded;
};

/*
 * 'Nicer' way to get the current graph instance.
 */
Graph2D* pyGraph();

#endif // PYPLOT_HPP
