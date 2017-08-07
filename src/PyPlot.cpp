#include "PyPlot.hpp"
#include <iostream>

Graph2D* Graph2D::_instance = nullptr;

Graph2D* Graph2D::Instance() {
    if (_instance == nullptr) {
        _instance = new Graph2D;
    }
    return _instance;
}

Graph2D::Graph2D() : axes3dloaded(false) {
    initPython();
}

Graph2D::~Graph2D() {
    Py_Finalize();
}

void Graph2D::loadAxes3D() {
    if (!axes3dloaded) {
        std::string s = "from mpl_toolkits.mplot3d import Axes3D\nimport matplotlib.tri as mtri\nfig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\nfig.suptitle('Test ellipsoid example',fontsize=20)\nax.set_xlabel('X axis')\nax.set_ylabel('Y axis')\nax.set_zlabel('Z axis')";
        PyRun_SimpleString(s.c_str());
        axes3dloaded = true;
    }
}

std::vector<double> Graph2D::linspace(double s, double e, int n) {
    std::vector<double> result;
    for (int i = 0; i < n; ++i) {
        result.push_back(i*(e-s) / n + s);
    }
    return result;
}

void Graph2D::x_lim(double x_min, double x_max) {
    std::string s = "plt.xlim(" + STR(x_min) + "," + STR(x_max) + ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::y_lim(double y_min, double y_max) {
    std::string s = "plt.ylim(" + STR(y_min) + "," + STR(y_max) + ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::xlabel(std::string axisName) {
    std::string s = "plt.xlabel('" + axisName + "')";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::ylabel(std::string axisName) {
    std::string s = "plt.ylabel('" + axisName + "')";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::title(std::string t) {
    PyRun_SimpleString(("plt.title('" + t + "')").c_str());
}


void Graph2D::legend(std::vector<std::string> labels, int location) {
    std::string s = "plt.legend([";
    bool first = true;
    for (int i = 0; i < labels.size(); ++i) {
        if (first) first = false;
        else s += ", ";
        s += "'" + labels[i] + "'";
    }
    s += "],loc =" + STR(location) + ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::grid(bool on) {
    if (on) {
        PyRun_SimpleString("plt.grid(True)");
    } else {
        PyRun_SimpleString("plt.grid(False)");
    }
}

void Graph2D::initPython() {
    Py_Initialize();
    PyRun_SimpleString("import matplotlib.pyplot as plt\nimport numpy as np\n");
}

std::string Graph2D::pyBool(bool b) {
    if (b) {
        return "True";
    } else {
        return "False";
    }
}

void Graph2D::plot(std::vector<double> x_vals, std::vector<double> y_vals, std::string color, double lw) {
    setValues(x_vals, "xes");
    setValues(y_vals, "yes");
    std::string s = "plt.plot(xes, yes, color=" + color + ",linewidth=" + STR(lw) + ")\ndel xes\ndel yes\n";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::plot(Mat x_vals, Mat y_vals, std::string color, double lw) {
    setValues(x_vals, "xes");
    setValues(y_vals, "yes");
    std::string s = "plt.plot(xes, yes, color=" + color + ",linewidth=" + STR(lw) + ")\ndel xes\ndel yes\n";
    PyRun_SimpleString(s.c_str());
}


void Graph2D::step(std::vector<double> x_vals, std::vector<double> y_vals) {
    setValues(x_vals, "step_xes");
    setValues(y_vals, "step_yes");
    PyRun_SimpleString("plt.step(step_xes, step_yes, where='post')\ndel step_xes\ndel step_yes\n");
}

void Graph2D::hist(std::vector<double> values, int bins, bool normed, double alpha) {
    setValues(values, "hist_v");
    //std::string s = "n, bins, patches = plt.hist(hist_v, bins =" + STR(bins) + ", normed = " + pyBool(normed) + ")";
    //PyRun_SimpleString(s.c_str()); //[0.0, 0.09, 0.18, 0.27, 0.36, 0.45, 0.55, 0.64, 0.73, 0.82, 0.91, 1.0]
    std::string s = "x, bins, p = plt.hist(hist_v, bins =" + STR(bins) + ", normed = " + pyBool(normed) + ", alpha = " + STR(alpha) + ")";
    if (normed) {
        s += "\nfor item in p:\n    item.set_height(item.get_height()/sum(x))";
    }
    PyRun_SimpleString(s.c_str());
}

void Graph2D::hist(std::vector<double> values, std::string bins, bool normed, double alpha) {
    setValues(values, "hist_v");
    //std::string s = "n, bins, patches = plt.hist(hist_v, bins =" + STR(bins) + ", normed = " + pyBool(normed) + ")";
    //PyRun_SimpleString(s.c_str()); //[0.0, 0.09, 0.18, 0.27, 0.36, 0.45, 0.55, 0.64, 0.73, 0.82, 0.91, 1.0]
    std::string s = "x, bins, p = plt.hist(hist_v, bins =" + bins + ", normed = " + pyBool(normed) + ", alpha = " + STR(alpha) + ")";
    if (normed) {
        s += "\nfor item in p:\n    item.set_height(item.get_height()/sum(x))";
    }
    PyRun_SimpleString(s.c_str());
}

void Graph2D::scatter(std::vector<double> x_vals, std::vector<double> y_vals, std::string color, double size) {
    setValues(x_vals, "sca_xes");
    setValues(y_vals, "sca_yes");
    std::string s = "plt.scatter(sca_xes, sca_yes, s=" + STR(size) + ", c=" + color +  ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::scatter(Mat x_vals, Mat y_vals, std::string color, double size) {
    setValues(x_vals, "sca_xes");
    setValues(y_vals, "sca_yes");
    std::string s = "plt.scatter(sca_xes, sca_yes, s=" + STR(size) + ", c=" + color +  ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::imshow(Mat image, float vmin, float vmax) {
    setValues(image, "image", true);
    if (vmin == vmax) {
        PyRun_SimpleString("plt.imshow(image, cmap='gray')");
    } else {
        std::string s = "plt.imshow(image, cmap='gray', vmin=" + STR(vmin) + ", vmax=" + STR(vmax) + ")";
        PyRun_SimpleString(s.c_str());
    }

}

void Graph2D::imshow(Mat image, Mat contour, float vmin, float vmax) {
    setValues(image, "image", true);
    setValues(contour, "contour", true);
    if (vmin == vmax) {
        PyRun_SimpleString("plt.imshow(image, cmap='gray')");
    } else {
        std::string s = "plt.imshow(image, cmap='gray', vmin=" + STR(vmin) + ", vmax=" + STR(vmax) + ")";
        PyRun_SimpleString(s.c_str());
    }
    PyRun_SimpleString("plt.contour(contour,levels = [0.1],colors='red')");
}



void Graph2D::scatter3d(Mat X, std::string color) {
    int N = X.rows();
    loadAxes3D();
    Mat x = X.block(0,0,N,1);
    Mat y = X.block(0,1,N,1);
    Mat z = X.block(0,2,N,1);
    setValues(x, "sca3d_xes");
    setValues(y, "sca3d_yes");
    setValues(z, "sca3d_zes");
    std::string s = "ax.scatter(sca3d_xes, sca3d_yes, sca3d_zes, s=5, c='" + color +  "')";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::plotSphere3d() {
    loadAxes3D();
    PyRun_SimpleString("import numpy as np\nfrom matplotlib import cm, colors\nr = 1\npi = np.pi\ncos = np.cos\nsin = np.sin\nphi, theta = np.mgrid[0.0:pi:25j, 0.0:2.0*pi:25j]\nx = r*sin(phi)*cos(theta)\ny = r*sin(phi)*sin(theta)\nz = r*cos(phi)\nax.plot_surface(x, y, z,  rstride=1, cstride=1, color='c', alpha=0.4, linewidth=0)");
    PyRun_SimpleString("ax.set_xlim([-1.5,1.5])\nax.set_ylim([-1.5,1.5])\nax.set_zlim([-1.5,1.5])\nax.set_aspect('equal')");
}

void Graph2D::plotTriSurfVesicle(Mat X) {

    int N = X.rows();

    loadAxes3D();

    Mat x = X.block(0,0,N,1);
    Mat y = X.block(0,1,N,1);
    Mat z = X.block(0,2,N,1);

    setValues(x, "trives3d_xes");
    setValues(y, "trives3d_yes");
    setValues(z, "trives3d_zes");

    PyRun_SimpleString("trives3d_xes = np.array(trives3d_xes)");
    PyRun_SimpleString("trives3d_yes = np.array(trives3d_yes)");
    PyRun_SimpleString("trives3d_zes = np.array(trives3d_zes)");

    PyRun_SimpleString("xa = trives3d_xes[trives3d_zes >= 0].flatten()\nya = trives3d_yes[trives3d_zes >= 0].flatten()\nza = trives3d_zes[trives3d_zes >= 0].flatten()");
    PyRun_SimpleString("minz = np.min(za)");
    PyRun_SimpleString("xb = trives3d_xes[trives3d_zes <= minz + 0.01].flatten()\nyb = trives3d_yes[trives3d_zes <= minz + 0.01].flatten()\nzb = trives3d_zes[trives3d_zes <= minz + 0.01].flatten()");
    PyRun_SimpleString("f = np.vectorize(lambda xy,z: xy/(0.1 + abs(z-minz)))");
    PyRun_SimpleString("txa = f(xa,za)");
    PyRun_SimpleString("txb = f(xb,zb)");
    PyRun_SimpleString("tya = f(ya,za)");
    PyRun_SimpleString("tyb = f(yb,zb)");
    PyRun_SimpleString("tria = mtri.Triangulation(txa, tya)\ntrib = mtri.Triangulation(txb, tyb)");

    PyRun_SimpleString("ax.scatter(xa, ya, za, c='red', alpha = 0.5)");
    PyRun_SimpleString("ax.scatter(xb, yb, zb, c='red', alpha = 0.5)");
    PyRun_SimpleString("ax.plot_trisurf(xa, ya, za, triangles=tria.triangles, color='orange')");
    PyRun_SimpleString("ax.plot_trisurf(xb, yb, zb, triangles=trib.triangles, color='orange')");

}


void Graph2D::subplot(int numRows, int numCols, int plotID) {
    std::string s = "plt.subplot(" + STR(numRows) + "," + STR(numCols) + "," + STR(plotID) + ")";
    PyRun_SimpleString(s.c_str());
}

void Graph2D::show() {
    PyRun_SimpleString("plt.show()\n");
}

void Graph2D::setValues(std::vector<double> &v, std::string varName) {

    PyRun_SimpleString((varName + " = []").c_str());
    for (int i = 0; i < v.size(); ++i) {
        std::string s = varName + ".append(" + STR(v[i]) + ")\n";
        PyRun_SimpleString(s.c_str());
    }
}

void Graph2D::setValues(Mat &v, std::string varName, bool conserveSize) {

    if (conserveSize) {

        int M = v.rows();
        int N = v.cols();
        PyRun_SimpleString((varName + " = np.empty((" + STR(M) + "," + STR(N) + "))\n").c_str());

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                std::string s = varName + "[" + STR(j) + "," + STR(i) + "] = " + STR(v(j,i)) + " \n";
                PyRun_SimpleString(s.c_str());
            }
        }
    } else {

        v.conservativeResize(v.rows() * v.cols(), 1);
        PyRun_SimpleString((varName + " = []").c_str());

        for (int i = 0; i < v.rows(); ++i) {

            std::string s = varName + ".append(" + STR(v(i,0)) + ")\n";
            PyRun_SimpleString(s.c_str());

        }
    }
}

Graph2D* pyGraph() {
    return Graph2D::Instance();
}