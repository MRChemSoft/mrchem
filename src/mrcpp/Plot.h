#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "mrcpp_declarations.h"

template<int D>
class Plot {
public:
    Plot(int npts = 1000, const double *a = 0, const double *b = 0);
    virtual ~Plot() { }

    void setRange(const double *a, const double *b);
    void setNPoints(int npts);
    void setSuffix(int t, const std::string &s);

    void linePlot(const RepresentableFunction<D> &func, const std::string &fname);
    void surfPlot(const RepresentableFunction<D> &func, const std::string &fname);
    void cubePlot(const RepresentableFunction<D> &func, const std::string &fname);

    void linePlot(FunctionTree<D> &func, const std::string &fname);
    void surfPlot(FunctionTree<D> &func, const std::string &fname);
    void cubePlot(FunctionTree<D> &func, const std::string &fname);
    void gridPlot(const MWTree<D> &tree, const std::string &fname);

    Eigen::VectorXd &linePlot(RepresentableFunction<D> &func);
    Eigen::VectorXd &surfPlot(RepresentableFunction<D> &func);
    Eigen::VectorXd &cubePlot(RepresentableFunction<D> &func);

    enum type {
        Line, Surface, Cube, Grid
    };

protected:
    std::ofstream fstrm;
    std::ofstream *fout;
    int nPoints;
    double A[D]; ///< lower left corner
    double B[D]; ///< upper right corner
    std::map<int, std::string> suffix;
    Eigen::MatrixXd coords;
    Eigen::VectorXd values;

    void calcLineCoordinates();
    void calcSurfCoordinates();
    void calcCubeCoordinates();

    void evaluateFunction(const RepresentableFunction<D> &func);
    void evaluateFunction(FunctionTree<D> &tree);
    bool verifyRange();

    void writeLineData();
    void writeSurfData();
    virtual void writeCubeData();

    void writeGrid(const MWTree<D> &tree);
    void writeNodeGrid(const MWNode<D> &node, const std::string &color);

private:
    void openPlot(const std::string &fname);
    void closePlot();
};

