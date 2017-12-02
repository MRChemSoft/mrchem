#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TreeCalculator.h"

template<int D>
class MWFDCalculator : public TreeCalculator<D> {
public:
    MWFDCalculator(int dir, int k, int n, FunctionTree<D> &f);
protected:
    const int apply_dir;
    const int diff_order;
    const int approx_order;
    FunctionTree<D> *f_tree;

    virtual void calcNode(MWNode<D> &node);

    void computeBlock(Eigen::VectorXd &out_pts, Eigen::MatrixXd &out_vals,
                      Eigen::VectorXd &inp_pts, Eigen::MatrixXd &inp_vals);
    void computeLine(Eigen::MatrixXd &output_data, Eigen::MatrixXd &input_data);
    double computePoint(int M, int N, double x0, Eigen::MatrixXd &input_data);

    void fetchInputData(const NodeIndex<D> &idx,
                        Eigen::VectorXd &pts,
                        Eigen::MatrixXd &vals);
    void fetchOutputPoints(const MWNode<D> &node,
                           Eigen::VectorXd &pts);
    void pushOutputValues(FunctionNode<D> &node,
                          Eigen::MatrixXd &vals);

    void getNodeData(FunctionNode<D> &node,
                     Eigen::VectorXd &pts,
                     Eigen::MatrixXd &vals);
};

