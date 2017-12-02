#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <string>

typedef Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic>
FilterBlock;

class MWFilter {
public:
    MWFilter(int k, int t, const std::string &lib = "");
    MWFilter(int t, const Eigen::MatrixXd &data);
    virtual ~MWFilter() { }

    virtual void apply(Eigen::MatrixXd &data) const;
    virtual void apply(Eigen::VectorXd &data) const;
    virtual void applyInverse(Eigen::MatrixXd &data) const;
    virtual void applyInverse(Eigen::VectorXd &data) const;

    int getOrder() const { return this->order; }
    int getType() const { return this->type; }

    const Eigen::MatrixXd &getFilter() const { return this->filter; }
    const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const;
    const Eigen::MatrixXd &getCompressionSubFilter(int i) const;
    const Eigen::MatrixXd &getReconstructionSubFilter(int i) const;

    static void setDefaultLibrary(const std::string &dir);
    static const std::string &getDefaultLibrary() { return default_filter_lib; }

protected:
    int type;
    int order;
    int dim;

    Eigen::MatrixXd filter; ///< Full MW-transformation matrix
    Eigen::MatrixXd G0;
    Eigen::MatrixXd G1;
    Eigen::MatrixXd H0;
    Eigen::MatrixXd H1;
    // Transpose
    Eigen::MatrixXd G0t;
    Eigen::MatrixXd G1t;
    Eigen::MatrixXd H0t;
    Eigen::MatrixXd H1t;

    std::string H_path;
    std::string G_path;
    static std::string default_filter_lib;

    void setFilterPaths(const std::string &lib);
    void readFilterBin();
    void fillFilterBlocks();
};


