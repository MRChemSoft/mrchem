#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <string>

class CrossCorrelation {
public:
    CrossCorrelation(int k, int t, const std::string &lib = "");
    CrossCorrelation(int t, const Eigen::MatrixXd &ldata,
                            const Eigen::MatrixXd &rdata);
    virtual ~CrossCorrelation() { }

    int getType() const { return this->type; }
    int getOrder() const { return this->order; }
    const Eigen::MatrixXd &getLMatrix() const { return this->Left; }
    const Eigen::MatrixXd &getRMatrix() const { return this->Right; }

    static void setDefaultLibrary(const std::string &dir);
    static const std::string &getDefaultLibrary() { return default_ccc_lib; }

protected:
    int type;
    int order;

    Eigen::MatrixXd Left;
    Eigen::MatrixXd Right;

    std::string L_path;
    std::string R_path;
    static std::string default_ccc_lib;

    void setCCCPaths(const std::string &lib);
    void readCCCBin();
};

