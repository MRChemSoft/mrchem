#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TelePrompter.h"
#include "mrcpp_declarations.h"

template<int D>
class OperatorStatistics {
public:
    OperatorStatistics();
    virtual ~OperatorStatistics();

    void flushNodeCounters();
    void incrementFNodeCounters(const MWNode<D> &fNode, int ft, int gt);
    void incrementGNodeCounters(const MWNode<D> &gNode);

    friend std::ostream& operator<<(std::ostream &o, const OperatorStatistics &os) {
        o << std::setw(8);
        o << "*OperatorFunc statistics: " << std::endl << std::endl;
        o << "  Total calculated gNodes      : " << os.totGCount << std::endl;
        o << "  Total applied fNodes         : " << os.totFCount << std::endl;
        o << "  Total applied genNodes       : " << os.totGenCount << std::endl << std::endl;
        o << "  By components:" << std::endl << *os.totCompCount << std::endl;
        return o;
    }

protected:
    int nThreads;
    int totFCount;
    int totGCount;
    int totGenCount;
    int *fCount;
    int *gCount;
    int *genCount;
    Eigen::Matrix<int, 8, 8> *totCompCount;
    Eigen::Matrix<int, 8, 8> **compCount;
};

