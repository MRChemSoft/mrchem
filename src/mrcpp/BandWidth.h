/*
 * BandWidth.h
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <iostream>

class BandWidth {
public:
    BandWidth(int depth = 0)
            : widths(depth + 1, 5) {
        this->clear();
    }
    BandWidth &operator=(const BandWidth &bw) {
        this->widths = bw.widths;
        return *this;
    }
    virtual ~BandWidth() { }

    void clear() {
        this->widths.setConstant(-1);
    }

    int getDepth() const {
        return this->widths.rows() - 1;
    }

    bool isEmpty(int depth) const {
        if (depth > getDepth()) {
            return true;
        }
        if (this->widths(depth, 4) < 0) {
            return true;
        }
        return false;
    }

    int getWidth(int depth, int index) const {
        assert(depth >= 0);
        assert(index >= 0 and index < 4);
        if (depth > getDepth()) { // No bw at requested depth
            return -1;
        }
        return this->widths(depth, index);
    }

    int getMaxWidth(int depth) const {
        assert(depth >= 0);
        if (depth > getDepth()) { // No bw at requested depth
            return -1;
        }
        return this->widths(depth, 4);
    }

    void setWidth(int depth, int index, int wd) {
        assert(depth >= 0 and depth < getDepth());
        assert(index >= 0 and index < 4);
        assert(wd >= 0);
        this->widths(depth, index) = wd;
        if (wd > this->widths(depth, 4)) {
            this->widths(depth, 4) = wd;
        }
    }

    friend std::ostream& operator<<(std::ostream &o, const BandWidth &bw) {
        o << "  *BandWidths:" << std::endl;
        o << "   n      T   C   B   A  |  max " << std::endl;
        o << " -------------------------------" << std::endl;
        for (int depth = 0; depth <= bw.getDepth(); depth++) {
            o << std::setw(4) << depth << " | ";
            o << std::setw(4) << bw.widths(depth, 0);
            o << std::setw(4) << bw.widths(depth, 1);
            o << std::setw(4) << bw.widths(depth, 2);
            o << std::setw(4) << bw.widths(depth, 3) << "  | ";
            o << std::setw(4) << bw.widths(depth, 4) << std::endl;
        }
        o << std::endl;
        return o;
    }

private:
    Eigen::MatrixXi widths; /// column 5 stores max width at depth
};

