/**
 * @file historylist.hpp
 * @author Moritz Gubler (moritz.gubler@unibas.ch)
 * @brief Implementation of the SQNM method. More informations about the algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
 * @date 2022-07-13
 * 
 */

#ifndef HISTORYLIST_HPP
#define HISTORYLIST_HPP
#include <Eigen/Dense>

namespace hlist_space{
  class HistoryList {
    private:
    int i_count = 0;
    int nhistx;
    int ndim;
    Eigen::VectorXd oldElem;

    public: 

    Eigen::MatrixXd hlist;
    Eigen::MatrixXd difflist;
    Eigen::MatrixXd normalized_difflist;


    HistoryList(int ndim_, int nhistx_) : nhistx(nhistx_), ndim(ndim_)
    {
      hlist.resize(ndim, nhistx);
      difflist.resize(ndim, nhistx);
      normalized_difflist.resize(ndim, nhistx);
    }

    int add(const Eigen::VectorXd &x_in)
    {
      if (i_count < nhistx){
        // list not yet full, append at the end.
        hlist.col(i_count) = x_in;
        i_count += 1;
        // calculate difference list if we already have more than one element
        if (i_count > 1)
        {
            difflist.block(0, 0, ndim, i_count - 1) = hlist.block(0, 1, ndim, i_count - 1) - hlist.block(0, 0, ndim, i_count -1);
            normalized_difflist.block(0, 0, ndim, i_count - 1) = difflist.block(0, 0, ndim, i_count - 1).colwise().normalized();
        }
        return i_count - 1;
      } else {
        // list is full, remove oldest element and shift all to new place.
        oldElem = hlist.col(0);
        hlist.block(0, 0, ndim, nhistx - 1) = hlist.block(0, 1, ndim, nhistx - 1);
        hlist.col(nhistx - 1) = x_in;
        // calculate difference list
        difflist.col(0) = hlist.col(0) - oldElem;
        normalized_difflist.col(0) = difflist.col(0) / difflist.col(0).norm();
        difflist.block(0, 1, ndim, nhistx - 1) = hlist.block(0, 1, ndim, nhistx - 1) - hlist.block(0, 0, ndim, nhistx - 1);
        normalized_difflist.block(0, 1, ndim, nhistx - 1) = difflist.block(0, 1, ndim, nhistx - 1).colwise().normalized();
        return nhistx;
        }
      }
    };
  }
#endif