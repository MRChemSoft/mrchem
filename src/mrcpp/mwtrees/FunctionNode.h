#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "MWNode.h"
#include "FunctionTree.h"

template<int D>
class FunctionNode : public MWNode<D> {
public:
    FunctionTree<D> &getFuncTree() { return static_cast<FunctionTree<D> &>(*this->tree); }
    FunctionNode<D> &getFuncParent() { return static_cast<FunctionNode<D> &>(*this->parent); }
    FunctionNode<D> &getFuncChild(int i) { return static_cast<FunctionNode<D> &>(*this->children[i]); }

    const FunctionTree<D> &getFuncTree() const { return static_cast<const FunctionTree<D> &>(*this->tree); }
    const FunctionNode<D> &getFuncParent() const { return static_cast<const FunctionNode<D> &>(*this->parent); }
    const FunctionNode<D> &getFuncChild(int i) const { return static_cast<const FunctionNode<D> &>(*this->children[i]); }

    virtual void setValues(const Eigen::VectorXd &vec);
    virtual void getValues(Eigen::VectorXd &vec);

    friend class FunctionTree<D>;

protected:
    FunctionNode() : MWNode<D>() { }
    virtual ~FunctionNode() { assert(this->tree == 0); }

    double evalf(const double *r);

    double integrate() const;
    double dotScaling(const FunctionNode<D> &ket) const;
    double dotWavelet(const FunctionNode<D> &ket) const;

    double evalScaling(const double *r) const;
    double integrateLegendre() const;
    double integrateInterpolating() const;

};

