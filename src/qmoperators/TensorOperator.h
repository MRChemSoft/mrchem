#pragma once

namespace mrchem {

template<int I, class T>
class TensorOperator {
public:
    TensorOperator() { }
    virtual ~TensorOperator() { }

    void setup(double prec) { for (int i = 0; i < I; i++) this->oper[i].setup(prec); }
    void clear() { for (int i = 0; i < I; i++) this->oper[i].clear(); }

    T& operator[](int i) { return this->oper[i]; }
    const T& operator[](int i) const { return this->oper[i]; }

protected:
    T oper[I];
};

} //namespace mrchem
