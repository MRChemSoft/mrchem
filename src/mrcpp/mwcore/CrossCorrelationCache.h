#pragma once

#include "CrossCorrelation.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

#define getCrossCorrelationCache(T,X)\
    CrossCorrelationCache<T> &X = CrossCorrelationCache<T>::getInstance()

template <int T>
class CrossCorrelationCache: public ObjectCache<CrossCorrelation> {
public:
    static CrossCorrelationCache<T> &getInstance() {
        static CrossCorrelationCache<T> theCrossCorrelationCache;
        return theCrossCorrelationCache;
    }
    virtual void load(int order);
    virtual CrossCorrelation &get(int order);

    const Eigen::MatrixXd &getLMatrix(int order);
    const Eigen::MatrixXd &getRMatrix(int order);

    int getType() const { return this->type; 	}

    const std::string &getLibPath() { return this->libPath; }
    void setLibPath(const std::string &path) { this->libPath = path; }
protected:
    int type;
    std::string libPath; ///< Base path to filter library
private:
    CrossCorrelationCache();
    virtual ~CrossCorrelationCache() { }
    CrossCorrelationCache(CrossCorrelationCache<T> const &ccc) { }
    CrossCorrelationCache<T> &operator=(CrossCorrelationCache<T> const &ccc) { return *this; }
};

