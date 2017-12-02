/*
 *
 * \breif FilterCache is a static class taking care of loading and
 * unloading MultiWavelet filters, and their tensor counter parts.
 *
 * All data in FilterCache is static, and thus shared amongst all
 * instance objects. The type of filter, Legendre or Interpolating is
 * determined by a template variable so that both types of filters can
 * co-exist.
 *
 */

#pragma once

#include "MWFilter.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

#define getFilterCache(T,X)\
    FilterCache<T> &X = FilterCache<T>::getInstance()
#define getLegendreFilterCache(X)\
    FilterCache<Legendre> &X = FilterCache<Legendre>::getInstance()
#define getInterpolatingFilterCache(X)\
    FilterCache<Interpol> &X = FilterCache<Interpol>::getInstance()

/** This class is an abstract base class for the various filter caches.
 * It's needed in order to be able to use the actual filter caches
 * without reference to the actual filter types */
class BaseFilterCache: public ObjectCache<MWFilter> {
public:
    virtual void load(int order) = 0;
    virtual MWFilter &get(int order) = 0;
    virtual const Eigen::MatrixXd &getFilterMatrix(int order) = 0;
    virtual const std::string &getLibPath() = 0;
    virtual void setLibPath(const std::string &path) = 0;
};

template <int T>
class FilterCache: public BaseFilterCache {
public:
    static FilterCache &getInstance() {
        static FilterCache theFilterCache;
        return theFilterCache;
    }

    virtual void load(int order);
    virtual MWFilter &get(int order);
    virtual const Eigen::MatrixXd &getFilterMatrix(int order);

    virtual const std::string &getLibPath() { return this->libPath; }
    virtual void setLibPath(const std::string &path) { this->libPath = path; }
protected:
    int type;
    std::string libPath; ///< Base path to filter library
private:
    FilterCache();
    virtual ~FilterCache() { }
    FilterCache(FilterCache<T> const &fc) { }
    FilterCache &operator=(FilterCache<T> const &fc) { return *this; }
};

