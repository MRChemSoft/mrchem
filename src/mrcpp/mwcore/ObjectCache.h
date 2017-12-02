#pragma once

#include <string>
#include <vector>
#include "TelePrompter.h"
#include "parallel.h"
#include "macros.h"

#define getObjectCache(T,X) \
    ObjectCache<T> &X = ObjectCache<T>::getInstance();

#ifdef HAVE_OPENMP
#define SET_CACHE_LOCK() omp_set_lock(&this->cache_lock)
#define UNSET_CACHE_LOCK() omp_unset_lock(&this->cache_lock)
#define TEST_CACHE_LOCK() omp_test_lock(&this->cache_lock)
#else
#define SET_CACHE_LOCK()
#define UNSET_CACHE_LOCK()
#define TEST_CACHE_LOCK()
#endif

template<class T>
class ObjectCache {
public:
    static ObjectCache<T> &getInstance() {
        static ObjectCache<T> theObjectCache;
        return theObjectCache;
    }

    virtual void clear() {
        for (unsigned int i = 0; i < this->objs.size(); i++) {
            if (this->objs[i] != 0) {
                unload(i);
            }
        }
    }

    virtual void load(int id) {
        MSG_INFO("This routine does nothing in this class.");
    }

    void load(int id, T *new_o, int memory) {
        if (id >= this->highWaterMark) {
            for (int i = 0; i < id - this->highWaterMark + 1; i++) {
                this->objs.push_back(0);
                this->mem.push_back(0);
            }
            this->highWaterMark = id;
        }
        if (this->objs[id] != 0) {
            return;
        }
        this->mem[id] = memory;
        this->memLoaded += memory;
        this->objs[id] = new_o;

    }
    virtual void unload(int id) {
        if (id < 0 or id > this->highWaterMark) {
            MSG_ERROR("Id out of bounds:" << id);
        }
        if (this->objs[id] == 0) {
            MSG_WARN ("Object not loaded.");
            return;
        }
        this->memLoaded -= this->mem[id];
        this->mem[id] = 0;
        delete this->objs[id];
        this->objs[id] = 0;
    }

    virtual T &get(int id) {
        if (id < 0) {
            MSG_ERROR("Id out of bounds:" << id);
        }
        if (this->objs[id] == 0) {
            MSG_ERROR("Object not loaded!");
        }
        return *(this->objs[id]);
    }
    bool hasId(int id) {
        if (id > this->highWaterMark)
            return false;
        if (this->objs[id] == 0)
            return false;
        return true;
    }

    int getNObjs() { return this->objs.size(); }
    int getMem() { return this->memLoaded; }
    int getMem(int id) { return this->mem[id]; }

protected:
    ObjectCache() {
        this->highWaterMark = 0;
        this->memLoaded = 0;
        this->objs.push_back(0);
        this->mem.push_back(0);
#ifdef HAVE_OPENMP
        omp_init_lock(&cache_lock);
#endif
    }

    virtual ~ObjectCache() {
        SET_CACHE_LOCK();
        clear();
        UNSET_CACHE_LOCK();
#ifdef HAVE_OPENMP
        omp_destroy_lock(&cache_lock);
#endif
    }
    ObjectCache(ObjectCache<T> const &oc) { }
    ObjectCache<T> &operator=(ObjectCache<T> const &oc) { return *this; }
#ifdef HAVE_OPENMP
    omp_lock_t cache_lock;
#endif
private:
    int highWaterMark;
    int memLoaded; ///< memory occupied by loaded objects
    std::vector<T *> objs; ///< objects store
    std::vector<int> mem; ///< mem per object
};

