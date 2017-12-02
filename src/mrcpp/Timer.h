#pragma once

#include <chrono>

#include "TelePrompter.h"

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timeT;

class Timer {
public:
    Timer(bool start_timer = true) : running(false), time_used(0.0) {
        if (start_timer) {
            start();
        }
    }
    Timer& operator=(const Timer &timer) {
        if (this != &timer) {
            this->time_used = timer.time_used;
            this->clock_start = timer.clock_start;
            this->running = timer.running;
        }
        return *this;
    }
    void start() {
        this->clock_start = now();
        this->time_used = 0.0;
        this->running = true;
    }
    void resume() {
        if (this->running) MSG_WARN("Timer already running");
        this->clock_start = now();
        this->running = true;
    }
    void stop() {
        if (not this->running) MSG_WARN("Timer not running");
        this->time_used += diffTime(now(), this->clock_start);
        this->running = false;
    }

    double getWallTime() const {
        if (this->running) MSG_WARN("Timer still running");
        return this->time_used;
    }
    double getUserTime() const {
        if (this->running) MSG_WARN("Timer still running");
        return this->time_used;
    }
    double getSystemTime() const {
        if (this->running) MSG_WARN("Timer still running");
        return this->time_used;
    }

    friend std::ostream& operator<<(std::ostream &o, const Timer &timer) {
        int old_prec;
        GET_PRINT_PRECISION(old_prec);
        SET_PRINT_PRECISION(5);
	//uncommented until properly defined:        o << "    user   " << std::setw(6) <<  timer.getUserTime();
        o <<  timer.getWallTime();
        SET_PRINT_PRECISION(old_prec);
        return o;
	}
private:
    bool running;
    double time_used;
    timeT clock_start;

    timeT now(){
        return std::chrono::high_resolution_clock::now();
    }

    double diffTime(timeT t2, timeT t1) {
        std::chrono::duration<double> diff = t2-t1;
        return diff.count();
    }
};
