/** \file TelePrompter.h
 * Collection of assertions and a standard error/warn/info/debug
 * message interface.
 *
 */
#pragma once

#include <cassert>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>

class Timer;

class TelePrompter {
public:
    static void init(int level = 0, bool teletype = false, const char *fil=0);
    static void printSeparator(int level, const char &sep, int newlines = 0);
    static void printHeader(int level, const std::string &str, int newlines = 0);
    static void printFooter(int level, const Timer &t, int newlines = 0);
    static void printDouble(int level, const std::string &name, double d);
    static void printTree(int level, const std::string &name, int n, double t);
    static void setOutputStream(std::ostream &o) { out = &o; }
    static int setPrintLevel(int i) {
        int oldLevel = printLevel;
        printLevel = i;
        return oldLevel;
    }
    static int setPrecision(int i) {
        int oldPrec = precision;
        precision = i;
        *out << std::scientific << std::setprecision(i);
        return oldPrec;
    }
    static int  getPrecision() { return precision; }
    static int  getPrintLevel() { return printLevel; }
    static std::ostream *out;
private:
    static int printLevel;
    static int precision;
};

#define STR_DEBUG(S,X) {std::ostringstream _str;\
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_INFO(S,X) {std::ostringstream _str;\
    _str << "Info: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_WARN(S,X) {std::ostringstream _str;\
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_ERROR(S,X) {std::ostringstream _str;\
    _str << "Error: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}

#define println(level, STR) if (level <= TelePrompter::getPrintLevel()) \
*TelePrompter::out << STR << std::endl;
#define printout(level, STR) if (level <= TelePrompter::getPrintLevel()) \
*TelePrompter::out << STR;

#define GET_PRINT_PRECISION(a) a = TelePrompter::getPrecision();
#define SET_PRINT_PRECISION(a) TelePrompter::setPrecision(a);
#define SET_PRINT_LEVEL(a) TelePrompter::setPrintLevel(a);
#define SET_MESSAGE_STREAM(s) TelePrompter::setOutputStream(s);
#define PRINT_LEVEL TelePrompter::getPrintLevel()

#define MSG_DEBUG(X) { *TelePrompter::out << "Debug: " << \
    __func__ << "(), line " << __LINE__ << "in " << \
    __FILE__ << ": " << X << std::endl;}
#define MSG_INFO(X) { *TelePrompter::out << "Info: " << __FILE__ << ": " <<\
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_NOTE(X) { *TelePrompter::out << "Note: " << __FILE__ << ": " <<\
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_WARN(X) { *TelePrompter::out << "Warning: " \
    << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_ERROR(X) { *TelePrompter::out << "Error: " << \
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_FATAL(X) { *TelePrompter::out << "Error: " << __FILE__ << ": " << \
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl; abort();}

#define MSG_INVALID_ARG(X) { *TelePrompter::out << \
    "Error, invalid argument passed: " << __func__ << "(), line " <<\
    __LINE__ << ": " << X << std::endl;}
#define INVALID_ARG_ABORT { *TelePrompter::out << \
    "Error, invalid argument passed: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}
#define NOT_REACHED_ABORT { *TelePrompter::out << \
    "Error, should not be reached: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}
#define INTERNAL_INCONSISTENCY { *TelePrompter::out << \
    "Internal inconsistency! You have found a bug: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}

#define NEEDS_TESTING {static bool __once = true;\
    if (__once) { __once = false;\
    *TelePrompter::out << "NEEDS TESTING: " << __FILE__ << ", " << \
    __func__ << "(), line " << __LINE__ <<  std::endl;}}

#define ASSERT_FILE(A,B) {if (A == NULL) {\
    *TelePrompter::out << "Error: " << __func__ << "(), line " << __LINE__ << \
    ": No such file, " << B << std::endl; abort();}}

#define NOT_IMPLEMENTED_ABORT {*TelePrompter::out << \
    "Error: Not implemented, " << __FILE__ ", " << __func__ << "(), line " << __LINE__ << \
    std::endl; abort();}

#define NOTE(X) {static bool __once = true;\
    if (__once) { __once = false; *TelePrompter::out << \
    "NOTE: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; }}

#define NEEDS_FIX(X) {static bool __once = true;\
    if (__once) { __once = false; *TelePrompter::out << \
    "NEEDS FIX: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; }}

#define WRONG(X) {*TelePrompter::out << \
    "WRONG: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; abort();}

#define STR_DEBUG(S,X) {std::ostringstream _str;\
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_INFO(S,X) {std::ostringstream _str;\
    _str << "Info: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_WARN(S,X) {std::ostringstream _str;\
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_ERROR(S,X) {std::ostringstream _str;\
    _str << "Error: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}

/* The quiet versions...
 #define SET_PRINT_LEVEL(a)
 #define SET_MESSAGE_STREAM(s)
 #define PRINT_LEVEL

 #define MSG_DEBUG(X)
 #define MSG_INFO(X)
 #define MSG_WARN(X)
 #define MSG_ERROR(X)
 #define MSG_FATAL(X) abort();


 #define MSG_INVALID_ARG
 #define INVALID_ARG_ABORT abort();
 #define NOT_REACHED_ABORT abort();

 #define NEEDS_TESTING

 #define ASSERT_FILE(A,B)
 #define NOT_IMPLEMENTED_ABORT abort();
 */

