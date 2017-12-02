/*
 *  \date Jul 24, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include <fstream>
#include "TelePrompter.h"
#include "Timer.h"
#include "parallel.h"

using namespace std;

static ofstream tp_outfile;
int TelePrompter::printLevel = 0;
int TelePrompter::precision = 12;
std::ostream *TelePrompter::out = &std::cout;

void TelePrompter::printSeparator(int level, const char &sep, int newlines) {
    int N = 60;
    for (int i = 0; i < N; i++) {
        printout(level, sep);
    }
    for (int i = 0; i <= newlines; i++) {
        printout(level, endl);
    }
}
void TelePrompter::printHeader(int level, const string &str, int newlines) {
    int N = 60;
    int len = str.size();
    printSeparator(level, '=', 0);
    int spaces = (N - len)/2;
    for (int i = 0; i < spaces; i++) {
        printout(level, " ");
    }
    println(level, str);
    printSeparator(level, '-', newlines);
}

void TelePrompter::printFooter(int level, const Timer &t, int newlines) {
    printSeparator(level, '-');
    println(level, "                 Wall time: " << setw(11) << t << " sec ");
    printSeparator(level, '=', newlines);
}

void TelePrompter::printDouble(int level, const std::string &name, double d) {
    char cName[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < name.size()) {
            cName[i+1] = name[i];
        }
    }
    int oldPrec = TelePrompter::setPrecision(5);
    println(level, cName << setw(29) << d);
    TelePrompter::setPrecision(oldPrec);
}

void TelePrompter::printTree(int level, const std::string &name, int n, double t) {
    char cName[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < name.size()) {
            cName[i+1] = name[i];
        }
    }
    int oldPrec = TelePrompter::setPrecision(5);
    println(level, cName << setw(12) << n << setw(17) << t);
    TelePrompter::setPrecision(oldPrec);
}

void TelePrompter::init(int printLevel, bool teletype, const char *fil) {
    SET_PRINT_PRECISION(15);
    cout << scientific << setprecision(14);

    int rank=mpiOrbRank;
    int world_size = 1;
    SET_PRINT_LEVEL(printLevel);
    if (teletype and fil != 0) {
        stringstream fname;
        if (rank == 0) {
            if (printLevel < 0) {
                if (world_size > 1) {
                    fname << fil << "-" << rank << ".out";
                } else {
                    fname << fil << ".out";
                }
            }
        } else {
            fname << fil << "-" << rank << ".out";
        }

        if (printLevel < 0) {
            SET_PRINT_LEVEL(-printLevel);
        }
        if (not fname.str().empty()) {
            tp_outfile.open(fname.str().c_str());
            SET_MESSAGE_STREAM(tp_outfile);
        }
    } else {
        if (rank > 0) {
            SET_PRINT_LEVEL(-10);
        }
    }
}
