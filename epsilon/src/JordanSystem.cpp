// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/JordanSystem.cpp
 * 
 *  Copyright (C) 2016 Mario Prausa 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <sstream>
#include <set>
#include <ctime>
#include <iomanip>
#include <JordanSystem.h>
#include <FermatArray.h>
#include <Eigenvalues.h>
using namespace std;


static void gramSchmidt(FermatArray &U) {
    Fermat *fermat = U.fer();
    (*fermat)("GramSchm(["+U.name()+"])");
}

static int kern(const FermatArray &mat, FermatArray &U) {
    struct colprops {
        int col;
        int nonzero;
        int deg;

        bool operator<(const struct colprops &other) const {
            return (nonzero < other.nonzero || (nonzero == other.nonzero && deg < other.deg) || (nonzero == other.nonzero && deg == other.deg && col < other.col));
        }
    };

    set<colprops> cols;
    stringstream strm;
    int rk;

    Fermat *fermat = mat.fer();

    FermatArray M(mat);
    FermatArray A(fermat,M.rows(),M.rows());
    FermatArray B(fermat,M.cols(),M.cols());

    rk = M.colReduce(A,B);
    U = FermatArray(B,1,B.rows(),rk+1,B.cols());
    
    // optimization
    for (int n=1; n<=U.cols(); ++n) {
        colprops cp;

        cp.col = n;
        cp.nonzero = 0;
        cp.deg = 0;

        for (int j=1; j<=U.rows(); ++j) {
            if (U(j,n).str() != "0") {
                cp.nonzero++;

                int ndeg = U(j,n).numer().deg("ep");
                int ddeg = U(j,n).denom().deg("ep");

                cp.deg = max(cp.deg,ndeg+ddeg);
            }
        }

        cols.insert(cp);
    }

    FermatArray swap(U.fer(),U.cols(),U.cols());
    int c=1;

    swap.assign("0");

    for (auto it=cols.begin(); it != cols.end(); ++it) {
        swap.set(it->col,c++,FermatExpression(U.fer(),"1"));
    }

    U = U*swap;

    return M.rows()-rk;
}

static void jordanDecomposition(const FermatArray &array, eigen_t ev, JordanSystem &system) {
    FermatArray mat(array.fer(),array.rows(),array.cols());
    stringstream strm;
    std::vector<int> as;
    std::vector<int> blocks;
    FermatArray U(array.fer(),array.cols(),array.cols());
    FermatArray id(array.fer(),array.cols(),array.cols());
    FermatArray v(array.fer(),array.cols(),1);
    std::vector<FermatArray> Us;

    id.assign("[1] + 0");

    strm << "(" << ev.u << "+(" << ev.v << ")*ep)*[1] + 0";
    
    mat.assign(strm.str());
   
    mat = array - mat;

    FermatArray mat2(mat);

    as.push_back(0);

    U.assign("0");
    Us.push_back(U);

    for (int s=1;; ++s) {
        int a = kern(mat2,U);

        if (s > 1) {
            int b = 2*(*as.rbegin()) - *(++as.rbegin()) - a;

            blocks.push_back(b);
        }

        if (as.size() && a == *as.rbegin()) break;

        as.push_back(a);
        Us.push_back(U);

        mat2 = mat2*mat;
    }
    
    for (int s=blocks.size(); s>0; --s) {
        int b = blocks[s-1];
        FermatArray B = Us[s-1].transpose();
        int pos;

        for (auto it=system.begin(); it != system.end(); ++it) {
            if (it->ev.u != ev.u || it->ev.v != ev.v) continue;

            v = it->rootvectors[s-1];
            
            B = B.concatenate(v.transpose());
        }

        pos = B.rows()+1;

        B = B.concatenate(Us[s].transpose()).transpose();

        gramSchmidt(B);

        B = FermatArray(B,1,B.rows(),pos,B.cols());

        for (int n=0,c=1; n<b; ++n) {
            JordanBlock blck;
            do {
                v = FermatArray(B,1,B.rows(),c,c);
                ++c;
            } while (v.isZero() && c<=B.cols());

            if (v.isZero()) {
                throw invalid_argument("jordan decomposition failed.");
            }

            blck.ev = ev;
            blck.rootvectors.push_front(v);

            for (int j=s-1; j>=1; --j) {
                v = mat*v;
            
                blck.rootvectors.push_front(v);
            }

            system.insert(blck);
        }
    }
}

void jordanSystem(const FermatArray &A, const eigenvalues_t &evs, JordanSystem &system) {
    for (auto &ev : evs) {
        jordanDecomposition(A,ev.first,system);  
    }
}

void Eigenvectors(const FermatArray &A, eigen_t ev, vector<FermatArray> &vectors) {
    FermatArray mat(A.fer());
    FermatArray U;
    stringstream strm;

    strm << "[" << A.name() << "] - ((" << ev.u << ")+(" << ev.v << ")*ep)*[1]";
    mat.assign(strm.str());

    kern(mat, U);
    gramSchmidt(U);

    vectors.clear();
    for (int c=1; c<=U.cols(); ++c) {
        FermatArray v(U,1,U.rows(),c,c);
        if (!v.isZero()) vectors.push_back(v);
    }
}

