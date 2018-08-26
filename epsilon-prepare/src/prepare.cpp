// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/prepare.cpp
 *
 *  Copyright (C) 2016, 2018 Mario Prausa
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

#include <prepare.h>
#include <sstream>
using namespace std;
using namespace GiNaC;

static void check(const matrix &mat, const symbol &x, const map<apart_sing,matrix> &A, const map<int,matrix> &B) {
    ex check = 0;

    for (auto &a : A) {
        check = (check + a.second/pow(x-a.first.xj,a.first.k+1)).evalm();
    }

    for (auto &b : B) {
        check = (check + pow(x,b.first)*b.second).evalm();
    }

    check = check.normal();

    for (int n=0; n<check.nops(); ++n) {
        ex numden = check.op(n).numer_denom();
        ex num=numden.op(0);
        ex den=numden.op(1);

        check.let_op(n) = modout(num)/modout(den);
    } 

    check = (check - mat).evalm().normal();

    if (!check.is_zero_matrix()) {
        throw runtime_error("check failed.");
    }            
}

void prepare(const matrix &mat, const symbol &x, map<apart_sing,matrix> &A, map<int,matrix> &B) {
    int rows=mat.rows();
    int cols=mat.cols();

    for (int r=0; r<rows; ++r) {
        for (int c=0; c<cols; ++c) {
            map<apart_sing,ex> coeffsA;
            map<int,ex> coeffsB;

            try {
                apart(mat(r,c),x,coeffsA,coeffsB);
            } catch(const string &s) {
                stringstream strm;
                strm << s << " @ (" << r+1 << "," << c+1 << ")";
                throw string(strm.str());
            }
            
            for (auto &e : coeffsA) {
                if (!A.count(e.first)) A[e.first] = matrix(rows,cols);
                A[e.first](r,c) = e.second;
            }

            for (auto &e : coeffsB) {
                if (!B.count(e.first)) B[e.first] = matrix(rows,cols);
                B[e.first](r,c) = e.second;
            }                
        }
    }        

    check(mat,x,A,B);
}

