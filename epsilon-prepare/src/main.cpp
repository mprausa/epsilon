// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/main.cpp
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
#include <fstream>
#include <prepare.h>
#include <export.h>
using namespace std;
using namespace GiNaC;

int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " matrix" << endl;
        return 1;
    }

    string str="",line;
    ifstream file(argv[1]);

    if (!file.is_open()) {
        throw invalid_argument("unable to open file");
    }

    while (getline(file,line)) {
        str += line;
    }

    file.close();

    symbol x("x"),ep("ep");
    symtab table;
    table["x"] = x;
    table["ep"] = ep;
    parser reader(table);

    ex e = reader(str).normal();
    matrix mat = ex_to<matrix>(lst_to_matrix(ex_to<lst>(e)));
    map<apart_sing,matrix> A;
    map<int,matrix> B;

    try {
        prepare(mat,x,A,B);
        exprt(A,B);
    } catch (const string &s) {
        cerr << s << endl;
        return 1;
    }

    return 0;
}


