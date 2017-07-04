// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/export.cpp
 * 
 *  Copyright (C) 2016, 2017 Mario Prausa 
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

#include <export.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;
using namespace GiNaC;

static string convert(const ex &e) {
    stringstream strm;
    strm << e;
    string s = strm.str();
    string s1;

    for (auto it=s.begin(); it != s.end(); ++it) {
        switch(*it) {
            case '[':
                s1 += '{';
                break;
            case ']':
                s1 += '}';
                break;
            case 'I':
                s1 += 'i';
                break;
            case '(':
                if (s1.size() >= 4 && s1.substr(s1.size()-4) == "sqrt") {
                    s1.erase(s1.size()-4);
                    ++it;
                    bool imag = false;                       
                    string num="";
                    if (*it == '-') {
                        imag = true;
                        ++it;
                    }
                    for (;it != s.end() && *it != ')'; ++it) {
                        if (!isdigit(*it)) {
                            throw invalid_argument("parse error.");
                        }                            
                        num += *it;
                    }
                    if (*it != ')') {
                        throw invalid_argument("parse error.");
                    }
                    s1 += (imag?"isqrt":"sqrt") + num;
                } else {
                    s1 += *it;
                }
                break;
            default:
                s1 += *it;
        }
    }        
    return s1;
}    

void exprt(map<apart_sing,matrix> &A, map<int,matrix> &B) {
    int width=0;

    for (auto &a : A) {
        stringstream strm;
        strm << "A[" << convert(a.first.xj) << "," << a.first.k << "]:";
        int len = strm.str().size();
        if (len>width) width=len;
    }
    for (auto &b : B) {
        stringstream strm;
        strm << "B[" << convert(b.first) << "]:";
        int len = strm.str().size();
        if (len>width) width=len;
    }

    for (auto &a : A) {
        stringstream strm;
        strm << "A[" << convert(a.first.xj) << "," << a.first.k << "]:";

        cout << setw(width) << left << strm.str() << "\t" << convert(a.second) << endl;
    }

    for (auto &b : B) {
        stringstream strm;
        strm << "B[" << b.first << "]:";
        cout << setw(width) << left << strm.str() << "\t" << convert(b.second) << endl;
    }        
}

