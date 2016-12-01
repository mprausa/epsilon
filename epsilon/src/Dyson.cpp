// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/Dyson.cpp
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

#include <Dyson.h>
#include <fstream>
#include <iostream>
#include <boost/functional/hash.hpp>

using namespace std;
using namespace boost;

Dyson::GPL::GPL() {
    _arg = "";
    _indices.clear();
    isHPL = true;
}

Dyson::GPL::GPL(const GPL &other) {
    _arg = other._arg;
    _indices = other._indices;
    isHPL = other.isHPL;
}

void Dyson::GPL::addIndex(const FermatExpression &xj) {
    _indices.push_front(xj);
    string s = xj.str();

    if (s != "-1" && s != "0" && s != "1") {
        isHPL = false;
    }
}

void Dyson::GPL::setArg(string arg) {
    _arg = arg;
}

string Dyson::GPL::arg() const {
    return _arg;
}

deque<FermatExpression> Dyson::GPL::indices() const {
    return _indices;
}

bool Dyson::GPL::operator<(const Dyson::GPL &other) const {
    if (_indices.size() > other._indices.size()) {
        return true;
    } else if (_indices.size() < other._indices.size()) {
        return false;
    }

    for (int n=0; n<_indices.size(); ++n) {
        if (_indices[n] < other._indices[n]) {
            return true;
        } else if (other._indices[n] < _indices[n]) {
            return false;
        }
    }

    return false;
}

bool Dyson::GPL::operator==(const Dyson::GPL &other) const {
    if (_arg != other._arg) return false;
    if (_indices != other._indices) return false;

    return true;
}

string Dyson::GPL::str(pltype_t type, format_t format) const {
    stringstream strm;

    if (_indices.size() == 0) return "1";

    if (!isHPL) type = tGPL;

    strm << (type==tGPL?"GPL":"HPL") << (format==fMma?"[{":"(");
    
    deque<FermatExpression> inds = type == tHPLalt ? indsHPLalt() : _indices;

    for (int n = 0; n<inds.size(); ++n) {
        strm << inds[n].str();
        if (n<inds.size()-1) {
            strm << ",";
        }
    }

    strm << (format==fMma?"},":",") << _arg << (format==fMma?"]":")");

    return strm.str();
}

int Dyson::GPL::sign(pltype_t type) const {
    if (!isHPL) type = tGPL;

    if (type == tGPL) return 1;

    int sgn=1;

    for (auto it = _indices.begin(); it != _indices.end(); ++it) {
        if (it->str() == "1") sgn *= -1;
    }

    return sgn;
}

deque<FermatExpression> Dyson::GPL::indsHPLalt() const {
    if (_indices.empty()) return _indices;

    deque<FermatExpression> inds;
    Fermat *fermat = _indices.begin()->fer();
    FermatExpression zero(fermat,"0");
    FermatExpression one(fermat,"1");
    FermatExpression dir=zero;

    for (auto it = _indices.rbegin(); it != _indices.rend(); ++it) {
        if (it->str() == "0") {
            if (dir.str() != "0") {
                *inds.begin() = *inds.begin() + dir;
            } else {
                inds.push_front(zero);
            }
        } else {
            inds.push_front(*it);
            dir = (it->str() == "-1") ? -one : one;
        }
    }

    return inds;
}

Dyson::Term::Term(const GPL &xgpl, const map<GPL,int> &x0gpl) {
    _xGPL = xgpl;
    _x0GPL = x0gpl;

    calchash();
}
 
const Dyson::GPL &Dyson::Term::xGPL() const {
    return _xGPL;
}

const map<Dyson::GPL,int> &Dyson::Term::x0GPL() const {
    return _x0GPL;
}

int Dyson::Term::sign(pltype_t type) const {
    int sgn = 1;

    sgn *= _xGPL.sign(type);

    for (auto it=_x0GPL.begin(); it != _x0GPL.end(); ++it) {
        if ((it->second % 2) == 0) continue;
        sgn *= it->first.sign(type);
    }

    return sgn;
}

size_t Dyson::Term::hash() const {
    return _hash;
}

void Dyson::Term::calchash() {
    _hash = 0;
    deque<FermatExpression> indices;

    indices = _xGPL.indices();
    
    hash_combine(_hash,_xGPL.arg());
    hash_combine(_hash,indices.size());
    for (auto it = indices.begin(); it != indices.end(); ++it) {
        hash_combine(_hash,it->str());
    }
    hash_combine(_hash,10000);
    hash_combine(_hash,_x0GPL.size());
    for (auto it = _x0GPL.begin(); it != _x0GPL.end(); ++it) {
        hash_combine(_hash,it->second);
        hash_combine(_hash,it->first.arg());
        
        indices = it->first.indices();

        for (auto it2 = indices.begin(); it2 != indices.end(); ++it2) {
            hash_combine(_hash,it2->str());
        }

        hash_combine(_hash,20000);
    }
    hash_combine(_hash,30000);
}

bool Dyson::Term::operator==(const Term &other) const {
    return _hash == other._hash;
}

Dyson::Expression::Expression() {
}

Dyson::Expression::~Expression() {
}

Dyson::Expression &Dyson::Expression::operator+=(const Dyson::Expression &other) {
    for (auto it = other.begin(); it != other.end(); ++it) {
        if ((*this)[it->first].fer()) {
            (*this)[it->first] = (*this)[it->first] + it->second;
        } else {
            (*this)[it->first] = it->second;
        }

        if ((*this)[it->first].str() == "0") erase(it->first);
    }
    return *this;
}

Dyson::Expression Dyson::Expression::operator*(const FermatExpression &factor) const {
    if (factor.str() == "0") return Expression();

    Expression n = *this;

    for (auto it = n.begin(); it != n.end(); ++it) {
        it->second = it->second*factor;
    }

    return n;
}

Dyson::Expression Dyson::Expression::integrate(const FermatExpression &xj) {
    Expression n;

    for (auto it = begin(); it != end(); ++it) {
        GPL xgpl = it->first.xGPL();
        FermatExpression prefactor = it->second;
        
        xgpl.addIndex(xj);

        Term t1(xgpl,it->first.x0GPL());

        if (n[t1].fer()) {
            n[t1] = n[t1]+prefactor;
        } else {
            n[t1] = prefactor;
        }

        if (n[t1].str() == "0") n.erase(t1);

        map<GPL,int> x0gpl = it->first.x0GPL();

        xgpl.setArg("x0");
        x0gpl[xgpl]++;

        GPL g;
        g.setArg("x");

        Term t2(g,x0gpl);

        if (n[t2].fer()) {
            n[t2] = n[t2]-prefactor;
        } else {
            n[t2] = -prefactor;
        }

        if (n[t2].str() == "0") n.erase(t2);
    }

    return n;
}

string Dyson::Expression::str(Dyson::pltype_t type, Dyson::format_t format) {
    stringstream strm;

    if (size() == 0) return "0";
    
    for (auto it = begin(); it != end(); ++it) {
        Term t = it->first;
        FermatExpression prefactor = it->second*it->first.sign(type);
        string xGPL = t.xGPL().str(type,format);
    
        if (it != begin()) strm << "+";
 
        switch(format) {
            case fMma:
                strm << "(" << prefactor.str() << ")";
                break;
            case fForm:
                strm << "rat(" << prefactor.numer().str() << "," << prefactor.denom().str() << ")";
                break;
        }                

        if (xGPL != "1" || !t.x0GPL().empty()) {
            strm << "*";
        }

        if (xGPL != "1") {
            strm << xGPL;
        
            if (t.x0GPL().size()) {
                strm << "*";
            }
        }

        for (auto it2=t.x0GPL().begin(); it2 != t.x0GPL().end(); ++it2) {
            if (it2 != t.x0GPL().begin()) {
                strm << "*";
            }
            strm << it2->first.str(type,format);
            if (it2->second > 1) {
                strm << "^" << it2->second;
            }
        }
    }

    return strm.str();
}

Dyson::ExMatrix::ExMatrix(int dim) : multi_array(extents[dim][dim]) {
}

Dyson::ExMatrix::~ExMatrix() {
}

Dyson::ExMatrix &Dyson::ExMatrix::operator+=(const ExMatrix &other) {
    int rows = shape()[0];
    int cols = shape()[0];
   
    for (int r=0; r<rows; ++r) {
        for (int c=0; c<cols; ++c) {
            (*this)[r][c] += other[r][c];
        }
    }

    return *this;
}

Dyson::ExMatrix Dyson::ExMatrix::integrate(const FermatExpression &xj) {
    int dim = shape()[0];
    ExMatrix n(dim);

    for (int r=0; r<dim; ++r) {
        for (int c=0; c<dim; ++c) {
            n[r][c] = (*this)[r][c].integrate(xj);
        }
    }

    return n;
}      
 
Dyson::ExMatrix Dyson::ExMatrix::lmul(const FermatArray &left) const {
    int dim = shape()[0];
    ExMatrix nn(dim);

    for (int r=0; r<dim; ++r) {
        for (int c=0; c<dim; ++c) {
            for (int n=0; n<dim; ++n) {
                Expression exr = (*this)[n][c];
                if (exr.empty()) continue;                

                FermatExpression exl = left(r+1,n+1);
                if (exl.str() == "0") continue;

                nn[r][c] += exr * exl;
            }
        }
    }

    return nn;
}


string Dyson::ExMatrix::str(pltype_t type) {
    stringstream strm;
    int rows = shape()[0];
    int cols = shape()[1];
    
    strm << "{";

    for (int r=0; r<rows; ++r) {
        strm << "{";
        for (int c=0; c<cols; ++c) {
            strm << (*this)[r][c].str(type,fMma);
            if (c < cols-1) strm << ",";
        }
        strm << "}";
        if (r < cols-1) strm << ",";
    }
    strm << "}";

    return strm.str();
}

Dyson::Dyson(const System &system) {
    map<FermatExpression,FermatArray> fuchs = system.exportFuchs();
    Fermat *fermat = system.fer();
    FermatExpression ep(fermat,"ep");

    dim = system.dimC();

    for (auto it = fuchs.begin(); it != fuchs.end(); ++it) {
        FermatArray M = it->second.subst("ep",1);
        if ((M*ep).str() != it->second.str()) {
            throw invalid_argument("system not in ep-form.");
        }
    
        Mxj.push_back({it->first,M});
    }

    ExMatrix U0(dim);
    Expression one;

    GPL xgpl;
    xgpl.setArg("x");
    
    Term oneterm(xgpl,map<GPL,int>());

    one[oneterm] = FermatExpression(fermat,"1");

    for (int n=0; n<dim; ++n) {
        U0[n][n] = one;
    }

    Un.push_back(U0);
}

void Dyson::expand(int order) {
    if (order < Un.size()) return;
    if (Un.size() < order) expand(order-1);

    cout << "order: " << order << endl;

    ExMatrix U(dim);

    for (auto it = Mxj.begin(); it != Mxj.end(); ++it) {
        U += Un[order-1].lmul(it->second).integrate(it->first);
    }

    Un.push_back(U);
}

void Dyson::write(string filename, pltype_t type, format_t format) {
    switch(format) {
        case fMma:
            writeMma(filename,type);
            break;
        case fForm:
            writeForm(filename,type);
            break;
    }            
}

void Dyson::writeMma(string filename, pltype_t type) {
    ofstream file(filename);

    file << "{" << endl << "  ";

    for (int r=0; r<dim; ++r) {
        file << "{" << endl;
        for (int c=0; c<dim; ++c) {
            int first=-1;
            file << "    SeriesData[ep,0,{";

            for (int n=0; n<Un.size(); ++n) {
                if (!Un[n][r][c].size() && first<0) continue;
                if (first < 0) first = n;

                file << Un[n][r][c].str(type,fMma);
                if (n<Un.size()-1) file << ",";
            }
            if (first < 0) first = Un.size();
            file << "}," << first << "," << Un.size() << ",1]";
            if (c < dim-1) file << ",";
            file << endl;
        }
        file << "  }";
        if (r<dim-1) file << ",";
    }

    file << endl << "}" << endl;

    file.close();
}

void Dyson::writeForm(string filename, pltype_t type) {
    ofstream file(filename);
   
    file << "*--#[ setup :" << endl << endl;

    file << "#define DYSONMAXEP \"" << Un.size()-1 << "\"" << endl;

    file << endl << "*--#] setup :" << endl << endl;

    file << "*--#[ identities :" << endl;

    for (int n=0; n<Un.size(); ++n) {
        file << endl << "* order ep^" << n << endl;
        for (int r=0; r<dim; ++r) {
            for (int c=0; c<dim; ++c) {
                file << "id dyson(" << r+1 << "," << c+1 << "," << n << ",x?,x0?) = " << Un[n][r][c].str(type,fForm) << ";" << endl;
            }
        }
    }
    
    
    file << endl << "*--#] identities :" << endl;

    file.close();
}


