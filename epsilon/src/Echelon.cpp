// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/Echelon.cpp
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

#include <Echelon.h>
#include <climits>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;

EchelonBase::Row::Row(Fermat *fermat) {
    this->fermat = fermat;
}

EchelonBase::Row EchelonBase::Row::operator-(const EchelonBase::Row &other) const {
    auto it2 = other.data.begin();
    EchelonBase::Row res(fermat);

    for (auto it1 = data.begin(); it1 != data.end(); ++it1) {
        for (;it2 != other.data.end() && it2->first < it1->first; ++it2) {
            res.set(it2->first,-it2->second);
        }

        if (it2 != other.data.end() && it2->first == it1->first) {
            res.set(it2->first,it1->second-it2->second);
            ++it2;
        } else {
            res.set(it1->first,it1->second);
        }
    }        

    for (; it2 != other.data.end(); ++it2) {
        res.set(it2->first,-it2->second);
    }        

    return res;
}

EchelonBase::Row EchelonBase::Row::operator* (const FermatExpression &f) const {
    EchelonBase::Row res(fermat);

    if (f.str() == "0") return res;
    
    for (auto &e : data) {
        res.set(e.first,e.second*f);
    }
    return res;
}

FermatExpression EchelonBase::Row::operator[] (int c) const {
    if (empty()) return FermatExpression(fermat,"0");

    if (c < data.front().first) return FermatExpression(fermat,"0");
    if (c > data.back().first) return FermatExpression(fermat,"0");

    if (c == data.front().first) return data.front().second;
    if (c == data.back().first) return data.back().second;

    size_t p0 = 0;
    size_t p1 = data.size()-1;

    while (p1 > p0+1) {
        size_t next = (p0+p1)/2;
        int col = data[next].first;

        if (col == c) return data[next].second;

        if (col < c) {
            p0 = next;
        } else {
            p1 = next;
        }
    }

    return FermatExpression(fermat,"0");
}
 
bool EchelonBase::Row::empty() const {
    return data.empty();
}

size_t EchelonBase::Row::size() const {
    return data.size();
}

int EchelonBase::Row::col1() const {
    if (empty()) {
        throw invalid_argument("row is empty.");
    }

    return data.front().first;
}

void EchelonBase::Row::set(int n, const FermatExpression &ex) {
    if (!data.empty() && n <= data.back().first) {
        throw invalid_argument("wrong order.");
    }

    if (ex.str() == "0") return;

    data.push_back({n,ex});
}

void EchelonBase::Row::clear() {
    data.clear();
}

void EchelonBase::Row::normalize() {
    if (empty()) return;

    FermatExpression norm = data.front().second;

    for (auto &e : data) {
        e.second = e.second/norm;
    }
}

vector<pair<int,FermatExpression>>::const_iterator EchelonBase::Row::begin() const {
    return data.begin();
}

vector<pair<int,FermatExpression>>::const_iterator EchelonBase::Row::end() const {
    return data.end();
}

EchelonBase::Iterator::Iterator(Fermat *fermat, std::map<int,std::vector<Row>> *rows, const map<int,vector<Row>>::iterator &it1, const vector<Row>::iterator &it2) : row(fermat) {
    this->rows = rows;
	this->it1 = it1;
	this->it2 = it2;
    this->array = NULL;
}

EchelonBase::Iterator::Iterator(FermatArray *array, int r) : row(array->fer()) {
    this->rows = NULL;
    this->array = array;
    rownum = r;
    valid=false;
}

EchelonBase::Iterator::Iterator(const Iterator &it) : row(it.row) {
    rows = it.rows;
    it1 = it.it1;
    it2 = it.it2;

    array = it.array;
    rownum = it.rownum;
    valid = false;
}

EchelonBase::Iterator &EchelonBase::Iterator::operator++() {
    if (rows) {
        ++it2;
        if (it2 == it1->second.end()) {
            ++it1;
            if (it1 == rows->end()) return *this;
            it2 = it1->second.begin();
        }
    } else {
        ++rownum;
        valid = false;
    }
    return *this;
}    

EchelonBase::Iterator EchelonBase::Iterator::operator++(int) {
    Iterator tmp(*this);
    operator++();
    return tmp;
}    

bool EchelonBase::Iterator::operator==(const EchelonBase::Iterator &other) const {
    if (rows) {
        if (it1 != other.it1) return false;
        if (it1 == rows->end()) return true;
        return it2 == other.it2;
    } else {
        return (array == other.array) && (rownum == other.rownum);
    }
}

bool EchelonBase::Iterator::operator!=(const EchelonBase::Iterator &other) const {
    return !(*this == other);
}

EchelonBase::Row &EchelonBase::Iterator::operator*() {
    if (rows) {
        return *it2;
    } else {
        if (!valid) {
            row.clear();
            for (int c=1; c<=array->cols(); ++c) {  //TODO: this could be optimized with sparse access loops
                row.set(c,(*array)(rownum,c));
            } 
            valid = true;
        }
        return row;
    }        
}

Echelon::Echelon(Fermat *fermat) {
    this->fermat = fermat;
}

void Echelon::set(const Row &row) {
    if (row.empty()) return;
    rows[row.col1()].push_back(row);
}
        
void Echelon::set(const std::map<int,FermatExpression> &m) {
    Row row(fermat);
    for (auto &e : m) {
        row.set(e.first,e.second);
    }
    set(row);
}

int Echelon::run() {
    size_t rnum=0;
    std::set<int> cols;

    for (auto &r0 : rows) {
        rnum += r0.second.size();
        for (auto &r : r0.second) {
            for (auto &e : r) {
                cols.insert(e.first);
            }
        }
    }

    int maxrk = min(rnum,cols.size());
    int cnt=0;
    int prog=50;

    cout << "forward elimination:  ";

    for (auto it = rows.begin(); it != rows.end(); ++it) {
        int newprog = 50 - (++cnt)*50/maxrk;

        if (newprog < prog) {
            for (; prog>newprog; --prog) {
                if ((prog % 5) == 0) {
                    cout << prog/5;
                } else {
                    cout << ".";
                }
            }
            cout << flush;
            prog = newprog;
        }

        normalize(it->second);

        int pivotn = findPivot(it->second);
        Row pivot = it->second[pivotn];

        for (int n=0; n<it->second.size(); ++n) {
            if (n == pivotn) continue;
            set(it->second[n] - pivot);
        }
        it->second.clear();
        set(pivot);
    }

    for (; prog>=0; --prog) {
        if ((prog % 5) == 0) {
            cout << prog/5;
        } else {
            cout << ".";
        }
    }
    cout << endl;

    int rk = rows.size();
    cnt = 0;
    prog=50;
    
    cout << "back substitution:    ";

    for (auto it = rows.rbegin(); it != rows.rend(); ++it) {
        int newprog = 50 - (++cnt)*50/rk;

        if (newprog < prog) {
            for (; prog>newprog; --prog) {
                if ((prog % 5) == 0) {
                    cout << prog/5;
                } else {
                    cout << ".";
                }
            }
            cout << flush;
            prog = newprog;
        }

        for (auto it2 = next(it); it2 != rows.rend(); ++it2) {
            for (auto &r : it2->second) {
                r = r - it->second.front()*r[it->first];
            }
        }
    }

    for (; prog>=0; --prog) {
        if ((prog % 5) == 0) {
            cout << prog/5;
        } else {
            cout << ".";
        }
    }
    cout << endl;

    return rk;
}

int Echelon::findPivot(const vector<Row> &rows) const {
    size_t min=SIZE_MAX;
    int pivot=-1;
   
    if (rows.size() == 1) return 0;
 
    for (int n=0; n<rows.size(); ++n) {
        size_t s = 0;
        for (auto &e : rows[n]) {
            s += e.second.str().size();
        }

        if (s < min) {
            min = s;
            pivot = n;
        }
    }

    return pivot; 
}

void Echelon::normalize(vector<Row> &rows) {
    for (auto &r : rows) {
        r.normalize();
    }
}

EchelonBase::Iterator Echelon::begin() {
    if (rows.empty()) return end();
    return Iterator(fermat,&rows,rows.begin(),rows.begin()->second.begin());
}

EchelonBase::Iterator Echelon::end() {
    std::vector<Row>::iterator dummy;
    return Iterator(fermat,&rows,rows.end(),dummy);
}

void Echelon::print(ostream &os) const {
    int cnt=1;
    os << "[";

    for (auto &r0 : rows) {
        for (auto &r : r0.second) {
            if (cnt > 1) os << ",";
            os << "[" << (cnt++);
            for (auto &e : r) {
                os << ",[" << e.first << "," << e.second.str() << "]";
            }
            os << "]";
        }
    }
    os << "]";
}

EchelonFermat::EchelonFermat(Fermat *fermat, int rows, int cols) {
    this->fermat = fermat;

    array = FermatArray(fermat,rows,cols,true);
    pos = 1;
}

void EchelonFermat::set(const Row &row) {
    if (row.empty()) return;
    for (auto &e : row) {
        array.set(pos,e.first,e.second);
    }
    ++pos;
}

void EchelonFermat::set(const map<int,FermatExpression> &m) {
    Row row(fermat);
    for (auto &e : m) {
        row.set(e.first,e.second);
    }
    set(row);
}

int EchelonFermat::run() {
    int rk = array.rowEchelon();
    pos = rk+1;
    return rk;
}

EchelonBase::Iterator EchelonFermat::begin() {
    return Iterator(&array,1);
}

EchelonBase::Iterator EchelonFermat::end() {
    return Iterator(&array,pos);
}

void EchelonFermat::print(ostream &os) const {
    os << array.sstr();
}

ostream &operator<<(ostream &os, const EchelonBase &e) {
    e.print(os);
    return os;
}

