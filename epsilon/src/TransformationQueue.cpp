// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/TransformationQueue.cpp
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

#include <TransformationQueue.h>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <System.h>
using namespace std;

TransformationQueue::TransformationQueue(const TransformationQueue &other) {
    before = after = 0;
    replaying = false;
    fermat = other.fermat;
    queue = other.queue;
    _filename = "";
}

TransformationQueue::TransformationQueue(Fermat *fermat) {
	before = after = 0;
    replaying = false;
    this->fermat = fermat;
    _filename = "";
}

TransformationQueue::~TransformationQueue() {
	if (file.is_open()) file.close();
}

void TransformationQueue::setpadding(int before, int after) {
    this->before = before;
    this->after = after;
}

void TransformationQueue::setfile(string _filename, bool append) {
    if (_filename == "") return;
    
    if (file.is_open()) file.close();

    this->_filename = _filename;

    file.open(_filename,append?ios::app:ios::out);	

    if (!file.is_open()) {
        throw invalid_argument("unable to open file.");
    }
}

string TransformationQueue::filename() {
    return _filename;
}

void TransformationQueue::load(string _filename) {
    ifstream file(_filename);
	string str;
    
    Fermat *fermat = infinity.fer();
    FermatExpression zero(fermat,"0");

    if (!file.is_open()) {
        throw invalid_argument("unable to open file.");
    }
	
    while (getline(file,str)) {
        transformation_t trans;

        str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
        if (str == "") continue;

        size_t colon = str.find(':');
        if (colon == string::npos) {
            throw invalid_argument("parse error");
        }
        
        string str0(str,0,colon);
        str = str.erase(0,colon+1);

        if (str0.substr(0,2) == "B(" && str0.back() == ')') {
            size_t comma = str0.find(',');
            if (comma == string::npos) {
                throw invalid_argument("parse error");
            }
            
            trans.type = transformation_t::Balance;
            string xstr;

            xstr = str0.substr(2,comma-2);
            if (xstr == "inf") {
                trans.x1 = infinity;
            } else {
                trans.x1 = FermatExpression(fermat,xstr);
            }

            xstr = str0.substr(comma+1,str0.size()-comma-2);
            if (xstr == "inf") {
                trans.x2 = infinity;
            } else {
                trans.x2 = FermatExpression(fermat,xstr);
            }

            trans.k = 0;

            trans.T = FermatArray(fermat,str);
        } else if (str0 == "T") {
            trans.type = transformation_t::Transformation;
            trans.x1 = trans.x2 = zero;
            trans.k = 0;

            trans.T = FermatArray(fermat,str);
        } else if (str0.substr(0,2) == "L(" && str0.back() == ')') {
            size_t comma = str0.find(',');
            if (comma == string::npos) {
                throw invalid_argument("parse error");
            }

            trans.type = transformation_t::LeftTrans;
            string xstr(str0,2,comma-2);

            if (xstr == "inf") {
                trans.x1 = infinity;
            } else {
                trans.x1 = FermatExpression(fermat,xstr);
            }
            
            trans.k = stoi(str0.substr(comma+1,str0.size()-comma-2));

            trans.x2 = zero;

            trans.T = FermatArray(fermat,str);
        } else {
            throw invalid_argument("parse error.");
        }

        queue.push_back(trans);
    }
    file.close();
}

void TransformationQueue::replay(System &system) {
    replaying = true;

    if (before != 0 || after != 0) {
        throw invalid_argument("Replay not possible. Full system must be activated.");
    }

    for (auto it = queue.begin(); it != queue.end(); ++it) {
        switch (it->type) {
            case transformation_t::Balance:
                cout << "balance [" << pstr(it->x1) << "," << pstr(it->x2) << "]" << endl;
                system.balance(it->T,it->x1,it->x2);
                break;
            case transformation_t::Transformation:
                cout << "transformation" << endl;
                system.transform(it->T);
                break;
            case transformation_t::LeftTrans:
                cout << "left transformation [" << pstr(it->x1) << "," << it->k << "]" << endl;

                system.lefttransformFull(it->T,it->x1,it->k);
                break;
        }
    }

    replaying = false;
}

void TransformationQueue::exporttrans(string filename) {
    if (!queue.size()) {
        throw invalid_argument("transformation queue is empty.");
    }

    Fermat *fermat = queue.begin()->T.fer();
    FermatExpression one(fermat,"1");
    int size = queue.begin()->T.rows();
    FermatArray T(fermat,size,size);
    T.assign("[1] + 0");
    
    fermat->addSymbol("x");

    for (auto it = queue.begin(); it != queue.end(); ++it) {
        stringstream strm;
        FermatArray xT(fermat);

        switch(it->type) {
            case transformation_t::Balance:
                cout << "balance (" << pstr(it->x1) << "," << pstr(it->x2) << ")" << endl;
                if (it->x1 == infinity) {
                    strm << "[" << it->T.name() << "]*(x-(" << (it->x2+one).str() << ")) + [1]";
                } else if (it->x2 == infinity) {
                    strm << "[" << it->T.name() << "]*((" << (it->x1+one).str() << ")-x)/(x-(" << it->x1.str() << ")) + [1]"; 
                } else {
                    strm << "[" << it->T.name() << "]*(" << (it->x1-it->x2).str() << ")/(x-(" << it->x1.str() << ")) + [1]";
                }
                xT.assign(strm.str());
                break;
            case transformation_t::Transformation:
                cout << "transformation" << endl;
                xT = it->T;
                break;
            case transformation_t::LeftTrans:
                cout << "left transformation(" << pstr(it->x1) << "," << it->k << ")" << endl;
                if (it->x1 == infinity) {
                    strm << "[" << it->T.name() << "]*x^(" << it->k << ") + [1]";
                } else {
                    strm << "[" << it->T.name() << "]/(x-(" << it->x1.str() << "))^(" << it->k << ") + [1]";
                }
                xT.assign(strm.str());
                break;
        }

        T = T*xT;
    }

    ofstream file(filename);
    file << T.str() << endl;
    file.close();

    fermat->dropSymbol("x");
}

void TransformationQueue::balance(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2) {
    if (replaying) return;

    Fermat *fermat = P.fer();
    FermatArray xP(fermat,P.rows()+before+after,P.cols()+before+after);
    stringstream strm;

    xP.assign("0");
    
    strm << "[" << xP.name() << "[" << before+1 << "~" << before+P.rows() << "," << before+1 << "~" << before+P.cols() << "]] := [" << P.name() << "]";
    (*fermat)(strm.str());

    if (file.is_open()) {
        file << "B(" << pstr(x1) << "," << pstr(x2) << "):  \t" << xP.str() << endl; 
    }

    transformation_t trans;

    trans.type = transformation_t::Balance;
    trans.x1 = x1;
    trans.x2 = x2;
    trans.k = 0;
    trans.T = xP;

    queue.push_back(trans);
}

void TransformationQueue::transform(const FermatArray &T) {
    if (replaying) return;

    Fermat *fermat = infinity.fer();
    FermatExpression zero(fermat,"0");
    FermatArray xT(fermat,T.rows()+before+after,T.cols()+before+after);
    stringstream strm;

    xT.assign("[1] + 0");
    
    strm << "[" << xT.name() << "[" << before+1 << "~" << before+T.rows() << "," << before+1 << "~" << before+T.cols() << "]] := [" << T.name() << "]";
    (*fermat)(strm.str());

    if (file.is_open()) {
        file << "T:        \t" << xT.str() << endl; 
    }

    transformation_t trans;

    trans.type = transformation_t::Transformation;
    trans.x1 = zero;
    trans.x2 = zero;
    trans.k = 0;
    trans.T = xT;

    queue.push_back(trans);
}
    
void TransformationQueue::lefttransform(const FermatArray &G, const FermatExpression &x1, int k) {
    if (replaying) return;

    Fermat *fermat = G.fer();
    FermatExpression zero(fermat,"0");
    FermatArray xG(fermat,G.rows()+before+after,G.rows()+before+after);
    stringstream strm;

    xG.assign("0");

    strm << "[" << xG.name() << "[" << before+1 << "~" << before+G.rows() << ",1~" << G.cols() << "]] := [" << G.name() << "]";
    (*fermat)(strm.str());

    if (file.is_open()) {
        file << "L(" << pstr(x1) << "," << k << "):  \t" << xG.str() << endl;
    }

    transformation_t trans;

    trans.type = transformation_t::LeftTrans;
    trans.x1 = x1;
    trans.x2 = zero;
    trans.k = k;
    trans.T = xG;

    queue.push_back(trans);
}

string TransformationQueue::pstr(const FermatExpression &x) const {
    if (x == infinity) {
        return "inf";
    } else {
        return x.str();
    }
}

