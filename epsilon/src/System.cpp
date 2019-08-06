// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/System.cpp
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

#include <System.h>
#include <Eigenvalues.h>
#include <JordanSystem.h>
#include <TransformationQueue.h>
#include <FermatException.h>
#include <Echelon.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
using namespace std;

FermatExpression infinity;

bool System::singLess::operator() (const FermatExpression &a, const FermatExpression &b) const {
    if (a == infinity) return false;
    if (b == infinity) return true;

    string as = a.str();
    string bs = b.str();

    bool anum=true,bnum=true;

    for (auto &c : as) {
        if (c != '-' && !isdigit(c)) {
            anum=false;
            break;
        }
    }
    for (auto &c : bs) {
        if (c != '-' && !isdigit(c)) {
            bnum=false;
            break;
        }
    }

    if (anum != bnum) return anum;

    if (anum && bnum) {
        int an = atoi(as.c_str());
        int bn = atoi(bs.c_str());
        return an < bn;
    } 

    return as < bs;
}

System::System(Fermat *fermat, bool echfer) : tqueue(fermat) {
    this->fermat = fermat;
    this->echfer = echfer;
}

System::System(Fermat *fermat, string filename, int start, int end, bool echfer) : tqueue(fermat) {
	string str;
	ifstream file(filename);
    int r;
    
    this->fermat = fermat;
    this->echfer = echfer;

    if (!file.is_open()) {
        throw invalid_argument("unable to open file.");
    }

    kmaxC = kmax = -1;

	while (getline(file,str)) {
        str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
        if (str == "") continue;

        size_t colon = str.find(':');
        if (colon == string::npos) {
            throw invalid_argument("parse error");
        }

        string str0(str,0,colon);
        str = str.erase(0,colon+1);
       
        if (str0.substr(0,2) == "A[" && str0.back() == ']') {
            size_t comma = str0.find(',');
            if (comma == string::npos) {
                throw invalid_argument("parse error");
            }

            sing_t singularity;
            TriangleBlockMatrix mat;

            singularity.point = FermatExpression(fermat,str0.substr(2,comma-2));
            singularity.rank = stoi(str0.substr(comma+1,str0.size()-comma-2));

            FermatArray array(fermat,str);
            r = array.rows();
            if (end<0) end=r;

            mat.A = FermatArray(array,1,start-1,1,start-1);
            mat.B = FermatArray(array,start,end,1,start-1);
            mat.C = FermatArray(array,start,end,start,end);
            mat.D = FermatArray(array,end+1,r,1,start-1);
            mat.E = FermatArray(array,end+1,r,start,end);
            mat.F = FermatArray(array,end+1,r,end+1,r);

            if (!FermatArray(array,1,start-1,start,r).isZero() || !FermatArray(array,start,end,end+1,r).isZero()) {
                throw invalid_argument("upper right blocks are not zero.");
            }

            _A[singularity] = mat;

            if (!mat.C.isZero()) {
                if (singularities.count(singularity.point)) {
                    singularities[singularity.point].rankC = max(singularities[singularity.point].rankC,singularity.rank);
                } else {
                    singularities[singularity.point].rankC = singularity.rank;
                    singularities[singularity.point].rank = -1;
                }
            }
            
            if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
                if (singularities.count(singularity.point)) {
                    singularities[singularity.point].rank = max(singularities[singularity.point].rank,singularity.rank);
                } else {
                    singularities[singularity.point].rankC = -1;
                    singularities[singularity.point].rank = singularity.rank;
                }
            }
        } else if (str0.substr(0,2) == "B[" && str0.back() == ']') {
            int k = stoi(str0.substr(2,str0.size()-3));
            TriangleBlockMatrix mat;
            int r;
            
            FermatArray array(fermat,str);
            r = array.rows();
            if (end<0) end=r;

            mat.A = FermatArray(array,1,start-1,1,start-1);
            mat.B = FermatArray(array,start,end,1,start-1);
            mat.C = FermatArray(array,start,end,start,end);
            mat.D = FermatArray(array,end+1,r,1,start-1);
            mat.E = FermatArray(array,end+1,r,start,end);
            mat.F = FermatArray(array,end+1,r,end+1,r);

            if (!FermatArray(array,1,start-1,start,r).isZero() || !FermatArray(array,start,end,end+1,r).isZero()) {
                throw invalid_argument("upper right blocks are not zero.");
            }

            _B[k] = mat;

            if (!mat.C.isZero()) {
                if (k>kmaxC) kmaxC=k;
            }
            if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
                if (k>kmax) kmax=k;
            }
        } else {
            throw invalid_argument("parse error.");
        }
	}
	file.close();

    nullMatrix.A = FermatArray(fermat,start-1,start-1);
    nullMatrix.B = FermatArray(fermat,end-start+1,start-1);
    nullMatrix.C = FermatArray(fermat,end-start+1,end-start+1);
    nullMatrix.D = FermatArray(fermat,r-end,start-1);
    nullMatrix.E = FermatArray(fermat,r-end,end-start+1);
    nullMatrix.F = FermatArray(fermat,r-end,r-end);

    tqueue.setpadding(start-1,r-end);

    singularities[infinity].rankC = -1;
    singularities[infinity].rank = -1;

    if (kmaxC<0) {
        if (!A(infinity,0).C.isZero()) {
            singularities[infinity].rankC = 0;
        } else {
            singularities[infinity].rankC = -1;
        }
    } else {
        singularities[infinity].rankC = kmaxC+1;
    }
    
    if (kmax<0) {
        TriangleBlockMatrix mat = A(infinity,0);
        if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
            singularities[infinity].rank = 0;
        } else {
            singularities[infinity].rank = -1;
        }
    } else {
        singularities[infinity].rank = kmax+1;
    }

    if (singularities[infinity].rank < 0) {
        singularities.erase(infinity);
    }
}

System::System(const System &orig, int start, int end) : tqueue(orig.tqueue) {
    int r;

    fermat = orig.fermat;
    echfer = orig.echfer;
 
    kmaxC = kmax = -1;

    for (auto it = orig._A.begin(); it != orig._A.end(); ++it) {
        FermatArray array = putTogether(it->second);
        TriangleBlockMatrix mat;
        r = array.rows();

        if (end<0) end=r;

        mat.A = FermatArray(array,1,start-1,1,start-1);
        mat.B = FermatArray(array,start,end,1,start-1);
        mat.C = FermatArray(array,start,end,start,end);
        mat.D = FermatArray(array,end+1,r,1,start-1);
        mat.E = FermatArray(array,end+1,r,start,end);
        mat.F = FermatArray(array,end+1,r,end+1,r);

        if (!FermatArray(array,1,start-1,start,r).isZero() || !FermatArray(array,start,end,end+1,r).isZero()) {
            throw invalid_argument("upper right blocks are not zero.");
        }

        _A[it->first] = mat;
            
        if (!mat.C.isZero()) {
            if (singularities.count(it->first.point)) {
                singularities[it->first.point].rankC = max(singularities[it->first.point].rankC,it->first.rank);
            } else {
                singularities[it->first.point].rankC = it->first.rank;
                singularities[it->first.point].rank = -1;
            }
        }
            
        if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
            if (singularities.count(it->first.point)) {
                singularities[it->first.point].rank = max(singularities[it->first.point].rank,it->first.rank);
            } else {
                singularities[it->first.point].rankC = -1;
                singularities[it->first.point].rank = it->first.rank;
            }
        }
    }

    for (auto it = orig._B.begin(); it != orig._B.end(); ++it) {
        FermatArray array = putTogether(it->second);
        TriangleBlockMatrix mat;
        r = array.rows();

        if (end<0) end=r;

        mat.A = FermatArray(array,1,start-1,1,start-1);
        mat.B = FermatArray(array,start,end,1,start-1);
        mat.C = FermatArray(array,start,end,start,end);
        mat.D = FermatArray(array,end+1,r,1,start-1);
        mat.E = FermatArray(array,end+1,r,start,end);
        mat.F = FermatArray(array,end+1,r,end+1,r);

        if (!FermatArray(array,1,start-1,start,r).isZero() || !FermatArray(array,start,end,end+1,r).isZero()) {
            throw invalid_argument("upper right blocks are not zero.");
        }

        _B[it->first] = mat;
            
        if (!mat.C.isZero()) {
            if (it->first>kmaxC) kmaxC=it->first;
        }
        if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
            if (it->first>kmax) kmax=it->first;
        }
    }
    
    nullMatrix.A = FermatArray(fermat,start-1,start-1);
    nullMatrix.B = FermatArray(fermat,end-start+1,start-1);
    nullMatrix.C = FermatArray(fermat,end-start+1,end-start+1);
    nullMatrix.D = FermatArray(fermat,r-end,start-1);
    nullMatrix.E = FermatArray(fermat,r-end,end-start+1);
    nullMatrix.F = FermatArray(fermat,r-end,r-end);

    tqueue.setpadding(start-1,r-end);
    
    singularities[infinity].rankC = -1;
    singularities[infinity].rank = -1;

    if (kmaxC<0) {
        if (!A(infinity,0).C.isZero()) {
            singularities[infinity].rankC = 0;
        } else {
            singularities[infinity].rankC = -1;
        }
    } else {
        singularities[infinity].rankC = kmaxC+1;
    }
    
    if (kmax<0) {
        TriangleBlockMatrix mat = A(infinity,0);
        if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
            singularities[infinity].rank = 0;
        } else {
            singularities[infinity].rank = -1;
        }
    } else {
        singularities[infinity].rank = kmax+1;
    }
}

System::System(const System &orig, const FermatArray &left, const FermatArray &right) : tqueue(orig.fermat) {
    TriangleBlockMatrix mat;

    fermat = orig.fermat;
    echfer = orig.echfer;
    nullMatrix = orig.nullMatrix;
    singularities = orig.singularities;
    kmaxC = orig.kmaxC;
    kmax = orig.kmax;

    mat.A = FermatArray();
    mat.D = FermatArray();
    mat.F = FermatArray();

    for (auto it = orig._A.begin(); it != orig._A.end(); ++it) {
        mat.C = FermatArray(left,it->second.C,right);

        if (it->second.B.cols() > 0) {
            mat.B = FermatArray(left,it->second.B);
        } else {
            mat.B = FermatArray();
        }

        if (it->second.E.rows() > 0) {
            mat.E = FermatArray(it->second.E,right);
        } else {
            mat.E = FermatArray();
        }
        
        _A[it->first] = mat;
    }

    for (auto it = orig._B.begin(); it != orig._B.end(); ++it) {
        mat.C = FermatArray(left,it->second.C,right);
        
        if (it->second.B.cols() > 0) {
            mat.B = FermatArray(left,it->second.B);
        } else {
            mat.B = FermatArray();
        }

        if (it->second.E.rows() > 0) {
            mat.E = FermatArray(it->second.E,right);
        } else {
            mat.E = FermatArray();
        }


        _B[it->first] = mat;
    }
}

System::System(const System &orig, int ep) : tqueue(orig.fermat) {
    TriangleBlockMatrix mat;
    
    fermat = orig.fermat;
    echfer = orig.echfer;
    nullMatrix = orig.nullMatrix;
    singularities = orig.singularities;
    kmaxC = orig.kmaxC;
    kmax = orig.kmax;

    mat.A = FermatArray();
    mat.B = FermatArray();
    mat.D = FermatArray();
    mat.E = FermatArray();
    mat.F = FermatArray();
    
    for (auto it = orig._A.begin(); it != orig._A.end(); ++it) {
        mat.C = it->second.C.subst("ep",ep);

        _A[it->first] = mat;
    }
    
    for (auto it = orig._B.begin(); it != orig._B.end(); ++it) {
        mat.C = it->second.C.subst("ep",ep);

        _B[it->first] = mat;
    }
}

Fermat *System::fer() const {
    return fermat;
}

int System::dimC() const {
    return nullMatrix.C.rows();
}

void System::write(string filename) const {
    ofstream file(filename);

    for (auto it = _A.begin(); it != _A.end(); ++it) {
        FermatArray A = putTogether(it->second);
        if (A.isZero()) continue;
        file << "A[" << it->first.point.str() << "," << it->first.rank << "]:  \t" << A.str() << endl;
    }

    for (auto it = _B.begin(); it != _B.end(); ++it) {
        FermatArray B = putTogether(it->second);
        if (B.isZero()) continue;
        file << "B[" << it->first << "]:    \t" << B.str() << endl;
    }

    file.close();
}

map<FermatExpression,FermatArray> System::exportFuchs() const {
    map<FermatExpression,FermatArray> fuchs;

    for (auto it = _B.begin(); it != _B.end(); ++it) {
        if (!it->second.C.isZero()) {
            throw invalid_argument("system not in fuchs form.");
        }
    }

    for (auto it = _A.begin(); it != _A.end(); ++it) {
        if (it->second.C.isZero()) continue;
        if (it->first.rank != 0) {
            throw invalid_argument("system not in fuchs form.");
        }

        fuchs[it->first.point] = it->second.C;
    }

    return fuchs;
}

TransformationQueue *System::transformationQueue() {
    return &tqueue;
}

void System::analyze() {
    auto check = [](const FermatArray &mat, int start, int end) -> int {
        for (int c=end+1; c<=mat.cols(); ++c) {
            FermatArray right(mat,start,end,c,mat.cols());

            if (right.isZero()) return c-1;
        }
        
        return mat.cols();
    };

    int start=1;
    vector<pair<int,int>> blocks;

    while(start <= nullMatrix.C.rows()) {
        int size=1;
        for(;;) {
            int size1=size;

            for (auto &a : _A) {
                size1 = check(a.second.C,start,start+size-1) - start + 1;
                if (size1 != size) break;
            }

            if (size1 != size) {
                size = size1;
                continue;
            }

            for (auto &b : _B) {
                size1 = check(b.second.C,start,start+size-1) - start + 1;
                if (size1 != size) break;
            }
                
            if (size1 == size) break;
            size = size1;
        }
        
        blocks.push_back({start,start+size-1});

        start += size;                
    }

    int cnt=1;
    int offset=nullMatrix.A.rows();

    for (auto &b : blocks) {
        int start = b.first;
        int end = b.second;
 
        cout << "block " << (cnt++) << ": [" << start+offset << "," << end+offset << "]" << endl;

        stringstream strm;
        for (auto &s : singularities) {
            int rank;
            FermatArray mat = nullMatrix.C;

            for (rank = s.second.rankC; rank >= 0; --rank) {
                FermatArray mat(A(s.first,rank).C,start,end,start,end);
                if (!mat.isZero()) break;
            }

            if (rank < 0) continue;

            strm << "  " << pstr(s.first) << ":" << rank;
        }

        if (strm.str().size()) {
            cout << "  singularities:" << strm.str() << endl;
        } else {
            cout << "  no singularities." << endl;
        }                    

        if (!nullMatrix.A.rows()) {
            strm.str("");
            strm.clear();

            for (auto &s : singularities) {
                int rank;
                FermatArray mat = nullMatrix.C;

                for (rank = s.second.rank; rank >= 0; --rank) {
                    FermatArray mat(A(s.first,rank).C,start,end,1,start-1);
                    if (!mat.isZero()) break;
                }

                if (rank < 0) continue;

                strm << "  " << pstr(s.first) << ":" << rank;
            }

            if (strm.str().size()) {
                cout << "  left singularities:" << strm.str() << endl;
            }
        }
    }
}
 
void System::fuchsify() {
    FermatExpression x1,x2;
    FermatArray Q;

    for(;;) {
        bool success=false;
        bool finished=true;

        printSingularities();

        // TODO: parallelize
        for (auto it = singularities.begin(); it != singularities.end() && !success; ++it) {
            if (it->second.rankC > 0) {
                finished = false;
                x1 = it->first;
                for (auto it2 = singularities.begin(); it2 != singularities.end(); ++it2) {
                    x2 = it2->first;
                    if (x1 == x2) continue;
                    if (it2->second.rankC < 0) continue;

                    if (projectorQ(x1,x2,Q)) {
                        success = true;
                        break;
                    }
                }
            }
        }
                        
        if (finished) break;

        if (!success) {
            for (auto it = singularities.begin(); it != singularities.end(); ++it) {
                if (it->second.rankC > 0) {
                    x1 = it->first;
                    x2 = regularPoint();

                    projectorP(x1,Q);
                    break;
                }
            }
        }
       
        balance(Q,x1,x2);
    }
}

void System::fuchsify(const FermatExpression &x1) {
    FermatExpression x2;
    FermatArray Q;

    if (!singularities.count(x1)) {
        cout << "no singularity at " << pstr(x1) << endl;
        return;
    }

    printSingularities();

    while (singularities[x1].rankC > 0) {
        bool success = false;

        for (auto &s : singularities) {
            x2 = s.first;
            if (x1 == x2) continue;
            if (s.second.rankC < 0) continue;

            if (projectorQ(x1,x2,Q)) {
                success = true;
                break;
            }
        }

        if (!success) {
            x2 = regularPoint();
            projectorP(x1,Q);
        }

        balance(Q,x1,x2);

        printSingularities();
    }
}

void System::normalize() {
    bool found=false;

    for (auto it=singularities.begin(); it != singularities.end(); ++it) {
        if (it->second.rankC > 0) {
            throw invalid_argument("not in fuchsian form");
        }
    }
 
    // find non apparent singularity
    FermatExpression x0;
    
    for (auto it=singularities.begin(); !found && it != singularities.end(); ++it) {
        if (it->first == infinity) continue;
        eigen(it->first);

        for (auto &e : eigenvalues[it->first]) {
            if (e.first.v != 0) {
                x0 = it->first;
                found = true;
                break;
            }
        }
    }

    if (!found) {
        auto it=singularities.rbegin();
        if (it->first == infinity) ++it;

        x0 = it->first;

        cout << "WARNING: All singularities have integer eigenvalues." << endl
             << "         choosing x0 = " << pstr(x0) << " which might be an apparent singularity." << endl;
    }

    printEigenvalues();
    
    for(;;) {
        FermatExpression x1,x2;
        FermatArray P;

        if (findBalance(x1,x2,P,FermatExpression())) {
            cout << "mutual balance [" << pstr(x1) << "," << pstr(x2) << "]" << endl;
        } else if (findBalance(x1,x2,P,x0)) {
            cout << "balance [" << pstr(x1) << "," << pstr(x2) << "]" << endl;
        } else {
            bool normalized=true;
            eigen(x0);

            for (auto &e : eigenvalues[x0]) {
                if (e.first.u != 0) {
                    normalized = false;
                    break;
                }
            }

            if (normalized) break;

            FermatExpression xr = regularPoint();

            if (!findBalance(x1,x2,P,xr)) {
                throw invalid_argument("unable to normalize system.");
            }
            cout << "balance with regular point [" << pstr(x1) << "," << pstr(x2) << "]" << endl;
        }

        balance(P,x1,x2);
        printEigenvalues();
    }
}

void System::factorep() {
    FermatExpression ep(fermat,"ep");
    int N = nullMatrix.C.rows();
    // TODO: check eigenvalues

    set<string> symbols;
    int rk;

    fermat->addSymbol("mu");
    FermatExpression mu(fermat,"mu");

    symbols.insert("mu");

    EchelonBase *echelon;

    if (echfer) {
        echelon = new EchelonFermat(fermat,N*N*singularities.size(),N*N+1);
    } else {
        echelon = new Echelon(fermat);
    }

    for (auto it = singularities.begin(); it != singularities.end(); ++it) {
        FermatExpression xj = it->first;
        if (xj == infinity) continue;

        FermatArray Aep = A(xj,0).C;
        FermatArray Amu = Aep.subst("ep",mu);

        #define pos(i,j) (((i)-1)*N+(j))

        for (int i=1; i<=N; ++i) {
            for (int k=1; k<=N; ++k) {
                map<int,FermatExpression> row;
                for (int j=1; j<=N; ++j) {
                    if (row.count(pos(j,k))) {
                        row[pos(j,k)] = row[pos(j,k)] + Aep(i,j)/ep;
                    } else {
                        row[pos(j,k)] = Aep(i,j)/ep;
                    }

                    if (row.count(pos(i,j))) {
                        row[pos(i,j)] = row[pos(i,j)] - Amu(j,k)/mu;
                    } else {
                        row[pos(i,j)] = -Amu(j,k)/mu;
                    }
                }
                echelon->set(row);
            }
        }

        #undef pos
    }

    rk = echelon->run(); 
   
    FermatArray T(fermat,N,N);
    T.assign("0");

    int pos=0;

    #define rpos(p) (((p)-1)/N+1)
    #define cpos(p) ((((p)-1)%N)+1)

    for (auto &r : *echelon) {
        for (int c=pos+1; c<r.col1(); ++c) {
            stringstream strm;
            string sym;

            strm << "t" << rpos(c) << "x" << cpos(c);
            sym = strm.str();

            if (!symbols.count(sym)) {
                fermat->addSymbol(sym);
                symbols.insert(sym);
                cout << "adding symbol " << sym << "." << endl;
            }

            T.set(rpos(c),cpos(c),FermatExpression(fermat,sym));
        }
       
        pos = r.col1();
 
        for (auto &e : r) {
            if (e.first == pos) {
                if (e.second.str() != "1") {
                    throw invalid_argument("not in row reduced echelon form (this is a bug)");
                }
                continue;
            }
            
            stringstream strm;
            string sym;
               
            strm << "t" << rpos(e.first) << "x" << cpos(e.first);
            sym = strm.str();
            
            if (!symbols.count(sym)) {
                fermat->addSymbol(sym);
                symbols.insert(sym);
                cout << "adding symbol " << sym << "." << endl;
            }

            T.set(rpos(pos),cpos(pos),T(rpos(pos),cpos(pos))-e.second*FermatExpression(fermat,sym));
        }
    }

    for (int c=pos+1; c<=N*N; ++c) {
        stringstream strm;
        string sym;

        strm << "t" << rpos(c) << "x" << cpos(c);
        sym = strm.str();

        if (!symbols.count(sym)) {
            fermat->addSymbol(sym);
            symbols.insert(sym);
            cout << "adding symbol " << sym << "." << endl;
        }

        T.set(rpos(c),cpos(c),FermatExpression(fermat,sym));
    }

    #undef rpos
    #undef cpos

    delete echelon;

    FermatExpression det = T.det();

    if (det.str() == "0") {
        throw invalid_argument("transformation is singular.");
    }
   
    for (auto it=symbols.begin(); it != symbols.end(); ++it) {
        for (int _sym0=0; _sym0 <= 200; ++_sym0) {
            int sym0 = ((_sym0&1)?1:-1)*((_sym0+1)>>1);
            try {
                FermatExpression det0 = det.subst(*it,sym0);
                if (det0.str() != "0") {
                    cout << *it << " -> " << sym0 << endl;
                    det = det0;
                    T = T.subst(*it,sym0);
                    fermat->dropSymbol(*it);
                    break;
                }
            } catch (const FermatDivByZero &e) {
            }
        }
    }

    transform(T);
}

void System::factorep(int mu) {
    if (mu == 0) {
        throw invalid_argument("mu must be != 0");
    }

    FermatExpression ep(fermat,"ep");
    int N = nullMatrix.C.rows();
    // TODO: check eigenvalues

    set<string> symbols;
    int rk;
    
    EchelonBase *echelon;

    if (echfer) {
        echelon = new EchelonFermat(fermat,N*N*singularities.size(),N*N+1);
    } else {
        echelon = new Echelon(fermat);
    }

    try {
        for (auto it = singularities.begin(); it != singularities.end(); ++it) {
            FermatExpression xj = it->first;
            if (xj == infinity) continue;

            FermatArray Aep = A(xj,0).C;
            FermatArray Amu = Aep.subst("ep",mu);

            #define pos(i,j) (((i)-1)*N+(j))

            for (int i=1; i<=N; ++i) {
                for (int k=1; k<=N; ++k) {
                    map<int,FermatExpression> row;
                    for (int j=1; j<=N; ++j) {
                        if (row.count(pos(j,k))) {
                            row[pos(j,k)] = row[pos(j,k)] + Aep(i,j)/ep;
                        } else {
                            row[pos(j,k)] = Aep(i,j)/ep;
                        }

                        if (row.count(pos(i,j))) {
                            row[pos(i,j)] = row[pos(i,j)] - Amu(j,k)/mu;
                        } else {
                            row[pos(i,j)] = -Amu(j,k)/mu;
                        }
                    }
                    echelon->set(row);
                }
            }

            #undef pos
        }
    } catch(const FermatDivByZero &e) {
        stringstream strm;
        strm << "system is singular at mu=" << mu << ".";
        throw invalid_argument(strm.str());
    }

    rk = echelon->run();
 
    FermatArray T(fermat,N,N);
    T.assign("0");

    int pos=0;

    #define rpos(p) (((p)-1)/N+1)
    #define cpos(p) ((((p)-1)%N)+1)

    for (auto &r : *echelon) {
        for (int c=pos+1; c<r.col1(); ++c) {
            stringstream strm;
            string sym;

            strm << "t" << rpos(c) << "x" << cpos(c);
            sym = strm.str();

            if (!symbols.count(sym)) {
                fermat->addSymbol(sym);
                symbols.insert(sym);
                cout << "adding symbol " << sym << "." << endl;
            }

            T.set(rpos(c),cpos(c),FermatExpression(fermat,sym));
        }
       
        pos = r.col1();
 
        for (auto &e : r) {
            if (e.first == pos) {
                if (e.second.str() != "1") {
                    throw invalid_argument("not in row reduced echelon form (this is a bug)");
                }
                continue;
            }
            
            stringstream strm;
            string sym;
               
            strm << "t" << rpos(e.first) << "x" << cpos(e.first);
            sym = strm.str();
            
            if (!symbols.count(sym)) {
                fermat->addSymbol(sym);
                symbols.insert(sym);
                cout << "adding symbol " << sym << "." << endl;
            }

            T.set(rpos(pos),cpos(pos),T(rpos(pos),cpos(pos))-e.second*FermatExpression(fermat,sym));
        }
    }

    for (int c=pos+1; c<=N*N; ++c) {
        stringstream strm;
        string sym;

        strm << "t" << rpos(c) << "x" << cpos(c);
        sym = strm.str();

        if (!symbols.count(sym)) {
            fermat->addSymbol(sym);
            symbols.insert(sym);
            cout << "adding symbol " << sym << "." << endl;
        }

        T.set(rpos(c),cpos(c),FermatExpression(fermat,sym));
    }
    
    #undef rpos
    #undef cpos

    delete echelon;

    FermatExpression det = T.det();

    if (det.str() == "0") {
        throw invalid_argument("transformation is singular.");
    }
   
    for (auto it=symbols.begin(); it != symbols.end(); ++it) {
        for (int _sym0=0; _sym0 <= 200; ++_sym0) {
            int sym0 = ((_sym0&1)?1:-1)*((_sym0+1)>>1);
            try {
                FermatExpression det0 = det.subst(*it,sym0);
                if (det0.str() != "0") {
                    cout << *it << " -> " << sym0 << endl;
                    det = det0;
                    T = T.subst(*it,sym0);
                    fermat->dropSymbol(*it);
                    break;
                }
            } catch (const FermatDivByZero &e) {
            }
        }
    }

    transform(T);
}

void System::leftranks() {
    for (auto &s : singularities) {
        FermatExpression xj = s.first;
        int k;

        for (k=s.second.rank; k>=0 && A(xj,k).B.isZero(); --k);

        if (k>=0) cout << "rank:\t" << pstr(xj) << ":" << k << endl;
    }
}

int System::leftreduce(const FermatExpression &xj) {
    int k;
    for (k=singularities.at(xj).rank; k>=0 && A(xj,k).B.isZero(); --k);

    if (k<=0) {
        cout << pstr(xj) << " is already a fuchsian singularity." << endl;
        return k;
    }

    FermatArray B = A(xj,k).B;
    TriangleBlockMatrix A0 = A(xj,0);

    EchelonBase *echelon;

    if (echfer) {
        echelon = new EchelonFermat(fermat,B.rows()*B.cols(),B.rows()*B.cols()+2);
    } else {
        echelon = new Echelon(fermat);
    }

    #define pos(i,j) (((i)-1)*B.cols()+(j))

    for (int i=1; i<=B.rows(); ++i) {
        for (int j=1; j<=B.cols(); ++j) {
            map<int,FermatExpression> row;

            row[B.rows()*B.cols()+1] = B(i,j);
            row[pos(i,j)] = FermatExpression(fermat,-k);

            for (int n=1; n<=B.rows(); ++n) {
                if (row.count(pos(n,j))) {
                    row[pos(n,j)] = row[pos(n,j)]-A0.C(i,n);
                } else {
                    row[pos(n,j)] = -A0.C(i,n);
                }
            }
            for (int n=1; n<=B.cols(); ++n) {
                if (row.count(pos(i,n))) {
                    row[pos(i,n)] = row[pos(i,n)]+A0.A(n,j);
                } else {
                    row[pos(i,n)] = A0.A(n,j);
                }
            }

            echelon->set(row);
        }
    }
    
    #undef pos

    int rk = echelon->run();

    FermatArray G(fermat,B.rows(),B.cols());
    G.assign("0");

    for (auto &r : *echelon) {
        auto it = r.begin();
        int pos = it->first;

        if (pos == B.rows()*B.cols()+1) {
            throw invalid_argument("linear system has no solution.");
        }

        if (r.size() < 2) continue;

        for (auto it2=r.begin(); it2 != r.end(); ++it2) {
            it = it2;
        }
       
        if (it->first != B.rows()*B.cols()+1) continue; 

        int row = ((pos-1)/B.cols())+1;
        int col = ((pos-1)%B.cols())+1;

        G.set(row,col,it->second);
    }

    delete echelon;

    lefttransform(G,xj,k);

    if (!A(xj,k).B.isZero()) {
        throw invalid_argument("transformation failed.");
    }

    for (--k; k>=0 && A(xj,k).B.isZero(); --k);

    cout << "new rank:\t" << pstr(xj) <<":" << k << endl;

    return k;
}

void System::leftfuchsify() {
    auto sings = singularities;

    for (auto &s : sings) {
        FermatExpression xj = s.first;
        int k;

        for (k=s.second.rank; k>=0 && A(xj,k).B.isZero(); --k);

        if (k>=0) cout << "rank:    \t" << pstr(xj) << ":" << k << endl;

        while (k>0) {
            k = leftreduce(xj);
        }            
    }
}

void System::leftfuchsify(const FermatExpression &xj) {
    int k;

    if (!singularities.count(xj)) {
        cout << "no singularity at " << pstr(xj) << endl;
        return;
    }
    
    for (k=singularities[xj].rank; k>=0 && A(xj,k).B.isZero(); --k);

    if (k<0) {
        cout << "off-diagonal block is not singular at " << pstr(xj) << endl;
        return;
    }

    cout << "rank:    \t" << pstr(xj) << ":" << k << endl;

    while (k>0) {
        k = leftreduce(xj);
    }            
}

void System::tjordan(const FermatExpression &xj, bool divep) {
    if (!singularities.count(xj)) {
        cout << "no singularity at " << xj.str() << endl;
        return;
    }

    FermatArray C = A(xj,singularities.at(xj).rankC).C;
    eigenvalues_t evs;
    JordanSystem sys;
    FermatArray T(fermat,nullMatrix.C.rows(),nullMatrix.C.cols());

    eigen(xj);

    if (divep) {
        C = C/FermatExpression(fermat,"ep");
        
        for (auto &ev : eigenvalues[xj]) {
            if (ev.first.u != 0) {
                throw invalid_argument("eigenvalues have to be proportional to ep");
            }
            eigen_t ev0;

            ev0.u = ev.first.v;
            ev0.v = 0;

            evs[ev0] = ev.second;
        }
    } else {
        evs = eigenvalues[xj];
    }                

    jordanSystem(C,evs,sys);

    int i = 1;
    for (auto &b : sys) {
        for (auto &v : b.rootvectors) {
            T.setColumn(i++,v);
        }
    }
    
    if (i != T.rows()+1) throw invalid_argument("wrong number of root vectors");

    transform(T);
}


bool System::projectorQ(const FermatExpression &x1, const FermatExpression &x2, FermatArray &Q) {
    list<JordanBlock> inv;
    int i,k,k0;
    set<int> S; 

    TriangleBlockMatrix A1 = A(x1,singularities[x1].rankC-1);

    jordan(x1);

    inverseJordan(x1,inv);

    FermatArray U0(fermat,nullMatrix.C.cols(),jordans[x1].size());
    FermatArray V0(fermat,inv.size(),nullMatrix.C.rows());

    i=1;
    for (auto &b : jordans[x1]) {
        U0.setColumn(i++,*(b.rootvectors.begin()));
    }

    i=1;
    for (auto &b : inv) {
        V0.setRow(i++,*(b.rootvectors.begin()));
    }
 
    FermatArray L0(V0,A1.C,U0);
    FermatArray L1(V0,U0);
    FermatArray Delta;

    for (k=0; k<L1.rows() && L1(k+1,k+1).str() == "0" ; ++k);

    k0 = reduceL0(L0,k,x1,S,Delta);

    FermatArray id(fermat,Delta.rows(),Delta.cols());
    id.assign("[1] + 0");

    U0 = U0*(id + Delta);

    S.insert(k0);

    FermatArray Uk(fermat,U0.rows(),S.size());
    FermatArray Vk;

    i=1;
    for (auto it=S.begin(); it != S.end(); ++it) {
        Uk.setColumn(i++,FermatArray(U0,1,U0.rows(),*it,*it));
    }

    if (!invariantSubspace(x2,Uk,Vk)) {
        return false;
    }

    Q = Uk*Vk.transpose();

    return true;
}

void System::projectorP(const FermatExpression &x1, FermatArray &P) {
    list<JordanBlock> inv;
    int i,k,k0;
    set<int> S; 

    TriangleBlockMatrix A1 = A(x1,singularities[x1].rankC-1);
    
    jordan(x1);

    inverseJordan(x1,inv);
    
    FermatArray U0(fermat,nullMatrix.C.cols(),jordans[x1].size());
    FermatArray V0(fermat,inv.size(),nullMatrix.C.rows());
    FermatArray Vn(fermat,nullMatrix.C.cols(),inv.size());
    
    i=1;
    for (auto &b : jordans[x1]) {
        U0.setColumn(i++,*(b.rootvectors.begin()));
    }
    
    i=1;
    for (auto it = inv.begin(); it != inv.end(); ++it) {
        V0.setRow(i,*(it->rootvectors.begin()));
        Vn.setColumn(i,it->rootvectors.back());
        ++i;
    }
    
    FermatArray L0(V0,A1.C,U0);
    FermatArray L1(V0,U0);
    FermatArray Delta;

    for (k=0; k<L1.rows() && L1(k+1,k+1).str() == "0" ; ++k);

    k0 = reduceL0(L0,k,x1,S,Delta);
    
    FermatArray id(fermat,Delta.rows(),Delta.cols());
    id.assign("[1] + 0");

    FermatArray xDelta;

    for (FermatArray sDelta(id); !sDelta.isZero(); sDelta = sDelta*(-Delta)) {
        xDelta += sDelta;
    }

    xDelta = xDelta.transpose();

    U0 = U0*(id + Delta);
    Vn = Vn*xDelta;

    S.insert(k0);

    FermatArray Uk(fermat,U0.rows(),S.size());
    FermatArray Vk(fermat,Vn.rows(),S.size());
    
    i=1;
    for (auto it=S.begin(); it != S.end(); ++it) {
        Uk.setColumn(i,FermatArray(U0,1,U0.rows(),*it,*it));
        Vk.setColumn(i,FermatArray(Vn,1,Vn.rows(),*it,*it));
        ++i;
    }

    P = Uk*Vk.transpose();
}

int System::reduceL0(FermatArray L0, int k, const FermatExpression &x1, set<int> &S, FermatArray &Delta) {
    int i;
    FermatArray id(fermat,L0.rows(),L0.cols());

    id.assign("[1] + 0");

    jordan(x1);

    S.clear();
    Delta = FermatArray(fermat,L0.rows(),L0.cols());
    Delta.assign("0");

    do {
        int rank=0;
        FermatArray L0x;
        FermatArray cj;
        multiset<JordanBlock>::iterator iit;

        for (i=1; i<=L0.rows(); ++i) {
            if (S.count(i)) continue;
            L0x = L0x.concatenate(FermatArray(L0,i,i,1,L0.cols()));
        }

        for (i=1,iit = jordans[x1].begin(); i<=L0x.cols(); ++i,++iit) {
            FermatArray L0p(L0x,1,L0x.rows(),1,i);
            int newrank = L0p.rank();
           
            if (!S.count(i) && newrank == rank) {
                if (FermatArray(L0p,1,L0p.rows(),i,i).isZero()) break;

                int rk = L0p.rowEchelon();
                
                cj = FermatArray(fermat,i-1,1);
                cj.assign("0");

                for (int n=1; n<=rk; ++n) {
                    for (int j=n; j<=i-1; ++j) {
                        if (L0p(n,j).str() == "1") {
                            cj.set(j,1,L0p(n,i));
                        }
                    }
                }

                break;
            }
            rank = newrank;
        }

        if (cj.rows()) {
            FermatArray Delta0(fermat,L0.rows(),L0.cols());
            FermatArray Delta0x(fermat,L0.rows(),L0.cols());
            multiset<JordanBlock>::iterator nit;
            int n;

            Delta0.assign("[1] + 0");
            Delta0x.assign("[1] + 0");

            for (n=1, nit = jordans[x1].begin(); n<=i-1;++n,++nit) {
                FermatExpression c = cj(n,1);
                Delta0.set(n,i,-c);
                if (nit->rootvectors.size() == iit->rootvectors.size()) {
                    Delta0x.set(n,i,c);
                }
            }
            L0 = Delta0x * L0 * Delta0;
            Delta = Delta0 - id + Delta * Delta0;
        }

        S.insert(i);
    } while(i>k);

    S.erase(i);
   
    return i;
}

bool System::invariantSubspace(const FermatExpression &x2, const FermatArray &Uk, FermatArray &Vk) {
    list<JordanBlock> inv;
    set<int> found;

    inverseJordan(x2,inv);

    Vk = FermatArray(fermat,Uk.rows(),Uk.cols());

    while(inv.size()) {
        for (auto block = inv.begin(); block != inv.end();) {
            FermatArray v = block->rootvectors.front();
            block->rootvectors.pop_front();

            FermatArray test = v.transpose() * Uk;
            FermatExpression norm;
            int pos=-1;

            for (int n=1; n<=test.cols(); ++n) {
                FermatExpression sprod = test(1,n);
                if (sprod.str() != "0") {
                    if (pos<0) {
                        pos = n;
                        norm = sprod;
                    } else {
                        pos = -1;
                        break;
                    }
                }
            }

            if (block->rootvectors.empty() || pos<0 || found.count(pos)) {
                block = inv.erase(block);
            } else {
                ++block;
            }

            if (pos>0) {
                v = v/norm;

                Vk.setColumn(pos,v);
                found.insert(pos);
            
                if (found.size() == Uk.cols()) {
                    return true;
                }
            }

        }
    }

    return false;
}

bool System::findBalance(FermatExpression &x1, FermatExpression &x2, FermatArray &P, const FermatExpression &x0) {
    map<FermatExpression,poincareRank,singLess> left,right;
    bool second=false;
    size_t len=0;

    if (x0.fer()) {
        if (singularities.count(x0)) {
            left[x0] = singularities[x0];
        } else {
            left[x0].rank = left[x0].rankC = -1; 
        }
        
        right = singularities;
    } else {
        left = right = singularities;
    }

    //TODO: parallelize
    do {
        for (auto &l : left) {
            if (x0.fer() && second && l.first == x0) continue;

            FermatArray A0 = A(l.first,0).C;

            eigen(l.first);

            for (auto &e1 : eigenvalues[l.first]) {
                if ((!x0.fer() || second) && e1.first.u >= 0) continue;

                vector<FermatArray> vectors1;
                Eigenvectors(A0,e1.first,vectors1);
                for (auto &r : right) {
                    if (l.first == r.first) continue;

                    FermatArray B0 = A(r.first,0).C.transpose();
                    eigen(r.first);

                    for (auto &e2 : eigenvalues[r.first]) {
                        if ((!x0.fer() || !second) && e2.first.u <= 0) continue;

                        vector<FermatArray> vectors2;
                        Eigenvectors(B0,e2.first,vectors2);

                        for (auto &v1 : vectors1) {
                            for (auto &v2 : vectors2) {
                                FermatExpression expr = (v1.transpose() * v2)(1,1);

                                if (expr.str() == "0") continue;

                                FermatArray P0 = v1 * v2.transpose() / expr;
                                size_t len0 = P0.str().size();

                                if (!len || len0 < len) {
                                    len = len0;
                                    P = P0;
                                    x1 = l.first;
                                    x2 = r.first;
                                }                                    
                            }
                        }                            
                    }
                }
            }
        }

        if (x0.fer()) {
            if (!second) {
                swap(left,right);
                second = true;
            } else {
                second = false;
            }
        }
    } while(second);

    return len>0;
}

void System::balance(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2) {
    if (x1 == infinity) {
        balance_inf_x2(P,x2);
    } else if (x2 == infinity) {
        balance_x1_inf(P,x1);
    } else {
        balance_x1_x2(P,x1,x2);
    }

    jordans.clear();
    eigenvalues.erase(x1);
    eigenvalues.erase(x2);

    updatePoincareRanks();

    tqueue.balance(P,x1,x2);
}

void System::balance_x1_x2(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2) {
    FermatArray id(fermat,P.rows(),P.cols());
    id.assign("[1] + 0");
    sing_t sing;

    System PMPbar(*this,P,id-P);
    System PbarMP(*this,id-P,P);

    // A(x1,0)
    
    sing.point=x1;
    sing.rank=0;
   
    if (!singularities.count(x1)) {
        singularities[x1].rank = -1;
        singularities[x1].rankC = -1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = nullMatrix.B;
        _A[sing].C = nullMatrix.C;
        _A[sing].D = nullMatrix.D;
        _A[sing].E = nullMatrix.E;
        _A[sing].F = nullMatrix.F;
    }

    for (int n=0; n <= singularities[x1].rank; ++n) {
        _A[sing].C -= PMPbar.A(x1,n).C/pow(x2-x1,n);
        _A[sing].B -= PMPbar.A(x1,n).B/pow(x2-x1,n);
    }

    for (auto it = PbarMP._A.begin(); it != PbarMP._A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;
        if (xj == x1) continue;

        _A[sing].C += it->second.C*(x1-x2)/pow(x1-xj,n+1);
        _A[sing].E += it->second.E*(x1-x2)/pow(x1-xj,n+1);
    }

    for (int n=0; n<=kmax; ++n) {
        _A[sing].C += PbarMP.B(n).C*pow(x1,n)*(x1 - x2);
        _A[sing].E += PbarMP.B(n).E*pow(x1,n)*(x1 - x2);
    }

    _A[sing].C += P;

    // A(x1,k>0)
    
    for (int k=1; k <= singularities[x1].rank; ++k) {
        sing.point = x1;
        sing.rank = k;

        _A[sing].C += PbarMP.A(x1,k-1).C*(x1-x2);
        _A[sing].E += PbarMP.A(x1,k-1).E*(x1-x2);
        
        for (int n=0; n+k<=singularities[x1].rank; ++n) {
            _A[sing].C -= PMPbar.A(x1,n+k).C/pow(x2-x1,n);
            _A[sing].B -= PMPbar.A(x1,n+k).B/pow(x2-x1,n);
        }
    }

    if (!PbarMP.A(x1,singularities[x1].rank).C.isZero() || !PbarMP.A(x1,singularities[x1].rank).E.isZero()) {
        sing.point = x1;
        sing.rank = singularities[x1].rank+1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = nullMatrix.B;
        _A[sing].C = PbarMP.A(x1,singularities[x1].rank).C*(x1-x2);
        _A[sing].D = nullMatrix.D;
        _A[sing].E = PbarMP.A(x1,singularities[x1].rank).E*(x1-x2);
        _A[sing].F = nullMatrix.F;
    } 
    
    // A(x2,0)
    
    sing.point=x2;
    sing.rank=0;

    if (!singularities.count(x2)) {
        singularities[x2].rank = -1;
        singularities[x2].rankC = -1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = nullMatrix.B;
        _A[sing].C = nullMatrix.C;
        _A[sing].D = nullMatrix.D;
        _A[sing].E = nullMatrix.E;
        _A[sing].F = nullMatrix.F;
    }

    for (int n=0; n<=singularities[x2].rank; ++n) {
        _A[sing].C -= PbarMP.A(x2,n).C/pow(x1-x2,n);
        _A[sing].E -= PbarMP.A(x2,n).E/pow(x1-x2,n);
    }

    for (auto it=PMPbar._A.begin(); it != PMPbar._A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;
        if (xj == x2) continue;

        _A[sing].C += it->second.C*(x2-x1)/pow(x2-xj,n+1);
        _A[sing].B += it->second.B*(x2-x1)/pow(x2-xj,n+1);
    }

    for (int n=0; n<=kmax; ++n) {
        _A[sing].C += PMPbar.B(n).C*pow(x2,n)*(x2 - x1);
        _A[sing].B += PMPbar.B(n).B*pow(x2,n)*(x2 - x1);
    }

    _A[sing].C -= P;

    // A(x2,k>0)
    
    for (int k=1; k <= singularities[x2].rank; ++k) {
        sing.point = x2;
        sing.rank = k;

        _A[sing].C += PMPbar.A(x2,k-1).C*(x2-x1);
        _A[sing].B += PMPbar.A(x2,k-1).B*(x2-x1);

        for (int n=0; n<=singularities[x2].rank-k; ++n) {
            _A[sing].C -= PbarMP.A(x2,n+k).C/pow(x1-x2,n);
            _A[sing].E -= PbarMP.A(x2,n+k).E/pow(x1-x2,n);
        }
    }
    
    if (!PMPbar.A(x2,singularities[x2].rank).B.isZero() || !PMPbar.A(x2,singularities[x2].rank).C.isZero()) {
        sing.point = x2;
        sing.rank = singularities[x2].rank+1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = PMPbar.A(x2,singularities[x2].rank).B*(x2-x1);
        _A[sing].C = PMPbar.A(x2,singularities[x2].rank).C*(x2-x1);
        _A[sing].D = nullMatrix.D;
        _A[sing].E = nullMatrix.E;
        _A[sing].F = nullMatrix.F;
    } 

    // A(xj != x1 && xj != x2,k)

    for (auto it=_A.begin(); it != _A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int k = it->first.rank;

        if (xj == x1 || xj == x2) continue;

        for (int n=0; n+k<=singularities[xj].rank; ++n) {
            it->second.C += PbarMP.A(xj,n+k).C*(x2-x1)/pow(x1-xj,n+1);
            it->second.E += PbarMP.A(xj,n+k).E*(x2-x1)/pow(x1-xj,n+1);
            it->second.C += PMPbar.A(xj,n+k).C*(x1-x2)/pow(x2-xj,n+1);
            it->second.B += PMPbar.A(xj,n+k).B*(x1-x2)/pow(x2-xj,n+1);
        }
    }

    // B(k)

    for (int k=0; k<=kmax; ++k) {
        for (int n=0; k+n+1<=kmax; ++n) {
            _B[k].C += (PbarMP.B(k+n+1).C*pow(x1,n) - PMPbar.B(k+n+1).C*pow(x2,n))*(x1-x2);
            _B[k].B -= PMPbar.B(k+n+1).B*pow(x2,n)*(x1-x2);
            _B[k].E += PbarMP.B(k+n+1).E*pow(x1,n)*(x1-x2);
        }
    }
}

void System::balance_x1_inf(const FermatArray &P, const FermatExpression &x1) {
    FermatArray id(fermat,P.rows(),P.cols());
    id.assign("[1] + 0");
    sing_t sing;

    System PMPbar(*this,P,id-P);
    System PbarMP(*this,id-P,P);

    // A(x1,0)
    
    sing.point=x1;
    sing.rank=0;
    
    _A[sing].C += -PbarMP.A(x1,0).C - PMPbar.A(x1,0).C + PMPbar.A(x1,1).C;
    _A[sing].B += -PMPbar.A(x1,0).B + PMPbar.A(x1,1).B;
    _A[sing].E -= PbarMP.A(x1,0).E;
    
    for (auto it = PbarMP._A.begin(); it != PbarMP._A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;
        if (xj == x1) continue;

        _A[sing].C += it->second.C/pow(x1-xj,n+1);
        _A[sing].E += it->second.E/pow(x1-xj,n+1);
    }

    for (int n=0; n<=kmax; ++n) {
        _A[sing].C += PbarMP.B(n).C*pow(x1,n);
        _A[sing].E += PbarMP.B(n).E*pow(x1,n);
    }

    _A[sing].C += P;

    // A(x1,k>0)
    for (int k=1; k<=singularities[x1].rank; ++k) {
        sing.point = x1;
        sing.rank = k;

        _A[sing].C += -PbarMP.A(x1,k).C - PMPbar.A(x1,k).C + PMPbar.A(x1,k+1).C + PbarMP.A(x1,k-1).C;
        _A[sing].B += -PMPbar.A(x1,k).B + PMPbar.A(x1,k+1).B;
        _A[sing].E += -PbarMP.A(x1,k).E + PbarMP.A(x1,k-1).E;
    }

    if (!PbarMP.A(x1,singularities[x1].rank).C.isZero() || !PbarMP.A(x1,singularities[x1].rank).E.isZero()) {
        sing.point = x1;
        sing.rank = singularities[x1].rank+1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = nullMatrix.B;
        _A[sing].C += PbarMP.A(x1,singularities[x1].rank).C;
        _A[sing].D = nullMatrix.D;
        _A[sing].E += PbarMP.A(x1,singularities[x1].rank).E;
        _A[sing].F = nullMatrix.F;
    }

    // A(xj != x1, k)

    for (auto it=_A.begin(); it != _A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int k = it->first.rank;

        if (xj == x1) continue;

        it->second.C += -PbarMP.A(xj,k).C - PMPbar.A(xj,k).C + PMPbar.A(xj,k+1).C + PMPbar.A(xj,k).C*(xj-x1);
        it->second.B += -PMPbar.A(xj,k).B + PMPbar.A(xj,k+1).B + PMPbar.A(xj,k).B*(xj-x1);
        it->second.E -= PbarMP.A(xj,k).E;

        for (int n=0; n+k<=singularities[xj].rank; ++n) {
            it->second.C -= PbarMP.A(xj,n+k).C/pow(x1-xj,n+1);
            it->second.E -= PbarMP.A(xj,n+k).E/pow(x1-xj,n+1);
        }
    }

    // B(0)

    if (!_B.count(0)) _B[0] = nullMatrix;

    _B[0].C += -PbarMP.B(0).C - PMPbar.B(0).C - PMPbar.B(0).C*x1;
    _B[0].B += -PMPbar.B(0).B - PMPbar.B(0).B*x1;
    _B[0].E -= PbarMP.B(0).E;

    for (auto it=singularities.begin(); it != singularities.end(); ++it) {
        FermatExpression xj = it->first;
        if (xj == infinity) continue;

        _B[0].C += PMPbar.A(xj,0).C;
        _B[0].B += PMPbar.A(xj,0).B;
    }

    for (int n=0; n+1<=kmax; ++n) {
        _B[0].C += PbarMP.B(n+1).C*pow(x1,n);
        _B[0].E += PbarMP.B(n+1).E*pow(x1,n);
    }

    // B(k > 0)

    for (int k=1; k<=kmax; ++k) {
        _B[k].C += -PbarMP.B(k).C - PMPbar.B(k).C + PMPbar.B(k-1).C - PMPbar.B(k).C*x1;
        _B[k].B += -PMPbar.B(k).B + PMPbar.B(k-1).B - PMPbar.B(k).B*x1;
        _B[k].E += -PbarMP.B(k).E;
        for (int n=0; k+n+1 <= kmax; ++n) {
            _B[k].C += PbarMP.B(k+n+1).C*pow(x1,n);
            _B[k].E += PbarMP.B(k+n+1).E*pow(x1,n);
        }
    }

    if (!PMPbar.B(kmax).B.isZero() || !PMPbar.B(kmax).C.isZero()) {
        _B[kmax+1].A = nullMatrix.A;
        _B[kmax+1].B = PMPbar.B(kmax).B;
        _B[kmax+1].C = PMPbar.B(kmax).C;
        _B[kmax+1].D = nullMatrix.D;
        _B[kmax+1].E = nullMatrix.E;
        _B[kmax+1].F = nullMatrix.F;
    }
}

void System::balance_inf_x2(const FermatArray &P, const FermatExpression &x2) {
    FermatArray id(fermat,P.rows(),P.cols());
    id.assign("[1] + 0");
    sing_t sing;

    System PMPbar(*this,P,id-P);
    System PbarMP(*this,id-P,P);
    
    // A(x2,0)

    sing.point=x2;
    sing.rank=0;
    
    if (!singularities.count(x2)) {
        singularities[x2].rank = -1;
        singularities[x2].rankC = -1;

        _A[sing].A = nullMatrix.A;
        _A[sing].B = nullMatrix.B;
        _A[sing].C = nullMatrix.C;
        _A[sing].D = nullMatrix.D;
        _A[sing].E = nullMatrix.E;
        _A[sing].F = nullMatrix.F;
    }

    _A[sing].C += -PbarMP.A(x2,0).C - PMPbar.A(x2,0).C  + PbarMP.A(x2,1).C;
    _A[sing].B -= PMPbar.A(x2,0).B;
    _A[sing].E += -PbarMP.A(x2,0).E + PbarMP.A(x2,1).E;

    for (auto it = PMPbar._A.begin(); it != PMPbar._A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;
        if (xj == x2) continue;

        _A[sing].C += it->second.C/pow(x2-xj,n+1);
        _A[sing].B += it->second.B/pow(x2-xj,n+1);
    }

    for (int n=0; n<=kmax; ++n) {
        _A[sing].C += PMPbar.B(n).C*pow(x2,n);
        _A[sing].B += PMPbar.B(n).B*pow(x2,n);
    }

    _A[sing].C -= P;

    // A(x2,k>0)
    for (int k=1; k<=singularities[x2].rank; ++k) {
        sing.point = x2;
        sing.rank = k;

        _A[sing].C += -PbarMP.A(x2,k).C - PMPbar.A(x2,k).C + PMPbar.A(x2,k-1).C + PbarMP.A(x2,k+1).C;
        _A[sing].B += -PMPbar.A(x2,k).B + PMPbar.A(x2,k-1).B;
        _A[sing].E += -PbarMP.A(x2,k).E + PbarMP.A(x2,k+1).E;
    }
    
    if (!PMPbar.A(x2,singularities[x2].rank).B.isZero() || !PMPbar.A(x2,singularities[x2].rank).C.isZero()) {
        sing.point = x2;
        sing.rank = singularities[x2].rank+1;
        
        _A[sing].A = nullMatrix.A;
        _A[sing].B = PMPbar.A(x2,singularities[x2].rank).B;
        _A[sing].C = PMPbar.A(x2,singularities[x2].rank).C;
        _A[sing].D = nullMatrix.D;
        _A[sing].E = nullMatrix.E;
        _A[sing].F = nullMatrix.F;
    }

    // A(xj != x2, k)

    for (auto it=_A.begin(); it != _A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int k = it->first.rank;

        if (xj == x2) continue;

        it->second.C += -PbarMP.A(xj,k).C - PMPbar.A(xj,k).C + PbarMP.A(xj,k+1).C + PbarMP.A(xj,k).C*(xj-x2);
        it->second.B -= PMPbar.A(xj,k).B;
        it->second.E += -PbarMP.A(xj,k).E + PbarMP.A(xj,k+1).E + PbarMP.A(xj,k).E*(xj-x2);

        for (int n=0; n+k<=singularities[xj].rank; ++n) {
            it->second.C -= PMPbar.A(xj,n+k).C/pow(x2-xj,n+1);
            it->second.B -= PMPbar.A(xj,n+k).B/pow(x2-xj,n+1);
        }
    }

    // B(0)

    if (!_B.count(0)) _B[0] = nullMatrix;

    _B[0].C += -PbarMP.B(0).C - PMPbar.B(0).C - PbarMP.B(0).C*x2;
    _B[0].B -= PMPbar.B(0).B;
    _B[0].E += -PbarMP.B(0).E - PbarMP.B(0).E*x2;

    for (auto it=singularities.begin(); it != singularities.end(); ++it) {
        FermatExpression xj = it->first;
        if (xj == infinity) continue;

        _B[0].C += PbarMP.A(xj,0).C;
        _B[0].E += PbarMP.A(xj,0).E;
    }

    for (int n=0; n+1<=kmax; ++n) {
        _B[0].C += PMPbar.B(n+1).C*pow(x2,n);
        _B[0].B += PMPbar.B(n+1).B*pow(x2,n);
    }

    // B(k > 0)

    for (int k=1; k<=kmax; ++k) {
        _B[k].C += -PbarMP.B(k).C - PMPbar.B(k).C + PbarMP.B(k-1).C - PbarMP.B(k).C*x2;
        _B[k].B -= PMPbar.B(k).B;
        _B[k].E += -PbarMP.B(k).E + PbarMP.B(k-1).E - PbarMP.B(k).E*x2;

        for (int n=0; k+n+1 <= kmax; ++n) {
            _B[k].C += PMPbar.B(k+n+1).C*pow(x2,n);
            _B[k].B += PMPbar.B(k+n+1).B*pow(x2,n);
        }
    }

    if (!PbarMP.B(kmax).C.isZero() || !PbarMP.B(kmax).E.isZero()) {
        _B[kmax+1].A = nullMatrix.A;
        _B[kmax+1].B = nullMatrix.B;
        _B[kmax+1].C = PbarMP.B(kmax).C;
        _B[kmax+1].D = nullMatrix.D;
        _B[kmax+1].E = PbarMP.B(kmax).E;
        _B[kmax+1].F = nullMatrix.F;
    }
}

void System::transform(const FermatArray &T) {
    FermatArray Tinv = T.inverse();

    for (auto it = _A.begin(); it != _A.end(); ++it) {
        it->second.B = Tinv * it->second.B;
        it->second.C = Tinv * it->second.C * T;
        it->second.E = it->second.E * T;
    }
    
    for (auto it = _B.begin(); it != _B.end(); ++it) {
        it->second.B = Tinv * it->second.B;
        it->second.C = Tinv * it->second.C * T;
        it->second.E = it->second.E * T;
    }
    
    tqueue.transform(T);
}

void System::lefttransform(const FermatArray &G, const FermatExpression &x1, int k) {
    if (x1 == infinity) {
        lefttransform_inf(G,k);
        return;
    }

    sing_t sing;

    //B
   
    for (auto &s : singularities) {
        FermatExpression xj = s.first;
        if (xj == x1 || xj == infinity) continue;

        for (int n=0; n<=s.second.rank; ++n) {
            FermatArray mat = nullMatrix.B;

            for (int i=0; n+i <= s.second.rank; ++i) {
                mat += (A(xj,n+i).C*G - G*A(xj,n+i).A) * powi(-1,k)*binomi(k+i-1,i)/pow(x1-xj,k+i);
            }

            if (mat.isZero()) continue;

            sing.point = xj;
            sing.rank = n;

            if (!_A.count(sing)) _A[sing] = nullMatrix;

            _A[sing].B += mat;
        }
    }

    for (int n=0; n<k; ++n) {
        FermatArray mat = nullMatrix.B;

        for (auto &a : _A) {
            FermatExpression xj = a.first.point;
            int i = a.first.rank;

            if (xj == x1) continue;

            mat += (a.second.C*G - G*a.second.A) * powi(-1,i+1) * binomi(k+i-n-1,i) / pow(xj-x1,k+i-n);
        }

        for (int i=0; i+k-n-1 <= kmax; ++i) {
            mat += (B(i+k-n-1).C*G - G*B(i+k-n-1).A)*pow(x1,i)*binomi(i+k-n-1,i);
        }
        
        if (mat.isZero()) continue;

        sing.point = x1;
        sing.rank = n;

        if (!_A.count(sing)) _A[sing] = nullMatrix;

        _A[sing].B += mat;
    }

    sing.point = x1;
    sing.rank = k;

    _A[sing].B += A(x1,0).C*G - G*A(x1,0).A + G*k;

    for (int n=k+1; n-k <= singularities[x1].rank; ++n) {
        FermatArray mat = A(x1,n-k).C*G - G*A(x1,n-k).A;

        if (mat.isZero()) continue;

        sing.point = x1;
        sing.rank = n;

        if (!_A.count(sing)) _A[sing] = nullMatrix;

        _A[sing].B += mat;
    }

    for (int n=0; n<=kmax; ++n) {
        FermatArray mat = nullMatrix.B;

        for (int m=0; n+m+k <= kmax; ++m) {
            for (int i=0; i+n+m+k <= kmax; ++i) {
                mat += (B(i+n+m+k).C*G - G*B(i+n+m+k).A) * powi(-1,m) * pow(x1,m+i) * binomi(n+m,n) * binomi(i+n+m+k,i);
            }
        }            

        if (mat.isZero()) continue;

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].B += mat;
    }

    //D
    
    for (auto &s : singularities) {
        FermatExpression xj = s.first;
        if (xj == x1 || xj == infinity) continue;

        for (int n=0; n<=s.second.rank; ++n) {
            FermatArray mat = nullMatrix.D;

            for (int i=0; n+i <= s.second.rank; ++i) {
                mat += A(xj,n+i).E*G * powi(-1,k)*binomi(k+i-1,i)/pow(x1-xj,k+i);
            }

            if (mat.isZero()) continue;

            sing.point = xj;
            sing.rank = n;

            if (!_A.count(sing)) _A[sing] = nullMatrix;

            _A[sing].D += mat;
        }
    }
    
    for (int n=0; n<k; ++n) {
        FermatArray mat = nullMatrix.D;

        for (auto &a : _A) {
            FermatExpression xj = a.first.point;
            int i = a.first.rank;

            if (xj == x1) continue;

            mat += a.second.E*G * powi(-1,i+1) * binomi(k+i-n-1,i) / pow(xj-x1,k+i-n);
        }

        for (int i=0; i+k-n-1 <= kmax; ++i) {
            mat += B(i+k-n-1).E*G * pow(x1,i) * binomi(i+k-n-1,i);
        }
        
        if (mat.isZero()) continue;

        sing.point = x1;
        sing.rank = n;

        if (!_A.count(sing)) _A[sing] = nullMatrix;

        _A[sing].D += mat;
    }
    
    for (int n=k; n-k <= singularities[x1].rank; ++n) {
        FermatArray mat = A(x1,n-k).E*G;

        if (mat.isZero()) continue;

        sing.point = x1;
        sing.rank = n;

        if (!_A.count(sing)) _A[sing] = nullMatrix;

        _A[sing].D += mat;
    }
    
    for (int n=0; n<=kmax; ++n) {
        FermatArray mat = nullMatrix.D;

        for (int m=0; n+m+k <= kmax; ++m) {
            for (int i=0; i+n+m+k <= kmax; ++i) {
                mat += B(i+n+m+k).E*G * powi(-1,m) * pow(x1,m+i) * binomi(n+m,n) * binomi(i+n+m+k,i);
            }
        }            

        if (mat.isZero()) continue;

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].D += mat;
    }

    updatePoincareRanks();
    tqueue.lefttransform(G,x1,k);
}

void System::lefttransform_inf(const FermatArray &G, int k) {
    //B
    
    for (auto &s : singularities) {
        FermatExpression xj = s.first;

        if (xj == infinity) continue;

        for (int n=0; n<=s.second.rank; ++n) {
            FermatArray mat = nullMatrix.B;

            for (int i=0;  i<=k; ++i) {
                mat += (A(xj,n+k-i).C*G - G*A(xj,n+k-i).A) * pow(xj,i) * binomi(k,k-i);
            }
            if (mat.isZero()) continue;

            sing_t sing;

            sing.point = xj;
            sing.rank = n;

            if (!_A.count(sing)) _A[sing] = nullMatrix;

            _A[sing].B += mat;
        }        
    }

    for (int n=0; n<=k-1; ++n) {
        FermatArray mat = nullMatrix.B;

        for (auto &s : singularities) {
            FermatExpression xj = s.first;
            if (xj == infinity) continue;

            for (int m=0; m<=k-n-1; ++m) {
                for (int i=0; i<=m; ++i) {
                    mat += (A(xj,i).C*G - G*A(xj,i).A) * powi(-1,k-n-m-1) * pow(xj,k-n-i-1) * binomi(k-m-1,n) * binomi(k,k+i-m);
                }
            }
        }

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].B += mat;
    }

    if (!_B.count(k-1)) _B[k-1] = nullMatrix;
    _B[k-1].B -= G*k;

    for (int n=k; n-k <= kmax; ++n) {
        FermatArray mat = B(n-k).C*G - G*B(n-k).A;

        if (mat.isZero()) continue;

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].B += mat;
    }

    //D
    
    for (auto &s : singularities) {
        FermatExpression xj = s.first;

        if (xj == infinity) continue;

        for (int n=0; n<=s.second.rank; ++n) {
            FermatArray mat = nullMatrix.D;

            for (int i=0;  i<=k; ++i) {
                mat += A(xj,n+k-i).E*G * pow(xj,i) * binomi(k,k-i);
            }
            if (mat.isZero()) continue;

            sing_t sing;

            sing.point = xj;
            sing.rank = n;

            if (!_A.count(sing)) _A[sing] = nullMatrix;

            _A[sing].D += mat;
        }        
    }

    for (int n=0; n<=k-1; ++n) {
        FermatArray mat = nullMatrix.D;

        for (auto &s : singularities) {
            FermatExpression xj = s.first;

            if (xj == infinity) continue;

            for (int m=0; m<=k-n-1; ++m) {
                for (int i=0; i<=m; ++i) {
                    mat += A(xj,i).E*G * powi(-1,k-n-m-1) * pow(xj,k-n-i-1) * binomi(k-m-1,n) * binomi(k,k+i-m);
                }
            }
        }

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].D += mat;
    }

    for (int n=k; n-k <= kmax; ++n) {
        FermatArray mat = B(n-k).E*G;

        if (mat.isZero()) continue;

        if (!_B.count(n)) _B[n] = nullMatrix;

        _B[n].D += mat;
    }
    
    updatePoincareRanks();
    tqueue.lefttransform(G,infinity,k);
}

void System::lefttransformFull(const FermatArray &G, const FermatExpression &x1, int k) {
    sing_t sing;

    if (!(G*G).isZero()) {
        throw invalid_argument("G^2 must be zero.");
    }

    if (x1 == infinity) {
        lefttransformFull_inf(G,k);
        return;
    }

    for (auto it = _A.begin(); it != _A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;
        if (xj == x1) continue;

        for (int i=0; n+i <= singularities[xj].rank; ++i) {
            it->second.C += (A(xj,n+i).C*G - G*A(xj,n+i).C)*powi(-1,k)*binomi(k+i-1,i)/pow(x1-xj,k+i);
        }
    }

    for (int n=0; n<k; ++n) {
        sing.point = x1;
        sing.rank = n;

        for (auto it=_A.begin(); it != _A.end(); ++it) {
            FermatExpression xj = it->first.point;
            int i = it->first.rank;
            if (xj == x1) continue;

            _A[sing].C += (it->second.C*G - G*it->second.C)*powi(-1,i+1)*binomi(k+i-n-1,i)/pow(xj-x1,k+i-n);
        }

        for (int i=0; i+k-n-1<=kmax; ++i) {
            _A[sing].C += (B(i+k-n-1).C*G - G*B(i+k-n-1).C)*pow(x1,i)*binomi(i+k-n-1,i);
        }
    }

    sing.point = x1;
    sing.rank = k;
    _A[sing].C += G*k;

    for (int n=k; n-k <= singularities[x1].rank; ++n) {
        sing.point = x1;
        sing.rank = n;

        if (!_A.count(sing)) _A[sing] = nullMatrix;

        _A[sing].C += A(x1,n-k).C*G - G*A(x1,n-k).C;
    }

    for (auto it=_B.begin(); it != _B.end(); ++it) {
        int n=it->first;

        for (int m=0; n+m+k<=kmax; ++m) {
            for (int i=0; i+n+m+k<=kmax; ++i) {
                _B[n].C += (B(i+n+m+k).C*G - G*B(i+n+m+k).C)*powi(-1,m)*pow(x1,m+i)*binomi(n+m,n)*binomi(i+n+m+k,i);
            }
        }
    }

    updatePoincareRanks();
}

void System::lefttransformFull_inf(const FermatArray &G, int k) {
    for (auto it = _A.begin(); it != _A.end(); ++it) {
        FermatExpression xj = it->first.point;
        int n = it->first.rank;

        for (int i=0; i<=k; ++i) {
            it->second.C += (A(xj,n+k-i).C*G - G*A(xj,n+k-i).C)*pow(xj,i)*binomi(k,k-i);
        }
    }

    for (int n=0; n<k-1; ++n) {
        for (auto it=singularities.begin(); it != singularities.end(); ++it) {
            FermatExpression xj = it->first;
            if (xj == infinity) continue;

            for (int m=0; m<=k-n-1; ++m) {
                for (int i=0; i<=m; ++i) {
                    _B[n].C += (A(xj,i).C*G - G*A(xj,i).C)*powi(-1,k-n-m-1)*pow(xj,k-n-i-1)*binomi(k-m-1,n)*binomi(k,k+i-m);
                }
            }
        }
    }

    _B[k-1].C -= G*k;
    for (auto it=singularities.begin(); it != singularities.end(); ++it) {
        FermatExpression xj = it->first;
        if (xj == infinity) continue;

        _B[k-1].C += A(xj,0).C*G - G*A(xj,0).C;
    }

    for (int n=k; n-k<=kmax; ++n) {
        if (!_B.count(n)) _B[n] = nullMatrix;
        _B[n].C += B(n-k).C*G - G*B(n-k).C;
    }
    
    updatePoincareRanks();
}

void System::updatePoincareRanks() {
    singularities.clear();
    kmaxC = kmax = -1;

    for (auto it = _A.begin(); it != _A.end(); ++it) {
        if (!it->second.C.isZero()) {
            if (singularities.count(it->first.point)) {
                singularities[it->first.point].rankC = max(singularities[it->first.point].rankC,it->first.rank);
            } else {
                singularities[it->first.point].rankC = it->first.rank;
                singularities[it->first.point].rank = -1;
            }
        }
        
        if (!it->second.A.isZero() || !it->second.B.isZero() || !it->second.C.isZero() || !it->second.D.isZero() || ! it->second.E.isZero() || !it->second.F.isZero()) {
            if (singularities.count(it->first.point)) {
                singularities[it->first.point].rank = max(singularities[it->first.point].rank,it->first.rank);
            } else {
                singularities[it->first.point].rankC = -1;
                singularities[it->first.point].rank = it->first.rank;
            }
        }
    }

    for (auto it = _B.begin(); it != _B.end(); ++it) {
        if (!it->second.C.isZero()) {
            if (it->first>kmaxC) kmaxC=it->first;
        }
        if (!it->second.A.isZero() || !it->second.B.isZero() || !it->second.C.isZero() || !it->second.D.isZero() || ! it->second.E.isZero() || !it->second.F.isZero()) {
            if (it->first>kmax) kmax=it->first;
        }
    }
    
    singularities[infinity].rankC = -1;
    singularities[infinity].rank = -1;

    if (kmaxC<0) {
        if (!A(infinity,0).C.isZero()) {
            singularities[infinity].rankC = 0;
        } else {
            singularities[infinity].rankC = -1;
        }
    } else {
        singularities[infinity].rankC = kmaxC+1;
    }
    
    if (kmax<0) {
        TriangleBlockMatrix mat = A(infinity,0);
        if (!mat.A.isZero() || !mat.B.isZero() || !mat.C.isZero() || !mat.D.isZero() || ! mat.E.isZero() || !mat.F.isZero()) {
            singularities[infinity].rank = 0;
        } else {
            singularities[infinity].rank = -1;
        }
    } else {
        singularities[infinity].rank = kmax+1;
    }
    
    if (singularities[infinity].rank < 0) {
        singularities.erase(infinity);
    }
}

FermatExpression System::regularPoint() {
    if (!singularities.count(infinity)) {
        return infinity;
    }

    if (!singularities.count(FermatExpression(fermat,"0"))) {
        return FermatExpression(fermat,"0");
    }

    for (int n=1; n<100; ++n) {
        stringstream strm;
        strm << n;
        FermatExpression fn(fermat,strm.str());

        if (!singularities.count(fn)) {
            return fn;
        }
        if (!singularities.count(-fn)) {
            return -fn;
        }
    }        

    throw invalid_argument("no regular point found");
}

void System::jordan(const FermatExpression &xj) {
    if (jordans.count(xj)) return;
    eigen(xj);
 
    FermatArray C = singularities.count(xj) ? A(xj,singularities.at(xj).rankC).C : nullMatrix.C;
    jordanSystem(C,eigenvalues[xj],jordans[xj]);
}

void System::eigen(const FermatExpression &xj) {
    if (eigenvalues.count(xj)) return;
    FermatArray C = singularities.count(xj) ? A(xj,singularities.at(xj).rankC).C : nullMatrix.C;
    eigenvalues[xj] = findEigenvalues(C,100);
}

  
void System::inverseJordan(const FermatExpression &xj, list<JordanBlock> &inv) {
    FermatArray U(fermat,nullMatrix.C.rows(),nullMatrix.C.cols());
    FermatArray V(fermat);
    int i;

    inv.clear();

    jordan(xj);

    i = 1;
    for (auto &b : jordans[xj]) {
        for (auto &v : b.rootvectors) {
            U.setRow(i++,v);
        }
    }

    if (i != U.rows()+1) throw invalid_argument("wrong number of root vectors");

    V.assign("1/["+U.name()+"]");

    i = 1;
    for (auto &b : jordans[xj]) {
        JordanBlock block;
        block.ev = b.ev;

        for (int j=0; j<b.rootvectors.size(); ++j) {
            FermatArray v(V,1,V.rows(),i,i);
            i++;

            block.rootvectors.push_front(v);
        }

        inv.push_back(block);
    }
}

System::TriangleBlockMatrix System::A(const FermatExpression &xj, int k) const {
    if (xj == infinity) return Ainf(k);

    sing_t sing;

    sing.point = xj;
    sing.rank = k;

    if (!_A.count(sing)) {
        return nullMatrix;
    }

    return _A.at(sing);
}

System::TriangleBlockMatrix System::B(int k) const {
    if (!_B.count(k)) {
        return nullMatrix;
    }

    return _B.at(k);
}

System::TriangleBlockMatrix System::Ainf(int k) const {
    if (k == 0) {
        TriangleBlockMatrix mat=nullMatrix;

        for (auto it = singularities.begin(); it != singularities.end(); ++it) {
            if (it->first == infinity) continue;

            TriangleBlockMatrix Aj0 = A(it->first,0);

            mat.A -= Aj0.A;
            mat.B -= Aj0.B;
            mat.C -= Aj0.C;
            mat.D -= Aj0.D;
            mat.E -= Aj0.E;
            mat.F -= Aj0.F;
        }

        return mat;
    } else {
        TriangleBlockMatrix mat=B(k-1);

        mat.A *= -1;
        mat.B *= -1;
        mat.C *= -1;
        mat.D *= -1;
        mat.E *= -1;
        mat.F *= -1;

        return mat;
    }
}

FermatArray System::putTogether(const TriangleBlockMatrix &A) const {
    FermatArray B(fermat,A.A.rows()+A.C.rows()+A.F.rows(),A.A.cols()+A.C.cols()+A.F.cols());
    stringstream strm;

    B.assign("0");

    if (A.A.cols() > 0) {
        strm.str("");
        strm.clear();
        strm << "[" << B.name() << "[1~" << A.A.rows() << ",1~" << A.A.cols() << "]] := [" << A.A.name() << "]";
    
        (*fermat)(strm.str());
        
        strm.str("");
        strm.clear();
        strm << "[" << B.name() << "[" << A.A.rows()+1 << "~" << A.A.rows()+A.B.rows() << ",1~" << A.B.cols() << "]] := [" << A.B.name() << "]";
    
        (*fermat)(strm.str());

        if (A.D.rows() > 0) {
            strm.str("");
            strm.clear();
            strm << "[" << B.name() << "[" << A.A.rows()+A.B.rows()+1 << "~" << A.A.rows()+A.B.rows()+A.D.rows() << ",1~" << A.D.cols() << "]] := [" << A.D.name() << "]";
        
            (*fermat)(strm.str());
        }
    }
        
    strm.str("");
    strm.clear();
    strm << "[" << B.name() << "[" << A.A.rows()+1 << "~" << A.A.rows()+A.C.rows() << "," << A.B.cols()+1 << "~" << A.B.cols()+A.C.cols() << "]] := [" << A.C.name() << "]";

    (*fermat)(strm.str());

    if (A.F.rows() > 0) {
        strm.str("");
        strm.clear();
        strm << "[" << B.name() << "[" << A.A.rows()+A.C.rows()+1 << "~" << A.A.rows()+A.C.rows()+A.E.rows() << "," << A.D.cols()+1 << "~" << A.D.cols()+A.E.cols() << "]] := [" << A.E.name() << "]";
    
        (*fermat)(strm.str());
        
        strm.str("");
        strm.clear();
        strm << "[" << B.name() << "[" << A.A.rows()+A.C.rows()+1 << "~" << A.A.rows()+A.C.rows()+A.F.rows() << "," << A.D.cols()+A.E.cols()+1 << "~" << A.D.cols()+A.E.cols()+A.F.cols() << "]] := [" << A.F.name() << "]";
    
        (*fermat)(strm.str());
    }

    return B;
}

void System::printSingularities() const {
    cout << "singularities:  ";
    for (auto it=singularities.begin(); it != singularities.end(); ++it) {
        if (it->second.rankC < 0) continue;

        cout << "  " << pstr(it->first) << ":" << it->second.rankC;
    }
    cout << endl;
}

void System::printEigenvalues() {
    int width=0;

    for (auto &s : singularities) {
        int len = pstr(s.first).size(); 
        if (len > width) width = len;
    }

    for (auto it = singularities.begin(); it != singularities.end(); ++it) {
        FermatExpression xj = it->first;

        if (it->second.rankC != 0) continue;

        if (it == singularities.begin()) {
            cout << "eigenvalues:  ";
        } else {
            cout << "              ";
        }
        cout << "[" << setw(width) << right << pstr(xj) << "]     ";

        eigen(xj);

        for (auto &ev0 : eigenvalues[xj]) {
            eigen_t ev = ev0.first;
            int mult = ev0.second;

            if (ev.u == 0) {
                stringstream strm;
                
                if (ev.v > 1 || ev.v < -1) {
                    strm << ev.v << "*ep";
                } else if (ev.v == 1) {
                    strm << "ep";
                } else if (ev.v == -1) {
                    strm << "-ep";
                } else {
                    strm << 0;
                }

                cout << setw(10) << right << strm.str() << ":" << mult << "\t";
            }
        }
        
        for (auto &ev0 : eigenvalues[xj]) {
            eigen_t ev = ev0.first;
            int mult = ev0.second;

            if (ev.u != 0) {
                stringstream strm;

                if (ev.u != 0 && ev.v > 1) {
                    strm << ev.u << "+" << ev.v << "*ep";
                } else if (ev.u != 0 && ev.v == 1) {
                    strm << ev.u << "+ep";
                } else if (ev.u != 0 && ev.v < -1) {
                    strm << ev.u << ev.v << "*ep";
                } else if (ev.u != 0 && ev.v == -1) {
                    strm << ev.u << "-ep";
                } else if (ev.u != 0 && ev.v == 0) {
                    strm << ev.u;
                } 

                cout << setw(10) << right << strm.str() << ":" << mult << "\t";
            }                
        }
        cout << endl;
    }
}

string System::pstr(const FermatExpression &x) const {
    if (x == infinity) {
        return "inf";
    } else {
        return x.str();
    }
}

FermatExpression System::powi(int b, int e) const {
    stringstream strm;
    strm << "(" << b << ")^(" << e << ")";
    
    return FermatExpression(fermat,strm.str());
}

FermatExpression System::pow(const FermatExpression &b, int e) const {
    stringstream strm;
    strm << "(" << b.name() << ")^(" << e << ")";
    
    return FermatExpression(fermat,strm.str());
}

FermatExpression System::binomi(int n, int k) const {
    stringstream strm;
    strm << "Bin(" << n << "," << k << ")";
    
    return FermatExpression(fermat,strm.str());
}



