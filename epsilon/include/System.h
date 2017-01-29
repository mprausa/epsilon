// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/System.h
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

#ifndef __SYSTEM_H
#define __SYSTEM_H

#include <string>
#include <map>
#include <set>
#include <list>
#include <climits>
#include <JordanSystem.h>
#include <FermatArray.h>
#include <TransformationQueue.h>

extern FermatExpression infinity;

class System {
    protected:
        typedef struct _sing {
            FermatExpression point;
            int rank;
            
            bool operator<(const _sing &other) const {
                return (point < other.point) || (point == other.point && rank < other.rank);
            }
        } sing_t; 

        typedef struct {
            // ( A 0 0 )
            // ( B C 0 )
            // ( D E F )
            FermatArray A,B,C,D,E,F;
        } TriangleBlockMatrix;

        typedef struct {
            int rank;
            int rankC;
        } poincareRank;

        Fermat *fermat;
        TriangleBlockMatrix nullMatrix;

        std::map<sing_t,TriangleBlockMatrix> _A;
        std::map<int,TriangleBlockMatrix> _B;
       
        struct singLess {
            bool operator() (const FermatExpression &a, const FermatExpression &b) const;
        };
 
        std::map<FermatExpression,poincareRank,singLess> singularities;
        int kmax;
        int kmaxC;

        std::map<FermatExpression,eigenvalues_t> eigenvalues;
        std::map<FermatExpression,JordanSystem> jordans;

        TransformationQueue tqueue;
        bool echfer;
    public:
        System(Fermat *fermat, bool echfer);
        System(Fermat *fermat, std::string filename, int start, int end, bool echfer); 
        System(const System &orig, int start, int end);
        System(const System &orig, const FermatArray &left, const FermatArray &right);
        System(const System &orig, int ep);

        Fermat *fer() const;
        int dimC() const;

        void write(std::string filename) const;
        TransformationQueue *transformationQueue();

        void analyze();
        void fuchsify();
        void normalize();
        void factorep();
        void factorep(int mu);
        void leftranks();
        int leftreduce(const FermatExpression &xj);
        void leftfuchsify();

        void balance(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2);
        void transform(const FermatArray &T); 
        void lefttransformFull(const FermatArray &G, const FermatExpression &x1, int k);

        std::map<FermatExpression,FermatArray> exportFuchs() const;
        
        void printEigenvalues();
    private:
        bool projectorQ(const FermatExpression &x1, const FermatExpression &x2, FermatArray &Q);
        void projectorP(const FermatExpression &x1, FermatArray &P);

        int reduceL0(FermatArray L0, int k, const FermatExpression &x1, std::set<int> &S, FermatArray &Delta);
        bool invariantSubspace(const FermatExpression &x2, const FermatArray &Uk, FermatArray &Vk);
        bool findBalance(FermatExpression &x1, FermatExpression &x2, FermatArray &P, const FermatExpression &x0);
        FermatExpression regularPoint();

        void balance_x1_x2(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2);
        void balance_inf_x2(const FermatArray &P, const FermatExpression &x2);
        void balance_x1_inf(const FermatArray &P, const FermatExpression &x1);

        void lefttransform(const FermatArray &G, const FermatExpression &x1, int k);
        void lefttransform_inf(const FermatArray &G, int k);

        void lefttransformFull_inf(const FermatArray &G, int k);

        void updatePoincareRanks();
        void jordan(const FermatExpression &xj);
        void eigen(const FermatExpression &xj);
        void inverseJordan(const FermatExpression &xj, std::list<JordanBlock> &inv);

        TriangleBlockMatrix A(const FermatExpression &xj, int k) const;
        TriangleBlockMatrix B(int k) const;
        TriangleBlockMatrix Ainf(int k) const;

        FermatArray putTogether(const TriangleBlockMatrix &A) const;

        std::string pstr(const FermatExpression &x) const;

        FermatExpression powi(int b, int e) const;
        FermatExpression pow(const FermatExpression &b, int e) const;
        FermatExpression binomi(int n, int k) const;

        void printSingularities() const;
};

#endif //__SYSTEM_H

