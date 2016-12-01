// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/Echelon.h
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

#ifndef __ECHELON_H
#define __ECHELON_H

#include <FermatExpression.h>
#include <FermatArray.h>
#include <vector>
#include <map>
#include <iterator>

class EchelonBase {
    public:
        class Row {
            protected:
                Fermat *fermat;
                std::vector<std::pair<int,FermatExpression>> data;
            public:
                Row(Fermat *fermat);

                Row operator-(const Row &other) const;
                Row operator*(const FermatExpression &f) const;
                FermatExpression operator[] (int c) const;

                bool empty() const;
                size_t size() const;
                int col1() const;
                int maxcol() const;

                void set(int n, const FermatExpression &ex);
                void clear();

                void normalize();

                std::vector<std::pair<int,FermatExpression>>::const_iterator begin() const;
                std::vector<std::pair<int,FermatExpression>>::const_iterator end() const;
        };

        class Iterator : public std::iterator<std::input_iterator_tag,int> {
            protected:
                //Echelon
                std::map<int,std::vector<Row>> *rows;
                std::map<int,std::vector<Row>>::iterator it1;
                std::vector<Row>::iterator it2;

                //EchelonFermat
                FermatArray *array;
                int rownum;
                bool valid;
                Row row;
            public:
                Iterator(Fermat *fermat, std::map<int,std::vector<Row>> *rows, const std::map<int,std::vector<Row>>::iterator &it1, const std::vector<Row>::iterator &it2);
                Iterator(FermatArray *array, int r);
                Iterator(const Iterator &it);
                Iterator &operator++();
                Iterator operator++(int);
                bool operator==(const Iterator &other) const;
                bool operator!=(const Iterator &other) const;
                Row &operator*();
        };

        virtual void set(const Row &row) = 0;
        virtual void set(const std::map<int,FermatExpression> &m) = 0;
        virtual int run() = 0; 

        virtual Iterator begin() = 0;
        virtual Iterator end() = 0;

        virtual void print(std::ostream &os) const = 0;
};

class Echelon : public EchelonBase {
    protected:
        Fermat *fermat;
        std::map<int,std::vector<Row>> rows;
    public:
        Echelon(Fermat *fermat);

        virtual void set(const Row &row);
        virtual void set(const std::map<int,FermatExpression> &m);
        virtual int run();

        virtual Iterator begin();
        virtual Iterator end();
        
        virtual void print(std::ostream &os) const;
    private:
        int findPivot(const std::vector<Row> &rows) const;
        void normalize(std::vector<Row> &rows);
};

class EchelonFermat : public EchelonBase {
    protected:
        Fermat *fermat;
        FermatArray array;
        int pos;
    public:
        EchelonFermat(Fermat *fermat, int rows, int cols);
    
        virtual void set(const Row &row);
        virtual void set(const std::map<int,FermatExpression> &m);
        virtual int run();

        virtual Iterator begin();
        virtual Iterator end();

        virtual void print(std::ostream &os) const;
};
    	
std::ostream &operator<<(std::ostream &os, const EchelonBase &e);

#endif //__ECHELON_H

