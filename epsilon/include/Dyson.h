// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/Dyson.h
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

#ifndef __DYSON_H
#define __DYSON_H

#include <System.h>
#include <unordered_map>

class Dyson {
    public:
        typedef enum {tGPL,tHPL,tHPLalt} pltype_t;
        typedef enum {fMma, fForm} format_t;
    protected:
        class GPL {
            protected:
                std::string _arg;
                std::deque<FermatExpression> _indices;  
                bool isHPL;
            public:
                GPL();
                GPL(const GPL &other);

                void addIndex(const FermatExpression &xj);
                void setArg(std::string arg);

                std::string arg() const;
                std::deque<FermatExpression> indices() const;

                bool operator<(const GPL &other) const;
                bool operator==(const GPL &other) const;

                std::string str(pltype_t type, format_t format) const;
                int sign(pltype_t type) const;
            private:
                std::deque<FermatExpression> indsHPLalt() const;
        };

        class Term {
            protected:
                GPL _xGPL;
                std::map<GPL,int> _x0GPL;
                std::size_t _hash;
            public:
                Term(const GPL &xgpl, const std::map<GPL,int> &x0gpl);

                const GPL &xGPL() const;
                const std::map<GPL,int> &x0GPL() const;

                int sign(pltype_t type) const;
                std::size_t hash() const;

                bool operator==(const Term &other) const;
            private:
                void calchash();
        };
        
        struct _exHash : std::unary_function<Term,size_t> {
            std::size_t operator()(const Term &t) const {
                return t.hash();
            }
        };

        class Expression : public std::unordered_map<Term,FermatExpression,_exHash> {
            public:
                Expression();
                virtual ~Expression();

                Expression &operator+=(const Expression &other);
                Expression operator*(const FermatExpression &factor) const;

                Expression integrate(const FermatExpression &xj);
                std::string str(pltype_t type, format_t format);
        };
       
        class ExMatrix {
            protected:
                int dim;
                std::vector<Expression> data;
            public:
                ExMatrix(int dim);
                virtual ~ExMatrix();

                ExMatrix &operator+=(const ExMatrix &other);
                ExMatrix integrate(const FermatExpression &xj);

                ExMatrix lmul(const FermatArray &left) const;

                Expression &operator() (int r, int c);
                const Expression &operator() (int r, int c) const;

                std::string str(pltype_t type);
        };
        
        int dim;        
        std::vector<std::pair<FermatExpression,FermatArray>> Mxj;
        std::vector<ExMatrix> Un;
    public:
        Dyson(const System &system);

        void expand(int order);
        void write(std::string filename, pltype_t type, format_t format);
    private:
        void writeMma(std::string filename, pltype_t type);
        void writeForm(std::string filename, pltype_t type);
};

#endif //__DYSON_H

