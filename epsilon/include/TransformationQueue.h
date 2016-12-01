// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/TransformationQueue.h
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

#ifndef __TRANSFORMATION_QUEUE_H
#define __TRANSFORMATION_QUEUE_H

#include <string>
#include <fstream>
#include <list>
#include <FermatArray.h>

class System;

class TransformationQueue {
    protected:
        std::ofstream file;
        std::string _filename;
        int before;
        int after;
        Fermat *fermat;
        bool replaying;

        typedef struct {
            enum {Balance, Transformation, LeftTrans} type;
            FermatExpression x1;
            FermatExpression x2;
            int k;
            FermatArray T;
        } transformation_t;

        std::list<transformation_t> queue;
    public:
        TransformationQueue(const TransformationQueue &other);
        TransformationQueue(Fermat *fermat);
        ~TransformationQueue();

        void setpadding(int before, int after);
        void setfile(std::string _filename, bool append=false);
        std::string filename();
        void load(std::string _filename);

        void replay(System &system);
        void exporttrans(std::string filename);

        void balance(const FermatArray &P, const FermatExpression &x1, const FermatExpression &x2);
        void transform(const FermatArray &T);
        void lefttransform(const FermatArray &G, const FermatExpression &x1, int k);
    private:
        std::string pstr(const FermatExpression &x) const;
};

#endif //__TRANSFORMATION_QUEUE_H

