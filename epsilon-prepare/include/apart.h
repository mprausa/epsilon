// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/apart.h
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

#ifndef __APART_H
#define __APART_H

#include <ginac/ginac.h>

struct apart_sing {
    GiNaC::ex xj;
    int k;

    bool operator<(const apart_sing &other) const {
        if (GiNaC::ex_is_less()(xj,other.xj)) return false;
        if (GiNaC::ex_is_less()(other.xj,xj)) return true;

        return k < other.k;
    }        
};                       

GiNaC::ex modout(GiNaC::ex expr);
void apart(const GiNaC::ex &expr, const GiNaC::symbol &x, std::map<apart_sing,GiNaC::ex> &coeffsA, std::map<int,GiNaC::ex> &coeffsB);

#endif //__APART_H


