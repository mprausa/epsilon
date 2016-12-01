// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/prepare.h
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

#ifndef __PREPARE_H
#define __PREPARE_H

#include <apart.h>

void prepare(const GiNaC::matrix &mat, const GiNaC::symbol &x, std::map<apart_sing,GiNaC::matrix> &A, std::map<int,GiNaC::matrix> &B);

#endif //__PREPARE_H


