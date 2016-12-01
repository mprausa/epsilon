// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/factor.h
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

#ifndef __FACTOR_H
#define __FACTOR_H

#include <ginac/ginac.h>

std::map<GiNaC::numeric,int> factor(GiNaC::numeric n);

#endif //__FACTOR_H

