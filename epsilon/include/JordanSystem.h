// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/JordanSystem.h
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

#ifndef __JORDAN_SYSTEM_H
#define __JORDAN_SYSTEM_H

#include <Eigenvalues.h>
#include <FermatArray.h>
#include <deque>
#include <set>

typedef struct _jordanBlock {
    eigen_t ev;
    std::deque<FermatArray> rootvectors;

    bool operator<(const struct _jordanBlock &other) const {
        return (rootvectors.size() > other.rootvectors.size()) || (rootvectors.size() == other.rootvectors.size() && ev < other.ev);
    }
} JordanBlock;

typedef std::multiset<JordanBlock> JordanSystem;

void jordanSystem(const FermatArray &A, const eigenvalues_t &evs, JordanSystem &system);
void Eigenvectors(const FermatArray &A, eigen_t ev, std::vector<FermatArray> &vectors);

#endif //__JORDAN_SYSTEM_H
