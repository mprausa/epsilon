// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/HalfInteger.h
 *
 *  Copyright (C) 2017, 2018 Mario Prausa
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

#ifndef __HALF_INTEGER_H
#define __HALF_INTEGER_H

#include <iostream>

class HalfInteger {
    protected:
    	int value;
    public:
        HalfInteger();
        HalfInteger(int i);
        HalfInteger(int i, int);

        bool operator<(const HalfInteger &other) const;
        bool operator<=(const HalfInteger &other) const;
        bool operator>(const HalfInteger &other) const;
        bool operator>=(const HalfInteger &other) const;
        bool operator==(const HalfInteger &other) const;
        bool operator!=(const HalfInteger &other) const;

        bool operator<(int other) const;
        bool operator<=(int other) const;
        bool operator>(int other) const;
        bool operator>=(int other) const;
        bool operator==(int other) const;
        bool operator!=(int other) const;

        HalfInteger operator-(const HalfInteger &other) const;
        HalfInteger operator-() const;

        HalfInteger &operator=(int i);

        void inc(bool half);
        void dec(bool half);

        bool isInteger() const;

        void print(std::ostream &os) const;
};

std::ostream &operator<<(std::ostream &os, const HalfInteger &hi);

extern const HalfInteger onehalf;

#endif //__HALF_INTEGER_H

