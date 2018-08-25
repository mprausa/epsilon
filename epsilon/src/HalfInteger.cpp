// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/HalfInteger.cpp
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

#include <HalfInteger.h>
#include <sstream>
using namespace std;

const HalfInteger onehalf(1,2);

HalfInteger::HalfInteger() {
    value = 0;
}

HalfInteger::HalfInteger(int i) {
	value = 2*i;
}

HalfInteger::HalfInteger(int i, int) {
	value = i;
}

bool HalfInteger::operator<(const HalfInteger &other) const {
	return value < other.value;
}

bool HalfInteger::operator<=(const HalfInteger &other) const {
	return value <= other.value;
}

bool HalfInteger::operator>(const HalfInteger &other) const {
	return value > other.value;
}

bool HalfInteger::operator>=(const HalfInteger &other) const {
	return value >= other.value;
}

bool HalfInteger::operator==(const HalfInteger &other) const {
	return value == other.value;
}

bool HalfInteger::operator!=(const HalfInteger &other) const {
	return value != other.value;
}

bool HalfInteger::operator<(int other) const {
	return value < 2*other;
}

bool HalfInteger::operator<=(int other) const {
	return value <= 2*other;
}

bool HalfInteger::operator>(int other) const {
	return value > 2*other;
}

bool HalfInteger::operator>=(int other) const {
	return value >= 2*other;
}

bool HalfInteger::operator==(int other) const {
	return value == 2*other;
}

bool HalfInteger::operator!=(int other) const {
	return value != 2*other;
}

HalfInteger HalfInteger::operator-(const HalfInteger &other) const {
    HalfInteger hi;

    hi.value = value - other.value;
    return hi;
}

HalfInteger HalfInteger::operator-() const {
    HalfInteger hi;

    hi.value = -value;
    return hi;
}

HalfInteger &HalfInteger::operator=(int i) {
    value = 2*i;
    return *this;
}

void HalfInteger::inc(bool half) {
    value += half?1:2;
}

void HalfInteger::dec(bool half) {
    value -= half?1:2;
}

bool HalfInteger::isInteger() const {
    return !(value & 1);
}

void HalfInteger::print(ostream &os) const {
    stringstream strm;
    if (value & 1) {
        strm << value << "/2";
    } else {
        strm << value/2;
    }
    os << strm.str();
}

ostream &operator<<(ostream &os, const HalfInteger &hi) {
    hi.print(os);
    return os;
}


