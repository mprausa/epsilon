// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  include/Fraction.h
 *
 *  Copyright (C) 2019 Mario Prausa
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

#ifndef __FRACTION_H
#define __FRACTION_H

#include <iostream>
#include <string>
#include <cassert>

class Fraction;
std::ostream &operator<<(std::ostream &os, const Fraction &f);

class Fraction {
    friend std::ostream &operator<<(std::ostream &os, const Fraction &f);
    protected:
        int _numer;
        int _denom;
    public:
        Fraction() : _numer(0), _denom(1) {}
        Fraction(int i) : _numer(i), _denom(1) {}
        Fraction(int n, int d) : _numer(n), _denom(d) {
            cancel();
        }

        Fraction(const std::string &str) {
            auto slash = str.find('/');

            if (slash == std::string::npos) {
                _numer = stoi(str);
                _denom = 1;
            } else {
                _numer = stoi(str.substr(0,slash));
                _denom = stoi(str.substr(slash+1));
                cancel();
            }
        }

        Fraction operator+(const Fraction &other) const {
            return Fraction(_numer*other._denom + other._numer*_denom,_denom*other._denom);
        }

        Fraction operator-(const Fraction &other) const {
            return Fraction(_numer*other._denom - other._numer*_denom,_denom*other._denom);
        }

        Fraction operator-() const {
            return Fraction(-_numer,_denom);
        }

        Fraction &operator=(int i) {
            _numer = i;
            _denom = 1;
            return *this;
        }

        Fraction &operator+=(const Fraction &other) {
            *this = *this + other;
            return *this;
        }

        Fraction &operator-=(const Fraction &other) {
            *this = *this - other;
            return *this;
        }

        bool operator<(const Fraction &other) const {
            return _numer*other._denom < other._numer*_denom;
        }

        bool operator<=(const Fraction &other) const {
            return _numer*other._denom <= other._numer*_denom;
        }

        bool operator>(const Fraction &other) const {
            return _numer*other._denom > other._numer*_denom;
        }

        bool operator>=(const Fraction &other) const {
            return _numer*other._denom >= other._numer*_denom;
        }

        bool operator==(const Fraction &other) const {
            return _numer == other._numer && _denom == other._denom;
        }

        bool operator!=(const Fraction &other) const {
            return _numer != other._numer || _denom != other._denom;
        }

        bool is_int() const {
            return _denom == 1;
        }

        int to_int() const {
            assert(is_int());
            return _numer;
        }

    private:
        void cancel() {
            assert(_denom != 0);

            if (_denom < 0) {
                _numer = -_numer;
                _denom = -_denom;
            }

            int a = std::abs(_numer);
            int b = std::abs(_denom);
            for (int c=a; b; c=a) {
               c = a;
               a = b;
               b = c%b;
            }

            _numer /= a;
            _denom /= a;
        }
};

inline std::ostream &operator<<(std::ostream &os, const Fraction &f) {
    if (f.is_int()) {
        os << f._numer;
    } else {
        os << f._numer << "/" << f._denom;
    }
    return os;
}

#endif //__FRACTION_H

