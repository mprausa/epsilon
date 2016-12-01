// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/apart.cpp
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

#include <apart.h>
#include <factor.h>
#include <iostream>
#include <sstream>
using namespace std;
using namespace GiNaC;

#define mod4(n) (((uint32_t)n)&3)

/*
 * polymods:
 *
 *   i:         i^2 + 1
 *   sqrtN:     sqrtN^2 - N
 *   isqrtN:    isqrtN^2 + N
 *   rN:        rN^2 - rN + (1+N)/4
 *   qN:        qN^2 - qN + (1-N)/4
 */

static map<int,ex> sqrts;
static exmap polymods;
 
static ex sqrtex(int n) {
    if (sqrts.count(n)) return sqrts[n];

    if (n == -1) {
        symbol i("i");
        sqrts[-1] = i;
        polymods[i] = i*i + 1;
    } else if (n < 0) {
        if (mod4(n) == 1) {
            stringstream strm;
            strm << "r" << (-n);
            symbol rN(strm.str());

            sqrts[n] = 2*rN - 1;
            polymods[rN] = rN*rN - rN + (1-n)/4; 
        } else {
            stringstream strm;
            strm << "isqrt" << (-n);
            symbol isqrtN(strm.str());

            sqrts[n] = isqrtN;
            polymods[isqrtN] = isqrtN*isqrtN - n; 
        }
    } else {
        if (mod4(n) == 1) {
            stringstream strm;
            strm << "q" << n;
            symbol qN(strm.str());

            sqrts[n] = 2*qN - 1;
            polymods[qN] = qN*qN - qN + (1-n)/4; 
        } else {
            stringstream strm;
            strm << "sqrt" << n;
            symbol sqrtN(strm.str());

            sqrts[n] = sqrtN;
            polymods[sqrtN] = sqrtN*sqrtN - n; 
        }
    }

    return sqrts[n];
}

static ex mysqrt(const ex &n) {
    if (!is_a<numeric>(n.expand())) return sqrt(n);
    numeric i = ex_to<numeric>(n.expand());

    if (i.to_int() == 0) return ex(0);

    bool neg = i.to_int() < 0; 
    if (neg) {
        i = -i;
    }

    map<numeric,int> factors = factor(i);

    ex a=1,b=1;

    for (auto &f : factors) {
        int exp = f.second/2;

        a *= pow(f.first,exp);
        b *= pow(f.first,f.second-exp*2);
    }

    return a*sqrtex(ex_to<numeric>(b).to_int()*(neg?-1:1));    
}

static void solve2(const ex &expr, int mul, const symbol &x, map<ex,int,ex_is_less> &zeros) {
    ex a = expr.coeff(x,2);
    ex b = expr.coeff(x,1);
    ex c = expr.coeff(x,0);
    
    ex p = (b/a).normal();
    ex q = (c/a).normal();

    if (!is_a<numeric>(p) || !is_a<numeric>(q)) {
        throw invalid_argument("polynomial has symbol coefficients.");
    }

    ex z1 = ((-b - mysqrt(b*b - 4*a*c))/2/a).normal();
    ex z2 = ((-b + mysqrt(b*b - 4*a*c))/2/a).normal();

    zeros[z1] += mul;
    zeros[z2] += mul;
}

static void solve1(const ex &expr, int mul, const symbol &x, map<ex,int,ex_is_less> &zeros) {
    ex a = expr.coeff(x,1);
    ex b = expr.coeff(x,0);
    ex z = (-b/a).normal();

    if (!is_a<numeric>(z.evalf())) {
        throw string("singularity is not a number.");
    }

    zeros[z] += mul;
}

static void solve(const ex &expr, const symbol &x, map<ex,int,ex_is_less> &zeros) {
    ex e = factor(expr);
    vector<pair<ex,int>> factors;

    if (is_a<mul>(e)) {
        for (auto f : e) {
            if (is_a<power>(f)) {
                factors.push_back({f.op(0),ex_to<numeric>(f.op(1)).to_int()});
            } else {
                factors.push_back({f,1});
            }
        }
    } else if (is_a<power>(e)) {
        factors.push_back({e.op(0),ex_to<numeric>(e.op(1)).to_int()});
    } else {
        factors.push_back({e,1});
    }

    for (auto &f : factors) {
        int deg = f.first.degree(x);
        switch(deg) {
            case 0:
                break;
            case 1:
                solve1(f.first,f.second,x,zeros);
                break;
            case 2:
                solve2(f.first,f.second,x,zeros);
                break;
            default: {
                stringstream strm;
                strm << "polynomials of degree " << deg << " unsupported.";
                throw string(strm.str());
            }                
        }            
    }
}

static void polydiv(const ex &a, const ex &b, const ex &x, ex &q, ex &r) {
    if (b.is_zero())
        throw(std::overflow_error("polydiv: division by zero"));

    if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)) {
        q = a/b;
        r = 0;
        return;
    }

    // Polynomial long division
    r = a.expand();
    if (r.is_zero()) {
        q = r;
        return;
    }

    int bdeg = b.degree(x);
    int rdeg = r.degree(x);
    ex blcoeff = b.expand().coeff(x, bdeg);
    exvector v; v.reserve(std::max(rdeg - bdeg + 1, 0));
    while (rdeg >= bdeg) {
        ex term, rcoeff = r.coeff(x, rdeg);
        term = (rcoeff / blcoeff).normal();
        term *= pow(x, rdeg - bdeg);
        v.push_back(term);
        r = (r - term * b).normal();
        r = r.numer().expand()/r.denom().expand();

        if (r.is_zero()) break;
        rdeg = r.degree(x);
    }

    q = dynallocate<add>(v);
}

ex modout(ex expr) {
    for (auto &pmod : polymods) {
        ex q,r;
        polydiv(expr,pmod.second,pmod.first,q,r);
        expr = r;
    }

    return expr;
}    

static ex limit(const ex &expr, const symbol &x, const ex &x0) {
    ex num = expr.numer();
    ex den = expr.denom();

    for(;;) {
        ex num0 = modout(num.subs(x == x0).expand()).normal();
        ex den0 = modout(den.subs(x == x0).expand()).normal();

        if (!den0.is_zero()) {
            exmap m;
            num0 = num0.to_rational(m);
            den0 = den0.to_rational(m);
            return (num0/den0).normal().subs(m);
        }

        if (!num0.is_zero()) {
            throw overflow_error("division by zero");
        }

        num = num.diff(x).normal();
        den = den.diff(x).normal();
    }        

    return expr.subs(x == x0);
}

void apart(const ex &expr, const symbol &x, map<apart_sing,ex> &coeffsA, map<int,ex> &coeffsB) {
    ex num = expr.numer().expand();
    ex den = expr.denom().expand();
    map<ex,int,ex_is_less> zeros;
    
    ex b = 0;

    if (num.degree(x) >= den.degree(x)) {
        polydiv(num,den,x,b,num);
    }

    solve(den,x,zeros);

    coeffsA.clear();
    coeffsB.clear();

    for (auto &xj : zeros) {
        ex f = (pow(x-xj.first,xj.second)*num/den).normal();

        for (int j=0; j<xj.second; ++j) {
            coeffsA[{xj.first,xj.second-j-1}] = (limit(f,x,xj.first)/factorial(j)).normal();
            f = f.diff(x).normal();
        }
    }

    if (!b.is_zero()) {
        for (int n=0,deg=b.degree(x); n<=deg; ++n) {
            coeffsB[n] = b.coeff(x,n);
        }        
    }        
}


