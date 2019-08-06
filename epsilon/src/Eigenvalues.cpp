// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/Eigenvalues.cpp
 *
 *  Copyright (C) 2016, 2019 Mario Prausa
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

#include <Eigenvalues.h>
#include <sstream>
#include <iostream>
using namespace std;

int EVdenom = 1;

static FermatExpression makeExpression(Fermat *fermat, Fraction u, Fraction v) {
    stringstream strm;

    strm << "t-(" << u << "+" << v << "*ep)";

    return FermatExpression(fermat,strm.str());
}


static bool checkEV(FermatExpression poly, Fraction u, Fraction v) {
    stringstream strm;
    strm << u << "+" << v << "*ep";

    FermatExpression expr(poly.fer(),strm.str());

    return poly.subst("t",expr).str() == "0";
}

eigenvalues_t findEigenvalues(const FermatArray &array, int max) {
	FermatExpression poly = array.chPoly();
    Fermat *fermat = array.fer();
	eigenvalues_t values;
    Fraction inc(1,EVdenom);
    int ctr=0;

	for (Fraction i=0; i<=max; i += inc) {
		for (Fraction j=-i; j<=i; j += inc) {
			eigen_t ev;
            if (ctr == array.rows()) return values;

			if (checkEV(poly,i,j)) {
				ev.u = i;
				ev.v = j;

				values[ev]++;
                ++ctr;

				poly = poly/makeExpression(fermat,i,j);

				j -= inc;
				continue;
			}
			if (i == 0) break;
			if (checkEV(poly,j,i)) {
				ev.u = j;
				ev.v = i;

				values[ev]++;
                ++ctr;

				poly = poly/makeExpression(fermat,j,i);

                j -= inc;
				continue;
			}
			if (checkEV(poly,-i,j)) {
				ev.u = -i;
				ev.v = j;

				values[ev]++;
                ++ctr;

				poly = poly/makeExpression(fermat,-i,j);

                j -= inc;
				continue;
			}
			if (checkEV(poly,j,-i)) {
				ev.u = j;
				ev.v = -i;

				values[ev]++;
                ++ctr;

				poly = poly/makeExpression(fermat,j,-i);

                j -= inc;
				continue;
			}

		}
	}

    throw runtime_error("unable to find all eigenvalues.\nmatrix was: "+array.str());
}
