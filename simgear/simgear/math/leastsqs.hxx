/**
 * \file leastsqs.hxx
 * Implements a simple linear least squares best fit routine.
 */

// Written by Curtis Olson, started September 1997.
//
// Copyright (C) 1997  Curtis L. Olson  - http://www.flightgear.org/~curt
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// $Id$


#ifndef _LEASTSQS_H
#define _LEASTSQS_H


#ifndef __cplusplus
# error This library requires C++
#endif


/**
Classical least squares fit:

\f[
    y = b_0 + b_1 * x
\f]

\f[
    b_1 = \frac{n * \sum_0^{i-1} (x_i*y_i) - \sum_0^{i-1} x_i* \sum_0^{i-1} y_i}
          {n*\sum_0^{i-1} x_i^2 - (\sum_0^{i-1} x_i)^2}
\f]

\f[
    b_0 = \frac{\sum_0^{i-1} y_i}{n} - b_1 * \frac{\sum_0^{i-1} x_i}{n}
\f]
*/
void least_squares(double *x, double *y, int n, double *m, double *b);


/**
 * Incrimentally update existing values with a new data point.
 */
void least_squares_update(double x, double y, double *m, double *b);


/**
  @return the least squares error:.
\f[

    \frac{(y_i - \hat{y}_i)^2}{n}
\f]
*/
double least_squares_error(double *x, double *y, int n, double m, double b);


/**
  @return the maximum least squares error.

\f[
    (y_i - \hat{y}_i)^2
\f]
*/
double least_squares_max_error(double *x, double *y, int n, double m, double b);


#endif // _LEASTSQS_H


