/**************************************************************************
 * venus.cxx
 * Written by Durk Talsma. Originally started October 1997, for distribution  
 * with the FlightGear project. Version 2 was written in August and 
 * September 1998. This code is based upon algorithms and data kindly 
 * provided by Mr. Paul Schlyter. (pausch@saaf.se). 
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 * $Id$
 **************************************************************************/

#include <cmath>

#include "venus.hxx"

/*************************************************************************
 * Venus::Venus(double mjd)
 * Public constructor for class Venus
 * Argument: The current time.
 * the hard coded orbital elements for Venus are passed to 
 * CelestialBody::CelestialBody();
 ************************************************************************/
Venus::Venus(double mjd) :
  CelestialBody(76.67990,  2.4659000E-5, 
		3.3946,    2.75E-8,
		54.89100,  1.3837400E-5,
		0.7233300, 0.000000,
		0.006773, -1.302E-9,
		48.00520,  1.60213022440, mjd)
{
}
Venus::Venus() :
  CelestialBody(76.67990,  2.4659000E-5, 
		3.3946,    2.75E-8,
		54.89100,  1.3837400E-5,
		0.7233300, 0.000000,
		0.006773, -1.302E-9,
		48.00520,  1.60213022440)
{
}

/*************************************************************************
 * void Venus::updatePosition(double mjd, Star *ourSun)
 * 
 * calculates the current position of Venus, by calling the base class,
 * CelestialBody::updatePosition(); The current magnitude is calculated using 
 * a Venus specific equation
 *************************************************************************/
void Venus::updatePosition(double mjd, Star *ourSun)
{
  CelestialBody::updatePosition(mjd, ourSun);
  magnitude = -4.34 + 5*log10( r*R ) + 0.013 * FV + 4.2E-07 * pow(FV,3);
}
