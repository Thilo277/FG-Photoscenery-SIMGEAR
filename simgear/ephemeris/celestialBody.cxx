/**************************************************************************
 * celestialBody.cxx
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

#include <simgear_config.h>
#include <simgear/debug/logstream.hxx>

#include <math.h>

#include "celestialBody.hxx"
#include "star.hxx"


/**************************************************************************
 * void CelestialBody::updatePosition(double mjd, Star *ourSun)
 *
 * Basically, this member function provides a general interface for 
 * calculating the right ascension and declinaion. This function is 
 * used for calculating the planetary positions. For the planets, an 
 * overloaded member function is provided to additionally calculate the
 * planet's magnitude. 
 * The sun and moon have their own overloaded updatePosition member, as their
 * position is calculated an a slightly different manner.  
 *
 * arguments:
 * double mjd: provides the modified julian date.
 * Star *ourSun: the sun's position is needed to convert heliocentric 
 *               coordinates into geocentric coordinates.
 *
 * return value: none
 *
 *************************************************************************/
void CelestialBody::updatePosition(double mjd, Star *ourSun)
{
  double eccAnom, v, ecl, actTime, 
    xv, yv, xh, yh, zh, xg, yg, zg, xe, ye, ze,
    cosN, sinN, cosvw, sinvw, sinvw_cosi, cosecl, sinecl;

  updateOrbElements(mjd);
  actTime = sgCalcActTime(mjd);

  // calcualate the angle bewteen ecliptic and equatorial coordinate system
  ecl = SGD_DEGREES_TO_RADIANS * (23.4393 - 3.563E-7 *actTime);
  
  eccAnom = sgCalcEccAnom(M, e);  //calculate the eccentric anomaly
  xv = a * (cos(eccAnom) - e);
  yv = a * (sqrt (1.0 - e*e) * sin(eccAnom));
  v = atan2(yv, xv);           // the planet's true anomaly
  r = sqrt (xv*xv + yv*yv);    // the planet's distance
  
  // repetitive calculations, minimised for speed
  cosN = cos(N);
  sinN = sin(N);
  cosvw = cos(v+w);
  sinvw = sin(v+w);
  sinvw_cosi = sinvw * cos(i);
  cosecl = cos(ecl);
  sinecl = sin(ecl);

  // calculate the planet's position in 3D space
  xh = r * (cosN * cosvw - sinN * sinvw_cosi);
  yh = r * (sinN * cosvw + cosN * sinvw_cosi);
  zh = r * (sinvw * sin(i));

  // calculate the ecliptic longitude and latitude
  xg = xh + ourSun->getxs();
  yg = yh + ourSun->getys();
  zg = zh;

  lonEcl = atan2(yh, xh);
  latEcl = atan2(zh, sqrt(xh*xh+yh*yh));

  xe = xg;
  ye = yg * cosecl - zg * sinecl;
  ze = yg * sinecl + zg * cosecl;
  rightAscension = atan2(ye, xe);
  declination = atan2(ze, sqrt(xe*xe + ye*ye));
  /* SG_LOG(SG_GENERAL, SG_INFO, "Planet found at : " 
	 << rightAscension << " (ra), " << declination << " (dec)" ); */

  //calculate some variables specific to calculating the magnitude 
  //of the planet
  R = sqrt (xg*xg + yg*yg + zg*zg);
  s = ourSun->getDistance();

  // It is possible from these calculations for the argument to acos
  // to exceed the valid range for acos(). So we do a little extra
  // checking.

  double tmp = (r*r + R*R - s*s) / (2*r*R);
  if ( tmp > 1.0) { 
      tmp = 1.0;
  } else if ( tmp < -1.0) {
      tmp = -1.0;
  }

  FV = SGD_RADIANS_TO_DEGREES * acos( tmp );
}

/****************************************************************************
 * double CelestialBody::sgCalcEccAnom(double M, double e)
 * this private member calculates the eccentric anomaly of a celestial body, 
 * given its mean anomaly and eccentricity.
 * 
 * -Mean anomaly: the approximate angle between the perihelion and the current
 *  position. this angle increases uniformly with time.
 *
 * True anomaly: the actual angle between perihelion and current position.
 *
 * Eccentric anomaly: this is an auxilary angle, used in calculating the true
 * anomaly from the mean anomaly.
 * 
 * -eccentricity. Indicates the amount in which the orbit deviates from a 
 *  circle (0 = circle, 0-1, is ellipse, 1 = parabola, > 1 = hyperbola).
 *
 * This function is also known as solveKeplersEquation()
 *
 * arguments: 
 * M: the mean anomaly
 * e: the eccentricity
 *
 * return value:
 * the eccentric anomaly
 *
 ****************************************************************************/
double CelestialBody::sgCalcEccAnom(double M, double e)
{
  double 
    eccAnom, E0, E1, diff;
  
  eccAnom = M + e * sin(M) * (1.0 + e * cos (M));
  // iterate to achieve a greater precision for larger eccentricities 
  if (e > 0.05)
    {
      E0 = eccAnom;
      do
	{
	  E1 = E0 - (E0 - e * sin(E0) - M) / (1 - e *cos(E0));
	  diff = fabs(E0 - E1);
	  E0 = E1;
	}
      while (diff > (SGD_DEGREES_TO_RADIANS * 0.001));
      return E0;
    }
  return eccAnom;
}

/*****************************************************************************
 * inline CelestialBody::CelestialBody
 * public constructor for a generic celestialBody object.
 * initializes the 6 primary orbital elements. The elements are:
 * N: longitude of the ascending node
 * i: inclination to the ecliptic
 * w: argument of perihelion
 * a: semi-major axis, or mean distance from the sun
 * e: eccenticity
 * M: mean anomaly
 * Each orbital element consists of a constant part and a variable part that 
 * gradually changes over time. 
 *
 * Argumetns:
 * the 13 arguments to the constructor constitute the first, constant 
 * ([NiwaeM]f) and the second variable ([NiwaeM]s) part of the orbital 
 * elements. The 13th argument is the current time. Note that the inclination
 * is written with a capital (If, Is), because 'if' is a reserved word in the 
 * C/C++ programming language.
 ***************************************************************************/ 
CelestialBody::CelestialBody(double Nf, double Ns,
				    double If, double Is,
				    double wf, double ws,
				    double af, double as,
				    double ef, double es,
				    double Mf, double Ms, double mjd)
{
  NFirst = Nf;     NSec = Ns;
  iFirst = If;     iSec = Is;
  wFirst = wf;     wSec = ws;
  aFirst = af;     aSec = as;
  eFirst = ef;     eSec = es;
  MFirst = Mf;     MSec = Ms;
  updateOrbElements(mjd);
}

CelestialBody::CelestialBody(double Nf, double Ns,
				    double If, double Is,
				    double wf, double ws,
				    double af, double as,
				    double ef, double es,
				    double Mf, double Ms)
{
  NFirst = Nf;     NSec = Ns;
  iFirst = If;     iSec = Is;
  wFirst = wf;     wSec = ws;
  aFirst = af;     aSec = as;
  eFirst = ef;     eSec = es;
  MFirst = Mf;     MSec = Ms;
}

/****************************************************************************
 * inline void CelestialBody::updateOrbElements(double mjd)
 * given the current time, this private member calculates the actual 
 * orbital elements
 *
 * Arguments: double mjd: the current modified julian date:
 *
 * return value: none
 ***************************************************************************/
void CelestialBody::updateOrbElements(double mjd)
{
  double actTime = sgCalcActTime(mjd);
   M = SGD_DEGREES_TO_RADIANS * (MFirst + (MSec * actTime));
   w = SGD_DEGREES_TO_RADIANS * (wFirst + (wSec * actTime));
   N = SGD_DEGREES_TO_RADIANS * (NFirst + (NSec * actTime));
   i = SGD_DEGREES_TO_RADIANS * (iFirst + (iSec * actTime));
   e = eFirst + (eSec * actTime);
   a = aFirst + (aSec * actTime);
}

/*****************************************************************************
 * inline double CelestialBody::sgCalcActTime(double mjd)
 * this private member function returns the offset in days from the epoch for
 * wich the orbital elements are calculated (Jan, 1st, 2000).
 * 
 * Argument: the current time
 *
 * return value: the (fractional) number of days until Jan 1, 2000.
 ****************************************************************************/
double CelestialBody::sgCalcActTime(double mjd)
{
  return (mjd - 36523.5);
}

/*****************************************************************************
 * inline void CelestialBody::getPos(double* ra, double* dec)
 * gives public access to Right Ascension and declination
 *
 ****************************************************************************/
void CelestialBody::getPos(double* ra, double* dec)
{
  *ra  = rightAscension;
  *dec = declination;
}

/*****************************************************************************
 * inline void CelestialBody::getPos(double* ra, double* dec, double* magnitude
 * gives public acces to the current Right ascension, declination, and 
 * magnitude
 ****************************************************************************/
void CelestialBody::getPos(double* ra, double* dec, double* magn)
{
  *ra = rightAscension;
  *dec = declination;
  *magn = magnitude;
}


