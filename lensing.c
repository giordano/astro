/* lensing.h
 *
 * Copyright (C) 2014 Mosè Giordano
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later *
 * version.  This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * References:
 *   + Agol 2002: "Occultation and microlensing", ApJ, 579:430-436, 2002.
 *     DOI:10.1086/342880.
 *   + Witt & Mao 1994: "Can lensed stars be regarded as pointlike for
 *     microlensing by MACHOs?", ApJ, 430:505-510, 1994.  DOI:10.1086/174426.
 */
#include "lensing.h"

/* Define the precision of elliptic integrals */
#define ELLINT_PRECISION GSL_PREC_DOUBLE

/* Compute the G function, defined in equation (18) of Agol 2002, for the
 * special case phi = pi/2 so we can directly use the complete elliptic
 * integrals.  Arguments:
 *   u  = distance between source and lens;
 *   rs = source radius.
 */
double agol_Gcomp(double u, double rs)
{
  double rsrs, u1, u2, u3, u1u2, n, k1;
  rsrs=rs*rs;
  u1=pow(u - rs,2.);
  u2=pow(u + rs,2.);
  u1u2=u1*u2;
  u3=u*u - rsrs;
  n=1. - u1/u2;
  k1=2.*sqrt((u2 - u1)/((4.*u2 + u1u2)));
  return ((4.*u2 + u1u2)*ellint_Ecomp(k1) - (u1u2 + 8.*u3)*ellint_Kcomp(k1) + 4.*u1*(1. + rsrs)*ellint_Pcomp(k1, n))/sqrt((4.*u2 + u1u2));
}

/* Compute the amplification for an extended source with uniform brightness and
 * unocculted images, using equations (14) and (25) in Agol 2002.  This is the
 * same as equation (9) of Witt & Mao 1994.  Arguments:
 *   u  = distance between source and lens;
 *   rs = source radius.
 */
double extended_uniform_source_amp(double u, double rs)
{
  return agol_Gcomp(u, rs)/(2.*M_PI*rs*rs);
}
