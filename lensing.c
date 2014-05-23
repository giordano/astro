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

#define PRECISION 2e-16

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
 *   uu = distance between source and lens, in units of Einstein radii;
 *   rs = source radius, in units of Einstein radii.
 */
double extended_uniform_source_amp(double uu, double rs, double rl)
{
  double muplus, muminus, UU, mu, v1, v2, u0, u1, u2, u3, phi0, phi1, phi2;
  UU=uu*uu;
  rlrl=rl*rl;
  rsrs=rs*rs;
  /* Conversion table
   * rs:      rs
   * rl:      rl
   * u:       be = bvec*rs
   * muplus:  muop
   * muminus: muom
   */

  /* Point-like source, use Paczinsky amplification. */
  if (fabs(rs) < PRECISION)
    { /* rs = 0 */
      if (fabs(rl) < PRECISION)
	{ /* rl = 0 */
	  /* Outer image, case I */
	  muplus=0.5*((2. + uu)/(u*sqrt(4. + uu)) + 1.);
	  if (u*rl - (1. - rlrl) < PRECISION)
	    /* Inner image, case I */
	    muminus=0.5*((2. + uu)/(u*sqrt(4. + uu)) - 1.);
	  else
	    /* Inner image, case II */
	    muminus=0.;
	}
      else
	{ /* rl != 0 */
	  if(u*rl - (rlrl - 1.) > PRECISION)
	    /* Outer image, case I */
	    muplus=0.5*((2. + uu)/(u*sqrt(4. + uu)) + 1.);
	  else
	    /* Outer image, case II */
	    muplus=0.;
	  muminus=0.;
	}
      return muplus + muminus;
    }

  /* Extended source, use amplification formulae by Witt & Mao 1994. */
  if (fabs(uu - rs) > PRECISION)
    /* Case III (u != rs) */
    mu=agol_Gcomp(u, rs)/(2.*M_PI*rs*rs);
  else
    /* Case IV */
    mu=2.*(rs + (1 + 2.*rs + rsrs)*atan(rs))/(M_PI*rsrs);

  muplus =0.5*(1. + mu);
  muminus=0.5*(1. - mu);

  if (fabs(rl) < PRECISION)
    /* rl = 0 */
    return mu;
  else if (fabs(rl - 1.) < PRECISION)
    /* rl = 1 */
    return muplus;
  else if (rl < 1. && fabs(rl) > PRECISION)
    { /* 0 < rl <1 */
      /* Project lens radius onto source plane */
      bl=1/rl - rl;
      if (u > bl+rs)
	/* Inner image case II, outer image case III */
	return muplus;
      else if (fabs(u - rs) < PRECISION && rs > 0.5*bl)
	{ /* Inner image case VI */
	  v1=sqrt(4. + bl*bl);
	  v2=sqrt(4.*rsrs - bl*bl);
	  a1=(0.25*v2*(bl-v1) + (1. + rsrs)*(acos(0.5*bl/rs) - atan(v2/v1)))/(M_PI*rsrs);
	  a2=rlrl/M_PI/rsrs*acos(0.5*bl/rs);
	  muminus+=a1 - a2;
	  return muplus + muminus;
	}
      else if (u < rs - bl)
	{ /* Inner image, case VII */
	  muminus=(1. - rlrl)/rsrs;
	}
      else
	{
	  u0=bl*bl;
	  u1=pow(uu - rs,2);
	  u2=pow(uu + rs,2);
	  u3=UU - rsrs;
	  phi0=acos(sqrt(u1*(u2 - u0)/(u0*(u2-u1))));
	  phi1=acos((u1 + u2 - 2.*u0)/(u2-u1));
	  phi2=acos((u3 + u0)/(2.*u*bl));
	  if (be < rs)
	    muminus=1./rsrs;
	  else
	    muminus=0;
	  muminus-=(-4.*u3/fabs(u3)*phi0 + 2.*(1. + rsrs)*phi1 + 4.*rlrl*phi2 + sqrt((u2 - u0)*(u0 - u1))*(sqrt(1. + 4./u0) - 1.) - agol_G(phi0, uu, rs))/4./M_PI/rsrs;
	}
    }
  bl=rl - 1./rl;
  if (uu > bl + rs)
    return muplus;
  else if (rs <= bl && uu <= bl - rs)
    return 0;
  else if (rs < bl && fabs(uu) < rs - bl)
    {
      muplus=mu + 0.5 - (rlrl - 1.)/rsrs;
      return muplus;
    }
  else
    {
      u1=pow(uu - rs,2);
      u2=pow(uu + rs,2);
      u3=UU - rsrs;
      phi0=acos(sqrt(u1*(u2 - u0)/(u0*(u2-u1))));
      phi1=acos((u1 + u2 - 2.*u0)/(u2-u1));
      phi2=acos((u3 + u0)/(2.*u*bl));
      psi1=M_PI_2 - phi0;
      psi2=M_PI - phi1 + 2.*acos(sqrt(u0*(4. + u0)/(u0*(4. + u1 + u2) - u1*u2)));
      if (fabs(u3) > PRECISION)
	muplus=(-4.*u3/fabs(u3)*psi1 + 2.*(1. + rsrs)*psi2 - 4.*rlrl*phi2 + sqrt((u2 - u0)*(u0 - u1))*(sqrt(u0/(4. + u0)) + 1.) + agol_G(psi0, uu, rs))/4./M_PI/rsrs;
      else
	muplus=(-4.*psi1 + 2.*(1. + rsrs)*psi2 - 4.*rlrl*phi2 + sqrt((u2 - u0)*(u0 - u1))*(sqrt(u0/(4. + u0)) + 1.) + agol_G(psi0, uu, rs))/4./M_PI/rsrs;
      muminus=0.;
      return muplus;
    }
}
