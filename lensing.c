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
/* Define precision accuracy. */
#define EPSILON GSL_DBL_EPSILON

/* Compute the G function, defined in equation (18) of Agol 2002.  Arguments:
 *   phi = variable of integration in equation (15) of Agol 2002.  See also
 *         equation (18);
 *   u   = distance between source and lens;
 *   rs  = source radius.
 */
double agol_G(double phi, double uu, double rs)
{
  double rsrs, u1, u2, u3, u1u2, nn, k1;
  rsrs=rs*rs;
  u1=pow(uu - rs,2.);
  u2=pow(uu + rs,2.);
  u1u2=u1*u2;
  u3=uu*uu - rsrs;
  nn=1. - u1/u2;
  k1=2.*sqrt((u2 - u1)/((4.*u2 + u1u2)));
  if (fabs(phi - M_PI_2) < EPSILON)
    /* For the special case phi = pi/2 use the dedicated complete elliptic
     * integral functions, which are slighly faster than the incomplete ones.
     */
    return ((4.*u2 + u1u2)*ellint_Ecomp(k1) - (u1u2 + 8.*u3)*ellint_Kcomp(k1)
	    + 4.*u1*(1. + rsrs)*ellint_Pcomp(k1, nn))/sqrt((4.*u2 + u1u2));
  else
    return ((4.*u2 + u1u2)*ellint_E(phi, k1) - (u1u2 + 8.*u3)*ellint_F(phi, k1)
	    + 4.*u1*(1. + rsrs)*ellint_P(phi, k1, nn))/sqrt((4.*u2 + u1u2));
}

/* Compute the amplification for an extended source with uniform brightness,
 * using equations (16)-(25) in Agol 2002.  Arguments:
 *   uu = distance between source and lens, in units of Einstein radii;
 *   rs = source radius, in units of Einstein radii.
 */
double extended_uniform_source_amp(double uu, double rs, double rl)
{
  double UU, rlrl, rsrs, muplus, muminus,  mu, v1, v2, u0, u1, u2, u3,
    phi0, phi1, phi2, psi0, psi1, psi2, a1, a2, bl;
  UU=uu*uu;
  rlrl=rl*rl;

  if (fabs(rs) < EPSILON)
    { /* rs = 0 */
      /* Point-like source, use Paczinsky amplification. */
      if (rl - 1. < EPSILON)
	{ /* rl < 1 */
	  /* Outer image, case I */
	  if (uu*rl - (1. - rlrl) < EPSILON)
	    /* zeta_0 < beta_l */
	    /* Inner image, case I */
	    return (2. + UU)/(uu*sqrt(4. + UU));
	  else
	    /* zeta_0 >= beta_l + rs  (remember: rs = 0) */
	    /* Inner image, case II */
	    return 0.5*((2. + UU)/(uu*sqrt(4. + UU)) + 1.);
	}
      else
	{ /* rl >= 1 */
	  /* Inner image, case II (mu_{-} = 0) */
	  if(uu*rl - (rlrl - 1.) > -EPSILON)
	    /* zeta_0 > beta_l */
	    /* Outer image, case I */
	    return 0.5*((2. + UU)/(uu*sqrt(4. + UU)) + 1.);
	  else
	    /* zeta_0 <= beta_l - rs  (remember: rs = 0) */
	    /* Outer image, case II */
	    return 0.;
	}
    }
  else
    { /* rs != 0 */
      /* Extended source, use amplification formulae by Witt & Mao 1994. */

      rsrs=rs*rs;
      /* Calculate amplification which will be used as base for successive
       * calculations.
       */
      if (fabs(uu - rs) > EPSILON)
	/* zeta_0 != rs */
	/* Both images, case III */
	mu=agol_G(M_PI_2, uu, rs)/(2.*M_PI*rsrs);
      else
	/* zeta_0 = rs */
	/* Both images, case IV */
	mu=2.*(rs + (1 + 2.*rs + rsrs)*atan(rs))/(M_PI*rsrs);

      /* Calculate amplifications of inner and outer images which will be used
       * in successive calculations.  In cases III and IV, mu_{±} = (···) ± 1/2
       * (see equations (22) and (25) of Agol 2002), so mu_{±} = (1 ± mu)/2.
       */
      muplus =0.5*(1. + mu);
      muminus=0.5*(1. - mu);

      if (fabs(rl) < EPSILON)
	{ /* rl = 0 */
	  /* Both images, cases III and IV */
	  return mu;
	}
      else if (fabs(rl - 1.) < EPSILON)
	{ /* rl = 1 */
	  /* Inner image, case II (mu_{-} = 0) */
	  return muplus;
	}
      else if (rl - 1. < EPSILON && rl > EPSILON)
	{ /* 0 < rl < 1 */
	  /* Project lens radius onto source plane */
	  bl=1/rl - rl;
	  if (uu - (bl + rs) > EPSILON)
	    { /* zeta_0 > beta_l + rs */
	      /* Inner image case II (mu_{-}=0), outer image case III */
	      return muplus;
	    }
	  else if (fabs(uu - rs) < EPSILON && rs - 0.5*bl > EPSILON)
	    { /* zeta_0 = rs && rs > beta_l/2  */
	      /* Outer image case IV (muplus already calculated), inner image
	       * case VI.
	       */
	      v1=sqrt(4. + bl*bl);
	      v2=sqrt(4.*rsrs - bl*bl);
	      /* TODO: check this formula */
	      a1=(0.25*v2*(bl - v1)
		  + (1. + rsrs)*(acos(0.5*bl/rs)
				 - atan(v2/v1)))/(M_PI*rsrs);
	      a2=rlrl/M_PI/rsrs*acos(0.5*bl/rs);
	      return muplus + muminus + a1 - a2;
	    }
	  else if (uu - (bl - rs) <= EPSILON)
	    { /* zeta_0 <= beta_l - rs */
	      /* Both images, case III */
	      return mu;
	    }
	  else if (uu - (rs - bl) <= EPSILON)
	    { /* zeta_0 <= rs - beta_l */
	      /* Inner image case VII, outer image case III */
	      return muplus + (1. - rlrl)/rsrs;
	    }
	  else
	    { /* zeta_0 > rs - beta_l */
	      /* Inner image, case V? */
	      u0=bl*bl;
	      u1=pow(uu - rs,2);
	      u2=pow(uu + rs,2);
	      u3=UU - rsrs;
	      phi0=acos(sqrt(u1*(u2 - u0)/(u0*(u2-u1))));
	      phi1=acos((u1 + u2 - 2.*u0)/(u2-u1));
	      phi2=acos((u3 + u0)/(2.*uu*bl));
	      return muplus + HEAVISIDE(uu - rs)/rsrs -
		(-4.*SGN(u3)*phi0 + 2.*(1. + rsrs)*phi1 + 4.*rlrl*phi2 +
		 sqrt((u2 - u0)*(u0 - u1))*(sqrt(1. + 4./u0) - 1.)
		 - agol_G(phi0, uu, rs))/4./M_PI/rsrs;
	    }
	}
      else
	{ /* rl > 1 */
	  /* Inner image, case II (mu_{-} = 0) */
	  bl=rl - 1./rl;
	  if (uu - (bl + rs) > EPSILON)
	    { /* zeta_0 > beta_l + rs */
	      /* Outer image, case III */
	      return muplus;
	    }
	  else if (rs - bl <= EPSILON && uu - (bl - rs) <= EPSILON)
	    { /* rs <= beta_l && zeta_0 <= beta_l - rs */
	      /* Outer image, case II */
	      return 0;
	    }
	  else if (rs - bl > EPSILON && uu - (rs - bl) <= EPSILON)
	    { /* rs > beta_l && zeta_0 <= rs - bl */
	      /* Outer image, case VI */
	      return muplus + (1. - rlrl)/rsrs;
	    }
	  else
	    { /* Outer image, case V? */
	      u0=bl*bl;
	      u1=pow(uu - rs,2);
	      u2=pow(uu + rs,2);
	      u3=UU - rsrs;
	      psi0=acos(sqrt((u0 - u1)*(4. + u2)/(4. + u0)/(u2 - u1)));
	      phi0=acos(sqrt(u1*(u2 - u0)/(u0*(u2-u1))));
	      phi1=acos((u1 + u2 - 2.*u0)/(u2-u1));
	      phi2=acos((u3 + u0)/(2.*uu*bl));
	      psi1=M_PI_2 - phi0;
	      psi2=M_PI - phi1
		+ 2.*acos(sqrt(u0*(4. + u0)/(u0*(4. + u1 + u2) - u1*u2)));
	      /* XXX: this is slighly different from Agol code: here we use
	       * SGN(u3) also for the case u3 = 0, Agol replaces SGN(u3) with 1
	       * in that case.
	       */
	      return (-4.*SGN(u3)*psi1 + 2.*(1. + rsrs)*psi2 - 4.*rlrl*phi2
		      + sqrt((u2 - u0)*(u0 - u1))*(sqrt(u0/(4. + u0)) + 1.)
		      + agol_G(psi0, uu, rs))/4./M_PI/rsrs;
	    }
	}
    }
}
