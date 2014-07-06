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
 *   + Bulirsch 1965: "Numerical Calculation of Elliptic Integrals and Elliptic
 *     Functions", Numerische Mathematik 7, 78-90 and 353-354, 1965.  DOIs:
 *     10.1007/BF01397975 and 10.1007/BF01436529.
 */
#include "lensing.h"

/* Define the precision of elliptic integrals */
#define ELLINT_PRECISION GSL_PREC_DOUBLE
/* Define precision accuracy. */
#define EPSILON GSL_DBL_EPSILON

/* Complete elliptic integrals, with Bulirsch 1965 algorithms. */
double bulirsch_ellint_Kcomp (double kk)
{
  double kc, mm, hh;
  kc = sqrt(1. - kk*kk);
  mm = 1.;
  while (1)
    {
      hh = mm;
      mm += kc;
      if (fabs(1. - kc/hh) - 1e-13 > EPSILON)
	{
	  kc = sqrt(hh*kc);
	  mm = mm/2.;
	}
      else
	break;
    }
  return M_PI/mm;
}

double bulirsch_ellint_Ecomp (double kk)
{
  double hh, mm, aa, bb, kc, cc;
  mm = 1.;
  aa = 1.;
  bb = 1. - kk*kk;
  kc = sqrt(bb);
  cc = aa;
  aa += bb;
  while (1)
    {
      bb = 2.*(cc*kc + bb);
      cc = aa;
      hh = mm;
      mm += kc;
      aa += bb/mm;
      if (fabs(1. - kc/hh) - 1e-13 > EPSILON)
	kc = 2.*sqrt(hh*kc);
      else
	break;
    }
  return M_PI_4*aa/mm;
}

double bulirsch_ellint_Pcomp (double kk, double nn)
{
  double kc, pp, m0, cc, ff, gg, dd, ee;
  kc = sqrt(1. - kk*kk);
  pp = nn + 1.;
  if(fabs(kc*pp) < -EPSILON)
    GSL_ERROR("kc*p = 0", GSL_EINVAL);
  else
    {
      ee = kc;
      m0 = 1.;
      if (pp > EPSILON)
	{
	  cc = 1.;
	  pp = sqrt(pp);
	  dd = 1./pp;
	}
      else
	{
	  gg = 1. - pp;
	  ff = kc*kc - pp;
	  pp = sqrt(ff/gg);
	}
      while (1)
	{
	  ff = cc;
	  cc += dd/pp;
	  gg = ee/pp;
	  dd = 2.*(ff*gg + dd);
	  pp += gg;
	  gg = m0;
	  m0 += kc;
	  if(fabs(1. - kc/gg) - 1e-13 > EPSILON)
	    {
	      kc = 2.*sqrt(ee);
	      ee = kc*m0;
	    }
	  else
	    break;
	}
    }
  return M_PI_2*(cc*m0 + dd)/(m0*(m0 + pp));
}

/* COMPUTE THE G function, defined in equation (18) of Agol 2002.  Arguments:
 *   phi = variable of integration in equation (15) of Agol 2002.  See also
 *         equation (18);
 *   uu  = distance between source and lens centers;
 *   rs  = source radius.
 */
double agol_G(double phi, double uu, double rs)
{
  double UU, rsrs, u1, u2, u3, u1u2, nn, k1;
  UU=uu*uu;
  rsrs=rs*rs;
  /* Using the explicit extended expression for `u1' and `u2' gives more
   * reliable results than using the `pow' function.
   */
  u1=UU - 2.*uu*rs + rsrs;
  u2=UU + 2.*uu*rs + rsrs;
  u1u2=u1*u2;
  u3=UU - rsrs;
  /* When `uu' and `rs' are different but almost equal, u1 is very small (of
   * order 10^-15 or less) and 1.-u1/u2 could be equal to 1., but last argument
   * of P elliptic integral must be less than 1.  Thus we use the following
   * equivalent expression for `nn' which in some cases avoids `nn' being
   * exactly equal to 1. */
  nn=4.*uu*rs/u2;
  /* This prevents `nn' from being exactly 1 in any case. */
  if (nn >= 1)
    nn=1.-EPSILON;
  k1=2.*sqrt(nn/(4. + u1));
  if (fabs(phi - M_PI_2) < EPSILON)
    /* For the special case phi = pi/2 use the dedicated complete elliptic
     * integral functions, which are slighly faster than the incomplete ones.
     */
    return ((4.*u2 + u1u2)*ellint_Ecomp(k1) - (u1u2 + 8.*u3)*ellint_Kcomp(k1)
	    + 4.*u1*(1. + rsrs)*ellint_Pcomp(k1, nn))/sqrt(4.*u2 + u1u2);
  else
    return ((4.*u2 + u1u2)*ellint_E(phi, k1) - (u1u2 + 8.*u3)*ellint_F(phi, k1)
	    + 4.*u1*(1. + rsrs)*ellint_P(phi, k1, nn))/sqrt(4.*u2 + u1u2);
}

/* Compute the amplification for an extended source with uniform brightness,
 * using equations (16)-(25) in Agol 2002.  Arguments:
 *   uu = distance between source and lens centers, in units of Einstein radii;
 *   rs = source radius, in units of Einstein radii;
 *   rl = lens radius, in units of Einstein radii.
 *
 * NOTA BENE: this function may not give reliable results when the distance
 * between source and lens is much larger than the source radius **and** the
 * source radius is small (e.g., `uu' is about ten order of magnitudes larger
 * than `rs' and `rs' is of the order of unity or less).  You can workaround
 * this issue by approximating the amplification with the amplification of a
 * point-like source.
 */
double extended_uniform_source_amp(double uu, double rs, double rl)
{
  double UU, rlrl, rsrs, muplus, muminus, mu, v1, v2, u0, u1, u2, u3, phi0,
    phi1, phi2, psi0, psi1, psi2, bl;
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
	/* There is a misprint in the ApJ paper: the factor in front of the atan
	 * is (1 + rs^2) instead of (1 + rs)^2.  It has already been fixed in
	 * the FORTRAN code provided by Agol.
	 */
	mu=2.*(rs + (1. + rsrs)*atan(rs))/(M_PI*rsrs);

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
	      return muplus + muminus
		+ (0.25*v2*(bl - v1) + (1. + rsrs)*(acos(0.5*bl/rs)
						    - atan(v2/v1)))/(M_PI*rsrs)
		- rlrl/M_PI/rsrs*acos(0.5*bl/rs);
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
	      u1=UU - 2.*uu*rs + rsrs;
	      u2=UU + 2.*uu*rs + rsrs;
	      u3=UU - rsrs;
	      phi0=acos(sqrt(u1*(u2 - u0)/(u0*(u2-u1))));
	      phi1=acos((u1 + u2 - 2.*u0)/(u2-u1));
	      phi2=acos((u3 + u0)/(2.*uu*bl));
	      return muplus + HEAVISIDE(rs - uu)/rsrs -
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
	      u1=UU - 2.*uu*rs + rsrs;
	      u2=UU + 2.*uu*rs + rsrs;
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
