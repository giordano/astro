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
 */
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>

#define ellint_Kcomp(k)   gsl_sf_ellint_Kcomp(k,ELLINT_PRECISION)
#define ellint_Ecomp(k)   gsl_sf_ellint_Ecomp(k,ELLINT_PRECISION)
#define ellint_Pcomp(k,n) gsl_sf_ellint_Pcomp(k,-n,ELLINT_PRECISION)

double agol_Gcomp(double, double);
double extended_uniform_source_amp(double, double);
