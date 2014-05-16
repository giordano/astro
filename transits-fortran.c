/* transits-fortran.c -- Fortran wrappers for functions in transits.c
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
#include "transits.h"

/* r1   (in)  = radius of the background source;
 * r2   (in)  = radius of the transiting object;
 * d    (in)  = distance between the centers of the objects;
 * x1   (in)  = x coordinate of the background source;
 * x2   (in)  = x coordinate of the transiting object;
 * f    (in)  = geometric factor (f = 4 normalizes the flux to `lum');
 * lum  (in)  = intrinsic luminosity of the source;
 * flux (out) = flux of the source.
 */
void transit_flux_(double *r1, double *r2, double *d, double *x1,
		   double *x2, double *f, double *lum, double *flux)
{
  (*flux)=flusso((*f), (*lum), (*r1),
		 area_coperta((*r1), (*r2), (*d), (*x1), (*x2)));
}
