/* lensing-fortran.c
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
#include "lensing-fortran.h"

/* Compute the amplification for an extended source with uniform brightness and
 * unocculted images, using equations (14) and (25) in Agol 2002.  This is the
 * same as equation (9) of Witt & Mao 1994.  Arguments:
 *   u   (in)  = distance between source and lens;
 *   rs  (in)  = source radius;
 *   amp (out) = amplification of the source.
 */
void extended_uniform_source_amp_(double *u, double *rs, double *amp)
{
  (*amp)=Gcomp((*u), (*rs))/(2.*M_PI*(*rs)*(*rs));
}
