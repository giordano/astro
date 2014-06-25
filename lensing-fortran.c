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

/* Compute the amplification for an extended source with uniform brightness,
 * using equations (16)-(25) in Agol 2002.  Arguments:
 *   uu  (in)  = distance between source and lens centers, in units of Einstein
 *     radii;
 *   rs  (in)  = source radius, in units of Einstein radii;
 *   rl  (in)  = lens radius, in units of Einstein radii;
 *   amp (out) = amplification of the source.
 *
 * NOTA BENE: this function may not give reliable results when the distance
 * between source and lens is much larger than the source radius **and** the
 * source radius is small (e.g., `uu' is about ten order of magnitudes larger
 * than `rs' and `rs' is of the order of unity or less).  You can workaround
 * this issue by approximating the amplification with the amplification by a
 * point-like lens.
 */
void extended_uniform_source_amp_(double *uu, double *rs, double *rl,
				  double *amp)
{
  (*amp)=extended_uniform_source_amp((*uu), (*rs), (*rl));
}

/* Generate pseudo-random numbers uniformly distributed in an interval.
 * Arguments:
 *   aa    (in) = infimum of the interval in which to generate random numbers;
 *   bb    (in) = extremum of the interval in which to generate random numbers;
 *   nn    (in) = number of random numbers to be generated;
 *   rand (out) = vector of length `n' of the generated random numbers.
 */
void rng_uniform_(double *aa, double *bb, int *nn, double rand[(*nn)])
{
  gsl_rng *rng;
  int ii;

  /* Allocate generator. */
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  /* Seed the generator. */
  gsl_rng_set(rng,(unsigned int)time(NULL));
  /* Fill `rand' with generated random numbers. */
  for (ii=0; ii<(*nn); ii++)
    rand[ii]=(*aa) + ((*bb) - (*aa))*gsl_rng_uniform(rng);
  /* Delete generator. */
  gsl_rng_free(rng);
}
