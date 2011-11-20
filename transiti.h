/* transiti.h -- Header della libreria per la simulazioni di transiti di pianeti
 * extrasolari davanti alle stelle compagne.
 *
 * Copyright (C) 2011 Mosè Giordano
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
#include "keplero.h"
#include <gsl/gsl_sf_ellint.h>

/* Definisco macro per le funzione degli integrali ellittici completi. */
#define ellint_K(k) gsl_sf_ellint_Kcomp(sqrt(k),GSL_PREC_SINGLE)
#define ellint_E(k) gsl_sf_ellint_Ecomp(sqrt(k),GSL_PREC_SINGLE)
#define ellint_Pi(n,k) gsl_sf_ellint_Pcomp(sqrt(k),-n,GSL_PREC_SINGLE)

/* Macro che fornisce la funzione gradino di Heaviside. */
#define HEAVISIDE(z) (z > 0) ? 1 : ((z == 0) ? 1/2 : 0)

double area_coperta(double, double, double, double, double);
double flusso(double, double, double, double);
double lambdae(double, double);
double flusso_ma(double, double, double, double);
