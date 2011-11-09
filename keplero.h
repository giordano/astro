/* keplero.h -- Header della libreria per la risoluzione dell'equazione di
 * Keplero.
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
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

/* Definisco funzioni in modo che richiamino le corrispondenti funzioni della
 * libreria GSL. Se si vuole utilizzare una libreria diversa dalla GSL sarà
 * sufficiente cambiare le seguenti macro lasciando invariato il resto.
 */
#define bessel_Jn gsl_sf_bessel_Jn

/* Numero di iterazioni da compiere con il metodo di Bessel prima di fermarsi */
#define MAX_BESSEL 200
/* Specifico il valore della precisione desiderata */
#define PRECISIONE 1e-12

double f(double, double, double);
double diff_besselj(int, double);
double psi_newton(double, double);
double psi_bessel(double, double);
double r_bessel(double, double, double);
double r(double, double, double);
double anomvera(double, double);
void pianodelcielo(double Qa[3], double, double, double Qb[3]);
void vettore_scalare(int, double a[], double b[], double);
