/* kepler-fortran.c -- Fortran wrappers for functions in kepler.c
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
#include "kepler.h"

/* Funzione che restituisce la posizione dei due corpi nella direzione y'',
 * rispetto al centro di massa.  Argomenti:
 *   anommedia (in) = anomalia media,  prodotto fra velocità angolare media e
 *                    tempo;
 *   phi  (in)  = angolo di rotazione fra sistema Ox'y'z' e Oxyz;
 *   i    (in)  = angolo di inclinazione, cioè angolo di rotazione fra sistema
 *              = Oxyz e Ox''y''z'';
 *   e    (in)  = eccentricità;
 *   a    (in)  = semiasse maggiore;
 *   m1   (in)  = massa del primo corpo;
 *   m2   (in)  = massa del secondo corpo;
 *   p1pc (out) = vettore delle coordinate del primo corpo;
 *   p2pc (out) = vettore delle coordinate del secondo corpo;
 */
void coordinates_(double *anommedia, double *phi, double *i, double *e,
		  double *a, double *m1, double *m2,
		  double p1pc[3], double p2pc[3])
{
  double ecc_anom, rr, theta, ppf[3], ppfpc[3], mu;
  mu=(*m1)*(*m2)/((*m1)+(*m2));
  /* calolo l'anomalia eccentrica */
  ecc_anom=psi_newton((*anommedia), (*e));
  /* calcolo la distanza dal fuoco */
  rr=r((*a), (*e), ecc_anom);
  /* calcolo l'anomalia vera */
  theta=anomvera((*e),ecc_anom);
  /* coordinata x della particella fittizia */
  ppf[0]=rr*cos(theta);
  /* coordinata y della particella fittizia */
  ppf[1]=rr*sin(theta);
  /* la coordinata z della particella fittizia è sempre nulla */
  ppf[2]=0.;
  /* Coordinate della particella fittizia e dei due corpi nel piano del cielo */
  pianodelcielo(ppf, (*phi), (*i), ppfpc);
  vettore_scalare(3,ppfpc,p1pc,-mu/(*m1));
  vettore_scalare(3,ppfpc,p2pc,mu/(*m2));
}
