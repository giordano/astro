/* transits.c -- Librery for the simulation of transits.
 *
 * Copyright (C) 2011, 2014 Mosè Giordano
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

/* Funzione che restituisce l'area di di sovrapposizione fra i due corpi.  `r1'
 * è il raggio della stella, `r2' del pianeta, `d' è la loro distanza
 * proiettata, `x1' e `x2' sono le coordinate, rispettivamente della stella e
 * del pianeta, nel piano del cielo dell'osservatore.
 */
double area_coperta(double r1, double r2, double d, double x1,
		    double x2){
  double theta1, theta2; /* vedi Figura 3.5 della tesi */
  double dA; /* area coperta */
  theta1=2*acos((r1*r1-r2*r2+d*d)/(2*r1*d));
  theta2=2*acos((r2*r2-r1*r1+d*d)/(2*r2*d));

  /* se per l'osservatore la stella è davanti al pianeta... */
  if(x1 >= x2)
    /* ...non si verifica l'eclissi (l'area coperta è nulla) */
    dA=0;
  else
    {
      if(d > r1+r2)
	dA=0;
      else if(r1-r2 <= d)
	dA=r1*r1/2*(theta1-sin(theta1))+r2*r2/2*(theta2-sin(theta2));
      else
	dA=M_PI*r2*r2;
    }
  return dA;
}

/* Funzione che restituisce il valore del flusso di una stella.  Vedi equazione
 * (3.27), pag. 46.  `f' è il fattore geometrico, `lum' è la luminosità
 * intrinseca della stella, `r' è il raggio della stella, `dA' è l'area coperta
 * dal pianeta.  Se si passa alla funzione f=4, il risultato sarà normalizzato a
 * `lum'.
 */
double flusso(double f, double lum, double r, double dA)
{
  return f*lum*(1-dA/(M_PI*r*r))/4;
}
