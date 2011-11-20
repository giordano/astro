/* transiti.c -- Libreria per la simulazioni di transiti di pianeti extrasolari
 * davanti alle stelle compagne.
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
 *
 * Riferimenti bibliografici usati in questa libreria:
 * [MA02] Mandel, K., Agol, E., 2002, ApJ, 580, L171. doi:10.1086/345520. Note:
 * consulta anche l'elenco delle errata corrige qui:
 * http://www.astro.washington.edu/users/agol/mandel_agol_errata.pdf
 */
#include "transiti.h"

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

/* Vedi Mandel e Agol [MA02].  p=raggio_pianeta/raggio_stella,
 * z=distanza_pianeta_stella/raggio_pianeta.  TODO: scoprire da dove esce il
 * termine sotto radice e provare a unificare questa funzione con `flusso'.
 */
double lambdae(double p, double z)
{
  double kappa0=acos((1-p*p+z*z)/(2*z)), kappa1=acos((p*p+z*z-1)/(2*p*z));
  if (1+p < z)
    return 0;
  else if (fabs(1-p) < z && z <= 1+p)
    return (p*p*kappa0+kappa1-sqrt((4*z*z-pow((1+z*z-p*p),2))/4))/M_PI;
  else if (z <= 1-p)
    return p*p;
  /* Questo ultimo `else' dovrebbe contemplare il caso (z <= p-1).  Non uso
   * `else if (z <= p-1)' perché altrimenti riceverei un warning in fase di
   * complilazione. */
  else
    return 1;
}

/* Funzione che restituisce il flusso luminoso di una stella tenendo in
 * considerazione l'effetto del limb darkening.  In questa funzione si fa uso
 * del limb darkening quadratico così come descritto in Mandel e Agol [MA02].
 * p=raggio_pianeta/raggio_stella, z=distanza_pianeta_stella/raggio_pianeta.
 * TODO: rivedere questa funzione perché non l'ho mai testata!
 */
double flusso_ma(double p, double z, double gamma1, double gamma2)
{
  double a=(z-p)*(z-p), b=(z+p)*(z+p);
  double c[5];
  c[2]=gamma1+2*gamma2;
  c[4]=-gamma2;
  c[0]=1-c[1]-c[2]-c[3]-c[4];
  double Omega=c[0]/(0+4)+c[2]/(2+4)+c[4]/(4+4); /* c[1] = c[3] = 0 */
  double k=sqrt((1-a)/(4*z*p)), q=p*p-z*z;
  double kappa0=acos((1-p*p+z*z)/(2*z)), kappa1=acos((p*p+z*z-1)/(2*p*z));
  double lambdad, etad;
  double eta1, eta2, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6;

  /* Definisco le funzioni ausiliarie utilizzate da Mandel e Agol nella
   * tabella 1 di [MA02]. */
  lambda1=(((1-b)*(2*b+a-3)-3*q*(b-2))*ellint_K(k) + 4*p*z*(z*z+7*p*p-4)*ellint_E(k)-3*q*ellint_Pi((a-1)/a,k)/a)/(9*M_PI*sqrt(p*z));
  lambda2=((1-5*z*z+p*p+q*q)*ellint_K(k)+(1-a)*(z*z+7*p*p-4)*ellint_E(k)-3*q*ellint_Pi((a-b)/a,k)/a)*2/(9*M_PI*sqrt(1-a));
  /* Vedi la definizione corretta di lambda3 nell'errata corrige di [MA02].
   * Nell'articolo originale c'era ellint_E(1/(2*k)) ed ellint_K(1/(2*k)) invece
   * dei corretti ellint_E(1/(2*p)) ed ellint_K(1/(2*p)) */
  lambda3=1/3+16*p*(2*p*p-1)*ellint_E(1/(2*p))/(9*M_PI)-(1-4*p*p)*(3-8*p*p)*ellint_K(1/(2*p))/(9*M_PI*p);
  /* Vedi la definizione corretta di lambda4 nell'errata corrige di [MA02].
   * Nell'articolo originale c'era ellint_E(2*k) ed ellint_K(2*k) invece dei
   * corretti ellint_E(2*p) ed ellint_K(2*p) */
  lambda4=1/3+2*((8*p*p-4)*ellint_E(2*p)+(1-4*p*p)*ellint_K(2*p))/(9*M_PI);
    /* Vedi la definizione corretta di lambda5 nell'errata corrige di [MA02]. */
  lambda5=2*acos(1-2*p)/(3*M_PI)-(12+8*p-32*p*p)*sqrt(p-p*p)/(9*M_PI)-2*HEAVISIDE(p-1/2)/3;
  lambda6=-2*pow((1-p*p),3/2)/3;
  eta2=p*p*(p*p+2*z*z)/2;
  eta1=(kappa1+2*eta2*kappa0-(1+5*p*p+z*z)*sqrt((1-a)*(b-1))/4)/(2*M_PI);

  if ((p > 0 && z >= 1+p) || (p == 0 && z >= 0)) /* Caso 1 */
    {
      lambdad=0;
      etad=0;
    }
  else if((p > 0 && z > 1/2. + fabs(p-1/2.) && z < 1+p) || (p > 1/2. && z>= fabs(1-p) && z < p)) /* Casi 2, 8 */
    {
      lambdad=lambda1;
      etad=eta1;
    }
  else if((p > 0 && p < 1/2. && z > p && z < 1-p) || (p > 0 && p < 1 && z> 0 && z < 1/2.-fabs(p-1/2.))) /* Casi 3, 9 */
    {
      lambdad=lambda2;
      etad=eta2;
    }
  else if(p > 0 && p < 1/2. && z == 1-p) /* Caso 4 */
    {
      lambdad=lambda5;
      etad=eta2;
    }
  else if (p > 0 && p < 1/2. && z == p) /* Caso 5 */
    {
      lambdad=lambda4;
      etad=eta2;
    }
  else if (p == 1/2. && z == 1/2.) /* Caso 6 */
    {
      lambdad=1/3.-4/(9.*M_PI);
      etad=3./32.;
    }
  else if (p > 1/2. && z == p) /* Caso 7 */
    {
      lambdad=lambda3;
      etad=eta1;
    }
  else if (p > 0 && p < 1 && z == 0) /* Caso 10 */
    {
      lambdad=lambda6;
      etad=eta2;
    }
  else if (p > 1 && z >= 0 && z < p-1) /* Caso 11 */
    {
      /* Vedi l'errata corrige di [MA02].  Nell'articolo originale c'era
       * lambdad=1 invece di lambdad=0 ed etad=1 invece di etad=1/2 */
      lambdad=0;
      etad=1/2;
    }
  /* restituisco il flusso */
  return 1-((1-c[2])*lambdae(p,z) + c[2]*(lambdad + 2*HEAVISIDE(p-z)/3)-c[4]*etad)/(4*Omega);
}
