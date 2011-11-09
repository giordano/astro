/* keplero.c -- Libreria per la risoluzione dell'equazione di Keplero.
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

/* Funzione di cui vogliamo trovare le radici. Il primo
 * argomento è l'anomalia eccentrica, il secondo argomento
 * è l'eccentricità, il terzo è l'anomalia media, cioè il
 * prodotto della velocità angolare media e del tempo.
 */
double f(double psi, double e, double phi)
{
  /* phi deve trovarsi nell'intervallo [0,2pi] */
  phi=fmod(phi,2*M_PI);
  return psi-e*sin(psi)-phi;
}

/* Funzione che restituisce la derivata rispetto a x del
 * coefficiente di Bessel J_n(x) sfruttando le proprietà dei
 * coefficienti di Bessel. Il primo argomento è l'ordine del
 * coefficiente di Bessel, il secondo è il suo argomento
 */
double diff_besselj(int n, double x)
{
  return (bessel_Jn(n-1,x)-bessel_Jn(n+1,x))/2.;
}

/* Funzione che restituisce il valore dell'anomalia eccentrica
 * psi in corrispondenza dell'anomalia media phi calcolata
 * usando il metodo di Newton. Il primo argomento è l'anomalia
 * media, il secondo è l'eccentricità
 */
double psi_newton(double phi, double e)
{
  double psi; /* anomalia eccentrica. */
  /* phi deve trovarsi nell'intervallo [0,2pi] */
  phi=fmod(phi,2*M_PI);
  /* punto iniziale per anomalia eccentrica = anomalia media */
  psi=phi;
  /* Se il valore assoluto della funzione valutata nel punto
   * iniziale è maggiore della precisione desiderata utilizzo
   * il metodo di Newton per cercare un nuovo punto.
   */
  while(fabs(f(psi,e,phi))>PRECISIONE)
    psi-=(psi-e*sin(psi)-phi)/(1-e*cos(psi));
  return psi;
}

/* Funzione che restituisce il valore dell'anomalia eccentrica
 * psi in corrispondenza dell'anomalia media phi calcolata
 * usando il metodo dei coefficienti di Bessel. Il primo
 * argomento è l'anomalia media, il secondo l'eccentricità
 */
double psi_bessel(double phi, double e)
{
  int n;
  double psi;
  /* phi deve trovarsi nell'intervallo [0,2pi] */
  phi=fmod(phi,2*M_PI);
  psi=phi;
  for(n=1;n<=MAX_BESSEL;n++)
    psi+=2*bessel_Jn(n,n*e)*sin(n*phi)/n;
  return psi;
}

/* Funzione che restituisce il valore della distanza r da
 * fuoco al tempo t, calcolata con il metodo dei coefficienti
 * di Bessel. Il primo argomento è l'anomalia media, il
 * secondo il semiasse maggiore, il terzo l'eccentricità
 */
double r_bessel(double phi, double semiasse, double e)
{
  int n;
  double distanza=semiasse*(1+e*e/2.);
  /* phi deve trovarsi nell'intervallo [0,2pi] */
  phi=fmod(phi,2*M_PI);
  for(n=1;n<=MAX_BESSEL;n++)
    distanza-=2*semiasse*e*diff_besselj(n,n*e)*cos(n*phi)/n;
  return distanza;
}

/* Funzione che restituisce il valore della distanza dal
 * fuoco. Il primo argomento è la lunghezza del semiasse
 * maggiore, il secondo è l'eccentricità dell'orbita, il
 * terzo è il valore dell'anomalia eccentrica.
 */
double r(double semiasse, double e, double psi)
{
  return semiasse*(1-e*cos(psi));
}

/* Funzione che restituisce il valore dell'anomalia vera. Il
 * primo argomento è l'eccentricità, il secondo è
 * l'anomalia eccentrica.
 */
double anomvera(double e, double psi)
{
  if(psi<=M_PI)
    return 2*atan(sqrt((1+e)/(1-e))*tan(psi/2));
  else
    return 2*(atan(sqrt((1+e)/(1-e))*tan(psi/2))+M_PI);
}

/* Funzione per calcolare l'anomalia vera facendo uso delle funzioni di Bessel
 * (vedi pagina 554 del libro di G. N. Watson sulle funzioni di Bessel). Questa
 * funzione non è molto efficiente (converge *molto* lentamente e i tempi di
 * esecuzione sono elevati) e non è stata testata a sufficienza. La metto qui
 * solo per curiosità ma il suo utilizzo è sconsigliato. */
double anomvera_bessel(double e, double phi)
{
  int n, m;
  double theta=fmod(phi,2*M_PI), C, D=0;
  double f;
  f=e/(1+sqrt(1-e*e));
  for(n=1;n<=MAX_BESSEL;n++)
    {
      for(m=1;m<=MAX_BESSEL;m++)
	{
	  /* esco dal ciclo in caso di numeri troppo piccoli per evitare underflow */
	  if(bessel_Jn(n+m,n*e)<1e-305 || bessel_Jn(n-m,n*e)<1e-305)
	    /* TODO: provare a salvare i risultati di queste operazioni in
	       variabili, di modo che se l'if da esito negativo non sia
	       necessario ripetere nuovamente il calcolo.*/
	    break ;
	  D+=pow(f,m)*(bessel_Jn(n-m,n*e) + bessel_Jn(n+m,n*e));
	}
      C=2./n*(bessel_Jn(n,n*e) + D);
      theta+=C*sin(n*phi);
    }
  return theta;
}

/* Funzione che trasforma le coordinate del punto Qa del sistema di riferimento
 * intrinseco al sistema binario nelle coordinate del punto Qb visto da un
 * osservatore nel proprio piano del cielo. Abbiamo sfruttato la formula (1.97)
 * di pagina 24 della tesi e le trasformazioni discusse nel capitolo 4 del
 * Goldstein, Poole e Safko (vedi la bibliografia). `phi' è l'angolo fra l'asse
 * x' e l'asse x, `i' è l'angolo fra l'asse z e x'', come spiegato nella tesi.
 */
void pianodelcielo(double Qa[3], double phi, double i, double Qb[3])
{
  Qb[0]=sin(i)*(Qa[0]*cos(phi)-Qa[1]*sin(phi))+Qa[2]*cos(i);
  Qb[1]=Qa[0]*sin(phi)+Qa[1]*cos(phi);
  Qb[2]=cos(i)*(Qa[1]*sin(phi)-Qa[0]*cos(phi))+Qa[2]*sin(i);
}

/* Funzione che calcola il prodotto fra il vettore `a' e lo scalare `c',
 * salvando il risultato nel vettore `b'. Il primo argomento è la lunghezza dei
 * due vettori `a' e `b'.
 */
void vettore_scalare(int n, double a[], double b[], double c)
{
  int i;
  for (i=0; i<n; i++)
    b[i]=c*a[i];
}
