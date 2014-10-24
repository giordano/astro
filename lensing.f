c     Deprecated routine, kept for backward compatibility.
c     ==================================================================
      subroutine Position (t,r,phi,xi,eta)
      Implicit double precision (a-h,o-z)
      common /Par_Evento/ t0, u0, tE, theta, rho, vv

      print *, 'POSITION: this is a deprecated routine!'
      x0  = (t-t0)*vv           ! posizione del centro sorgente
      y0  = u0
      xs  = x0 + r*cos(phi)
      ys  = y0 + r*sin(phi)
      xi=xs*cos(theta)-ys*sin(theta)
      eta=xs*sin(theta)+ys*cos(theta)
      end subroutine Position
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Returns at time `tt' the (xi, eta) coordinates, in frame of
c     reference of the center of mass of the lens(es), of a point placed
c     at (rr, phi) polar cooordinates referred to the center of mass of
c     the source.  Thus, in order to get the position of the center of
c     the source with respect to the center of mass of the lens(es) at
c     time `tt', set rr = phi = 0.
c     Input: t0    (double precision): time of closest approach
c            u0    (double precision): distance of closest approach
c            theta (double precision): angle between the apparent
c                                      trajectory of the source in the
c                                      lens plane and the xi axis
c            vv    (double precision): apparent speed of the source,
c                                      Einstein radius/Einstein time
c            tt    (double precision): current time
c            rr    (double precision): radial polar coordinate from source
c                                    center
c            phi (double precision): angular polar coordinate from
c                                    source center
c     Output: xi  (double precision): xi coordinate of the given point
c             eta (double precision): eta coordinate of the given point
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine source_position (t0, u0, theta, vv, tt, rr, phi, xi,
     &    eta)
      implicit none
      double precision t0, u0, theta, vv, tt, rr, phi, xi, eta, xs, ys

c     Position of the source center in the non rotated frame of
c     reference.
      x0  = (tt-t0)*vv
      y0  = u0
c     Position of the point at (rr, phi) from the source center.
      xs  = x0 + rr*cos(phi)
      ys  = y0 + rr*sin(phi)
c     Rotate coordinates by theta.
      xi  = xs*cos(theta) - ys*sin(theta)
      eta = xs*sin(theta) + ys*cos(theta)
      end subroutine source_position
c     ==================================================================

c     ==================================================================
      subroutine paczynski_amplification (u, ampl)
      implicit none
      double precision u, ampl
      ampl = (u**2 + 2d0)/(u*sqrt(u**2 + 4d0))
      end subroutine paczynski_amplification
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Routine per il calcolo dell'amplificazione di una sorgente estesa
c     a un fissato istante di tempo.
c     Input: tempo (double precision): istante di tempo a cui calcolare
c            l'amplificazione.
c     Output: hex_amplification (double precision): amplificazione di
c             una sorgente estesa calcolata con il metodo
c             dell'esadecapolo.
c             A0 (double precision): amplificazione di una sorgente
c             puntiforme, calcolata con il metodo di Witt & Mao.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine hexadecapole_amplification (b, q, xi, eta, chi,
     &    hex_amplification, A0)
      implicit none
      common /Par_Evento/ t0, u0, tE, theta, rho, vv
      common /constants/ pi
      common /limb_darkening/ Gamma
      common /Par_ampli/ tfirst, tlast, NStepTime

      double precision hex_amplification, q, b, t0, u0, te, theta, rho,
     &    xi, eta, pi, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10,a11,
     &    a12, Ahalfrho, Arhoc, Arhop, A2rho2, A4rho4, xipm, etapm, chi,
     &    tfirst, tlast,  Gamma, phi, vv
      integer NStepTime

c     phi: angolo rispetto alle direzioni cardinali.
      phi = 0d0

c     La coordinata `xi' che riceviamo da fuori è la posizione della
c     sorgente nel sistema di riferimento del centro di massa delle
c     lenti.  Il metodo di Witt e Mao ha bisogno delle coordinate della
c     sorgente rispetto al punto medio delle lenti, quindi a Witt e Mao
c     passiamo le coordinate `xipm' e `etapm' nel sistema di riferimento
c     del punto medio delle lenti.  Le operazioni da fare sono una
c     rotazione dell'angolo `-chi' e una traslazione lungo `xi' pari
c     alla `xi' del centro di massa.
      xipm  = xi*cos(-chi) - eta*sin(-chi) + (b/2d0)*(q-1d0)/(1d0+q)
      etapm = xi*sin(-chi) + eta*cos(-chi)

      call witt_mao(xipm, etapm, b, q, A0)
      call witt_mao(xipm + rho/2d0*cos(phi), etapm + rho*sin(phi), b, q,
     &    A1)
      call witt_mao(xipm + rho/2d0*cos(phi + pi/2d0), etapm + rho /2d0
     &    *sin(phi+ pi/2d0), b, q, A2)
      call witt_mao(xipm + rho/2d0*cos(phi + pi), etapm + rho/2d0
     &    *sin(phi +pi), b, q, A3)
      call witt_mao(xipm + rho/2.*cos(phi + 3d0*pi/2d0), etapm + rho/2d0
     &    *sin(phi + 3d0*pi/2d0), b, q, A4)
      call witt_mao(xipm + rho*cos(phi), etapm + rho*sin(phi), b, q, A5)
      call witt_mao(xipm + rho*cos(phi + pi/2d0), etapm + rho*sin(phi +
     &    pi/2d0),b, q, A6)
      call witt_mao(xipm + rho*cos(phi + pi), etapm + rho*sin(phi + pi),
     &    b,q, A7)
      call witt_mao(xipm + rho*cos(phi + 3d0*pi/2d0), etapm + rho
     &    *sin(phi +3d0*pi/2d0), b, q, A8)
      call witt_mao(xipm + rho*cos(phi + pi/4d0), etapm + rho*sin(phi +
     &    pi/4d0),b, q, A9)
      call witt_mao(xipm + rho*cos(phi + 3d0*pi/4d0), etapm + rho
     &    *sin(phi +3d0*pi/4d0), b, q, A10)
      call witt_mao(xipm + rho*cos(phi + 5d0*pi/4d0), etapm + rho
     &    *sin(phi +5d0*pi/4d0), b, q, A11)
      call witt_mao(xipm + rho*cos(phi + 7d0*pi/4d0), etapm + rho
     &    *sin(phi +7d0*pi/4d0), b, q, A12)

c     A_{rho/2, +} in Gould, vedi equazione (7)
      Ahalfrho = (A1 + A2 + A3 + A4)/4d0 - A0
c     A_{rho, +} in Gould, vedi equazione (7)
      Arhop = (A5 + A6 + A7 + A8)/4d0 - A0
c     A_{rho, ×} in Gould, vedi equazione (7)
      Arhoc = (A9 + A10 + A11 + A12)/4d0 - A0
c     A2*rho**2, vedi equazione (9) di Gould
      A2rho2 = (16d0*Ahalfrho - Arhop)/3d0
c     A4*rho**4, vedi equazione (9) di Gould
      A4rho4 = (Arhop + Arhoc)/2d0 - A2rho2
c     amplificazione, vedi equazione (6) di Gould
      hex_amplification = A0 + A2rho2*(1d0 - Gamma/5d0)/2d0 + A4rho4
     &    *(1d0 -11d0*Gamma/35d0)/3d0

      end subroutine hexadecapole_amplification
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Scrive su file le coordinate dei punti che formano le caustiche,
c     calcolati con la routine `CalCaustiche'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine  write_caustic(b, q, chi, m)
      implicit double precision (a-h,o-y)
      parameter(N=2)
      parameter(N_theta=3600)
      common/parca/theta(N_theta),xi(2*N,N_theta),eta(2*N,N_theta)
      dimension xitemp(2*N,N_theta), etatemp(2*N,N_theta)
      character(50) filename
      write(filename,'(a,i0.3,a)') 'curves/caustics', m, '.res'

c     Le caustiche sono calcolate con il metodo di Witt e Mao che è
c     centrato nel punto medio delle lenti, però vogliamo scrivere il
c     risultato nel sistema di riferimento del centro di massa delle
c     lenti.  Allora nel file scriviamo le coordinate `xitemp' ed
c     `etatemp' ottenute da `xi' ed `eta' mediate una traslazione lungo
c     l'asse `chi' pari alla posizione del centro di massa seguita da
c     una rotazione di `chi'
      do j=1,N_theta
        do i=1,2*N
          xitemp(i,j)=(xi(i,j)-(b/2d0)*(q-1)/(q+1))*cos(chi)-eta(i,j)
     &        *sin(chi)
          etatemp(i,j)=(xi(i,j)-(b/2d0)*(q-1)/(q+1))*sin(chi)+eta(i,j)
     &        *cos(chi)
        enddo
      enddo

      open(unit=10,file=filename)
      icount=1
      do j=1,N_theta
        do i=1,2*N
          write(10,'(3f17.8)') theta(j), xitemp(i,j), etatemp(i,j)
          icount=icount+1
        enddo
      enddo
      close(10)
      end subroutine  write_caustic
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Scrive su file le coordinate dei punti che formano le curve
c     critiche, calcolati da `CalCaustiche'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine write_criticcurves(m)
      implicit double precision (a-h,o-y)
      parameter(N=2)
      parameter(N_theta=3600)
      common/parcc/xcr(2*N,N_theta),ycr(2*N,N_theta)
      character(50) filename
      write(filename,'(a,i0.3,a)') 'curves/criticcurves', m, '.res'
      open(unit=10,file=filename)
      icount=1
      do j=1,N_theta
        do i=1,2*N
          write(10,'(2f17.8)') xcr(i,j),ycr(i,j)
          icount=icount+1
        enddo
      enddo
      close(10)
      end subroutine write_criticcurves
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Risolve l'equazione delle curve critiche, mappando poi le
c     soluzioni anche nel piano della sorgente, per ottenere le
c     soluzioni che descrivono le caustiche.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine calCaustiche(b0,q)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter(N=2)
      parameter(N_theta=3600)
      complex*16 Zsol(2*N),Zcoeff(2*N+1),r
c     Variabili complesse in singola precisione da usare nel caso in cui
c     si utilizzi zroots.for per la ricerca degli zeri.
      complex zsol_single(2*N), zcoeff_single(2*N+1)
      double precision m1,m2
      logical polish
      common/parca/theta(N_theta),xi(2*N,N_theta),eta(2*N,N_theta)
      common/parcc/xcr(2*N,N_theta),ycr(2*N,N_theta)
      common /constants/ pi
      common /mode_zroots/ mode_zroots

      polish=.true.

      grado=2d0*N
      Theta_step=360d0/float(N_theta)
c     Parametri del sistema: masse e distanze
      m1 =   q/(1d0+q)          !   pianeta
      m2 = 1d0/(1d0+q)          !   STELLA
      z1=cmplx(b0/2d0,0d0)
      z2=cmplx(-b0/2d0,0d0)
      z1q=z1*z1
      z2q=z2*z2

c     Costruzione del vettore di angoli
      do i=1,N_theta
        theta(i)=(i-1d0)*Theta_step*(pi/180d0)
      enddo

c     Coefficienti polinomio complesso
      do j=1,N_theta
        r=cmplx(cos(theta(j)),-sin(theta(j)))
        Zcoeff(1)= m1*z2q + m2*z1q - r*z1q*z2q
        Zcoeff(2)= 2d0*(-m2*z1 - m1*z2 + r*z1*z2q + r*z1q*z2)
        Zcoeff(3)= m2 + m1 - r*z2q - r*z1q - 4d0*r*z1*z2
        Zcoeff(4)= 2d0*r*(z1 + z2)
        Zcoeff(5)= -r
        if (mode_zroots.eq.0) then
          call cmplx_roots_gen(zsol, zcoeff, 2*N, polish, .false.)
        else if (mode_zroots.eq.1) then
          zcoeff_single = zcoeff
          call zroots(zcoeff_single, 2*N, zsol_single, polish)
          zsol = zsol_single
        end if

        do i=1,2*N
c     Curve critiche
          xcr(i,j)=real(Zsol(i))
          ycr(i,j)=imag(Zsol(i))

          a=xcr(i,j)-b0/2d0
          b=xcr(i,j)+b0/2d0
          y2=ycr(i,j)*ycr(i,j)
          den1 = a*a+y2
          den2 = b*b+y2

c     Curve caustiche
          xi(i,j)=xcr(i,j)-m1*a/den1-m2*b/den2
          eta(i,j)=ycr(i,j)-m1*ycr(i,j)/den1-m2*ycr(i,j)/den2

        enddo
      enddo
      end subroutine calCaustiche
c     ==================================================================

cccccccccccccc ROUTINE PER IL CALCOLO DELLE SOLE CURVE CRITICHE cccccccc
c     Serve per un calcolo più veloce delle curve critiche, che fungono
c     poi da input alla routine `distance'.  L'algoritmo è identico a
c     `CalCaustiche', senza però mappare nel piano della sorgente.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ==================================================================
      subroutine calCritCur(b0,q)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter(N=2)
      parameter(N_theta=3600)
      complex*16 Zsol(2*N),Zcoeff(2*N+1),r
c     Variabili complesse in singola precisione da usare nel caso in cui
c     si utilizzi zroots.for per la ricerca degli zeri.
      complex zsol_single(2*N), zcoeff_single(2*N+1)
      double precision m1,m2
      logical polish
      common/parca/theta(N_theta),xi(2*N,N_theta),eta(2*N,N_theta)
      common/parcc/xcr(2*N,N_theta),ycr(2*N,N_theta)
      common /constants/ pi
      common /mode_zroots/ mode_zroots

      polish=.true.

      grado=2d0*N
      Theta_step=360d0/float(N_theta)
c     Parametri del sistema: masse e distanze
      m1 =   q/(1d0+q)          !   pianeta
      m2 = 1d0/(1d0+q)          !   STELLA
      z1=cmplx(b0/2d0,0)
      z2=cmplx(-b0/2d0,0)
      z1q=z1*z1
      z2q=z2*z2

c     Costruzione del vettore di angoli
      do i=1,N_theta
        theta(i)=(i-1)*Theta_step*(pi/180d0)
      enddo

c     Coefficienti polinomio complesso
      do j=1,N_theta
        r=cmplx(cos(theta(j)),-sin(theta(j)))
        Zcoeff(1)= m1*z2q + m2*z1q - r*z1q*z2q
        Zcoeff(2)= 2d0*(-m2*z1 - m1*z2 + r*z1*z2q + r*z1q*z2)
        Zcoeff(3)= m2 + m1 - r*z2q - r*z1q - 4d0*r*z1*z2
        Zcoeff(4)= 2d0*r*(z1 + z2)
        Zcoeff(5)= -r
        if (mode_zroots.eq.0) then
          call cmplx_roots_gen(zsol, zcoeff, 2*N, polish, .false.)
        else if (mode_zroots.eq.1) then
          zcoeff_single = zcoeff
          call zroots(zcoeff_single, 2*N, zsol_single, polish)
          zsol = zsol_single
        end if
        do i=1,2*N
c     Curve critiche
          xcr(i,j)=real(Zsol(i))
          ycr(i,j)=imag(Zsol(i))
        enddo
      enddo
      end subroutine calCritCur
c     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Witt_Mao
c     Input (NOTA: le coordinate devono essere date nel sistema di
c     riferimento del punto medio delle lenti):
c       xi  = coordinata xi attuale della sorgente
c       eta = coordinata eta attuale della sorgente
c       b   = separazione tra le lenti
c       q   = rapporto di massa tra le lenti
c     Output: Ampli(xi, eta)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Witt_Mao (xi, eta, b, q, Ampli)
      implicit double precision (a-h,o-v)
      implicit complex*16 (z)
      common /mode_zroots/ mode_zroots
      double precision xi, eta
      double precision gm1, gm2, dm, sm
      logical first_3_roots_order_changed
c     Variabili complesse in singola precisione da usare nel caso in cui
c     si utilizzi zroots.for per la ricerca degli zeri.
      complex zi_single(5), zc_single(6)

c     Complex coefficients of 5th order polynomial
      dimension zc(6)
c     Image positions
      dimension zi(5)

c     Highly optimized code (e.g., code compiled with gfortran and
c     `-Ofast' flag) may return wrong amplification values.  The
c     following dummy `write' prints nothing and should work around this
c     bug.  Please, DO NOT remove it!
      write (*, fmt='(a)', advance="no") ''

      gm1 =   q/(1d0+q)         ! pianeta
      gm2 = 1d0/(1d0+q)         ! STELLA
      dm  = (gm2-gm1)           !/2.0
      sm  = (gm1+gm2)           !/2.0

c     Criterion for checking the validity of solutions: Map back the
c     images to source position - if solution within EP from source
c     position then image is declared valid.
      ep = 1d-3                ! 0.0001 ! 0.001 !0.1!5.e-2
c     half binary separation
      d = b*0.5d0
c     lens configuration
      Z1 = CMPLX( d, 0d0)
      Z2 = CMPLX(-d, 0d0)
c     source position
      ZS = CMPLX(xi, eta)

c     complex conjugates, recall that `Z1' and `Z2' are real.
      Z1C = Z1
      Z2C = Z2
      ZSC = CONJG(zs)

c     complex coefficients of 5th order polynomial
      Z1_Z1   = Z1*Z1
      ZSC_ZSC = ZSC*ZSC
      ZC(6)   = Z1_Z1 - ZSC_ZSC

      ZF     = ZS*ZSC_ZSC - Z1*(DM + Z1*ZS)
      ZSC_SM = ZSC*SM
      ZC(5)  = ZF - ZSC_SM

      ZA     = 2d0*SM*ZS*ZSC
      Z1_Z6  = Z1*ZC(6)
      ZB     = 2d0*DM*ZSC
      ZC(4)  = ZA + Z1*(ZB - 2d0*Z1_Z6)

      SDM   = SM*DM
      DDM   = DM*DM
      SSM   = SM*SM
      ZS_ZB = ZS*ZB
      ZC(3) = SSM*ZS + Z1*(SDM - ZS_ZB - 2d0*Z1*ZF)

      ZC(2) = -Z1*(2d0*SDM*ZS + Z1*(DDM + SSM + ZA + Z1*(ZB - Z1_Z6)))

      ZC(1) = Z1_Z1*(DDM*ZS + Z1*(SDM + ZS_ZB + Z1*(ZF + ZSC_SM)))

c     Find roots.  `mode_zroots' selects the polynomial solver routine.
      if (mode_zroots.eq.0) then
        call cmplx_roots_5(zi, first_3_roots_order_changed, zc, .true.)
      else if (mode_zroots.eq.1) then
        zc_single = zc
        call zroots(zc_single, 5, zi_single, .true.)
        zi = zi_single
      end if
      eps2=1d-6
c     Initialize amplification.
      Ampli = 0d0
c     Sum amplification over images.
      DO I = 1, 5
        ZIC  = CONJG(ZI(I))
        ZDC1 = ZIC - Z1C
        ZDC2 = ZIC - Z2C
c     ZE map back position of source
        ZE   = ZI(I) - GM1/ZDC1 - GM2/ZDC2

c     reject invalid image
        IF(ZABS(ZE - ZS) .LE. EP) THEN
c     Determinant of the amplification matrix.
          detJ = 1d0 - ZABS(GM1/ZDC1**2 + GM2/ZDC2**2)**2
c     sum amplification
          if(abs(detJ) .eq. eps2) then
            Ampli = Ampli + 1d30
          else
            Ampli = Ampli + abs(1d0/detJ)
          end if
        END IF
      END DO
      END subroutine Witt_Mao
c     ==================================================================

c     ==================================================================
c     Routine to calculate the amplification using the approximation by
c     Lee et al (doi:10.1088/0004-637X/695/1/200) for a finite source,
c     without limb darkening.  Arguments:
c       u   (in) = distance between source and lens
c       rho (in) = source radius projected onto the lens projected onto
c                  the lens plane, in Einstein radii
c       n   (in) = number of points to be used in the Cavalieri-Simpson
c                  rule
c       a  (out) = amplification
      subroutine lee_uniform(u, rho, n, a)
      implicit none
      double precision rho, a, f_lee, u, pi_2, pi, sum1, sum2
      integer n, k
      external f_lee
      pi = 3.14159265358979323846264338328d0
      pi_2 = 1.57079632679489661923132169164d0
      sum1 = 0d0
      sum2 = 0d0
      if (u.le.rho) then
        do k = 1, n-1
          sum1 = sum1 + f_lee(k*pi/n, u, rho)
          sum2 = sum2 + f_lee((2d0*k-1d0)*pi_2/n, u, rho)
        end do
        sum2 = sum2 + f_lee((2d0*n-1d0)*pi_2/n, u, rho)
        a = (((u + rho)*sqrt((u + rho)**2 + 4d0) - (u - rho)*sqrt((u -
     &      rho)**2 + 4d0))/3d0 + 2d0*sum1/3d0 + 4d0*sum2/3d0)/(rho**2
     &      *2d0*n)
      else
        do k = 1, n/2-1
          sum1 = sum1 + f_lee(2d0*k*asin(rho/u)/n, u, rho)
          sum2 = sum2 + f_lee((2d0*k - 1d0)*asin(rho/u)/n, u, rho)
        end do
        sum2 = sum2 + f_lee((n - 1d0)*asin(rho/u)/n, u, rho)
        a = asin(rho/u)/(pi*rho**2*n)*(((u + rho)*sqrt((u + rho)**2 +
     &      4d0) - (u - rho)*sqrt((u -rho)**2 + 4d0))/3d0 + 2d0*sum1/3d0
     &      + 4d0*sum2/3d0)
      end if
      end subroutine lee_uniform
c     ==================================================================

c     ==================================================================
      double precision function u1(theta, u, rho)
      implicit none
      double precision theta, u, rho
      if ((u.gt.rho) .and. (theta.le.asin(rho/u))) then
        u1 = u*cos(theta) - sqrt(rho**2 - u**2*sin(theta)**2)
      else
        u1 = 0
      end if
      return
      end function u1
c     ==================================================================

c     ==================================================================
      double precision function u2(theta, u, rho)
      implicit none
      double precision theta, u, rho
      if ((u.gt.rho) .and. (theta.gt.asin(rho/u))) then
        u2 = 0
      else
        u2 = u*cos(theta) + sqrt(rho**2 - u**2*sin(theta)**2)
      end if
      return
      end function u2
c     ==================================================================

c     ==================================================================
      double precision function f_lee(theta, u, rho)
      implicit none
      double precision theta, u, rho, u1, u2
      external u1, u2
      f_lee = u2(theta, u, rho)*sqrt(u2(theta, u, rho)**2 + 4d0)
     &    -u1(theta, u, rho)*sqrt(u1(theta, u, rho)**2 +4d0)
      return
      end function f_lee
c     ==================================================================

c     Routine di J. Skowron & A. Gould per la ricerca degli zeri di
c     polinomi a coefficienti complessi.  Vedi
c     http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
c     http://arxiv.org/abs/1203.1034
      include 'cmplx_roots_sg_77.f'

c     Routine per la ricerca degli zeri di polinomi a coefficienti
c     complessi incluse in Numerical Recipes.
      include 'zroots.for'
      include 'laguer.for'
