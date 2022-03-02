! =============== ALNORM : calcul de proba normale ==========




      function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
      implicit none

      real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
      real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
      real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
      real ( kind = 8 ) alnorm
      real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
      real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
      real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
      real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
      real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
      real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
      real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
      real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
      real ( kind = 8 ), parameter :: con = 1.28D+00
      real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
      real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
      real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
      real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
      real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
      real ( kind = 8 ), parameter :: ltone = 7.0D+00
      real ( kind = 8 ), parameter :: p = 0.398942280444D+00
      real ( kind = 8 ), parameter :: q = 0.39990348504D+00
      real ( kind = 8 ), parameter :: r = 0.398942280385D+00
      logical up
      logical upper
      real ( kind = 8 ), parameter :: utzero = 18.66D+00
      real ( kind = 8 ) x
      real ( kind = 8 ) y
      real ( kind = 8 ) z

      up = upper
      z = x

      if ( z < 0.0D+00 ) then
         up = .not. up
         z = - z
      end if

      if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

         if ( up ) then
            alnorm = 0.0D+00
         else
            alnorm = 1.0D+00
         end if

         return

      end if

      y = 0.5D+00 * z * z

      if ( z <= con ) then

         alnorm = 0.5D+00 - z * ( p - q * y  &
         / ( y + a1 + b1/ ( y + a2 + b2 / ( y + a3 ))))

      else

         alnorm = r * exp ( - y )   &
          / ( z + c1 + d1 / ( z + c2 + d2 / ( z + c3 + d3  &
          / ( z + c4 + d4 / ( z + c5 + d5 / ( z + c6 ))))))

      end if

      if ( .not. up ) then
         alnorm = 1.0D+00 - alnorm
      end if

      return
      end function alnorm



! =================================================================
!              Densite d'une beta
!=================================================================

      double precision Function beta_densite(X,a,b)

      implicit none

      double precision :: beta,a,b,gammln,X

      beta=exp(gammln(a+b)-gammln(a)-gammln(b))

      beta_densite=((X)**(a-1))*((1-X)**(b-1))*beta

      return

      end Function beta_densite


! Calcul de gamma (gammln(a))

      double precision Function gammln(xx)

! retourne la valeur ln(gamma(xx)) pour xx>0
      implicit none

      integer::j
      double precision:: ser,stp,tmp,x,y,xx
      double precision,dimension(6)::cof

      save cof,stp

      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
      24.01409824083091d0, -1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      end do
      gammln=tmp+log(stp*ser/x)

      return

      end Function gammln







      

      ! =================================================================
!              CDF incomplete d'une beta
!=================================================================


      double precision Function betai(a,b,x)

      implicit none

      double precision :: a,b,x,bt,betaCF,gammln,temp

      if(x.lt.0.d0.or.x.gt.1.d0) then
            betai=999.d0
            return
      end if
      if (x.eq.0..or.X.eq.1.) then
         bt=0.
      else
         bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      end if

      if (x.lt.(a+1.)/(a+B+2.)) then
         temp=betaCF(a,b,x)
         if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=bt*temp/A
         return
      else
        temp=betaCF(b,a,1.-x)
        if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=1.-bt*temp/b
         return
      end if

      end Function betai

! betaCF est utilisï¿½ par betai

      double precision Function betaCF(a,b,x)

      implicit none

      integer ::m,m2
      integer,parameter ::maxit=100
      double precision ::a,b,x,aa,c,d,del,h,qab,qam,qap
      double precision,parameter::eps=3.e-7,fpmin=1.e-30


      qab=a+b
      qap=a+1
      qam=a-1
      c=1.       ! first step of Lentz's method
      d=1.-qab*x/qap
      if (abs(d).lt.fpmin) d=fpmin
      d=1./d
      h=d
      do m=1,maxit
         m2=2*m
         aa=m*(b-m)*x/((qam+m2)*(a+m2))
         d=1.+aa*d        ! one step (the even one) of the recurrence
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         h=h*d*c
         aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
         d=1.+aa*d        ! next step of the recurrence (the odd one)
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h=h*del
         if (abs(del-1.).lt.eps) goto 1
      end do

!  pause 'a or b too big, or maxit too small in betaCF'
            betaCF=999.d0
!            write(*,*)'problem'
            return

 1    betaCF=h
      return
      end Function betaCF




      !C ******************** BGOS ********************************


      SUBROUTINE BGOS(SX,ID,X1,X2,RO)


!C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)

      implicit none
      double precision ::RO,SX
      integer ::ID
      double precision ::F,V1,V2,S,DLS,RO2
      double precision ::X1,X2,UNIRAN
!C     write(*,*)'dans bgos'


 5    CONTINUE

!C     write(*,*)'avant rand :'

!C     X1=RAND()
!C     X2=RAND()

      X1=UNIRAN()
      X2=UNIRAN()

      IF(ID.NE.1) GO TO 10
      F=2.*SQRT(3.)
      X1=(X1-0.5)*F
      X2=(X2-0.5)*F
      GO TO 20
10    CONTINUE
      V1=2.*X1-1
      V2=2.*X2-1
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 5
      DLS=SQRT(-2.*LOG(S)/S)
      X1=V1*DLS
      X2=V2*DLS
20    CONTINUE
      RO2=RO*RO
      IF(ABS(RO).GT.1.E-10) X2=(X1+X2*SQRT(1./RO2-1.))*RO
      X1=X1*SX
      X2=X2*SX

!C      write(*,*) 'X1 ',X1,' X2 ',X2
!C OK, X1 et X2 sont créés

!C      write(*,*)'fin bgos'

      RETURN
      END subroutine bgos


!C ------------------- FIN SUBROUTINE BGOS -----------------


!C ------------------------------------------------------

      DOUBLE PRECISION FUNCTION UNIRAN()
!C
!C     Random number generator(RCARRY), adapted from F. James
!C     "A Review of Random Number Generators"
!C      Comp. Phys. Comm. 60(1990), pp. 329-344.
!C
      implicit none
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /      &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,      &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )

      return
      END function uniran





!!!!!!!!!!!!!!!!!!!  DMFSD et DSINV !!!!!!!!!!!!!!!!
      
      
      subroutine dmfsd(a,n,eps,ier)
!
!   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
!   MATRICE = TRANSPOSEE(T)*T
!   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
!            PAR COLONNE DE LA METRICE A FACTORISER
!   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
!
!   SUBROUTINE APPELE PAR DSINV
!
!   N : DIM. MATRICE
!   EPS : SEUIL DE TOLERANCE
!   IER = 0 PAS D'ERREUR
!   IER = -1 ERREUR
!   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(out)::ier
      double precision,intent(in)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

!
!   TEST ON WRONG INPUT PARAMETER N
!
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
!
!   INITIALIZE DIAGONAL-LOOP
!
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
!
!   CALCULATE TOLERANCE
!
          tol=dabs(eps*sngl(A(kpiv)))
!
!   START FACTORIZATION-LOOP OVER K-TH ROW
!
         do i=k,n
            dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
!
!   START INNER LOOP
!
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
               dsum=dsum+A(lanf)*A(lind)
            end do

!
!   END OF INNEF LOOP
!
!   TRANSFORM ELEMENT A(IND)
!
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
!   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
!


5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
!
!   COMPUTE PIVOT ELEMENT
!
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
!
!   CALCULATE TERMS IN ROW
!
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
      end do

!
!   END OF DIAGONAL-LOOP
!
      return
12    ier=-1
      return

      end subroutine dmfsd


!------------------------------------------------------------
!                            DSINV
!------------------------------------------------------------


      subroutine dsinv(A,N,EPS,IER,DET)

!
!     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
!
!     MATRICE = TRANSPOSEE(T)*T
!     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
!
!     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
!         STOCKEE COLONNE PAR COLONNE
!     DIM. MATRICE A INVERSER = N
!     DIM. TABLEAU A = N*(N+1)/2
!
!     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
!           COMME NUL
!
!     IER : CODE D'ERREUR
!         IER=0 PAS D'ERREUR
!         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
!         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(out)::ier
      double precision,intent(in)::eps
      double precision,intent(out),optional::det
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision::din,work
      integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf

!
!     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
!     A=TRANSPOSE(T) * T
!

      call dmfsd(A,n,eps,ier)

      det=0.d0

      if (ier.lt.0) goto 9
      if (ier.ge.0) det=0.d0
!
!     INVERT UPPER TRIANGULAR MATRIX T
!     PREPARE INVERSION-LOOP
!
!
! calcul du log du determinant

      do i=1,n
         det=det+dlog(A(i*(i+1)/2))
      end do
      det=2*det
      ipiv=n*(n+1)/2
      ind=ipiv
!
!     INITIALIZE INVERSION-LOOP
!
      do i=1,n
         din=1.d0/A(ipiv)
         A(ipiv)=din
         min=n
         kend=i-1
         lanf=n-kend
         if (kend.le.0) goto 5
         if (kend.gt.0) j=ind
!
!     INITIALIZE ROW-LOOP
!
         do k=1,kend
            work=0.d0
            min=min-1
            lhor=ipiv
            lver=j
!
!     START INNER LOOP
!
            do l=lanf,min
                lver=lver+1
                lhor=lhor+l
                work=work+A(lver)*A(lhor)
            end do
!
!     END OF INNER LOOP
!
            A(j)=-work*din
            j=j-min
         end do

!
!     END OF ROW-LOOP
!
5        ipiv=ipiv-min
         ind=ind-1
      end do

!
!     END OF INVERSION-LOOP
!
!     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
!     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
!     INITIALIZE MULTIPLICATION-LOOP
!
      do i=1,n
         ipiv=ipiv+i
         j=ipiv
!
!     INITIALIZE ROW-LOOP
!
         do k=i,n
            work=0.d0
            lhor=j
!
!     START INNER LOOP
!
            do l=k,n
                lver=lhor+k-i
                work=work+A(lhor)*A(lver)
                lhor=lhor+l
            end do
!
!     END OF INNER LOOP
!
            A(j)=work
            j=j+k
         end do
      end do

!
!     END OF ROW-AND MULTIPLICATION-LOOP
!
9     return
      end subroutine dsinv
