!===============================================================================================================================
!
! This subroutine calculates quantities based on linear wave theory
!
!===============================================================================================================================
subroutine LWAVE(J, L, WD, QDISP)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   integer, intent(in) :: J, L
   double precision, intent(in) :: WD, QDISP
   
   ! LINEAR WAVE parameters
   ! WD     = mean water depth from Main Program
   ! QDISP  = water flux affecting wave period
   ! TP     = representative wave period specified as input
   ! WKP    = wave number at node J
   ! WT(J)  = wave period at node J
   ! CP(J)  = phase velocity of peak frequency at node J
   ! WN(J)  = ratio of group velocity to phase velocity at node J

   ! Solve linear wave dispersion relation with no current to find WKP
   if (IWCINT == 0 .OR. QDISP == 0.D0) then
      D = WD*WKPO
      if (J == 1) then
         X = D/DSQRT(DTANH(D))
      else
         X = WKP*WD
      endif
      
10    COTH = 1.D0/DTANH(X)
      XNEW = X - (X-D*COTH)/(1.D0+D*(COTH**2.D0-1.D0))
      
      IF (DABS(XNEW - X) > 1.D-7) then
         X = XNEW
         goto 10
      endif
      
      AF = TWOPI/TP

      ! Solve linear wave dispersion relation with current to find WKP
   else
      B = TP*QDISP/TWOPI/WD/WD
      D = WD*WKPO
      
      if (J == 1) then
         X = D/DSQRT(DTANH(D))
      else
         X = WKP*WD
      endif
      
11    COTH = 1.D0/DTANH(X)
      C = 1.D0 - B*X
      F = X - D*C*C*COTH
      FD = 1.D0 + D*C*(2.D0*B*COTH + C*(COTH*COTH - 1.D0))
      XNEW = X - F/FD
      
      IF (DABS(XNEW - X) > 1.D-7) then
         X = XNEW
         goto 11
      endif
      AF = DSQRT(GRAV*XNEW*DTANH(XNEW)/WD)
   endif

   X = XNEW
   WKP = X/WD
   WN(J) = 0.5D0*(1.D0 + 2.D0*X/DSINH(2.D0*X))
   WT(J) = TWOPI/AF
   CP(J) = AF/WKP
   FSX = 2.D0*WN(J) - 0.5D0
   FSY = 0.D0
   FE  = WN(J)*CP(J)*WT(J)

   ! If IANGLE=0, normally incident waves
   if (IANGLE == 0) then
      STHETA(J) = 0.D0
      CTHETA(J) = 1.D0
      goto 100
   endif

   ! Otherwise, compute wave angle THETA in radians at node J using Snell's Law where ANGLE = incident wave angle in degrees at node J= 1 and WKPSIN = constant
   if (J == 1) then
      THETA = ANGLE*PI/180.D0
      STHETA(1) = DSIN(THETA)
      CTHETA(1) = DCOS(THETA)
      WKPSIN = WKP*STHETA(1)
   else
      STHETA(J) = WKPSIN/WKP
      
      ! DFM safety check bodge ==============================
      STHETA(J) = MIN(1.D0, STHETA(J))
      STHETA(J) = MAX(-1.D0, STHETA(J))
      ! DFM safety check bodge END ==========================

      THETA = DASIN(STHETA(J))
      CTHETA(J) = DCOS(THETA)
   endif

   FSX = FSX - WN(J)*STHETA(J)*STHETA(J)
   FSY = WN(J)*STHETA(J)*CTHETA(J)
   FE = FE*CTHETA(J)

100 if (IWCINT == 1) FE=FE+WT(J)*QWX/WD

   ! Compute RX, RY and RE related to roller momentum and energy fluxes as well as RBETA =wave-front slope
   if (IROLL == 1) then
      if (IANGLE == 0) then
         RX(J)=CP(J)/GRAV
         RE(J)=RX(J)*CP(J)
      else
         DUM=CP(J)*CTHETA(J)/GRAV
         RX(J)=DUM*CTHETA(J)
         RY(J)=DUM*STHETA(J)
         RE(J)=DUM*CP(J)
      endif
      RBETA(J)=RBZERO
      if (BSLOPE(J,L) > 0.D0) RBETA(J)=RBETA(J)+BSLOPE(J,L)*CTHETA(J)
   endif

   return
end subroutine LWAVE


