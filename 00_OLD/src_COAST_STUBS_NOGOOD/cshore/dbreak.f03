!===============================================================================================================================
!
! This subroutine calculates QBREAK and DBSTA for wave breaking
!
!===============================================================================================================================
subroutine DBREAK(J, L, WHRMS, D)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   integer, intent(in) :: J, L
   double precision, intent(in) :: WHRMS, D
   
   ! Calculate energy dissipation factor ABREAK(J) for steep slope where D = mean water depth from main program
   ABREAK(J) = (TWOPI/WKP/D)*BSLOPE(J,L)*CTHETA(J)/3.D0
   if (ABREAK(J) < 1.D0) ABREAK(J) = 1.D0
   
   ! mg  Allow for variable gamma
   if (GAMMA < 0) then
      ! mg  Compute deep water wave height
      CO = GRAV*TP/TWOPI
      THETAO=DASIN(CO/CP(1)*STHETA(1))
      HRMSO = HRMS(1)*DSQRT((CP(1)*WN(1)*CTHETA(1))/ (0.5D0*CO*DCOS(THETAO)))
      ! mg  Alex Apotsos et al. 2008, Coastal Engineering 55 (2008) 224-235.  (Eq 23)
      GAMMA_TEMP = 0.3D0 + 0.45D0 * TANH(0.9D0 * HRMSO)
   else
      GAMMA_TEMP = GAMMA
   endif

   ! FRACTION OF BREAKING WAVES AND ASSOCIATED DISSIPATION
   ! QBREAK(J) = Fraction of breaking waves at node J
   ! DBSTA(J)  = Time averaged normalized energy dissipation due to
   ! wave breaking at node J
   ! mg
   ! mg  HM = 0.88D0/WKP*DTANH(GAMMA*WKP*D/0.88D0)
   HM = 0.88D0/WKP*DTANH(GAMMA_TEMP*WKP*D/0.88D0)
   
   ! Compute QBREAK = fraction of breaking waves
   B = (WHRMS/HM)**2.D0
   
   ! bdj if (B.LT.0.99999D0) then
   if (B < 0.99999D0 .and. WHRMS > 1D-10) then ! BDJ
      QBOLD = B/2.D0
10    QBREAK(J) = QBOLD - (1.D0-QBOLD + B*DLOG(QBOLD))/(B/QBOLD-1.D0)
      if (QBREAK(J) <= 0.D0) QBREAK(J) = QBOLD/2.D0
      if (DABS(QBREAK(J)-QBOLD) > 1.D-6) then
         QBOLD = QBREAK(J)
         goto 10
      endif
   else
      QBREAK(J) = 1.D0
      HM=WHRMS
   endif

   DBSTA(J) = 0.25D0*ABREAK(J)*QBREAK(J)*HM*HM/WT(J)

   ! Reduce SIGSTA if WHRMS is larger than GAMMA*D
   ! (used only for wave transmission over submerged breakwater)
   ! GAMD = GAMMA*D
   ! if (WHRMS.LE.GAMD) then
   SISMAX = 1.D0
   ! else
   ! SISMAX = DSQRT(GAMMA*WHRMS/D/8.0D0)
   ! endif

   return
end subroutine DBREAK

