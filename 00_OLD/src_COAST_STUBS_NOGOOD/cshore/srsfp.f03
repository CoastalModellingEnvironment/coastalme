!===============================================================================================================================
!
! This subroutine computes the surf similarity parameter
!
!===============================================================================================================================
subroutine SRFSP(L)
   use CShoreShared, only : TP, GRAV, SQR2, TWOPI, JR, WKPO, STHETA, CP, CTHETA, WN, HRMS, XB, ZB
   
   implicit none
   
   integer, intent(in) :: L
   double precision :: ARG, CO, HRMSO, SSP, TANB, THETAO
   
   CO = GRAV * TP / TWOPI
   ARG = ABS(CO / CP(1) * STHETA(1))
   
   ARG = MIN(ARG, SQR2 / 2.D0)         ! Arbitrary max deep water angle of 45 deg
   
   THETAO = DASIN(ARG)
   HRMSO = HRMS(1) * DSQRT((CP(1) * WN(1) * CTHETA(1))/ (0.5D0 * CO * DCOS(THETAO)))

   ! First guess at slope uses SWS slope
   if (JR == 1) then
      ! DFM Handle beginning of profile differently, to avoid invalid array index
      TANB = (ZB(JR+1, L) - ZB(JR, L)) / (XB(JR+1) - XB(JR))
   else if (JR == size(XB, 1)) then
      ! DFM Handle end of profile differently, to avoid invalid array index
      TANB = (ZB(JR, L) - ZB(JR-1, L)) / (XB(JR) - XB(JR-1))
   else
      ! DFM Use original approach for every other profile point   
      TANB = (ZB(JR+1, L) - ZB(JR-1, L)) / (XB(JR+1) - XB(JR-1))
   endif
   
   SSP = TANB / DSQRT(SQR2 * HRMSO / (TWOPI / WKPO))
   
   ! Just to improve slope estimate, estimate Runup with Mase 1989:
   ! R2P = SQR2 * HRMSO * 1.86D0 * SSP ** 0.71D0
   
   return
end subroutine SRFSP
