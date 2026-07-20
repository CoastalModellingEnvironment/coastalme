!===============================================================================================================================
!
! This subroutine computes GBX_IN and CF_IN for specified CTHETA_IN, USIGT, STHETA_IN and VSIGT for Gaussian variable R
!
!===============================================================================================================================
subroutine GBXAGF(CTHETA_IN, USIGT, STHETA_IN, VSIGT, GBX_IN, CF_IN)
   use CShoreShared, only : IANGLE, SQRG1, SQRG2, SQR2

   implicit none
   
   double precision, intent(in) :: CTHETA_IN, USIGT, STHETA_IN, VSIGT
   double precision, intent(out) :: GBX_IN, CF_IN
   double precision :: RM, AFM, DUM, C1, C2, C3, ERFCC
   
   interface
      function ERFCC(X)
         double precision, intent(in) :: X
      end function ERFCC
   end interface
   
   ! For obliquely incident waves, use approximate equations
   if (IANGLE == 1) then
      RM  = -USIGT * CTHETA_IN - VSIGT * STHETA_IN
      AFM = DABS(VSIGT * CTHETA_IN - USIGT * STHETA_IN)
      DUM = USIGT * USIGT + VSIGT * VSIGT
      GBX_IN = SQRG1 * (USIGT - RM * CTHETA_IN) + USIGT * AFM
      CF_IN  = SQRG2 + (1.D0 + DUM) * AFM + SQRG1 * (DUM + 2.D0 * RM * RM)
   endif

   ! For normally incident waves, use analytical expressions involving complementary error function ERFCC
   if (IANGLE == 0) then
      C1 = 1.D0 - ERFCC(USIGT / SQR2)
      C2 = SQRG1 * DEXP(-USIGT * USIGT / 2.D0)
      C3 = 1.D0 + USIGT * USIGT
      GBX_IN = C3 * C1 + C2 * USIGT
      CF_IN = USIGT * (C3 + 2.D0) * C1 + (C3 + 1.D0) * C2
   endif

   return
end subroutine GBXAGF
