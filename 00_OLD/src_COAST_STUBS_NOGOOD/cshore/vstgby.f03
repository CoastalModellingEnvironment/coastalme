!===============================================================================================================================
!
! This subroutine computes VSIGT = VMEAN / SIGT for specified GBY_IN, CTHETA_IN, USIGT = UMEAN / SIGT, and STHETA_IN
! DFM This subroutine computes VSIGT = VMEAN / SIGT for specified GBY_IN, CTHETA_IN, and STHETA_IN
!
!===============================================================================================================================
! DFM In original code, input parameter USIGT is passed but not used
!subroutine VSTGBY(CTHETA_IN, USIGT, STHETA_IN, VSIGT, GBY_IN)

subroutine VSTGBY(CTHETA_IN, STHETA_IN, VSIGT, GBY_IN)

   use CShoreShared, only : SQRG1

   implicit none
   ! DFM USIGT not used
   ! double precision, intent(in) :: CTHETA_IN, USIGT, STHETA_IN, GBY_IN
   double precision, intent(in) :: CTHETA_IN, STHETA_IN, GBY_IN
   double precision, intent(out) :: VSIGT
   double precision :: B, C, D

   VSIGT = 0.D0
   if (GBY_IN == 0.D0) goto 100
   
   B = SQRG1 * (1.D0 + STHETA_IN * STHETA_IN)
   C = GBY_IN
   
   if (GBY_IN > 0.D0) then
      D = B * B + 4.D0 * CTHETA_IN * C
      if (D >= 0.D0) VSIGT = 0.5D0 * (DSQRT(D) - B) / CTHETA_IN
      if (VSIGT < 0.D0) VSIGT = 0.D0
   else
      D = B * B - 4.0D0 * CTHETA_IN * C
      if (D >= 0.D0) VSIGT = 0.5D0 * (B - DSQRT(D)) / CTHETA_IN
      if (VSIGT > 0.D0) VSIGT = 0.D0
   endif

100 continue

   return
end subroutine VSTGBY
