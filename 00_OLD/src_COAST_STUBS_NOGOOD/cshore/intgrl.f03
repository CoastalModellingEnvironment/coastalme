!===============================================================================================================================
!
! This subroutine computes an integral
!
!===============================================================================================================================
subroutine INTGRL(NUM, DEL, F, G)
   implicit none
   
   integer, intent(in) :: NUM
   double precision, intent(in) :: DEL
   double precision, intent(in), dimension(:) ::  F
   double precision, intent(out) :: G
   
   integer :: I, IMAX, NEND, NEVEN
   double precision :: DUM, SE, SO

   ! NUM can be even or odd integer
   IMAX = (NUM-1) / 2
   DUM = DBLE(NUM-1) / 2.D0
   
   if (DBLE(IMAX) < DUM) then
      NEND = NUM - 1
      NEVEN = 1
   else
      NEND = NUM
   endif
   
   SE = F(2)
   SO = 0.D0
   
   do 500 I = 2, IMAX
      SE = SE + F(I * 2)
      SO = SO + F(I * 2 - 1)
500 end do

   G = DEL / 3.D0 * (F(1) + 4.D0 * SE + 2.D0 * SO + F(NEND))
   if (NEVEN == 1) G = G + (F(NEND) + F(NUM)) * DEL / 2.D0
   
   return
end subroutine INTGRL
