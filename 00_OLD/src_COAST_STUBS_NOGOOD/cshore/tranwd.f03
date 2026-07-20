!===============================================================================================================================
!
! This subroutine connects vector F1(J) with J= 1 to JR with vector F2(J) with J=JS to JE where JE is not less than JR
!
!===============================================================================================================================
subroutine TRANWD(F1, JR, F2, JS, JE)
   implicit none

   integer, intent(in) :: JR, JS, JE
   double precision, intent(in), dimension(:) :: F2
   double precision, intent(out), dimension(:) :: F1
   
   integer :: J
   double precision :: DSR, DUM

   if (JR > JS) then
      do 100 J = JS, JR
         ! F1(J) = 0.5D0 * (F1(J) + F2(J))   ! A mean value
         ! F1(J) = F2(J)                     ! No transition, just accept F2
         F1(J) = (F1(J) * DBLE(JR - J) + F2(J) * DBLE(J - JS)) / DBLE(JR - JS)
100   end do

      do 105 J = (JR+1), JE
         F1(J) = F2(J)
105   end do
   endif
   
   if (JR == JS) then
      F1(JR) = 0.5D0 * (F1(JR) + F2(JR))
   endif
   
   if (JR < JS) then
      DSR = DBLE(JS - JR)
      
      do 200 J = (JR+1), JS
         DUM = DBLE(JS - J) / DSR
         F1(J) = F1(JR) * DUM + F2(JS) * (1.D0 - DUM)
200   end do

      do 205 J = (JS+1), JE
         F1(J) = F2(J)
205   end do
   endif

   return
end subroutine TRANWD
