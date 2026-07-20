!===============================================================================================================================
!
! This subroutine smooths the vector RAW using NPT computed in Subroutine BOTTOM where NPFS = (2NPT+1) = number of points for smoothing
!
!===============================================================================================================================
subroutine SMOOTH(NUM, RAW, F)
   use CShoreShared, only : NPT

   implicit none
   
   integer, intent(in) :: NUM
   double precision, intent(in), dimension(:) :: RAW
   double precision, intent(out), dimension(:) :: F
   
   integer :: J, JEND, JJ, JSTA, NPFS
   double precision :: SUM1, TOTJ
   
   do 201 J = 1, NUM
      JSTA = J - NPT
      JEND = J + NPT
      
      if (JSTA < 1) then
         JSTA = 1
         JEND = 2 * J - 1
      endif
      
      if (JEND > NUM) then
         JSTA = 2 * J - NUM
         JEND = NUM
      endif 
      
      NPFS = JEND - JSTA + 1
      TOTJ = DBLE(NPFS)
      SUM1 = 0.D0
      
      do 202 JJ = JSTA, JEND
         SUM1 = SUM1 + RAW(JJ)
202   end do

      F(J) = SUM1 / TOTJ
      
      !$$$  BDJ 2012-10-03 Remember to remove!
      !$$$        F(J) = RAW(J) !BDJ remove!
      !$$$        write (*,*) 'Smoothing is turned off!'
      !$$$  BDJ 2012-10-03
201 end do

   return
end subroutine SMOOTH
