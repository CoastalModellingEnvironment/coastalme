!===============================================================================================================================
!
! This subroutine extrapolates vector F from node J1 to node J2 where values of F at nodes J= 1 to (J1-1) are computed
!
!===============================================================================================================================
subroutine EXTRAPO(J1, J2, F)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   integer, intent(in) :: J1, J2
   double precision, intent(out), dimension(:) ::  F

   ! NPE = number of points (computed in Subroutine BOTTOM) extrapolated to avoid a sudden jump from computed F(J1-1) to zero
   JJ=J1+NPE
   if (JJ <= J2) then
      JM = J1-1
      Y= F(JM)
      DELY = Y/DBLE(NPE+1)
      do 100 J= 1,NPE
         F(JM+J) = Y-DELY*DBLE(J)
100   end do
      F(JJ:J2) = 0.D0
   else
      if (J1 <= J2) F(J1:J2)=0.D0
   endif

   return
end subroutine EXTRAPO
