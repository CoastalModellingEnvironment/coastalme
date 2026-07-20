!===============================================================================================================================
!
! This function relates to bottom shear stress in the wet and dry zone, it is called from Subroutine WETDRY
!
!===============================================================================================================================
function GBWD(R)
   double precision, intent(in) :: R
   double precision :: R2, GBWD, ERFCC
   
   interface
      function ERFCC(X)
         double precision, intent(in) :: X
      end function ERFCC
   end interface
   
   if (R >= 0.D0) then
      GBWD = 1.D0 + 1.77245D0 * R + R * R
   else
      R2 = R * R
      GBWD = 2.D0 * DEXP(-R2) - R2 - 1.D0 + 1.77245D0 * R * (3.D0 - 2.D0 * ERFCC(R))
   endif

   return
end function GBWD
