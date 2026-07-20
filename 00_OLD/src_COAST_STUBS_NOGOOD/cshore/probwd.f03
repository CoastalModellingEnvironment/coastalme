!===============================================================================================================================
!
! This subroutine computes bedload probability PBWD(J) and suspended load probability PSWD(J) in the wet and dry zone for Subroutine SEDTRA
!
!===============================================================================================================================
subroutine PROBWD(PW, A, US, UC, P)
   implicit none
   
   double precision, intent(in) :: PW, A, US, UC
   double precision, intent(out) :: P

   if (DABS(US) <= UC) then
      P = PW * DEXP(-A * (UC - US)**2)
   else
      if (US > UC) then
         P = PW
      else
         P = PW * (1.D0 - DEXP(-A * (UC + US)**2) + DEXP(-A * (UC - US)**2))
      endif
   endif

   return
end subroutine PROBWD
