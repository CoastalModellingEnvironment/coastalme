!===============================================================================================================================
!
! This subroutine interpolates time series W1(I) specified at time levels T1(I) with I = 1, 2, ..., (N1+1) to obtain time series F(K) at time levels T2(K) with K = 1, 2, ..., (N2+1) and computes average value W2(K) with K = 1, 2, ..., N2 between time levels T2(K) and T2(K+1)
!
!===============================================================================================================================
subroutine TSINTP(N1, T1, W1, N2, T2, W2)
   use CShoreShared, only : NWAVE

   implicit none
   
   integer, intent(in) :: N1, N2
   double precision, intent(in), dimension(:) :: T1, W1, T2
   double precision, intent(out), dimension(:) :: W2
   
   double precision :: DUM
   double precision, allocatable, dimension(:) :: F   
   integer :: I, I1, K

   allocate(F(NWAVE+1))
   
   F(1) = W1(1)
   F(N2+1) = W1(N1+1)
   
   if (N2 >= 2) then
      do 100 K = 2, N2
         do 200 I = 1, N1
            I1 = I+1
            if (T1(I) <= T2(K) .and. T2(K) < T1(I1)) then
               DUM = (T2(K) - T1(I)) / (T1(I1) - T1(I))
               F(K) = (1.D0 - DUM) * W1(I) + DUM * W1(I1)
               goto 100
            endif
200      end do
100   end do
   endif

   do 300 K = 1, N2
      W2(K) = 0.5D0 * (F(K) + F(K+1))
300 end do

   deallocate(F)
   
   return
end subroutine TSINTP
