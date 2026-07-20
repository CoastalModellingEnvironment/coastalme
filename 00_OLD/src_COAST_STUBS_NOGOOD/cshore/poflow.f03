!===============================================================================================================================
!
! This subroutine computes mean and standard deviation of porous flow velocity and wave energy dissipation rate DPSTA for given PKHSIG and DEDX at node J in the wet zone
!
!===============================================================================================================================
subroutine POFLOW(J, L, PKHSIG, DEDX)
   use CShoreShared, only : UPMEAN, UPSTD, HP, DPSTA, WT, BESTA1, BESTA2, CTHETA, QP, ALSTA, SQRG1

   implicit integer (I-N)
   implicit double precision(A-H, O-Z)
   
   integer, intent(in) :: J, L
   double precision, intent(in) :: PKHSIG, DEDX
   
   ! For porous layer thickness HP(J,L)=0.0, no velocity and dissipation
   if (HP(J,L) == 0.D0) then
      UPMEAN(J) = 0.D0
      UPSTD(J) = 0.D0
      DPSTA(J) = 0.D0
   endif

   if (HP(J,L) > 0.D0) then
      A = 1.9D0*BESTA1
      B2 = BESTA2/WT(J)
      B = ALSTA + 1.9D0*B2
      UPSTD(J) = 0.5D0*(DSQRT(B*B+4.D0*A*PKHSIG)-B)/A
      A = SQRG1*(B2+BESTA1*UPSTD(J))
      C = CTHETA(J)*CTHETA(J)
      UPMEAN(J) = -DEDX/(ALSTA+A*(1.D0+C))

      ! To reduce numerical oscillations of UPMEAN(J), adjust
      RATIO = UPMEAN(J)/UPSTD(J)
      if (RATIO > 0.5D0) UPMEAN(J)=0.5D0*UPSTD(J)
      if (RATIO < -0.5D0) UPMEAN(J)=-0.5D0*UPSTD(J)
      QP(J)=UPMEAN(J)*HP(J,L)

      A2 = UPMEAN(J)*UPMEAN(J)
      B2 = UPSTD(J)*UPSTD(J)
      DPSTA(J) = HP(J,L)*(ALSTA*(A2+B2)+A*(2.D0*B2+A2*(1.D0+2.D0*C)))

   endif

   return
end subroutine POFLOW
