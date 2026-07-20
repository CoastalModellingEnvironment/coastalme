!===============================================================================================================================
!
! This subroutine computes the bottom elevation change using the volume conservation of bottom sediment
!
!===============================================================================================================================
subroutine CHANGE(ITIME, L, IEND, ICALL)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   interface 
      subroutine SMOOTH(NUM, RAW, F)
         integer, intent(in) :: NUM
         double precision, intent(in), dimension(:) :: RAW
         double precision, intent(out), dimension(:) :: F
      end subroutine SMOOTH
      
      subroutine INTGRL(NUM, DEL, F, G)
         integer, intent(in) :: NUM
         double precision, intent(in) :: DEL
         double precision, intent(in), dimension(:) ::  F
         double precision, intent(out) :: G
      end subroutine INTGRL
   end interface   

   integer, intent(in) :: ITIME, L, ICALL
   integer, intent(out) :: IEND
   
   double precision, allocatable, dimension(:) :: DZBDT, CB, R, DELZBRW, DELZBJ, V, AVY
   double precision, allocatable, dimension(:) :: VDY, ADZX
   
   allocate(DZBDT(JMAXMAX), CB(JMAXMAX), R(JMAXMAX), DELZBRW(JMAXMAX), DELZBJ(JMAXMAX), V(JMAXMAX), AVY(JMAXMAX))
   allocate(VDY(ILINE), ADZX(ILINE))

   ! If ICALL == 1, alongshore uniformity is assumed for profile change. Compute the first-order rate of the bottom elevation change where sediment transport rate Q(J) computed in Subroutine SEDTRA where the seaward boundary location at node 1 is chosen such that bottom change is negligible seaward of node 1. Also at node JMAX
   if (ICALL == 1) then
      JMAXM1 = JMAX(L) - 1
      DZBDT(1) = 0.D0
      DZBDT(JMAX(L)) = 0.D0
      do 100 J = 2, JMAXM1
         DZBDT(J) = (Q(J-1) - Q(J+1)) / DX2
100   end do
      if (IVWALL(L) == 2) DZBDT(JMAXM1) = (Q(JMAX(L) - 2) - Q(JMAXM1)) / DX
      ! Where backward finite difference is used at the node next to the wall if the vertical wall is exposed to wave action.

      ! Find the time step DELT using the numerical stability criterion but the value of DELT is limited by the end time TIMEBC(ITIME+1) for given TIME
      ! Compute CB(J)=bottom profile phase velocity
      ! CBMAX = 0.001D0
      ! Increase of CBMAX tends to reduce DELT and improve numerical stability
      CBMAX = 0.004D0
      DZBMAX = 0.1D0 * DX
      do 115 J = 1, JMAX(L)
         if (J == 1) then
            DELQ = Q(2) - Q(1)
            DELZB1 = ZB(2,L) - ZB(1,L)
         elseif (J == JMAX(L)) then
            J1 = JMAXM1
            DELQ = Q(J) - Q(J1)
            DELZB1 = ZB(J,L) - ZB(J1,L)
         else
            JP1 = J + 1
            JM1 = J - 1
            DELQ = Q(JP1) - Q(JM1)
            DELZB1 = ZB(JP1,L) - ZB(JM1,L)
         endif
         if (DABS(DELZB1) > DZBMAX) then
            CB(J) = DELQ/DELZB1
         else
            CB(J) = 0.D0
         endif
         DUMC = DABS(CB(J))
         if (DUMC > CBMAX) CBMAX = DUMC
115   end do
      DELT = DX / CBMAX
      IDUM = ITIME + 1
      DUM = (TIMEBC(IDUM) - TIMEBC(ITIME)) / 2.D0
      if (DELT > DUM) DELT=DUM

      DUM = TIME+DELT
      if (DUM > TIMEBC(IDUM)) then
         DELT = TIMEBC(IDUM) - TIME
         IEND = 1
      endif

      ! Compute DELZBRW(J)=first-order bottom elevation change before smoothing
      do 120 J = 1, JMAX(L)
         DELZBRW(J) = DELT * DZBDT(J)
120   end do

      ! Add second-order correction to DELZBRW(J)
      DTDT = DELT * DELT
      do 121 J= 1, JMAX(L)
         R(J) = DTDT * CB(J) * CB(J) / DXDX
121   end do
      do 122 J = 2, JMAXM1
         JP1 = J + 1
         JM1 = J - 1
         DUM = ZB(JP1,L) * (R(JP1) + R(J)) / 4.D0 - ZB(J,L) * (R(J) / 2.D0 + (R(JP1) + R(JM1)) / 4.D0) + ZB(JM1,L) * (R(J) + R(JM1)) / 4.D0
         DELZBRW(J) = DELZBRW(J) + DUM

         ! if (IPOND.EQ.0) then
         !   if (IVWALL(L).EQ.2.and.J.GT.JWD) then
         !     DUM = SWLBC(ITIME) - ZB(J,L)
         !     if (DELZBRW(J).LT.DUM) DELZBRW(J) = DUM
         !   endif
         ! endif
122   end do

      JMAXL=JMAX(L)
      
      ! Smoothing DELZBRW using Subroutine SMOOTH      
      call SMOOTH(JMAXL, DELZBRW, DELZBJ)

      ! Adjust smoothed bottom elevation change DELZB to satisfy the volume conservation between J = 1 to JMAX. If ISEDAV == 1 (hard bottom), sediment layer thickness HP(J,L) must be positive or zero
      DUM = DELT * (Q(1) - Q(JMAX(L)))
      call INTGRL(JMAXL, DX, DELZBJ, AREA)
      ADJUST = (DUM - AREA) / (XB(JMAX(L)) - XB(1))
      do 130 J= 1,JMAX(L)
         DELZB(J,L) = ADJUST+DELZBJ(J)

         if (ISEDAV == 1) then
            if ((DELZB(J,L)+HP(J,L)) < 0.D0) DELZB(J,L)=-HP(J,L)
         endif
130   end do
   endif
   ! End of ICALL = 1

   ! If ICALL == 2 from main program, the profile change due to alongshore gradient of longshore sediment transport is included if IQYDY == 1 and computed when IEND == 1 and L == ILINE
   if (ICALL == 2) then
      do 200 LL= 1, ILINE
         JMAXL=JMAX(LL)
         do 210 J= 1,JMAXL
            R(J)=VY(J,LL)
            DELZBRW(J)=DABS(ZB(J,LL)-DZX(J,LL))
210      end do
         call SMOOTH(JMAXL,R,CB)
         call SMOOTH(JMAXL,DELZBRW,DELZBJ)
         call INTGRL(JMAXL, DX, CB, AREA)
         AVY(LL)=AREA
         call INTGRL(JMAXL, DX, DELZBJ, AREA)
         ADZX(LL)=AREA
         do 211 J= 1,JMAXL
            VY(J,LL)=CB(J)
            DZX(J,LL)=DELZBJ(J)
211      end do
200   end do
      call SMOOTH(ILINE,AVY,V)
      ILINE1=ILINE-1
      do 220 LL= 1, ILINE1
         VDY(LL)=(V(LL+1)-V(LL))/DYLINE(LL)
220   end do
      AVY(1)=VDY(1)
      AVY(ILINE)=VDY(ILINE1)
      do 230 LL= 2, ILINE1
         ! Use upstream finite difference method
         DUM=V(LL)*DYLINE(LL)
         if (DUM >= 0.D0) then
            AVY(LL)=VDY(LL-1)
         else
            AVY(LL)=VDY(LL)
         endif
230   end do
      do 240 LL= 1, ILINE
         if (ADZX(LL) < 1.D-6) then
            AVY(LL)=AVY(LL)/1.D-6/SPORO1
         else
            AVY(LL)=AVY(LL)/ADZX(LL)/SPORO1
         endif
240   end do
      do 250 LL= 1, ILINE
         JMAXL=JMAX(LL)
!     DUM=XB(JMAXL)*SPORO1
         do 260 J= 1, JMAXL
            DELZBRW(J)=-DZX(J,LL)*AVY(LL)
!     DELZBRW(J)=-AVY(LL)/DUM
260      end do
         call SMOOTH(JMAXL,DELZBRW,DELZBJ)
         do 270 J= 1, JMAXL
            ZB(J,LL)=DELZBJ(J)+ZB(J,LL)
270      end do
250   end do
   endif
   ! End of ICALL = 2
   
   deallocate(DZBDT, CB, R, DELZBRW, DELZBJ, V, AVY, VDY, ADZX)

   return
end subroutine CHANGE

