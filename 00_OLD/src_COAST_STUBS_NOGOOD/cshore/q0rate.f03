!===============================================================================================================================
!
! This subroutine computes overtopping, overflow and seepage rates
!
!===============================================================================================================================
subroutine QORATE(ITIME, L, ITEQO, ICONV, ICALL)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision(A-H, O-Z)
   
   integer, intent(in) :: ITIME, L, ITEQO, ICALL
   integer, intent(inout) :: ICONV
   
   double precision, allocatable, dimension(:) :: WSET, ZRW, SIGPT
   
   interface
      subroutine WETDRY(ITIME, L, ITEQO, SPRATE)
         integer, intent(in) :: ITIME, L, ITEQO
         double precision, intent(out) :: SPRATE
      end subroutine WETDRY
   end interface
   
   allocate(WSET(JMAXMAX), ZRW(JMAXMAX), SIGPT(JMAXMAX))

   ! Find overtopping, overflow and seepage rates during ICALL=0
   ! Start of ICALL = 0
   if (ICALL == 0) then
      ! Predict combined wave overtopping and overflow rate QOTF by calling Subroutine16 WETDRY for wet and dry zone
      call WETDRY(ITIME, L, ITEQO, SPRATE)

      ! Compute new combined rate QONEW and check convergency (one percent). Allowable error is increased for QONEW less than 1.D-2(m*m/s)
      QONEW = QOTF
      if (IPERM == 1) QONEW = QONEW + SPRATE
      if (QONEW < 1.D-5) then
         ICONV = 0
         QO(L) = QONEW
         goto 99
      endif
      
      DUM = DABS(QONEW-QO(L))/QONEW
      AER= 1.D-4/QONEW
      if (AER < 1.D-2) AER= 1.D-2
      if (DUM < AER) then
         ICONV = 0
         QO(L)=QONEW
         goto 99
      else
         ICONV = 1
         ! To avoid numerical oscillation and accelerate convergance use FRACTN of previous value and (1.0-FRACTN) of new value
         FRACTN = 0.5D0 + 0.1D0 * DBLE(ITEQO)
         if (FRACTN > 0.9D0) FRACTN=0.9D0
         if (ITEQO == 10) FRACTN=0.5D0
         if (ITEQO == 20) FRACTN=0.D0
         QO(L) = FRACTN*QO(L) + (1.D0-FRACTN)*QONEW
      endif
99    continue

      if (IPOND == 1) then
         if (NOPOND == 1) QM=QO(L)
         if (JCREST(L) == JXW) QM=QO(L)
         if (ZW >= ZB(JMAX(L),L)) QM=QO(L)
      endif
   endif
   ! End of ICALL = 0

   ! Start of ICALL = 1
   ! Output computed values in file 20 'ODOC' if ICALL= 1 in Subroutine8 OUTPUT
   if (ICALL == 1) then
      ! Mean (ERMEAN) above datum Z=0 and standard deviation (SIGRUN) of runup
      ! WSET(J)=wave setup above datum Z=0.0 during wet+dry duration
      ! SIGPT(J)=standard deviation during wet+dry duration
      ! ZRW(J)=runup wire elevation (RWH above bottom) at node J above Z=0.0

      do 170 J= 1,JCREST(L)
         if (J <= JDRY) then
            WSET(J) = H(J)*PWET(J)+ZB(J,L)
         else
            WSET(J)=ZB(J,L)
         endif
         SIGPT(J)=SIGMA(J)*PWET(J)
         ZRW(J)=ZB(J,L)+RWH
170   end do

      ! K= 1, 2 and 3 correspond to intersections of ZRW with WSET, (WSET-SIGPT) and (WSET+SIGPT), respectively
      do 100 K= 1, 3
         J=JDRY
         if (JDRY > JCREST(L)) J=JCREST(L)
         DUM1=ZRW(J)-WSET(J)
         if (K == 2) DUM1=DUM1+SIGPT(J)
         if (K == 3) DUM1=DUM1-SIGPT(J)
         if (DUM1 < 0.D0) then
            if (K == 1) then
               ETARUN=WSET(J)
               goto 100
            endif
            if (K == 2) then
               Z1RUN=WSET(J)-SIGPT(J)
               ! start bdj 2012-05-17
               X1RUN=XB(J) ! BDJ
               ! end bdj 2012-05-17
               goto 100
            endif
            if (K == 3) then
               Z2RUN=WSET(J)+SIGPT(J)
               ! start bdj 2012-05-17
               X2RUN=XB(J) ! BDJ
               if (X2RUN <= X1RUN) X2RUN=X1RUN+DX ! BDJ
               ! end bdj 2012-05-17
               goto 100
            endif
         endif
105      J=J-1
         DUM2=ZRW(J)-WSET(J)
         if (K == 2) DUM2=DUM2+SIGPT(J)
         if (K == 3) DUM2=DUM2-SIGPT(J)
         if (DUM2 > 0.D0) then
            DUM1=DUM2
            goto 105
         else
            DUM3=DUM1-DUM2
            DUMJ1=-DUM2/DUM3
            DUMJ=DUM1/DUM3
            DUMETA=DUMJ*WSET(J)+DUMJ1*WSET(J+1)
            if (K == 1) ETARUN=DUMETA
            if (K == 2) then
               Z1RUN=DUMETA-DUMJ*SIGPT(J)-DUMJ1*SIGPT(J+1)
               X1RUN=DUMJ*XB(J)+DUMJ1*XB(J+1)
            endif
            if (K == 3) then
               Z2RUN=DUMETA+DUMJ*SIGPT(J)+DUMJ1*SIGPT(J+1)
               X2RUN=DUMJ*XB(J)+DUMJ1*XB(J+1)
               
               ! bdj 2013-03-04
               if ((WSET(J+1) - WSET(J)) / SIGPT(J) > (10.0D0 * DX)) then
                  DUMETA=WSET(J)  ! BDJ
                  Z2RUN=DUMETA+SIGPT(J) ! BDJ
                  X2RUN=XB(J)  ! BDJ
                  if (x2run-x1run <= 01D0*DX) then
                     Z2RUN=Z1RUN  + .01D0*DX*BSLOPE(J,L)
                     X2RUN=X1RUN + .01D0*DX
                  endif
               endif
               ! end bdj 2013-03-04
            endif
         endif
100   end do
      SIGRUN=(Z2RUN-Z1RUN)/2.D0
      ERMEAN=(Z1RUN+ETARUN+Z2RUN)/3.D0
      SLPRUN=(Z2RUN-Z1RUN)/(X2RUN-X1RUN)
      
      ! bdj 2014-04-29 added catch for negative slopes
      SIGRUN=max(0.D0,SIGRUN)
      ERMEAN=max(z1run,ERMEAN)
      SLPRUN=max(0.D0,SLPRUN)
      ! end bdj 2014-04-29 added catch for negative slopes
      
      ! bdj 2015-05-06  added catch for cases where waves are very small
      if (J < NINT(DBLE(JSWL(L)) / 2.0D0)) then
         SIGRUN=0.D0
         ERMEAN=SWLBC(ITIME)
      endif
      ! end bdj 2015-05-06  added catch for cases where waves are very small

      ! R13=significant runup height above Z=0.0
      ! R2P=two percent runup height above Z=0.0
      ! RKAPPA=Kappa for runup probability distribution
      if (IPERM == 1) then
         R13=ERMEAN+(2.D0+SLPRUN)*SIGRUN
         RSC=(RCREST(L)-ERMEAN)/(R13-ERMEAN)
         RKAPPA= 2.0D0+0.5D0/RSC**3.D0
      else
         DUM= 4.D0*SLPRUN
         ! bdj 2013-03-04
         ! if (DUM.GT.2.D0) DUM= 2.D0
         if (DUM > 1.D0) DUM= 1.D0
         ! end bdj 2013-03-04
         
         R13=(ERMEAN-SWLBC(ITIME)+2.D0*SIGRUN)*(1.D0+DUM)+SWLBC(ITIME)
         RKAPPA= 2.0D0
      endif
      
      if (RCREST(L) > ERMEAN) then
         R2P=ERMEAN+(R13-ERMEAN)*1.4D0**(2.D0/RKAPPA)
         R1P=ERMEAN+(R13-ERMEAN)*1.52D0**(2.D0/RKAPPA)
      else
         RKAPPA= 1000.D0
         R2P=R13
      endif

      ! Output swash hydrodynamics computed in Subroutine16 WETDRY..........
      if (JDRY >= JCREST(L)) then
         POTF=(DTANH(5.D0*PWET(JCREST(L))))**0.8D0
      else
         POTF=0.D0
      endif
      
      ! Depth H, velocity U and discharge Q corresponding to exceedance probability EWD specified in Subroutine04 PARAM
      if (JWD <= JDRY) then
         do 300 J=JWD, JDRY
            DUM = PWET(J)/EWD
            if (DUM < 1.1D0) DUM= 1.1D0
            HEWD(J)=(H(J)/PWET(J))*DLOG(DUM)
            DUM=USWD(J)
            if (DUM < 0.D0) DUM=0.D0
            UEWD(J) = AWD*DSQRT(GRAV*HEWD(J))+DUM
            QEWD(J) = HEWD(J)*UEWD(J)
300      end do
      endif      
      ! Where computed HEWD(J), UEWD(J) and QEWD(J) are stored in Subroutine OUTPUT
      
#if defined FILEINOUT || ARGINBOTHOUT
      ! DFM Got FP exception with undefined SPRATE when ICALL == 1
      SPRATE = 0.0D0
      write (20, 920) SWLBC(ITIME), L, RCREST(L), JSWL(L), JWD, H1, JDRY, POTF, (QO(L)-SPRATE), SPRATE, QO(L), ITEQO
      
920   format('COMBINED WAVE OVERTOPPING AND OVERFLOW'/ &
         'Still water level above Z == 0 (m)           SWL =',F13.6/ &
         'Cross-shore line number                        L =',I3/ &
         'Structure or dune crest elevation (m)     RCREST =',F13.6/ &
         'Node number at SWL                          JSWL =',I6/ &
         'Wet and dry transition node                  JWD =',I6/ &
         'Mean water depth H1(m) at node JWD            H1 =',F13.6/ &
         'End node of wet and dry zone                JDRY =',I6/ &
         'Wave overtopping probability at JCREST      POTF =',F13.6/ &
         'Comb. overtopping and overflow rate(m*m/s)  QOTF =',F13.9/ &
         'Seepage rate(m*m/s) at JCREST                 QP =',F13.9/ &
         'Total rate (QOTF+QP)(m*m/s)                      =',F13.9/ &
         'QO iteration number                        ITEQO =',I3/)

      ! Output empirical runup
      write (20, 900) L, SLPRUN, X1RUN, X2RUN, Z1RUN, Z2RUN, ERMEAN, SIGRUN, R13, R2P, R1P

900   format('EMPIRICAL WAVE RUNUP'/ &
         'Cross-shore line number                        L =',I3/ &
         'Swash zone bottom slope                   SLPRUN =',F13.6/ &
         '   computed from                           X1RUN =',F13.6/ &
         '   to                                      X2RUN =',F13.6/ &
         '   with                                    Z1RUN =',F13.6/ &
         '   and                                     Z2RUN =',F13.6/ &
         'Mean runup elevation above Z=0 (m)        ERMEAN =',F13.6/ &
         'Runup standard deviation (m)              SIGRUN =',F13.6/ &
         'Significant runup height above Z=0 (m)       R13 =',F13.6/ &
         '2 percent runup height above Z=0 (m)         R2P =',F13.6/ &
         '1 Percent runup height above z=0 (m)         R1P =',F13.6/)

      if (IWTRAN == 1) then
         if (JDRY == JSL1) then
            write (20, 940)L, JSL, XB(JSL), WSETUP(JSL), SIGMA(JSL), XB(JMAX(L)), WSETUP(JMAX(L)), SIGMA(JMAX(L)), (SIGMA(JMAX(L)) / SIGMA(1))
         else
            write (20, 941) JDRY, JSL1
         endif
      endif
      
940   format('WAVE TRANSMISSION DUE TO OVERTOPPING'/ &
         'Cross-shore line number                        L =', I3/ &
         'Starting node for wave transmission          JSL =', I6/ &
         'X-coordinate (m) at node JSL                  XB =', F13.6/ &
         'Wave setup (m) at node JSL                WSETUP =', F13.6/ &
         'Standard deviation (m) at node JSL         SIGMA =', F13.6/ &
         'X-coordinate (m) at landward end node JMAX       =', F13.6/ &
         'Wave setup (m) at landward end node JMAX         =', F13.6/ &
         'Standard dev. (m) at landward end node JMAX      =', F13.6/ &
         'Wave transmission coefficient at JMAX            =', F13.6/)
         
941   format(/'IWTRAN= 1 BUT NO WAVE TRANSMISSION'/'JDRY=', I6, '  is less than (JSL-1) =', I6/)

      if (IPOND == 1 .and. NOPOND == 0) then
         write (20, 960) L, JCREST(L), JXW, JX2, ZW, QD, QM
      endif

960   format('PONDED WATER IN RUNNEL'/ &
         'Cross-shore line number                        L =', I3/ &
         'Ridge crest node                          JCREST =', I6/ &
         'Ponded water nodes from                      JXW =', I6/ &
         '                     to                      JX2 =', I6/ &
         'Ponded water level (m)                        ZW =', F13.6/ &
         'Wave-induced volume flux (m*m/s)              QD =', F13.6/ &
         'Wave overtopping rate (m*m/s) at JMAX         QM =', F13.6/)

#endif

      !.................................End of ICALL= 1.............................
   endif

   deallocate(WSET, ZRW, SIGPT)
   
   return
end subroutine QORATE
