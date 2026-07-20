!===============================================================================================================================
!
! This subroutine computes swash hydrodynamics in the wet/dry zones (possibly with multiple bottom peaks), and combined wave overtopping and overflow rate QOTF
!
!===============================================================================================================================
subroutine WETDRY(ITIME, L, ITEQO, SPRATE)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   integer, intent(in) :: ITIME, L, ITEQO
   double precision, intent(out) :: SPRATE
   
   integer :: ITEH, IUPSLP, J, JC, JDUM, JEND, JP1, JSEP, LSTART

   double precision, allocatable, dimension(:) :: G, DG, ETA, ETAP
   
   interface
      function GBWD(R)
         double precision, intent(in) :: R
      end function GBWD
   end interface
   
   allocate(G(JMAXMAX), DG(JMAXMAX), ETA(JMAXMAX), ETAP(JMAXMAX))
   
   ! Compute swash variables for node J=JWD to JDRY
   ! JWD = wet and dry transition node
   ! JR = landward limit of wet computation in Main Program
   ! JDRY = landward limit of wet and dry zone
   ! PWET(J) = wet probability at node J where PWET(JWD)= 1.0
   ! HWD(J) = mean water depth where HWD(JWD)=H1
   ! USWD(J) = steady swash velocity component
   ! UMEAWD(J) = mean velocity in wet and dry zone
   ! QP(J) = water flux in permeable layer
   ! If INFILT= 1, QP(J) =infiltration rate between dune crest and landward node J where QP(J) = 0.0 assumed seaward of dune crest
   ! UPMWD(J) = mean discharge velocity in permeable layer
   
   if (ITEQO <= 2) then
      JWD = JSWL(L)
      if (JWD > JR) JWD=JR
         ! If IPROFL= 1 and JWD=JR, no wet and dry zone to avoid possible overwash under small waves
         ! bdj 2013-03-04
         ! if (IPROFL.EQ.1.and.JSWL(L).GT.JR) then
         if (JSWL(L) > JR) then
         ! bdj end 2013-03-04
         if (IPOND == 0) then
            JDRY=JR
            goto 110
         endif
      endif
      H1 = H(JWD)
   endif

   HWD(JWD) = H1
   BGH3=BWD*GRAV*H1*H1*H1
   
   ! BDJ added 2012-10-10
   SSP_50 = 1.D0     ! The value of SSP such that the CORRECT is .5*(1+GAMMA/SQR8)
   A = 1.D0          ! Dictates the steepness of blending curve near SSP_50
   CORRECT = GAMMA/SQR8
   CORRECT = 0.5D0*(1.D0+CORRECT)+ 0.5D0*(1.D0-CORRECT)*tanh(a*(SSP-SSP_50))
   ! CORRECT = 1.D0  ! Comment this to get correction
   SIGWD(JWD) = CORRECT*H1
   ! end BDJ added 2012-10-10

   PWET(JWD) = 1.D0
   if (IPERM == 0) then
      QWX=QO(L)
      if (INFILT == 1) QP(JWD)=0.D0
   else
      QWX=QO(L)-QP(JWD)
      UPMWD(JWD)=UPMEAN(JWD)
      ETA(JWD) = H1 + ZB(JWD,L)
      ETAP(JWD)=ZB(JWD,L)
      PMG1=AWD/DSQRT(GRAV)
      PMGH1=PMG1/H1
   endif
   QS = QWX-AQWD*H1*DSQRT(GRAV*H1)
   if (QS > 0.D0) QS=0.D0
   USWD(JWD) = QS/H1
   UMEAWD(JWD) = AUWD*DSQRT(GRAV*H1) + USWD(JWD)
   DUM=AGWD*GRAV*H1 - (UMEAWD(JWD)-USWD(JWD))**2.D0
   USTDWD(JWD)=DSQRT(DUM)
   A = QWX*QWX/BGH3
   A1=A

   ! Empirical formula for wet probability parameter n=WDN
   WDN= 1.01D0+0.98D0*DTANH(QO(L)*QO(L)/BGH3)**0.3D0
   W1=WDN-1.D0
   BNWD=BWD*(1.D0+A1)*(2.D0-WDN)/(WDN-1.D0)

   ! LANDWARD MARCHING COMPUTATION
   ! If IWTRAN= 1, the landward wet zone starts from node J=JSL
   JEND=JMAX(L)-1
   if (IWTRAN == 1) then
      JEND=JSL1-1
   endif
   
   ! LSTART= 1 indicates beginning of upslope computation
   LSTART = 1
   if (JWD > JEND) then
      JDRY = JMAX(L)
      goto 110
   endif
   do 100 J = JWD, JEND
      JP1 = J + 1

      ! BOTTOM ELEVATION INCREASING LANDWARD 
      ! On the seaward upslope and crest (J < JCREST but J < JMAX if IPOND == 1) use an empirical formula for wet probability PWET and compute mean depth and return flow velocity USWD for IUPSLP == 1
      IUPSLP = 0
      JDUM = JCREST(L)
      if (IPOND == 1) JDUM = JMAX(L)
      if (J < JDUM .and. ZB(JP1,L) >= ZB(J,L)) IUPSLP= 1
      ! For J=JWD, IUPSLP = 1 for any slope
      if (J == JWD) IUPSLP= 1
      ! If IPOND == 1 and NOPOND == 0, ponded water zone is treated like downslope zone with IUPSLP == 0
      if (IPOND == 1 .and. NOPOND == 0) then
         if (JP1 >= JXW .and. JP1 <= JX2) IUPSLP = 0
      endif
      if (IUPSLP == 1) then
         if (LSTART == 1) then
            H2 = HWD(J)
            BN12 = BNWD * (H1 / H2)**W1
            D = BN12 - ZB(J,L) / H1
            AH = AGWD / H1
            G(J) = 0.D0
            DUM = QWX - QS
            
            if (DABS(DUM) < 1.D-3) then
               R = 0.D0
            else
               R = CWD * QS / DUM
            endif
            
            DG(J) = AH * DX * FB2(J,L) * GBWD(R)
            
            !  BDJ 2012-10-26 added to kill friction in wetdry
            DG(J) = 0.D0
            !  end BDJ 2012-10-26

            LSTART=0
         endif
         CX = D + ZB(JP1,L)/H1
         if (IPERM == 0) then
            WPGH=0.D0
         else
            if (HP(JP1,L) <= 0.D0) then
               WPGH=0.D0
            else
               DUM=0.5D0*(HP(J,L)+HP(JP1,L))/SDP
               ! if (DUM.GT.10.D0) DUM= 10.D0
               PMGH=PMGH1*DUM**0.3D0
               WPGH=PMGH*WPM*PWET(J)*DX/DSQRT(HWD(J))
            endif
         endif
         DGJP1 = DG(J)+WPGH
         do 200 ITEH= 1, 20
            G(JP1) = G(J) + DGJP1
            C = CX + G(JP1)
            if (C <= 0.D0) then
               JDRY = J
               goto 110
            else
               Y = (C/BN12)**(1.D0/W1)
            endif
            HWD(JP1)=H2/Y
            Y=H1/HWD(JP1)
            DUM = (1.D0 + A1)*Y**WDN - A*Y*Y*Y
            if (DUM < 1.D0) then
               PWET(JP1) = PWET(J)
            else
               PWET(JP1) = 1.D0/DUM
               if (PWET(JP1) > PWET(J)) PWET(JP1)=PWET(J)
            endif
            QWAVE=AQWD*HWD(JP1)*DSQRT(GRAV*HWD(JP1)/PWET(JP1))

            ! Compute QP and UPMWD in permeable layer if IPERM= 1
            ! ETAP(JP1)=mean water level above datum inside permeable layer where ETA(JP1) and ETAP(JP1) are mean water levels above datum
            if (IPERM == 1) then
               if (HP(JP1,L) <= 0.D0) then
                  UPMWD(JP1)=0.D0
                  QP(JP1)=0.D0
                  WPGH=0.D0
               else
                  ETA(JP1)=HWD(JP1)+ZB(JP1,L)
                  DUM=ZP(JP1,L)
                  if (DUM < SWLBC(ITIME) .and. ZP(JP1,L) >= ZP(J,L)) DUM=SWLBC(ITIME)
                  ETAP(JP1)=ZB(JP1,L)*PWET(JP1)+DUM*(1.D0-PWET(JP1))
                  if (ETAP(JP1) < ZP(JP1,L)) ETAP(JP1)=ZP(JP1,L)
                  C=(ETA(JP1)-ETA(J))/DX
                  DUM=DSQRT(ALSTA2+BE4*DABS(C))
                  UPMWD(JP1)=(DUM-ALSTA)/BE2
                  if (C > 0.D0) UPMWD(JP1)=-UPMWD(JP1)
                  QP(JP1)=UPMWD(JP1)*(ETAP(JP1)-ZP(JP1,L))*PWET(JP1)
                  DUM=QO(L)-QWAVE
                  if (QP(JP1) < DUM .and. DABS(QP(JP1)) > 1.D-5) then
                     UPMWD(JP1)=UPMWD(JP1)*DUM/QP(JP1)
                     QP(JP1)=DUM
                  endif
                  DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                  WPGH=PMGH*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
               endif
               QWX=QO(L)-QP(JP1)
               A=QWX*QWX/BGH3
            endif
            
            if (INFILT == 1) QP(JP1)=0.D0

            QS = QWX-QWAVE
            ! QS = return flux must be zero or negative for J < JCREST
            if (QS > 0.D0) QS=0.D0
            USWD(JP1) = QS/HWD(JP1)
            DUM=QWX-QS
            if (DABS(DUM) < 1.D-3) then
               R=0.D0
            else
               R = CWD*QS/DUM
            endif
            DG(JP1) = AH*DX*FB2(JP1,L)*GBWD(R)
            
            ! BDJ 2012-10-26 added to kill friction in wetdry
            DG(JP1) = 0.D0
            ! end BDJ 2012-10-26
            DUM=0.5D0*(DG(J)+DG(JP1))+WPGH
            if (DABS(DUM-DGJP1) > 1.D-5) then
               DGJP1 = DUM
               goto 200
            else
               G(JP1)=G(J)+DUM
               
               ! LSTART= 2 indicates that bottom elevation is peaked at node JC
               if (JP1 < JMAX(L)) then
                  if (ZB(J+2,L) < ZB(JP1,L)) then
                     JC=JP1
                     HC=HWD(JC)
                     PC = PWET(JC)
                     QWC=QWX
                     LSTART= 2
                  else
                     if (J == JWD) LSTART= 1
                  endif
               endif
               goto 220
            endif
200      end do
      else

         ! BOTTOM ELEVATION DECREASING LANDWARD OR J > JCREST --------------
         ! On the landward slope (J>JCREST) or downslope zone for J<JCREST or ponded water zone, PWET=constant on impermeable bottom and compute HWD and USWD for IUPSLP=0
         if (LSTART == 2) then
            PCI= 1.D0/PC
            QWC2=QWC*QWC
            BG=BWD*GRAV
            CPC = 0.5D0*PC/BWD/HC
            AB = 0.25D0*PC*QWC2/(BG*HC*HC*HC)
            G(J) = 0.D0
            QS=USWD(J)*HWD(J)
            DUM=QWC-QS
            if (DABS(DUM) < 1.D-3) then
               R=0.D0
            else
               R=CWD*QS/DUM
            endif
            DG(J)=AGWD*FB2(J,L)*DX*GBWD(R)
            LSTART=0
         endif
         
         DZB=ZB(JC,L)-ZB(JP1,L)
         if (IPOND == 1 .and. NOPOND == 0) then
            if (JP1 >= JXW .and. JP1 <= JX2) then
               if (JX2 < JMAX(L)) DZB=ZB(JC,L)-ZW
               QWX=QO(L)-(QO(L)-QM)*(XB(JP1)-XB(JXW))/(XB(JX2)-XB(JXW))
               A=QWX*QWX/BGH3
            endif
         endif
         
         if (IPERM == 0) then
            WPGH=0.D0
            if (INFILT == 1) then
               WPGH=PMG1*WPM*PWET(J)*DX/DSQRT(HWD(J))
            endif
         else
            if (HP(JP1,L) <= 0.D0) then
               WPGH=0.D0
            else
               DUM=0.5D0*(HP(J,L)+HP(JP1,L))/SDP
               ! if (DUM.GT.10.D0) DUM= 10.D0
               PMG=PMG1*DUM**0.3D0
               WPGH=PMG*WPM*PWET(J)*DX/DSQRT(HWD(J))
            endif
         endif
         
         DGJP1 = DG(J)+WPGH
         
         do 210 ITEH= 1, 20
            G(JP1)= G(J)+DGJP1+WPGH
            C=CPC*(DZB-G(JP1))
            if (C < 0.D0) C=0.D0
            Y= HWD(J)/HC
205         DUM= 1.D0/Y/Y
            F=Y-1.D0+AB*(DUM-1.D0)-C
            DF= 1.D0-2.D0*AB*DUM/Y
            if (DABS(DF) < 1.D-6) then
               JDRY=J
               goto 110
            endif
            YNEW=Y-F/DF
            if (DABS(YNEW-Y) > 1.D-6) then
               Y=YNEW
               goto 205
            endif
            if (YNEW <= 0.D0) then
               JDRY=J
               goto 110
            endif
            HWD(JP1)=YNEW*HC
            if (HWD(JP1) > HWD(J)) HWD(JP1)=HWD(J)
            if (IPERM == 0 .and. INFILT == 0) then
               PWET(JP1)=PC
            else
               DUM=PCI+(QWC2-QWX*QWX)/BG/HWD(JP1)**3.D0
               if (DUM < PCI) DUM=PCI
               PWET(JP1)= 1.D0/DUM
               if (PWET(JP1) > PWET(J)) PWET(JP1)=PWET(J)
            endif
            QWAVE=AQWD*HWD(JP1)*DSQRT(GRAV*HWD(JP1)/PWET(JP1))

            ! Compute QP and UPMWD in permeable layer if IPERM= 1 as above
            if (IPERM == 1) then
               if (HP(JP1,L) <= 0.D0) then
                  UPMWD(JP1)=0.D0
                  QP(JP1)=0.D0
                  WPGH=0.D0
               else
                  ETA(JP1)=HWD(JP1)+ZB(JP1,L)
                  DUM=ZP(JP1,L)
                  if (DUM < SWLBC(ITIME) .and. ZP(JP1,L) >= ZP(J,L)) DUM=SWLBC(ITIME)
                  ETAP(JP1)=ZB(JP1,L)*PWET(JP1)+DUM*(1.D0-PWET(JP1))
                  if (ETAP(JP1) < ZP(JP1,L)) ETAP(JP1)=ZP(JP1,L)
                  C=(ETA(JP1)-ETA(J))/DX
                  DUM=DSQRT(ALSTA2+BE4*DABS(C))
                  UPMWD(JP1)=(DUM-ALSTA)/BE2
                  if (C > 0.D0) UPMWD(JP1)=-UPMWD(JP1)
                  QP(JP1)=UPMWD(JP1)*(ETAP(JP1)-ZP(JP1,L))*PWET(JP1)
                  
                  ! DUM=QO(L)-QWAVE
                  ! if (J.GE.JCREST(L).and.QP(JP1).GT.DUM) then
                  ! UPMWD(JP1)=UPMWD(JP1)*DUM/QP(JP1)
                  ! QP(JP1)=DUM
                  ! endif
                  
                  DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                  WPGH=PMG*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
               endif
               QWX=QO(L)-QP(JP1)
               A=QWX*QWX/BGH3
            endif

            ! Compute QP(J)=infiltration rate landward of dune crest if INFILT= 1
            if (INFILT == 1) then
               if (JP1 > JCREST(L)) then
                  QP(JP1)=QP(J)+0.5D0*DX*WPM*(PWET(J)+PWET(JP1))
                  DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                  WPGH=PMG1*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
                  QWX=QO(L)-QP(JP1)
                  A=QWX*QWX/BGH3
               else
                  QP(JP1)=0.D0
                  WPGH=0.D0
               endif
            endif
            QS = QWX - QWAVE
            
            ! QS = steady flux on landward slope (J>JCREST) must be zero or positive
            if (IPOND == 0) then
               if (J >= JCREST(L) .and. QS < 0.D0) QS=0.D0
            endif
            
            if (HWD(JP1) < 1.D-3 .and. QS > 1.D-3) QS= 1.D-3
            
            USWD(JP1) = QS/HWD(JP1)
            DUM=QWX-QS
            if (DABS(DUM) < 1.D-3) then
               R=0.D0
            else
               R=CWD*QS/DUM
            endif
            
            DG(JP1) = AGWD*FB2(JP1,L)*DX*GBWD(R)
            DUM=0.5D0*(DG(J)+DG(JP1))+WPGH
            if (DABS(DUM-DGJP1) > 1.D-5) then
               DGJP1 = DUM
               goto 210
            else
               G(JP1)=G(J)+DUM
               
               ! LSTART= 1 indicates beginning of upslope computation
               if (IPOND == 0 .OR. NOPOND == 1) then
                  if (JP1 < JCREST(L)) then
                     if (ZB(J+2,L) >= ZB(JP1,L)) LSTART= 1
                  endif
               else
                  if (JP1 == JX2) then
                     LSTART= 1
                     QWX=QM
                     A=QWX*QWX/BGH3
                  endif
               endif
               goto 220
            endif
210      end do
      endif
      ! END OF BOTTOM ELEVATION INCREASING OR DECREASING

      ! Compute mean velocity UMEAWD and standard deviations SIGWD and USTDWD in wet and dry zone
220   UMEAWD(JP1)=AUWD*DSQRT(GRAV*PWET(JP1)*HWD(JP1)) +PWET(JP1)*USWD(JP1)
      ! BDJ added 2012-10-10
      ! SIGWD(JP1)=HWD(JP1)*DSQRT(2.D0/PWET(JP1)-2.D0+PWET(JP1))
      SIGWD(JP1)=CORRECT*HWD(JP1)*DSQRT(2.D0/PWET(JP1)-2.D0+PWET(JP1))
      ! BDJ added 2012-10-10
      DUM = UMEAWD(JP1) - USWD(JP1)
      DUM1 = PWET(JP1)*DUM**2.D0-2.D0*DUM*(UMEAWD(JP1)- PWET(JP1)*USWD(JP1))
      DUM = AGWD*GRAV*HWD(JP1)+DUM1
      if (DUM > 0.D0) then
         USTDWD(JP1) = DSQRT(DUM)
      else
         USTDWD(JP1) = 0.D0
      endif
      if (IANGLE == 1) then
         STHETA(JP1)=STHETA(JWD)
         VMEAWD(JP1)=AUWD*DSQRT(GRAV*PWET(JP1)*HWD(JP1))*STHETA(JP1)
         DUM= 1.D0-0.25D0*PI*PWET(JP1)*(2.D0-PWET(JP1))
         VSTDWD(JP1)=AWD*DSQRT(GRAV*HWD(JP1)*DUM)*DABS(STHETA(JP1))
      endif
      
      ! Mean water depth limited by HWDMIN specified in Subroutine2 INPUT
      ! Horizontal distance of wet and dry zone is limited because of assumed alongshore uniformity
      ! DUM=(XB(JP1)-XB(JWD))/H1
      ! if (HWD(JP1).LT.HWDMIN.OR.DUM.GT.1000D0) then
      if (HWD(JP1) < HWDMIN) then
         JDRY = JP1
         ! if (DUM.GT.1000D0.and.JP1.GT.JCREST(L)) JMAX(L)=JP1
         goto 110
      endif

      if (J == JEND) JDRY=JP1
100 end do
   ! END OF LANDWARD MARCHING

110 continue

   ! QOTF = Combined overtopping and overflow rate QOTF
   ! SPRATE = seepage rate through permeable layer predicted by modified formula of Kobayashi and de los Santos(2007) for no overtopping, where USWD = 0.0 (unidirectional flow) at JCREST is assumed
   QOTF=0.D0
   SPRATE=0.D0
   if (JDRY >= JCREST(L)) then
      J=JCREST(L)
      if (JWD == JMAX(L)) J=JMAX(L)
      QOTF = AQWD*HWD(J)*DSQRT(GRAV*HWD(J)/PWET(J))
      if (IPERM == 1) then
         if (JDRY == JMAX(L) .OR. JDRY == JSL1) then
            SPRATE=QP(JCREST(L))
            if (SPRATE < 0.D0) SPRATE=0.D0
         else
            QOTF=0.D0
         endif
      endif
   endif
   
   if (IPOND == 1 .and. NOPOND == 0) then
      QD=QOTF
      if (JDRY == JMAX(L)) then
         if (ZW < ZB(JMAX(L),L)) then
            QM=AQWD*HWD(JMAX(L))*DSQRT(GRAV*HWD(JMAX(L))/PWET(JMAX(L)))
            if (QM > QOTF) QM=QOTF
         else
            QM=QOTF
         endif
      else
         QM=0.D0
      endif
      if (JCREST(L) == JXW) QOTF=QM
   endif
   if (IPERM == 1 .and. QOTF == 0.D0) then
      ! Find node JDUM for highest and most landward point of ZP(J,L)
      JSEP=JR
      JDUM=JSEP
      DUM=ZP(JSEP,L)
      do 300 J=(JSEP+1),JMAX(L)
         if (ZP(J,L) >= DUM) then
            DUM=ZP(J,L)
            JDUM=J
         endif
300   end do
      DETA=ZB(JSEP,L)-ZP(JDUM,L)
      if (DETA > 0.D0) then
         DUM=XB(JDUM)-XB(JSEP)
         if (DUM <= 0.D0) DUM=DX
         SPRATE=0.2D0*DETA**1.5D0/DSQRT(BESTA1*DUM)
      endif
   endif

   deallocate(G, DG, ETA, ETAP)

   return
end subroutine WETDRY
