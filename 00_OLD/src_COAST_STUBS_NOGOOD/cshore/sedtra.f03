!===============================================================================================================================
!
! This subroutine calculates cross-shore and longshore sediment transport
!
!===============================================================================================================================
subroutine SEDTRA(L)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   interface 
      subroutine SMOOTH(NUM, RAW, F)
         integer, intent(in) :: NUM
         double precision, intent(in), dimension(:) :: RAW
         double precision, intent(out), dimension(:) :: F
      end subroutine SMOOTH
      
      subroutine TRANWD(F1, JR, F2, JS, JE)
         integer, intent(in) :: JR, JS, JE
         double precision, intent(in), dimension(:) :: F2
         double precision, intent(out), dimension(:) :: F1
      end subroutine TRANWD
      
      subroutine EXTRAPO(J1, J2, F)
         integer, intent(in) :: J1, J2
         double precision, intent(out), dimension(:) ::  F
      end subroutine EXTRAPO
      
      function ERFCC(X)
         double precision, intent(in) :: X
      end function ERFCC
      
      subroutine PROBWD(PW, A, US, UC, P)
         double precision, intent(in) :: PW, A, US, UC
         double precision, intent(out) :: P
      end subroutine PROBWD
   end interface

   integer, intent(in) :: L         
   double precision :: ERFCC
   double precision, allocatable, dimension(:) :: QRAW, GSLRAW, ASLRAW, ASLOPE, RS, RB, PBWD, PSWD, VSWD, QSXWD, QBXWD, QRAWD, HDIP, QSYWD, QBYWD
   
   ! GSLMAX = Maximum absolute value of GSLOPE function
   data GSLMAX,BQCOEFF  /10.D0, 8.D0/

   allocate(QRAW(JMAXMAX), GSLRAW(JMAXMAX), ASLRAW(JMAXMAX), ASLOPE(JMAXMAX), RS(JMAXMAX), RB(JMAXMAX), PBWD(JMAXMAX), PSWD(JMAXMAX), VSWD(JMAXMAX), QSXWD(JMAXMAX), QBXWD(JMAXMAX), QRAWD(JMAXMAX), HDIP(JMAXMAX), QSYWD(JMAXMAX), QBYWD(JMAXMAX))

   ! Cross-Shore and Longshore Sediment Transport at Node J
   ! RB(J) = Sediment movement initiation parameter
   ! RS(J) = Sediment suspension initiation parameter
   ! PB(J) = bedload probability
   ! PS(J) = suspended load probability
   ! VS(J) = suspended sediment volume per unit area (m)
   ! GSLOPE(J) = bed slope correction for QBX(J)
   ! ASLOPE(J) = bed slope correction for suspended load parameter SLP
   ! QBX(J)= Cross-shore  bedload transport rate per unit width (m*m/s)
   ! QBY(J)= Longshore  bedload transport rate per unit width (m*m/s)
   ! BRF   = Bedload reduction factor for hard bottom (ISEDAV= 1)
   ! QSX(J)= Cross-shore suspended sediment transport rate (m*m/s)
   ! QSY(J)= Longshore suspended sediment transport rate (m*m/s)
   ! Q(J)  = total cross-shore sedimet transport rate including void
   ! (m*m/s) used for beach profile change computation

   if (TIME == 0.D0) then
      BSLOP1 = -TANPHI*(GSLMAX-1.D0)/GSLMAX
      BSLOP2 =  TANPHI*(GSLMAX+1.D0)/(GSLMAX+2.D0)
      IGMILD=0
      if (IPERM == 1) then
         DUM=0.5D0*TANPHI
         if (BSLOPE(JSWL(L),L) < DUM) IGMILD= 1
      endif
      
      ! Where input bedload parameter is increased in the surf zone in the following if IGMILD = 1 (based on two gravel tests MH and MB only)
   endif
   
   if (IVWALL(L) == 2) then
      JDUM=JMAX(L)
      BSLOPE(JDUM,L)=0.D0
      BSLOPE(JDUM-1,L)=(ZB(JDUM-1,L)-ZB(JDUM-2,L))/DX
   endif
   
   do 100 J= 1,JMAX(L)
      if (BSLOPE(J,L) < 0.D0) then
         if (BSLOPE(J,L) > BSLOP1) then
            GSLRAW(J) = TANPHI/(TANPHI + BSLOPE(J,L))
         else
            GSLRAW(J) = GSLMAX
         endif
      else
         if (BSLOPE(J,L) < BSLOP2) then
            GSLRAW(J) = (TANPHI - 2.D0*BSLOPE(J,L))/(TANPHI-BSLOPE(J,L))
         else
            GSLRAW(J) = -GSLMAX
         endif
         if (IGMILD == 1) then
            if (GSLRAW(J) < 0.D0) GSLRAW(J)=0.D0
         endif
      endif
      ASLRAW(J) = SLP
      if (TIME > 0.D0) then
         SUM=0.D0
         do 300 K=JSWL(L),JMAX(L)
            SUM=DELZB(K,L)
300      end do
         if (SUM < 0.D0) then
            if (BSLOPE(J,L) > 0.D0) ASLRAW(J)=SLP+(BSLOPE(J,L)/TANPHI) **0.3D0
         endif
      endif
100 end do

   ! Smoothing GSLOPE before Q is computed in Subroutine SMOOTH
   JMAXL=JMAX(L)
   call SMOOTH(JMAXL, GSLRAW, GSLOPE)
   call SMOOTH(JMAXL, ASLRAW, ASLOPE)

   ! Sediment transport rates are computed for normally incident waves in wet zone (IANGLE=0); wet and dry zone (IOVER= 1) for IANGLE=0 and 1; and obliquelly incident waves in wet zone (IANGLE= 1)
   if (IANGLE == 1) goto 888
   ! Normally Incident Waves in wet zone
   if (IANGLE == 0) then
      LMAX=LANCOM+1
      do 111 K= 1,LMAX
         if (K == 1) then
            JSTART= 1
            JEND=JR
         else
            JSTART=JSL
            JEND=JMAX(L)
         endif
         do 110 J =JSTART,JEND
            if (D50 < CSEDIA) then
               RB(J) = DSQRT(GSD50S/FB2(J,L))/USTD(J)
            else
               RB(J)=GSD50S/USTD(J)
            endif
            RS(J) = WF/USTD(J)/FB2(J,L)**0.3333D0
            US = USTA(J)
            PB(J)=0.5D0*(ERFCC((RB(J)+US)/SQR2 )+ERFCC((RB(J)-US)/SQR2))
            PS(J)=0.5D0*(ERFCC((RS(J)+US)/SQR2 )+ERFCC((RS(J)-US)/SQR2))
            if (PS(J) > PB(J)) PS(J) = PB(J)
            if (IROLL == 0) then
               VS(J) = PS(J)*(EFFF*DFSTA(J) + EFFB*DBSTA(J))/WFSGM1
            else
               VS(J) = PS(J)*(EFFF*DFSTA(J) + EFFB*RBETA(J)*RQ(J))/WFSGM1
            endif
            VS(J) = VS(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
            BQ=BLD
            
            ! Input bedload parameter in Subroutine2 INPUT is adjusted to account for QBREAK=fraction(0.0 to 1.0) of breaking waves. Comment out the next line for no adjustment
            BQ=BQ*(0.5D0+QBREAK(J))
            if (IGMILD == 1) then
               BQ=BLD*(1.D0+BQCOEFF*QBREAK(J))
            endif
            
            ! BDJ added 2012-10-23
            DECAYL = MIN(XB(JSWL(L))/4.D0, 2.D0*TP*CP(1)) ! The decay length
            JDECAY = NINT(DECAYL/DX)! index of decay intrusion length
            ! JDECAY = 1
            ! end BDJ added 2012-10-23
            
            QBX(J) = BQ*PB(J)*GSLOPE(J)*USTD(J)**3
            if (ISEDAV == 1) then
               if (HP(J,L) >= D50) then
                  BRF= 1.D0
               else
                  BRF=(HP(J,L)/D50)**BEDLM
               endif
               VS(J)=BRF*VS(J)
               QBX(J)=BRF*QBX(J)
            endif
            QSX(J) = ASLOPE(J)*UMEAN(J)*VS(J)
            
            ! Add onshore suspended sediment transport due to wave overtopping
            if (IOVER == 1) then
               DUM = H(J)
               if (DUM < HWDMIN) DUM = HWDMIN
               AO=SLPOT
               DUMQ=QO(L)
               QSX(J)=QSX(J)+AO*VS(J)*DUMQ/DUM
            endif
            QRAW(J) = (QBX(J) + QSX(J))/SPORO1
110      end do
111   end do

      ! BDJ added on 2012-10-24
      QSX(1:JDECAY) = QSX(JDECAY)
      QBX(1:JDECAY) = QBX(JDECAY)
      QRAW(1:JDECAY) = QRAW(JDECAY)
      ! end BDJ added on 2012-10-24

      ! If IOVER=0 or JDRY.LE.JR, no wet and dry zone and use scarping formula
      ! If IOVER= 1, compute sediment transport in wet and dry zone
      if (IOVER == 0 .OR. JDRY <= JR) then
         ! Linear extrapolation for scarped slope exceeding TANPHI only if QRAW(JR) is offshore
         JR1 = JR+1
         JE = JR1
         if (QRAW(JR) < 0.D0) then
102         if (BSLOPE(JE,L) > TANPHI) then
               JE = JE+1
               if (JE >= JMAX(L)) goto 103
               goto 102
            endif
         endif
103      JD = JE-JR
         if (JD >= 2) then
            do 104 J=JR1,JE-1
               DUM=DBLE(JE-J)/DBLE(JD)
               QRAW(J)=DUM*QRAW(JR)
104         end do
         endif
         
         ! Subroutine EXTRAPO extrapolates for nodes from J1 to J2
         call EXTRAPO(JR1, JMAXL, QBX)
         call EXTRAPO(JR1, JMAXL, QSX)
         call EXTRAPO(JE, JMAXL, QRAW)
         goto 900
      endif
   endif
   ! End of IANGLE = 0 in wet zone

   ! Wet and Dry Zone for IANGLE = 0 and 1
   ! For node J = JWD to JDRY in wet and dry (WD) zone
   ! PBWD(J) = bedload probability computed in Subroutine PROBWD
   ! PSWD(J) = suspended load probability computed in Subroutine PROBWD
   ! VSWD(J) = suspended sediment volume per unit area(m)
   ! QSXWD(J) = cross-shore suspended sediment transport rate(m*m/s)
   ! QBXWD(J) = cross-shore bedload sediment transport rate(m*m/s)
   ! where hydrodynamic variables in WD zone are computed in Subroutine WETDRY
   !   HDIP(J) = mean water depth adjusted for dip in wet and dry zone used for suspended sediment transport rate
   !   if IVWALL(L) = 0 (no vertical wall at landward end)

999 if (IOVER == 1 .and. JDRY > JR) then
      if (IVWALL(L) == 0) then
         J=JWD
         HDIP(J)=H(J)
         ZBPEAK=ZB(J,L)
140      J=J+1
         if (J == JDRY) goto 142
         if (J > JDRY) goto 145
         if (ZB(J-1,L) < ZB(J,L) .and. ZB(J,L) >= ZB(J+1,L)) ZBPEAK=ZB(J,L)
         DUM=ZBPEAK-ZB(J,L)
         if (DUM > H(J)) then
            HDIP(J)=DUM
         else
            HDIP(J)=H(J)
         endif
         if (J < JCREST(L)) goto 140
142      J=JDRY
         HDIP(J)=H(J)
         ZBPEAK=ZB(J,L)
141      J=J-1
         if (J <= JCREST(L)) goto 145
         if (ZB(J-1,L) < ZB(J,L) .and. ZB(J,L) >= ZB(J+1,L)) ZBPEAK=ZB(J,L)
         DUM=ZBPEAK-ZB(J,L)
         if (DUM > H(J)) then
            HDIP(J)=DUM
         else
            HDIP(J)=H(J)
         endif
         goto 141
      endif
145   continue

      ! For gravel tests MH and MB (IGMILD = 1), landward extension of bedload was necessary
      if (IGMILD == 1) JEXT=JWD+NINT(4.2D0*HRMS(1)/DX)

      do 150 J=JWD,JDRY
         if (IPERM == 0 .and. INFILT == 0) then
            QWX=QO(L)
            if (IPOND == 1 .and. NOPOND == 0) then
               if (J >= JX2) QWX=QM
               if (J > JXW .and. J < JX2) then
                  QWX=QO(L)-(QO(L)-QM)*(XB(J)-XB(JXW))/(XB(JX2)-XB(JXW))
               endif
            endif
         else
            QWX=QO(L)-QP(J)
         endif
         USWD(J)=QWX/H(J)-AQWD*DSQRT(GRAV*H(J)/PWET(J))
         if (D50 < CSEDIA) then
            UCB=DSQRT(GSD50S/FB2(J,L))
         else
            UCB=GSD50S
         endif
         PWAGH=PWET(J)/AGWD/GRAV/H(J)
         call PROBWD(PWET(J),PWAGH,USWD(J),UCB,PBWD(J))
         UCS=WF/FB2(J,L)**0.333333D0
         call PROBWD(PWET(J),PWAGH,USWD(J),UCS,PSWD(J))
         if (PSWD(J) > PBWD(J)) PSWD(J)=PBWD(J)

         ! Suspended load VBF and bedload factor BLDS in wet and dry zone are adjusted so that VS(J)=VSWD(J) and QBX(J)=QBXWD(J) at J=JWD
         if (J == JWD) then
            VBF= 1.D0
            BLDS= 1.D0
         endif
         VSWD(J)=VBF*PSWD(J)
         VSWD(J)=VSWD(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
         QBXWD(J)=BLDS*PBWD(J)*GSLOPE(J)*USTD(J)**3.D0
         if (J == JWD) then
            VBF=VS(J)/VSWD(J)
            VSWD(J)=VS(J)
            BLDS=QBX(J)/QBXWD(J)
            QBXWD(J)=QBX(J)
         endif
         if (ISEDAV == 1) then
            if (HP(J,L) >= D50) then
               BRF= 1.D0
            else
               BRF=(HP(J,L)/D50)**BEDLM
            endif
            if (IVWALL(L) == 0) VSWD(J)=BRF*VSWD(J)
            QBXWD(J)=BRF*QBXWD(J)
         endif
         QSXWD(J)=ASLOPE(J)*UMEAN(J)*VSWD(J)
         if (IOVER == 1) then
            DUM = H(J)
            if (IVWALL(L) == 0) DUM=HDIP(J)
            if (DUM < HWDMIN) DUM = HWDMIN
            AO=SLPOT
            DUMQ=QO(L)
            if (IPOND == 1 .and. NOPOND == 0) then
               if (J >= JCREST(L) .and. J < JX2) DUMQ=QD
               if (J >= JX2) DUMQ=QM
            endif
            QSXWD(J)=QSXWD(J)+AO*VSWD(J)*DUMQ/DUM
         endif

         ! If IGMILD= 1, adjust QBXWD as follows
         if (IGMILD == 1) then
            if (J <= JEXT) QBXWD(J)=QBXWD(JWD)
            if (J == JWD) then
               JWD1=JWD+1
               if (JWD1 < JR) then
                  do 149 JJ=JWD1,JR
                     QBX(JJ)=QBX(JWD)
149               end do
               endif
            endif
         endif
         QRAWD(J)=(QBXWD(J)+QSXWD(J))/SPORO1
         if (IANGLE == 1) then
            US=UMEAN(J)/USTD(J)
            DUM=(1.D0+US*US)*VMEAN(J)/USTD(J)+2.D0*US*STHETA(J)
            QBYWD(J)=QBXWD(J)*DUM/GSLOPE(J)
            QSYWD(J)=VMEAN(J)*VSWD(J)
         endif
150   end do

      ! If IPOND= 1 and NOPOND=0, ponded water exists between nodes J=JXW and JX2. Ponded water is assumed to cause sedimentation where WF=sediment fall velocity and QD=wave-induced onshore volume flux at ridge crest node JCREST for deposition
      if (IPOND == 1 .and. NOPOND == 0) then
         JDUM=JDRY
         DLEN=(XB(JX2)-XB(JXW))/SLPOT
         ! DUM=QD/WF
         ! if (DLEN.LT.DUM) DLEN=DUM
         if (JDUM > JXW) then
            JXW1=JXW+1
            do 151 J=JXW1, JDUM
               DUM=DEXP(-(XB(J)-XB(JXW))/DLEN)
               PBWD(J)=PBWD(J)*DUM
               VSWD(J)=VSWD(J)*DUM
               PSWD(J)=PSWD(J)*DUM
               QBXWD(J)=QBXWD(J)*DUM
               QSXWD(J)=QSXWD(J)*DUM
               QRAWD(J)=(QBXWD(J)+QSXWD(J))/SPORO1
               if (IANGLE == 1) then
                  QBYWD(J)=QBYWD(J)*DUM
                  QSYWD(J)=QSYWD(J)*DUM
               endif
151         end do
         endif
      endif

      ! Connect wet variables (J= 1 to JR) with WD variables (J = JWD to JDRY) using Subroutine TRANWD
      if (JDRY > JR) then
         call TRANWD(PB,JR,PBWD,JWD,JDRY)
         call TRANWD(PS,JR,PSWD,JWD,JDRY)
         call TRANWD(VS,JR,VSWD,JWD,JDRY)
         call TRANWD(QSX,JR,QSXWD,JWD,JDRY)
         call TRANWD(QBX,JR,QBXWD,JWD,JDRY)
         call TRANWD(QRAW,JR,QRAWD,JWD,JDRY)
         if (IANGLE == 1) then
            call TRANWD(QSY,JR,QSYWD,JWD,JDRY)
            call TRANWD(QBY,JR,QBYWD,JWD,JDRY)
         endif
      endif
      
      ! Connect QSX(J), QBX(J) and QRAW(J)=0.0 landward of JDRY
      JDUM=JDRY
      if (IWTRAN == 1 .and. JDRY == JSL1) JDUM=JMAX(L)
      if (JDUM < JMAX(L)) then
         JDRY1=JDRY+1
         call EXTRAPO(JDRY1,JMAXL,QSX)
         call EXTRAPO(JDRY1,JMAXL,QBX)
         call EXTRAPO(JDRY1,JMAXL,QRAW)
         if (IANGLE == 1) then
            call EXTRAPO(JDRY1,JMAXL,QSY)
            call EXTRAPO(JDRY1,JMAXL,QBY)
         endif
      endif
      
      ! Smooth connected PB(J) and PS(J) using Subroutine14 SMOOTH where QRAW(J) is smoothed at end of this subroutine and VS(J),QSX(J),QBX(J),QSY(J) and QBY(J) are smoothed in Subroutine OUTPUT
      PBWD=PB
      call SMOOTH(JDUM,PBWD,PB)
      PSWD=PS
      call SMOOTH(JDUM,PSWD,PS)
      goto 900

   endif
   ! End of Wet and Dry Zone for IANGLE = 0 and 1

   ! Obliquely Incident Waves in wet zone
888 if (IANGLE == 1) then
      LMAX=LANCOM+1
      do 191 K= 1,LMAX
         if (K == 1) then
            JSTART= 1
            JEND=JR
         else
            JSTART=JSL
            JEND=JMAX(L)
         endif
         do 190 J=JSTART,JEND
            SIGT = USTD(J)/CTHETA(J)
            if (D50 < CSEDIA) then
               RB(J)= DSQRT(GSD50S/FB2(J,L))/SIGT
            else
               RB(J)=GSD50S/SIGT
            endif
            RS(J)= WF/SIGT/FB2(J,L)**0.3333D0
            WSTA = USTA(J)*CTHETA(J) + VSTA(J)*STHETA(J)
            VCUS = VSTA(J)*CTHETA(J) -USTA(J)*STHETA(J)
            FS = RS(J)*RS(J) - VCUS*VCUS
            if (FS < 0.D0) then
               PS(J) = 1.D0
            else
               FS = DSQRT(FS)
               PS(J)= 0.5D0*(ERFCC((FS+WSTA)/SQR2)+ERFCC((FS-WSTA)/SQR2))
            endif
            FB = RB(J)*RB(J) - VCUS*VCUS
            if (FB < 0.D0) then
               PB(J) = 1.D0
            else
               FB = DSQRT(FB)
               PB(J)= 0.5D0*(ERFCC((FB+WSTA)/SQR2)+ERFCC((FB-WSTA)/SQR2))
            endif
            if (PS(J) > PB(J)) PS(J)=PB(J)
            if (IROLL == 0) then
               VS(J) = PS(J)*(EFFF*DFSTA(J)+EFFB*DBSTA(J))/WFSGM1
            else
               VS(J) = PS(J)*(EFFF*DFSTA(J)+EFFB*RBETA(J)*RQ(J))/WFSGM1
            endif
            VS(J) = VS(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
            VSTA2 = VSTA(J)*VSTA(J)
            TWOS = 2.D0*STHETA(J)
            BQ=BLD
            
            ! Input bedload parameter in Subroutine2 INPUT is adjusted to account for QBREAK=fraction(0.0 to 1.0) of breaking waves. Comment out the next line for no adjustment
            BQ=BQ*(0.5D0+QBREAK(J))
            if (IGMILD == 1) then
               BQ=BLD*(1.D0+BQCOEFF*QBREAK(J))
            endif

            ! BDJ added 2012-10-23
            DECAYL = MIN(XB(JSWL(L))/4.D0, 2.D0*TP*CP(1)) ! The decay length
            JDECAY = NINT(DECAYL/DX)! index of decay intrusion length
            ! end BDJ added 2012-10-23
            
            DUM = BQ*PB(J)*(USTD(J)*USTD(J) + VSTD(J)*VSTD(J))**1.5D0
            QBX(J) = DUM*GSLOPE(J)*(1.D0+USTA(J)*VSTA2+TWOS*VCUS)
            QBY(J) = DUM*(VSTA(J)*(1.D0 + USTA(J)*USTA(J)+ VSTA2)+ TWOS*WSTA)
            if (ISEDAV == 1) then
               if (HP(J,L) >= D50) then
                  BRF= 1.D0
               else
                  BRF=(HP(J,L)/D50)**BEDLM
               endif
               VS(J)=BRF*VS(J)
               QBX(J)=BRF*QBX(J)
               QBY(J)=BRF*QBY(J)
            endif
            QSX(J) = ASLOPE(J)*UMEAN(J)*VS(J)
            if (IOVER == 1) then
               DUM = H(J)
               if (DUM < HWDMIN) DUM = HWDMIN
               AO=SLPOT
               DUMQ=QO(L)
               QSX(J)=QSX(J)+AO*VS(J)*DUMQ/DUM
            endif
            QSY(J) = VMEAN(J)*VS(J)
            QRAW(J) = (QBX(J) + QSX(J))/SPORO1
190      end do
191   end do

      ! BDJ added on 2012-10-24
      QSX(1:JDECAY) = QSX(JDECAY)
      QBX(1:JDECAY) = QBX(JDECAY)
      QRAW(1:JDECAY) = QRAW(JDECAY)
      ! end BDJ added on 2012-10-24

      ! Scarping extrapolation is included for oblique waves as well
      if (IOVER == 0 .OR. JDRY <= JR) then
         JR1 = JR+1
         JE = JR1
         if (QRAW(JR) < 0.D0) then
202         if (BSLOPE(JE,L) > TANPHI) then
               JE = JE+1
               if (JE >= JMAX(L)) goto 203
               goto 202
            endif
         endif
203      JD = JE-JR
         if (JD >= 2) then
            do 204 J=JR1, JE-1
               DUM=DBLE(JE-J)/DBLE(JD)
               QRAW(J) =DUM*QRAW(JR)
204         end do
         endif

         ! Subroutine EXTRAPO extrapolates for nodes from J1 to J2
         call EXTRAPO(JR1, JMAXL, QBX)
         call EXTRAPO(JR1, JMAXL, QSX)
         call EXTRAPO(JR1, JMAXL, QBY)
         call EXTRAPO(JR1, JMAXL, QSY)
         call EXTRAPO(JE, JMAXL, QRAW)
         goto 900
      else
         goto 999
      endif

   endif
   ! End of IANGLE= 1 in wet zone

   ! Adjust computed QSX(1) and QBX(1) at node 1 to be consistent with the boundary condition used in Subroutine CHANGE
900 QSX(1)=QSX(2)
   QBX(1)=QBX(2)
   QRAW(1)=QRAW(2)
   
   ! Adjust sediment transport rates at node JMAX to be consitent with the boundary condition used in Subroutine CHANGE
   JMAX1=JMAX(L)-1
   QSX(JMAXL)=QSX(JMAX1)
   QBX(JMAXL)=QBX(JMAX1)
   QRAW(JMAXL)=QRAW(JMAX1)
   
   ! Smoothing QRAW (before DELZB is computed) using Subroutine SMOOTH
   call SMOOTH(JMAXL,QRAW,Q)
   
   deallocate(QRAW, GSLRAW, ASLRAW, ASLOPE, RS, RB, PBWD, PSWD, VSWD, QSXWD, QBXWD, QRAWD, HDIP, QSYWD, QBYWD)

   return
end subroutine SEDTRA
