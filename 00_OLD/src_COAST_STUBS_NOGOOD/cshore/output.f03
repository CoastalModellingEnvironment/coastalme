!===============================================================================================================================
!
! This subroutine stores computed and input quantities
!
!===============================================================================================================================
subroutine OUTPUT(ITIME, L, ITEQO, ICONV)
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
      
      subroutine QORATE(ITIME, L, ITEQO, ICONV, ICALL)
         integer, intent(in) :: ITIME, L, ITEQO, ICALL
         integer, intent(inout) :: ICONV
      end subroutine QORATE
   end interface
   
   integer, intent(in) :: ITIME, L, ITEQO
   integer, intent(out) :: ICONV
   
   double precision, allocatable, dimension(:) :: EDEPTH
   double precision, allocatable, dimension(:, :) :: ZB0
   
   allocate(EDEPTH(JMAXMAX))
   allocate(ZB0(JMAXMAX, ILINE))

   ! OUTPUT ONLY WHEN ITIME = 0
   if (ITIME == 0) then
      if (L > 1) goto 888

      ! COMPUTATIONAL OPTION
      ! ILINE = number of cross-shore lines
      
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 890) ILINE, IQYDY
890   format('ILINE (number of profiles) =', I3 / 'Alongshore gradient IQYDY =', I3 /)

      ! IPROFL == 0 for fixed bottom profile
      ! IPROFL == 1 for profile evolution computation
      if (IPROFL == 0) then
         write (20, 900) IPROFL
      endif
900   format('COMPUTATION OPTION IPROFL =', I3, ' (profile bottom is fixed, sediment transport not calculated)'/)

      if (IPROFL == 1) then
         write (20, 901) IPROFL, TIMEBC(NTIME+1), NTIME         
         if (ISEDAV == 1) write (20, 902) ISEDAV, BEDLM
      endif      
901   format('COMPUTATION OPTION IPROFL =', I3 / 'Profile evolution is computed from time = 0.0'/ 'to time =', F13.1, 'for NTIME =', I4 /)
902   format(/'ISEDAV =', I3, ' for hard bottom with bedload reduction factor BEDLM =', F4.1/)

      if (IROLL == 0) write (20, 910)
910   format('Rollers not included in computation'/)

      if (IROLL == 1) write (20, 911) RBZERO
911   format('Rollers included in computation' / 'Roller slope betazero =', F9.3 /)

      if (IWCINT == 0) write (20, 920)
920   format('Wave and current interaction not considered'/)

      if (IWCINT == 1) write (20, 921)
921   format('Wave and current interaction considered'/)

      if (IOVER == 0) write (20, 930)
930   format('No wave overtopping, overflow and seepage considered'/)
      
      if (IOVER == 1 .and. IPOND == 0) then
         write (20, 931) RWH, JCREST(L), RCREST(L), AWD, EWD
      endif
      
      if (IOVER == 1 .and. IPOND == 1) write (20, 932) RWH, AWD, EWD, ZW
931   format('Wwave overtopping, overflow and seepage' / &
         'Runup wire height (m)               RWH =', F9.3 / &
         'Initial crest location for L= 1      JCREST =', I6/ &
         'Initial crest height (m) for L= 1    RCREST =', F9.3/ &
         'Swash velocity parameter            AWD =', F9.3/ &
         'Output exceedance probability       EWD =', F9.3/)
932   format('PONDED WATER IN RUNNEL'/ &
         'Runup wire height (m)               RWH =', F9.3/ &
         'Swash velocity parameter            AWD =', F9.3/ &
         'Output exceedance probability       EWD =', F9.3/ &
         'Initial ponded water level (m)       ZW =', F9.3/)

      if (IPERM == 0) write (20, 940)
940   format('Impermeable bottom assumed'/)

      if (IPERM == 1) write (20, 941) SNP, SDP, WNU, WPM
941   format('Permeable bottom consisting of'/ &
         'Stone porosity                       NP =', F9.3/ &
         'Nominal stone diameter (m)         DN50 =', F9.4/ &
         'Water kinematic viscosity(m*m/s)        =', F9.7/ &
         'Maximum seepage velocity (m/s)      WPM =', F9.5/)

      if (IWIND == 0) write (20, 950)
950   format('Wind shear stresses not considered'/)

      if (IWIND == 1) write (20, 951) NWIND
951   format('Wind shear stresses in momentum equations'/ &
         'Number of wind speed and direction input =', I4/)

      write (20, 960) GAMMA
960   format('Breaker ratio parameter'/ 'Input gamma =', F5.2/)

      if (IPROFL == 1) write (20, 970) D50 * 1000.D0, WF, SG, EFFB, SLP, TANPHI, BLP
970   format('Sediment parameter if IPROFL= 1'/ &
         'Median diameter (mm)                D50 = ', F5.2/ &
         'Fall velocity (m/s)                  Wf = ', F6.4/ &
         'Specific gravity                     Sg = ', F5.2/ &
         'Suspension efficiency                eB = ', F6.3/ &
         'Suspended load parameter                = ', F5.2/ &
         'Tangent(friction angle)                 = ', F5.2/ &
         'Bedload parameter                     b = ', F6.4)      
      
      if (IPROFL == 1 .and. IOVER == 1) then
         write (20, 971) SLPOT
971   format('Suspended load parameter (IOVER == 1)              = ', F5.2/)
         
         if (INFILT == 1) write (20, 972) WPM
972   format('INFILT == 1 and infiltration velocity (m/s)   = ', F7.5/)

      endif      
#endif

      ! WAVE AND WATER LEVEL
      ! mg   nout = 10000
      nout = 1000
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1001) NTIME, nout
1001  format(/'Input wave and water level: NTIME =', I6,' from TIME = 0.0'/ 'Output lines are limited to first and last', I6, ' lines'/ &
         'End Time(sec) Tp (sec)  Hrms(m) Wave Setup(m) Storm tide(m) Angle(deg)'/)
#endif

      NTIM9 = nout + 1
      if (NTIME > (2*nout)) NTIM9=NTIME-(nout-1)
      
      do 130 I = 1,NTIME
         if (I <= nout .OR. I >= NTIM9) then
#if defined FILEINOUT || ARGINBOTHOUT
            write (20, 1002) TIMEBC(I+1), TPBC(I), HRMSBC(I), WSETBC(I), SWLBC(I), WANGBC(I)
1002  format(F11.1, 5F11.4)
#endif

         endif
130   end do

888   continue
      ! End of L= 1 output................................................

      ! OUTPUT BOTTOM GEOMETRY
      ! The bottom geometry is divided into segments of different inclination and roughness starting from the seaward boundary for ILINE cross-shore lines.
      ! NBINP(L)  = number of input points for L line
      ! XBINP(J,L)  = horizontal distance from seaward boundary to landward-end of segment (J-1) in meters
      ! ZBINP(J,L)  = dimensional vertical coordinate (+ above datum) of the landward end of segment (J-1) in meters
      ! FBINP(J,L) = bottom friction factor
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1100) L, YLINE(L), 0.D0-ZBINP(1,L), NBINP(L), DX, XS(L), JMAX(L)
1100  format (/'INPUT BEACH AND STRUCTURE GEOMETRY'/ &
         'Cross-shore line number               L =', I3/ &
         'Alongshore coordinate             YLINE =', F13.4/ &
         'Depth at seaward boundary (m)           =', F13.6/ &
         'Number of input points            NBINP =', I8/ &
         'Output lines are limited to first and last 5 lines'/ &
         'Node spacing (m)                     DX =', F13.6/ &
         'Shoreline location (m) of Zb=0       Xs =', F13.6/ &
         'Maximum landward node              JMAX =', I8// &
         'X  (m)        Zb  (m)    Fric.factor')
#endif

         
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1200) XBINP(1,L), ZBINP(1,L)
#endif

      NBINP4= 6
      if (NBINP(L) > 10) NBINP4=NBINP(L)-4
      
      do 140 J = 2,NBINP(L)
         if (J <= 5 .OR. J >= NBINP4) then
#if defined FILEINOUT || ARGINBOTHOUT
            write (20, 1200) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L)
1200  format(3F10.3)
#endif
         endif
140   end do

      if (IPERM == 1 .OR. ISEDAV == 1) then
#if defined FILEINOUT || ARGINBOTHOUT
         write (20, 1150) NPINP(L)
1150  format(/'INPUT IMPERMEABLE HARD BOTTOM GEOMETRY'/ 'Number of input points for line L =',I3, 'NPINP = ',I5 // 'X  (m)        ZP  (m)  ')
#endif
         NPINP4= 6
         if (NPINP(L) > 10) NPINP4=NPINP(L)-4
         do 141 J= 1,NPINP(L)
            if (J <= 5 .OR. J >= NPINP4) then
            
#if defined FILEINOUT || ARGINBOTHOUT
               write (20, 1201) XPINP(J,L), ZPINP(J,L)
1201  format(2F10.3)
#endif
            endif
141      end do
      endif
      
      ! Where the number of the output lines is limited to 10 or less to reduce the length of the output file ODOC.

      ! INPUT WIND SHEAR STRESSES FOR IWIND == 1
      if (L > 1) goto 889
      
      if (IWIND == 1) then
      
#if defined FILEINOUT || ARGINBOTHOUT
         write (20, 1370)
1370  format(/'INPUT WIND SPEED, DIRECTION AND STRESSES' / 'Start & End Time(s) Speed(m/s) Dir(deg) DragCoef'/)
#endif

         do 142 I= 1,NTIME
            if (I <= 10 .OR. I >= NTIM9) then
            
#if defined FILEINOUT || ARGINBOTHOUT
               write (20, 1371) TIMEBC(I),TIMEBC(I+1),W10(I),WANGLE(I), WINDCD(I)
1371  format(2F11.1, 2F11.2,E11.4)
#endif
            endif
142      end do
      endif

      ! INPUT LANDWARD STILL WATER LEVEL FOR IWTRAN == 1
      if (IWTRAN == 1) then
         if (ISWLSL == 0) then
#if defined FILEINOUT || ARGINBOTHOUT
            write (20, 1380)
1380  format(/'INPUT LANDWARD STILL WATER LEVEL'/ 'same as input seaward still water level'/)
#endif
         else
#if defined FILEINOUT || ARGINBOTHOUT
            write (20, 1381)
1381  format(/'INPUT LANDWARD STILL WATER LEVEL'/ 'Start & End Time(s) SWL(m) above datum'/)
#endif
            do 143 I= 1,NTIME
               if (I <= 10 .OR. I >= NTIM9) then
               
#if defined FILEINOUT || ARGINBOTHOUT
                  write (20, 1382) TIMEBC(I), TIMEBC(I+1), SWLAND(I)
1382  format(2F11.1,F11.4)
#endif

               endif
143         end do
         endif
      endif
      
      ! INPUT ALONGSHORE WATER LEVEL GRADIENT FOR ITIDE= 1
      if (ITIDE == 1) then
#if defined FILEINOUT || ARGINBOTHOUT
         write (20, 1390)
1390  format(/'INPUT ALONGSHORE WATER LEVEL GRADIENT'/ 'Start & End Time(s)      DETA/DY alongshore'/)
#endif

         do 144 I= 1,NTIME
            if (I <= 10 .OR. I >= NTIM9) then     
            
#if defined FILEINOUT || ARGINBOTHOUT
               write (20, 1383) TIMEBC(I), TIMEBC(I+1), DETADY(I)
1383     format(2F11.1,F11.7)
#endif
            endif
144      end do
      endif

889   continue
      ! End of L= 1 OUTPUT
   endif
   ! --------------------- END OF OUTPUT ONLY WHEN ITIME = 0 --------------

   ! ------------------- COMPUTED CROSS-SHORE VARIATIONS -----------------
   ! For each cross-shore line L of ILINE lines stored at Time = 0.0 and at the end of constant wave and water level at the seaward boundary if laboratory data (ILAB = 1).
   ! For field data (ILAB = 0), stored at the beginning, end, and every ten storage time levels (goto 200 goes to end of this subroutine)
   if (ITIME == 0) then
#if defined FILEINOUT || ARGINBOTHOUT
      ! write (21, 1490) L, JMAX(L), TIMEBC(ITIME)
      
      ! DFM Got run-time error accessing TIMEBC(0) since TIMEBC starts with element 1, so changed this
      write (21, 1490) L, JMAX(L), TIMEBC(ITIME+1)
#endif
      do 199 J= 1,JMAX(L)
         if (IPERM == 0 .and. ISEDAV == 0) then
#if defined FILEINOUT || ARGINBOTHOUT
            write (21, 1500) XB(J), ZB(J,L)
#endif
         else         
#if defined FILEINOUT || ARGINBOTHOUT
            write (21, 1500) XB(J), ZB(J,L), ZP(J,L)
#endif
            if (IPROFL == 1) ZB0(J,L)=ZB(J,L)
         endif
199   end do
      goto 200
   endif

   TIMOUT = TIME
   ! mg - explicit declaration of laboratory/field data sets
   if (ILAB == 0) then
      ! mg - ensure output at end of simulation for field data sets
      if (ITIME == NTIME) goto 201
      ! mg
      ! DUM=DBLE(ITIME)/10.D0-DBLE(ITIME/10)
      ! if (DUM.NE.0.D0) goto 200
   endif
201 continue

   if (IPROFL == 0) then
      TIMOUT = TIMEBC(ITIME+1)
      
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1440) TIMOUT, L
1440 format(/'**********COMPUTED CROSS-SHORE VARIATIONS**********'/ 'on input bottom profile at TIME =',F11.1, ' Line L =',I3/)

#endif

   else
   
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1441) TIMOUT, L
1441 format(/'**********COMPUTED CROSS-SHORE VARIATIONS**********'/ 'on bottom profile computed at TIME (s) = ', F11.1, ' Line L =',I3/)
#endif

   endif
   
#if defined FILEINOUT || ARGINBOTHOUT
   write (20, 1450) JR, XR, ZR, H(JR)
1450 format('LANDWARD WET COMPUTATION LIMIT'/ 'Most landward node of wet zone computation    JR =', I8/ 'X-coordinate at JR (m)                        XR =', F13.6/ 'Bottom elevation at JR (m)                    ZR =', F13.6/ 'Mean water depth at this node (m)          H(JR) =', F13.6/)
#endif


   ! Wave Reflection Coeffiecient at node J= 1 only for IOVER=0
   if (IOVER == 0) then
      if (JR > JSWL(L) .and. JSWL(L) < JMAX(L)) then
         DUM = SIGMA(JSWL(L))*SIGMA(JSWL(L))*CP(JSWL(L))*WN(JSWL(L))
         DUM = DUM/WN(1)/CP(1)
         SIGREF=DSQRT(DUM)
         if (IANGLE == 1) SIGREF=SIGREF/DSQRT(CTHETA(1)/CTHETA(JSWL(L)))
         REFCOF=SIGREF/SIGMA(1)
#if defined FILEINOUT || ARGINBOTHOUT
         write (20, 1460) REFCOF, JSWL(L)
1460 format('WAVE REFLECTION COEFFICIENT'/ 'Wave reflection coefficient (at x = 0) = ', F9.6/ 'Still water shoreline node location JSWL =', I5/)
#endif
      endif
   endif

   ! Output computed wave overtopping, overflow and seepage rates in Subroutine QORATE
   if (IOVER == 1) call QORATE(ITIME, L, ITEQO, ICONV, 1)

   ! Longshore (Suspended Sediment plus Bedload) Transport Rate
   if (IPROFL == 1 .and. IANGLE == 1) then
      DUM = 0.5D0*(QBY(1)+QSY(1))
      do 145 J = 2,JDRY-1
         DUM = DUM + (QBY(J)+QSY(J))
145   end do
      DUM = DUM + 0.5D0*(QBY(JDRY)+QSY(JDRY))
      QLONG = DUM * DX
      SIGMAX = SIGMA(1)
      JB= 1
      do 150 J = 2, JR
         if (SIGMA(J) > SIGMAX) then
            SIGMAX = SIGMA(J)
            JB = J
         endif
150   end do
      DUM = SIGMA(JB)**2.D0*CP(JB)*WN(JB)*CTHETA(JB)*STHETA(JB)
      CERCK = (SG-1.D0)*QLONG/DUM
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1470) QLONG, CERCK, STHETA(JB)
1470 format('LONGSHORE SUSPENDED AND BEDLOAD SAND TRANSPORT RATE'/ 'Transport Rate (m**3/s) =',E14.5/'CERC Formula K=',F11.3/ 'sin(breaker angle)=',F11.5/)
#endif

   endif

   ! Damage (normalized eroded area) of stone structure
   ! EDMAX = normalized maximum vertical erosion depth
   if (IPROFL == 1 .and. IPERM == 1) then
      EDMAX=0.D0
      do 300 J= 1,JMAX(L)
         EDEPTH(J)=ZB0(J,L)-ZB(J,L)
         if (EDEPTH(J) > EDMAX) EDMAX=EDEPTH(J)
         if (EDEPTH(J) < 0.D0) EDEPTH(J)=0.D0
300   end do

      EDMAX=EDMAX/SDP
      JMAXL=JMAX(L)
      call INTGRL(JMAXL, DX, EDEPTH, AREA)
      DAMAGE=AREA/SDP/SDP
      STABNO=SQR2*HRMS(1)/SDP/(SG-1.D0)
      
#if defined FILEINOUT || ARGINBOTHOUT
      write (20, 1480) DAMAGE, EDMAX, STABNO
1480 format('STONE STRUCTURE DAMAGE'/ 'Damage S=',F10.3/ 'Normalized erosion depth E=',F10.3/ 'Stability number Nmo=',F8.3/)
#endif

   endif

   ! COMPUTED CROSS-SHORE VARIATIONS
   ! Indicate the number of lines used for storage at specified time for each cross-shore line L= 1, 2,...,ILINE
   JSWASH = JDRY - JWD +1
   JDUM = JR
   
   if (IOVER == 1) then
      JDUM = JDRY
      if (IWTRAN == 1 .and. JDRY == JSL1) JDUM = JMAX(L)
   endif
   
#if defined FILEINOUT || ARGINBOTHOUT
   write (22, 1490) L, JDUM, TIMOUT
   write (23, 1490) L, JR, TIMOUT
   write (24, 1490) L, JR, TIMOUT
   write (25, 1490) L, JR, TIMOUT
   write (26, 1490) L, JR, TIMOUT
   write (27, 1490) L, JDUM, TIMOUT
   write (28, 1490) L, JDUM, TIMOUT
   write (29, 1490) L, JR, TIMOUT
   write (30, 1490) L, JDUM, TIMOUT
   write (31, 1490) L, JR, TIMOUT
   write (32, 1490) L, JMAX(L), TIMOUT
   write (33, 1490) L, JMAX(L), TIMOUT
   write (37, 1490) L, JMAX(L), TIMOUT
   write (38, 1490) L, JMAX(L), TIMOUT
#endif

   if (IOVER == 1) then
#if defined FILEINOUT || ARGINBOTHOUT
      write (34, 1490) L, JDUM, TIMOUT
      write (35, 1490) L, JSWASH, TIMOUT
1490 format(2I8,F11.1)
#endif

      TIMID=0.5D0*(TIMEBC(ITIME)+TIMEBC(ITIME+1))
      DUM=TIMEBC(ITIME+1)-TIMEBC(ITIME)
      
#if defined FILEINOUT || ARGINBOTHOUT
      write (36, 1491) L, TIMID, (TSQO(L) / DUM), (TSQBX(L) / DUM), (TSQSX(L) / DUM)
1491 format(I8, 4F17.9)
#endif

   endif
   
   if (IPROFL == 1 .and. L == ILINE) then
      do 181 LL= 1,ILINE
#if defined FILEINOUT || ARGINBOTHOUT
         write (21, 1490) LL, JMAX(LL), TIMOUT
#endif
         do 180 J= 1,JMAX(LL)
            if (IPERM == 0 .and. ISEDAV == 0) then
#if defined FILEINOUT || ARGINBOTHOUT
               write (21, 1500) XB(J), ZB(J,LL)
#endif
            else
#if defined FILEINOUT || ARGINBOTHOUT
               write (21, 1500) XB(J),ZB(J,LL),ZP(J,LL)
#endif
            endif
180      end do
181   end do
   endif

   ! Smooth computed VS(J) before storing and plotting
   DUMVEC=VS
   call SMOOTH(JDUM, DUMVEC, VS)

   do 160 J = 1, JR
#if defined FILEINOUT || ARGINBOTHOUT   
      write (22, 1500) XB(J), (H(J)+ZB(J,L)), H(J), SIGMA(J)
      
      write (23, 1500) XB(J), WT(J), QBREAK(J), SIGSTA(J)
      
      write (24, 1500) XB(J), SXXSTA(J), TBXSTA(J)
      
      if (IANGLE == 1) write (25, 1500) XB(J), SXYSTA(J), TBYSTA(J)
      
      write (26, 1500) XB(J), EFSTA(J) / WT(J), DBSTA(J), DFSTA(J)
      
      if (IPERM == 0) then
         write (27, 1500) XB(J), UMEAN(J), USTD(J)
      else
         write (27, 1500) XB(J), UMEAN(J), USTD(J), UPMEAN(J)
      endif
      
      if (IANGLE == 1) then
         write (28, 1500) XB(J), STHETA(J), VMEAN(J), VSTD(J)
      else
         write (28, 1500) XB(J), 0., 0., 0.
      endif
      
      if (IROLL == 1) write (29, 1500) XB(J), RQ(J)
      
      if (IPROFL == 1) write (30, 1500) XB(J), PB(J), PS(J), VS(J)
      
      if (IPERM == 1) write (31, 1500) XB(J), UPSTD(J), DPSTA(J)
#endif

#if defined ARGINOUT 
      ! DFM The second dimension of each output array is always 1 at present, since CShore currently calculates only one coastline-normal profile per CShore run

      ! For OSETUP, write XB(J), (H(J)+ZB(J,L)), H(J), SIGMA(J)   
      XYDist(J, 1) = XB(J)
      FreeSurfaceStd(J, 1) = SIGMA(J)
      WaveSetupSurge(J, 1) = H(J) + ZB(J,1) ! setup plus surge
      ! StormSurge(J, 1) = 0 ! maybe at future
      
      ! For OYVELO, write XB(J), STHETA(J), VMEAN(J), VSTD(J)
      if (IANGLE == 1) then
         SinWaveAngleRadians(J, 1) = STHETA(J)
      else
         SinWaveAngleRadians(J, 1) = 0.0d0
      endif
      
      ! For OPARAM, write XB(J), WT(J), QBREAK(J), SIGSTA(J)
      FractionBreakingWaves(J, 1) = QBREAK(J)    

#endif

160 end do

   if (IOVER == 1) then
      ! Store mean values over wet duration
      if (JR < JDRY) then
         do 170 J=(JR+1),JDUM
            DUM=H(J)+ZB(J,L)
            if (IPOND == 1 .and. NOPOND == 0) then
               if (JX2 < JMAX(L)) then
                  if (JXW <= J .and. J <= JX2) then
                     DUM=H(J)+ZW
                     PWET(J)= 1.D0
                  endif
               endif
            endif
            
#if defined FILEINOUT || ARGINBOTHOUT
            write (22, 1500) XB(J),DUM,H(J),SIGMA(J)
            if (IPERM == 0) then
               write (27, 1500) XB(J),UMEAN(J),USTD(J)
            else
               write (27, 1500) XB(J),UMEAN(J),USTD(J),UPMEAN(J)
            endif
            if (IANGLE == 1) then
               write (28, 1500) XB(J),STHETA(J),VMEAN(J),VSTD(J)
            else
               write (28, 1500) XB(J), 0., 0., 0.
            endif
            write (30, 1500) XB(J),PB(J),PS(J),VS(J)
#endif
170      end do
      endif
      
      ! Where UPMEAN, PB, PS, VS, and QP include effects of PWET
      do 171 J= 1,JDUM
#if defined FILEINOUT || ARGINBOTHOUT
         if (IPERM == 0) write (34, 1500) XB(J),PWET(J)
         if (IPERM == 1) write (34, 1500) XB(J),PWET(J),QP(J)
#endif
171   end do
      do 161 J=JWD,JDRY
#if defined FILEINOUT || ARGINBOTHOUT     
         write (35, 1500) XB(J),HEWD(J),UEWD(J),QEWD(J)
#endif
161   end do
   endif

   if (IPROFL == 1) then
      ! Smooth computed QBX(J) and QSX(J) before storing and plotting
      JMAXL=JMAX(L)
      DUMVEC = QBX
      call SMOOTH(JMAXL,DUMVEC,QBX)
      DUMVEC = QSX
      call SMOOTH(JMAXL,DUMVEC,QSX)
      
      ! Smooth computed QBY(J) and QSY(J) if IANGLE= 1
      if (IANGLE == 1) then
         DUMVEC=QBY
         call SMOOTH(JMAXL,DUMVEC,QBY)
         DUMVEC=QSY
         call SMOOTH(JMAXL,DUMVEC,QSY)
      endif
      do 162 J= 1,JMAX(L)
#if defined FILEINOUT || ARGINBOTHOUT
         write (32, 1500) XB(J),QBX(J),QSX(J),(QBX(J)+QSX(J))
         if (IANGLE == 1) then
            write (33, 1500) XB(J),QBY(J),QSY(J),(QBY(J) + QSY(J))
         else
            write (33, 1500) XB(J), 0., 0., 0.
         endif
#endif
162   end do

      ! Store sediment transport volume per unit width during Time = 0.0 to Time = TIMOUT
      do 163 J= 1,JMAX(L)
#if defined FILEINOUT || ARGINBOTHOUT
         write (37, 1500) XB(J),VBX(J,L),VSX(J,L),(VBX(J,L)+VSX(J,L))
         if (IANGLE == 1) write (38, 1500) XB(J),VBY(J,L),VSY(J,L), (VBY(J,L)+VSY(J,L))
1500 format(4F17.9)
#endif
163   end do
   endif

200 continue

   deallocate(EDEPTH, ZB0)

   return
end subroutine OUTPUT
