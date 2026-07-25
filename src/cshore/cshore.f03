!===============================================================================================================================
!
! The CShore subroutine initialises data then marches from the offshore boundary node toward the shoreline
!
!===============================================================================================================================
#if defined EXE
program CShore
#else
#if defined FILEINOUT
subroutine CShore(NRET) bind(c, name = "CShore")
#else
subroutine CShore(NRET)   
#endif
#endif

   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)   
  
!#if ! defined EXE || defined FILEINOUT
!   integer, parameter, intent(inout) :: NRET
   integer :: NRET
!#endif

   ! For iteration convergence, MAXITE is maximum number of iterations (DFM originally set to 20)
   integer, parameter :: MAXITE = 40
   
   ! For iteration convergence, EPS1 = 0.001 for depth (m), height (m) and velocity (m/s)
   ! For iteration convergence, EPS2 = 0.000001 for roller volume flux (m*m/s)
   double precision, parameter :: EPS1 = 1.D-3, EPS2 = 1.D-6
   
   interface 
      subroutine BOTTOM
      end subroutine BOTTOM
      
      subroutine CHANGE(ITIME, L, IEND, ICALL)
         integer, intent(in) :: ITIME, L, ICALL
         integer, intent(out) :: IEND
      end subroutine CHANGE
      
      subroutine DBREAK(J, L, WHRMS, D)
         integer, intent(in) :: J, L
         double precision, intent(in) :: WHRMS, D
      end subroutine DBREAK
      
      subroutine LWAVE(J, L, WD, QDISP)
         integer, intent(in) :: J, L
         double precision, intent(in) :: WD, QDISP
      end subroutine LWAVE
      
      subroutine PARAM
      end subroutine PARAM
      
      subroutine OUTPUT(ITIME, L, ITEQO, ICONV)
         integer, intent(in) :: ITIME, L, ITEQO
         integer, intent(out) :: ICONV
      end subroutine OUTPUT
      
      subroutine POFLOW(J, L, PKHSIG, DEDX)
         integer, intent(in) :: J, L
         double precision, intent(in) :: PKHSIG, DEDX
      end subroutine POFLOW   
      
      subroutine SEDTRA(L)
         integer, intent(in) :: L   
      end subroutine SEDTRA
      
      subroutine SMOOTH(NUM, RAW, F)
         integer, intent(in) :: NUM
         double precision, intent(in), dimension(:) :: RAW
         double precision, intent(out), dimension(:) :: F
      end subroutine SMOOTH
      
      subroutine SRFSP(L)
         integer, intent(in) :: L
      end subroutine SRFSP
   
      subroutine TRANWD(F1, JR, F2, JS, JE)
         integer, intent(in) :: JR, JS, JE
         double precision, intent(in), dimension(:) :: F2
         double precision, intent(out), dimension(:) :: F1
      end subroutine TRANWD
      
      subroutine GBXAGF(CTHETA_IN, USIGT, STHETA_IN, VSIGT, GBX_IN, CF_IN)
         double precision, intent(in) :: CTHETA_IN, USIGT, STHETA_IN, VSIGT
         double precision, intent(out) :: GBX_IN, CF_IN
      end subroutine GBXAGF
      
      subroutine QORATE(ITIME, L, ITEQO, ICONV, ICALL)
         integer, intent(in) :: ITIME, L, ITEQO, ICALL
         integer, intent(inout) :: ICONV
      end subroutine QORATE
      
      subroutine PONDED(L)
         integer, intent (in) :: L
      end subroutine PONDED

! DFM In original code, parameter USIGT is passed but not used
!      subroutine VSTGBY(CTHETA_IN, USIGT, STHETA_IN, VSIGT, GBY_IN)
!         double precision, intent(in) :: CTHETA_IN, USIGT, STHETA_IN, GBY_IN
!         double precision, intent(out) :: VSIGT
!      end subroutine VSTGBY
      subroutine VSTGBY(CTHETA_IN, STHETA_IN, VSIGT, GBY_IN)
         double precision, intent(in) :: CTHETA_IN, STHETA_IN, GBY_IN
         double precision, intent(out) :: VSIGT
      end subroutine VSTGBY

      
#if defined FILEINOUT
      subroutine INPUT_OPENER      
      end subroutine INPUT_OPENER
      
      subroutine READ_INPUT
      end subroutine READ_INPUT
#endif
      
#if defined FILEINOUT || ARGINBOTHOUT
      subroutine OUTPUT_OPENER      
      end subroutine OUTPUT_OPENER
#endif
   end interface   
   
   ! Start of executable section
   ! write (*, *) "In CShore"   
   
   ! Set flag for normal completion
   NRET = 0  
   
#if defined FILEINOUT
   ! We are reading input data from an ASCII file, and writing output data to other ASCII files
   call INPUT_OPENER
   call OUTPUT_OPENER

   write (20, *) VER
   write (20, *) "Input and output via ASCII files"
   write (20, *)
   
   call READ_INPUT
   
#elif defined ARGINBOTHOUT
   ! We are reading input data as arguments from the calling program, and returning output data the same way. But we also are outputting data to ASCII files for debugging purposes, so call OUTPUT_OPENER to open output files
   call OUTPUT_OPENER
      
   write (20, *) VER
   write (20, *) "Input and output via arguments, also outputting to ASCII files"
   write (20, *)
#endif

   ! Subroutine 3 BOTTOM computes initial bathymetry at each node
   call BOTTOM
   
   ! Subroutine 4 PARAM calculates constants
   call PARAM

   ! ************* START OF TIME MARCHING COMPUTATION ***********
   TIME = 0.D0
   ITIME = 0
   
   nDummy = 1
   do 1111 L = 1, ILINE
      call OUTPUT(ITIME, L, 0, nDummy)
1111 end do

   ! NTIME sets of constant wave and water level at the seaward boundary x = 0 for all ILINE cross-shore lines
   do 999 ITIME = 1, NTIME
      do 998 L = 1, ILINE
      
         ! IEND=0 during the constant wave and water level
         ! IEND= 1 at the end of each ITIME
         ! VY(J,L)=total longshore sediment transport rate integrated from
         ! If IPOND= 1, QO=wave overtopping rate at ridge crest and QM=wave overtopping rate at landward end node JMAX            
         QO(L) = 0.D0
         QWX = 0.D0
         
         if (IPOND == 1) QM = 0.D0
         
888      IEND = 0

         ! ..... PREPARE FOR LANDWARD MARCHING COMPUTATION
         ! SWLDEP(J,L) = still water depth at node J for the present landward marching computation along cross-shore line L            
         ICHECK = 0
         do 90 J = 1, JMAX(L)
            SWLDEP(J, L) = SWLBC(ITIME) - ZB(J, L)
            if (ICHECK == 0) then
               if (SWLDEP(J,L) < 0.D0) THEN
                  JSWL(L) = J
                  ICHECK = 1
               endif
            endif
90       end do

         if (ICHECK == 0) JSWL(L) = JMAX(L)
         
         ! If ITIDE= 1 and ILAB=0, compute cross-shore tidal water flux QTIDE at wet node J
         
         if (ITIDE == 1) then
            do 91 J = 1, JMAX(L)
               SMDEDY(J) = DETADY(ITIME)
               SMDEDY(J) = SMDEDY(J) * (0.5D0 + 0.5D0 * DTANH((XB(J) - 6.D0) / 1.D0))
               
               ! where the above transition function is specifICALLy for LSTF pumping system
               
91          end do

            if (ILAB == 0) then
               do 92 J = 1, JMAX(L)
                  if (J < JSWL(L)) then
                     QTIDE(J) = (XB(JSWL(L)) - XB(J)) * DSWLDT(ITIME)
                  else
                     QTIDE(J) = 0.D0
                  endif
92             end do

               QWX = QTIDE(1)
               
            endif
         endif

         ! If IWTRAN= 1 and JSWL(L) is less than JMAX(L), a bay exists landward of the emerged structure or dune crest. The landward still water level is SWLAND(ITIME) for nodes J=JSL,(JSL+1),...,JMAX. Integer LANCOM=0 implies no computation for transmitted waves in the bay. If this computation is needed, LANCOM= 1 is set after QO ITERATION            
         LANCOM = 0
         if (IWTRAN == 1 .and. JSWL(L) < JMAX(L)) then
            ICHECK = 0
            
            do 95 J = (JSWL(L) + 1), JMAX(L)
               DUM = SWLAND(ITIME) - ZB(J, L)
               if (DUM > 0.D0) then
                  SWLDEP(J, L) = DUM
                  if (ICHECK == 0) then
                     JSL = J
                     JSL1 = JSL - 1
                     ICHECK = 1
                  endif
               endif
95          end do
         else
            JSL = JMAX(L)
            JSL1 = JMAX(L) - 1
         endif

         ! If IPOND= 1, Subroutine20 PONDED finds ponded water zone
         
         if (IPOND == 1) then
            call PONDED(L)
         endif

         !.....ITERATION TO FIND QO(L) ............................................
         ! At beginning of each ITIME, QO(L) = 0.0 as specified above. During each ITIME for profile evolution computation with IPROFL = 1, converged QO(L) is used as an initial quess for the next profile change computation with ITEMAX = 4. If IOVER = 0, QO(L) = 0.0 always and no iteration.
         
         if (IOVER == 0) then
            ITEMAX = 1
         else
            ITEMAX = 20
            
            ! Computed wave overtopping rates QO(L) for ITEMAX= 10-30 changed very little for fixed coastal structures with IPROFL = 0
            
            if (IPROFL == 1) ITEMAX = 4
            
            ! Computed overwashed dune profile evolutions changed very little for ITEMAX == 3-4.
            
         endif

         ITEQO = 0
777      ITEQO = ITEQO + 1

         SIGMA(1) = HRMS(1) / SQR8
         H(1) = WSETUP(1) + SWLDEP(1,L)
         

         ! BDJ added on 2012-09-28
         if (H(1) <= 0) then
!             write (*,*) "CShore ERROR: ended with negative depth at the first node at time =", TIME

            NRET = -1

#if defined EXE
            stop 1
#else
            return
#endif
         endif
         ! end BDJ added on 2012-09-28

         SIGSTA(1) = SIGMA(1)/H(1)

         ! Subroutine LWAVE returns linear wave number WKP, phase velocity CP(J) ratio WN(J) of group velocity to phase velocity and sin STHETA(J) and cos CTHETA(J) of  wave angle for given QDISP = water flux in dispersion relation of linear waves. QDISP = 0.0 is assumed for J = 1 for simplicity
         
         QDISP = 0.D0
         call LWAVE(1, L, H(1), QDISP)

         ! Tentatively assume VMEAN(1) = 0.0         
         VMEAN(1) = 0.D0
         VSIGT = 0.D0
         
         ! At node J= 1, no porous layer         
         QWX = QO(L)
         if (IPERM == 1) then
            UPMEAN(1) = 0.D0
            QP(1) = 0.D0
            UPSTD(1) = 0.D0
            DPSTA(1) = 0.D0
         endif

         if (IROLL == 1) RQ(1) = 0.D0
         
         SIGMA2 =  SIGMA(1) ** 2.D0
         QWY = GRAV * SIGMA2 * STHETA(1) / CP(1)
         SXXSTA(1) = SIGMA2 * FSX
         EFSTA(1) = SIGMA2 * FE
         
         if (IANGLE == 1) SXYSTA(1) = SIGMA2 * FSY
         
         if (IWCINT == 1) then
            DUM = GRAV * H(1)
            SXXSTA(1) = SXXSTA(1) + QWX ** 2.D0 / DUM
            if (IANGLE == 1) SXYSTA(1) = SXYSTA(1) + QWX * QWY / DUM
         endif
         
         ! where roller volume flux RQ(1)=0 is assumed for SXXSTA(1), SXYSTA(1), USIGT=UMEAN/SIGT, and QWY        
         
         if (FB2(1, L) > 0.D0) then
            ! Bottom friction coefficient is positive,
            USIGT = -SIGSTA(1) * GRAV * H(1) / CP(1) / CP(1)
            
            if (IANGLE == 1) USIGT = USIGT * CTHETA(1)
            
            DUM = SIGSTA(1) * CP(1)
            ! bdj USIGT = USIGT + QWX / H(1) / DUM
            
            if (DUM > 1D-10) USIGT = USIGT + QWX / H(1) / DUM        ! BDJ
            
            ! Subroutine GBXAGF returns approximate analytical values for the Gbx and Gf factors used in calculating cross-shore bottom shear stress and energy dissipation.
            call GBXAGF(CTHETA(1), USIGT, STHETA(1), VSIGT, GBX(1), GF(1))
            
            TBXSTA(1) = FB2(1, L) * GBX(1) * DUM ** 2.D0 / GRAV
            DFSTA(1) = FB2(1,L) * GF(1) * DUM ** 3.D0 / GRAV
         else
            TBXSTA(1) = 0.D0
            DFSTA(1) = 0.D0
         endif

         ! Subroutine DBREAK computes the fraction of breaking waves and the associated wave energy dissipation and returns DBSTA(1)         
         call DBREAK(1, L, HRMS(1), H(1))

         ! ------------ LANDWARD MARCHING COMPUTATION -----------------------
         ! Computation marching landward from seaward boundary, J = 1. Compute unknown variables at node JP1 = (J+1) along line L
         J = 0
100      J = J + 1
         JP1 = J + 1
         ITE = 0

         ! DFM safety check ============================== NOTE could try goto 400 here
         if (J > 999) J = 999
         if (JP1 > 1000) JP1 = 1000
         ! DFM safety check END ==========================

         DUM = DFSTA(J) + DBSTA(J)
         
         if (IPERM == 1) DUM = DUM + DPSTA(J)
         
         DUM = DUM * WT(J)
         DUM = (EFSTA(J) - DX * DUM) / FE
         
         if (DUM <= 0.D0) then
            ! DUM (which is the square of sigma SIGTIE) is zero or negative
#if defined EXE
            write (*, 2902) JP1, L, TIME, DUM, ITEQO, ITE, QO(L)
2902        format('CShore WARNING 01: at end of landward marching computation, DUM (which is the square of sigma SIGTIE) <= 0 at node', I4, ' line', I3, ' time', F13.3, ', DUM =', F13.3, ' ITEQO =', I2, ' ITE =', I2, ' QO(L) =', F13.9)
#endif
            ! Set a warning flag
            NRET = 2

            ! Accept the computed results up to node JP1 - 1 and end landward marching computation
            JP1 = JP1 - 1

!           ***

            goto 400
         endif

         ! DFM moved, was at *** inside above if... endif ===========================
         ! BDJ added on 2012-09-28
         if (JP1 == 1 .and. EFSTA(1) > 1D-5) then
#if defined EXE
            write (*,*) 'CShore WARNING 02: large energy gradients at the first node at time =', TIME, ' (small waves with short period at sea boundary)'
#endif
            ! Set a warning flag
            NRET = 3

            ! STOP  %BDJ 2015-05-06
         endif

         if (JP1 == 1 .and. EFSTA(1) < 1D-5) then
#if defined EXE
            write (*,*) 'CShore WARNING 03: zero energy at the first node at time =', TIME
#endif

            ! Set a warning flag
            NRET = 4
         endif
         ! end BDJ added on 2012-09-28
         ! DFM moved END, was at *** inside above if... endif ===========================

         SIGITE = DSQRT(DUM)

         SXXSTA(JP1) = FSX * SIGITE ** 2.D0
         
         if (IROLL == 1) SXXSTA(JP1) = SXXSTA(JP1) + RX(J) * RQ(J)
         
         if (IWCINT == 1) SXXSTA(JP1) = SXXSTA(JP1) + QWX * QWX / GRAV / H(J)
         
         WSETUP(JP1) = WSETUP(J) - (SXXSTA(JP1) - SXXSTA(J) + (TBXSTA(J) - TWXSTA(ITIME)) * DX) / H(J)

         HITE = WSETUP(JP1) + SWLDEP(JP1, L)

         if (HITE < EPS1) then
            ! Water depth HITE is less than EPS1
#if defined EXE
            write (*, 2905) JP1, L, TIME, HITE, EPS1, ITEQO, QO(L)
2905        format('CShore WARNING 04: at end of landward marching computation, insufficient water depth HITE < EPS1 at node', I4, ' line', I3, ' time', F13.3, ', HITE =', F13.3, ' EPS1 =', F13.3, ' ITEQO =', I2, ' QO(L) =', F13.9)
#endif

            ! Set a warning flag 
            NRET = 5

            ! Accept the computed results up to node JP1 - 1 and end landward marching computation
            JP1 = JP1 -1
            
            goto 400
         endif

         QWX = QO(L)
         
         if (IPERM == 1) QWX = QO(L) - QP(J)
         
         if (ITIDE == 1 .and. ILAB == 0) QWX = QWX + QTIDE(JP1)
         
         if (IWCINT == 1) then
            if (IANGLE == 0) then
               QDISP = QWX
            else
               QWY = HITE * VMEAN(J) + GRAV * SIGITE ** 2.D0 * STHETA(J) / CP(J)
               
               if (IROLL == 1) QWY = QWY + RQ(J) * STHETA(J)
               
               QDISP = QWX * CTHETA(J) + QWY * STHETA(J)
            endif
         endif
         
         call LWAVE(JP1, L, HITE, QDISP)

         if (IANGLE == 1) then
            DUM1 = SIGITE ** 2.D0
            SXYSTA(JP1) = FSY * DUM1
            
            if (IROLL == 1) SXYSTA(JP1) = SXYSTA(JP1) + RY(J) * RQ(J)
            
            if (IWCINT == 1) SXYSTA(JP1) = SXYSTA(JP1) + QWX * QWY / GRAV / HITE
            
            DUM2 = SXYSTA(JP1) - SXYSTA(J)
            SIGN = STHETA(JP1) * DUM2
            
            if (SIGN > 0.D0) DUM2 = 0.D0
            
            TBYSTA(JP1) = -DUM2 / DX + TWYSTA(ITIME)
            
            if (ITIDE == 1) TBYSTA(JP1) = TBYSTA(JP1) - HITE * SMDEDY(JP1)
            
            DUM = SIGITE / HITE
            
            if (DUM > SISMAX) DUM = SISMAX
            
            DUM3 = CP(JP1)*CP(JP1)/GRAV

            GBY(JP1) = TBYSTA(JP1)/FB2(JP1,L)/DUM3/DUM/DUM
            
            ! Subroutine VSTGBY computes VSIGT for specified GBY, CTHETA, USIGT and STHETA
            ! DFM However in original code, USIGT is passed but not used
            USIGT = -CTHETA(J)*DUM*HITE/DUM3
            
            if (IROLL == 1) then
               USIGT = USIGT*(1.D0+ (CP(JP1)/GRAV)*RQ(J)/SIGITE**2.D0)
            endif
            
            SIGT = DUM*CP(JP1)
            
            if (IWCINT == 1) USIGT=USIGT+QWX/HITE/SIGT
            
            ! DFM In original code, parameter ISIGT is passed but not used
            ! call VSTGBY(CTHETA(J), USIGT, STHETA(J), VSIGT, GBY(JP1))
            call VSTGBY(CTHETA(J), STHETA(J), VSIGT, GBY(JP1))
            
            VITE = VSIGT*SIGT
         endif

         if (IROLL == 1) then
            RQITE = RQ(J) + DX*(DBSTA(J)-RBETA(J)*RQ(J))/RE(J)
            
            if (RQITE < 0.D0) RQITE=0.D0
         endif

         ! Begin iteration for improved Euler finite difference method
         do 200 ITE = 1, MAXITE
            HRMITE = SIGITE*SQR8

            call DBREAK(JP1, L, HRMITE, HITE)
            SIGSTA(JP1) = SIGITE / HITE
            if (SIGSTA(JP1) > SISMAX) SIGSTA(JP1) = SISMAX

            SIGT = CP(JP1)*SIGSTA(JP1)
            if (IANGLE == 0) then
               VSIGT = 0.D0
            else
               VSIGT = VITE/SIGT
            endif

            ! If IPERM= 1, Subroutine POFLOW computes porous flow variables.
            ! UPMEAN(J) = mean of horizontal discharge velocity UP
            ! UPSTD(J) = standard deviation of UP
            ! DPSTA(J) = energy dissipation rate of porous flow
            
            QWX=QO(L)
            
            if (IPERM == 1) then
               PKHSIG = WKP*HITE*SIGSTA(JP1)
               DEDX = (WSETUP(JP1) - WSETUP(J))/DX
               call POFLOW(JP1,L,PKHSIG,DEDX)
               QWX = QO(L) - QP(JP1)
            endif
            if (ITIDE == 1 .and. ILAB == 0) QWX=QWX+QTIDE(JP1)

            if (FB2(JP1,L) > 0.D0) then
               DUM = GRAV*HITE/CP(JP1)/CP(JP1)
               USIGT = -CTHETA(JP1)*SIGSTA(JP1)*DUM
               if (IROLL == 1) then
                  USIGT = USIGT*(1.D0+(CP(JP1)/GRAV)*RQITE/SIGITE**2.D0)
               endif
               
               if (IWCINT == 1) USIGT=USIGT+QWX/HITE/SIGT
               
               call GBXAGF(CTHETA(JP1),USIGT,STHETA(JP1),VSIGT, GBX(JP1), GF(JP1))
               
               TBXSTA(JP1) = FB2(JP1,L)*GBX(JP1)*SIGT**2.D0/GRAV
               DFSTA(JP1) = FB2(JP1,L)*GF(JP1)*SIGT**3.D0/GRAV
            else
               TBXSTA(JP1) = 0.D0
               DFSTA(JP1) = 0.D0
            endif

            DUMD = DFSTA(JP1) + DFSTA(J) + DBSTA(JP1) + DBSTA(J)
            
            if (IPERM == 1) DUMD = DPSTA(JP1) + DPSTA(J) + DUMD
            
            DUMD = DUMD * (WT(J) + WT(JP1)) / 2.D0
            DUM = (EFSTA(J) - DXD2 * DUMD) / FE
            
            if (DUM <= 0.D0) then
               ! DUM (which is the square of sigma SIGTIE) is zero or negative
#if defined EXE
               write (*, 2903) JP1, L, TIME, DUM, ITEQO, ITE, QO(L)
2903           format(/'CShore WARNING 05: at end of landward marching computation, DUM (which is the square of sigma SIGTIE) <= 0 at node', I4, ' line', I3, ' time', F13.3, ' DUM =', F13.3, ' ITEQO =', I2, ' ITE =', I2, ' QO(L) =', F13.9)
#endif

               ! Set a warning flag
               NRET = 2

               ! Accept the computed results up to node JP1-1 and end landward marching computation
               JP1 = JP1 - 1

               goto 400
            endif

            SIGMA(JP1) = DSQRT(DUM)            
            SXXSTA(JP1) = FSX*SIGMA(JP1)**2.D0
            
            if (IROLL == 1) SXXSTA(JP1)=SXXSTA(JP1)+RX(JP1)*RQITE
            
            if (IWCINT == 1) SXXSTA(JP1)=SXXSTA(JP1)+QWX*QWX/GRAV/HITE
            
            WSETUP(JP1) = WSETUP(J) - (2.D0* (SXXSTA(JP1)-SXXSTA(J)) + DX*(TBXSTA(JP1)+TBXSTA(J)-2.D0*TWXSTA(ITIME)))/ (HITE+H(J))

            H(JP1) = WSETUP(JP1) + SWLDEP(JP1,L)

            SIGSTA(JP1) = SIGMA(JP1)/H(JP1)
            
            if (SIGSTA(JP1) > SISMAX) SIGSTA(JP1)=SISMAX

            if (H(JP1) <= EPS1) then
               ! Water depth H(JP1) is less than EPS1
#if defined EXE
               write (*, 2904) JP1, L, TIME, H(JP1), EPS1, ITEQO, QO(L)
2904           format(/'CShore WARNING 06: at end of landward marching computation, insufficient water depth at node', I6, ' line', I3, ' time', F13.3, ', water depth H(JP1) =', F13.3, ' is less than EPS1 =', F13.3, ', ITEQO =', I2, ' QO(L) =', F13.9)
#endif
               ! Set a warning flag
               NRET = 5

               ! Accept the computed results up to node JP1-1 and end landward marching computation
               JP1 = JP1 - 1
               
               goto 400
            endif

            if (IWCINT == 1) then
               if (IANGLE == 0) then
                  QDISP = QWX
               else
                  QWY=H(JP1)*VITE+ GRAV*SIGMA(JP1)**2.D0*STHETA(JP1)/CP(JP1)
                  
                  if (IROLL == 1) QWY=QWY+RQITE*STHETA(JP1)
                  
                  QDISP = QWX*CTHETA(JP1) + QWY*STHETA(JP1)
               endif
            endif
            
            call LWAVE(JP1,L,H(JP1),QDISP)

            if (IANGLE == 1) then
               DUM1 = SIGMA(JP1)**2.D0
               SXYSTA(JP1) = FSY*DUM1
               if (IROLL == 1) SXYSTA(JP1)=SXYSTA(JP1)+RY(JP1)*RQITE
               if (IWCINT == 1) SXYSTA(JP1)=SXYSTA(JP1)+QWX*QWY/GRAV/ H(JP1)
               DUM2 = SXYSTA(JP1) - SXYSTA(J)
               SIGN=STHETA(JP1)*DUM2
               if (SIGN > 0.D0) DUM2=0.D0
               TBYSTA(JP1) = -DUM2/DX + TWYSTA(ITIME)
               if (ITIDE == 1) TBYSTA(JP1)=TBYSTA(JP1)-H(JP1)*SMDEDY(JP1)
               DUM3 = CP(JP1)*CP(JP1)/GRAV
               GBY(JP1)=TBYSTA(JP1)/FB2(JP1,L)/DUM3/ SIGSTA(JP1)/SIGSTA(JP1)
               USIGT = -CTHETA(JP1)*SIGSTA(JP1)*H(JP1)/DUM3
               if (IROLL == 1) then
                  USIGT=USIGT*(1.D0+(CP(JP1)/GRAV)*RQITE/SIGMA(JP1)**2.D0)
               endif
               SIGT = SIGSTA(JP1)*CP(JP1)
               if (IWCINT == 1) USIGT=USIGT+QWX/H(JP1)/SIGT
               
               ! DFM In original code, parameter ISIGT is passed but not used
               ! call VSTGBY(CTHETA(JP1), USIGT, STHETA(JP1), VSIGT, GBY(JP1))
               call VSTGBY(CTHETA(JP1), STHETA(JP1), VSIGT, GBY(JP1))

               VMEAN(JP1) = VSIGT*SIGT
            endif

            if (IROLL == 1) then
               DUM1 = RE(JP1) + DXD2*RBETA(JP1)
               DUM2 = (RE(J) - DXD2*RBETA(J))*RQ(J) + DXD2*(DBSTA(JP1) + DBSTA(J))
               RQ(JP1) = DUM2/DUM1
            endif

            ! Check for convergence
            ESIGMA = DABS(SIGMA(JP1) - SIGITE)
            EH = DABS(H(JP1) - HITE)
            if (IANGLE == 1) EV = DABS(VMEAN(JP1) - VITE)
            if (IROLL == 1) ERQ = DABS(RQ(JP1) - RQITE)
            if (ESIGMA < EPS1 .and. EH < EPS1) then
               if (IANGLE == 0) then
                  goto 199
               else
                  if (EV < EPS1) goto 199
                  goto 198
               endif
199            if (IROLL == 0) then
                  goto 210
               else
                  if (ERQ < EPS2) goto 210
                  goto 198
               endif
            endif

            ! Average new and previous values to accelerate convergence
198         SIGITE =  0.5D0 * (SIGMA(JP1) + SIGITE)
            HITE = 0.5D0 * (H(JP1) + HITE)
            if (IANGLE == 1) VITE = 0.5D0 * (VMEAN(JP1) + VITE)
            if (IROLL == 1) RQITE = 0.5D0 * (RQ(JP1)+RQITE)
200      end do
         ! End of iteration: do 200 ITE = 1 to MAXITE

         ! The iteration did not converge
#if defined EXE
         write (*, 29050) MAXITE, EPS1, JP1, L, TIME
29050    format('CShore WARNING 07: did not reach convergence after MAXITE =', I4, ' iterations with relative error EPS1 =', E17.5, ' at node JP1 =', I4, ' Line L =', I3, ' TIME =', F13.3)
#endif

         ! Set a warning flag 
         NRET = 1

         ! And adopt the last iteration values
         SIGMA(JP1) = SIGITE
         H(JP1) = HITE
         if (IANGLE == 1) VMEAN(JP1) = VITE
         if (IROLL == 1) RQ(JP1) = RQITE

         ! Wave transmission computation for LANCOM == 1 starts from here
210      HRMS(JP1) = SQR8*SIGMA(JP1)

         WSETUP(JP1) = H(JP1) - SWLDEP(JP1,L)

         if (IWCINT == 1) then
            if (IANGLE == 0 .OR. LANCOM == 1) then
               QDISP = QWX
            else
               QWY=H(JP1)*VMEAN(JP1) + GRAV*SIGMA(JP1)**2.D0*STHETA(JP1) / CP(JP1)
               if (IROLL == 1) QWY=QWY + RQ(JP1)*STHETA(JP1)
               QDISP = QWX*CTHETA(JP1) + QWY*STHETA(JP1)
            endif
         endif

         call LWAVE(JP1, L, H(JP1), QDISP)
         call DBREAK(JP1, L, HRMS(JP1), H(JP1))
         SIGSTA(JP1) = SIGMA(JP1)/H(JP1)
         if (SIGSTA(JP1) > SISMAX) SIGSTA(JP1) = SISMAX
         SIGT = SIGSTA(JP1)*CP(JP1)
         if (IANGLE == 0) then
            VSIGT = 0.D0
         else
            VSIGT = VMEAN(JP1)/SIGT
         endif

         QWX=QO(L)
         if (IPERM == 1) then
            PKHSIG = WKP*H(JP1)*SIGSTA(JP1)
            DEDX = (WSETUP(JP1) - WSETUP(J))/DX
            call POFLOW(JP1,L,PKHSIG,DEDX)
            QWX = QO(L) - QP(JP1)
         endif
         if (ITIDE == 1 .and. ILAB == 0) QWX=QWX+QTIDE(JP1)

         SIGMA2 = SIGMA(JP1)**2.D0
         SXXSTA(JP1) = SIGMA2*FSX
         if (IROLL == 1) SXXSTA(JP1)=SXXSTA(JP1)+RX(JP1)*RQ(JP1)
         if (IWCINT == 1) SXXSTA(JP1)=SXXSTA(JP1)+QWX*QWX/GRAV/H(JP1)
         EFSTA(JP1) = SIGMA2*FE
         DUM3 = CP(JP1)*CP(JP1)/GRAV
         if (FB2(JP1,L) > 0.D0) then
            USIGT = -CTHETA(JP1)*SIGSTA(JP1)*H(JP1)/DUM3
            if (IROLL == 1) then
               USIGT = USIGT*(1.D0+(CP(JP1)/GRAV)*RQ(JP1)/SIGMA2)
            endif
            if (IWCINT == 1) USIGT=USIGT+QWX/H(JP1)/SIGT
            
            call GBXAGF(CTHETA(JP1), USIGT, STHETA(JP1), VSIGT, GBX(JP1), GF(JP1))
            
            TBXSTA(JP1)=FB2(JP1,L)*GBX(JP1)*SIGT**2.D0/GRAV
            DFSTA(JP1)=FB2(JP1,L)*GF(JP1)*SIGT**3.D0/GRAV
         else
            TBXSTA(JP1) = 0.D0
            DFSTA(JP1) = 0.D0
         endif

         if (IANGLE == 1) then
            SXYSTA(JP1) = FSY*SIGMA2
            if (IROLL == 1) SXYSTA(JP1)=SXYSTA(JP1)+RY(JP1)*RQ(JP1)
            if (IWCINT == 1) SXYSTA(JP1)=SXYSTA(JP1)+QWX*QWY/GRAV/H(JP1)
            DUM2 = SXYSTA(JP1) - SXYSTA(J)
            SIGN=STHETA(JP1)*DUM2
            if (SIGN > 0.D0) DUM2=0.D0
            TBYSTA(JP1) = -DUM2/DX + TWYSTA(ITIME)
            if (ITIDE == 1) TBYSTA(JP1)=TBYSTA(JP1)-H(JP1)*SMDEDY(JP1)
            if (J == 1) then
               TBYSTA(J) = TBYSTA(JP1)
               VMEAN(J) = VMEAN(JP1)
            endif
            GBY(JP1)=TBYSTA(JP1)/FB2(JP1,L)/DUM3/SIGSTA(JP1)/SIGSTA(JP1)
         endif

         JDUM = JMAX(L)
         if (RCREST(L) > SWLBC(ITIME) .and. LANCOM == 0) JDUM = JCREST(L)
         if (H(JP1) < EPS1 .OR. JP1 == JDUM) goto 400
         if (JP1 == JMAX(L) .and. LANCOM == 1) goto 406

         goto 100
         !----------------End of LANDWARD MARCHING COMPUTATION -------------
         
400      continue

         JR = JP1
                     
!         ! BDJ added on 2012-09-27
!         if (JR.EQ.1) then
!            write (*,*) 'CShore ERROR: ended landward marching computation at the first node at time = ', TIME
!            NRET = -2
!            
!#if defined EXE
!            stop 1
!#else
!            return
!#endif
!         endif
!         ! end BDJ added on 2012-09-27

         XR = XB(JR)
         ZR = ZB(JR,L)

         ! BDJ added on 2012-10-09
         call SRFSP(L)
         ! end BDJ added on 2012-10-09


         ! If IOVER= 1, Subroutine10 QORATE computes for cross-shore line L
         ! QO(L) = sum of wave overtopping, overflow and seepage rates
         ! If IOVER=0, QO(L)=0.0, no iteration and no wet/dry zone
         ! bdj 2012-08-20 added HRMS(1)>0 switch to avoid odd behavior
         !      if (IOVER.EQ.1) then

         if (IOVER == 1 .and. HRMS(1) > 1D-10) then ! BDJ
            ICONV = 1
            QOUSED = QO(L)
            
            call QORATE(ITIME, L, ITEQO, ICONV, 0)
            
            if (ICONV == 0) goto 405
            
#if defined FILEINOUT || ARGINBOTHOUT
            if (ICONV == 1) write (40, 2906) JR, L, TIME, ITEQO, QOUSED, QO(L)
2906        format(/'NO CONVERGENCE OF QO ITERATION' / 'Landward end node JR =', I6, 'Line =', I3, 'TIME =', F13.3 / 'Iteration number ITEQO =', I3, 'assumed QO =', F13.9 /           'computed QO =', F13.9)
#endif
         else
            QO(L) = 0.D0
            JWD = JR
            JDRY = JR
            goto 405
         endif

         if (ITEQO < ITEMAX) goto 777

         !....................END OF QO ITERATION...........................
405      continue

         ! If IWTRAN= 1 and JDRY=JSL1=(JSL-1), compute transmitted waves for nodes J=JSL,(JSL+1),...,JMAX after QO ITERATION is completed because QO is not affected by transmitted waves. This landward marching computation is indicated by LANCOM= 1. Assumed boundary conditions at JSL are specified below where the wet probability PWET changes suddenly from PWET(JDRY) to PWET(JSL)= 1. Wave setup and standard deviation for the entire duration are matched.

         if (IWTRAN == 1 .and. JDRY == JSL1) then
            LANCOM = 1
            ! Wave setup is defined above the landward still water level SWLAND(ITIME)
            WSETUP(JDRY) = H(JDRY) * PWET(JDRY) + ZB(JDRY,L) - SWLAND(ITIME)
            WSETUP(JSL) = WSETUP(JDRY) * PWET(JDRY)
            H(JSL) = WSETUP(JSL) + SWLDEP(JSL,L)
            SIGMA(JSL) = SIGMA(JDRY) * PWET(JDRY)
            RQ(JSL) = 0.D0
            VMEAN(JSL) = 0.D0
            SXYSTA(JDRY) = 0.D0
            QWX = QO(L)
            QWY = 0.D0
            JP1 = JSL

            goto 210
         endif
         ! goto 406 after LANCOM = 1 computation

406      continue

         ! Calculate the standard deviation and mean of the horizontal velocities U and V
         LMAX=LANCOM+1
         do 411 K= 1,LMAX
            if (K == 1) then
               ISTART= 1
               IFIN=JR
            else
               ISTART=JSL
               IFIN=JMAX(L)
            endif
            do 410 I =ISTART,IFIN
               SIGT = CP(I)*SIGSTA(I)

               ! DFM safety check ==============================
               if (SIGT == 0.D0) SIGT = 1.0D-9
               ! DFM safety check END ==========================

               USTD(I) = SIGT*CTHETA(I)
               UMEAN(I)= -USTD(I)*SIGSTA(I)*GRAV*H(I)/CP(I)/CP(I)
               
               if (IROLL == 1) UMEAN(I)=UMEAN(I)*(1.D0+(CP(I)/GRAV) * RQ(I)/SIGMA(I)**2.D0)
               QWX = QO(L)
               
               if (IPERM == 1) QWX=QO(L)-HP(I,L)*UPMEAN(I)
               if (ITIDE == 1 .and. ILAB == 0) QWX=QWX+QTIDE(JP1)
               
               UMEAN(I) = UMEAN(I) + QWX/H(I)
               ! bdj USTA(I)=UMEAN(I)/SIGT
               
               if (SIGT > 1D-10) then ! BDJ
                  USTA(I)=UMEAN(I)/SIGT ! BDJ
               else ! BDJ
                  USTA(I) = 0.0D0   ! BDJ
               endif ! BDJ
               if (IANGLE == 1) then
                  VSTD(I) = SIGT*DABS(STHETA(I))
                  VSTA(I) = VMEAN(I)/SIGT
               endif
410         end do
411      end do

         ! If IOVER= 1, connect H(J) and UMEAN(J) with J= 1 to JR with wet/dry-zone HWD(J) and UMEAWD(J) with J=JWD to JDRY using Subroutine17 TRANWD also connect the corresponding standard deviations
         if (IOVER == 1) then
            PWET(1:JWD)= 1.D0
            if (IWTRAN == 1 .and. JDRY == JSL1) PWET(JSL:JMAX(L))= 1.D0
            if (JDRY > JR) then
               call TRANWD(H,JR,HWD,JWD,JDRY)
               call TRANWD(SIGMA,JR,SIGWD,JWD,JDRY)
               call TRANWD(UMEAN,JR,UMEAWD,JWD,JDRY)
               call TRANWD(USTD,JR,USTDWD,JWD,JDRY)
               if (IPERM == 1) call TRANWD(UPMEAN,JR,UPMWD,JWD,JDRY)
               if (IANGLE == 1) then
                  call TRANWD(VMEAN,JR,VMEAWD,JWD,JDRY)
                  call TRANWD(VSTD,JR,VSTDWD,JWD,JDRY)
               endif
            else
               JDRY=JR
               if (JWD < JR) then
                  JDUM=JWD+1
                  PWET(JDUM:JR)= 1.D0
               endif
            endif
         endif
         
         ! Smooth computed H(J), USTD(J), UMEAN(J), USTA(J), DFSTA(J), DBSTA(J),RQ(J), VMEAN(J), VSTD(J) and VSTA(J) using Subroutine 14 SMOOTH
         JDUM=JDRY
         if (IWTRAN == 1 .and. JDRY == JSL1) JDUM=JMAX(L)
         DUMVEC = H
         call SMOOTH(JDUM,DUMVEC,H)
         DUMVEC = SIGMA
         call SMOOTH(JDUM,DUMVEC,SIGMA)
         DUMVEC = USTD
         call SMOOTH(JDUM,DUMVEC,USTD)
         DUMVEC = UMEAN
         call SMOOTH(JDUM,DUMVEC,UMEAN)
         DUMVEC = USTA
         call SMOOTH(JR,DUMVEC,USTA)
         DUMVEC = DFSTA
         call SMOOTH(JR,DUMVEC,DFSTA)
         if (IPERM == 1) then
            DUMVEC=UPMEAN
            call SMOOTH(JDUM,DUMVEC,UPMEAN)
            if (IOVER == 1) then
               do 420 J= 2,JDUM
                  DUM=ZP(J,L)
                  if (DUM < SWLBC(ITIME) .and. ZP(J,L) >= ZP(J-1,L)) DUM=SWLBC(ITIME)
                  ETAPOR=ZB(J,L)*PWET(J)+DUM*(1.D0-PWET(J))
                  QP(J)=UPMEAN(J)*(ETAPOR-ZP(J,L))*PWET(J)
420            end do
            endif
         endif
         if (IROLL == 0) then
            DUMVEC=DBSTA
            call SMOOTH(JR,DUMVEC,DBSTA)
         else
            DUMVEC = RQ
            call SMOOTH(JR,DUMVEC,RQ)
         endif
         if (IANGLE == 1) then
            DUMVEC = VMEAN
            call SMOOTH(JDRY,DUMVEC,VMEAN)
            DUMVEC = VSTD
            call SMOOTH(JDRY, DUMVEC, VSTD)
            DUMVEC = VSTA
            call SMOOTH(JR,DUMVEC,VSTA)
         endif

         ! Subroutine11 SEDTRA computes the cross-shore and longshore sediment transport rates if IPROFL = 1. If a vertical wall exists, IVWALL= 2 indicates exposure to wave action along cross-shore line L
         if (IPROFL == 1) then
            if (IVWALL(L) >= 1) then
               if (ZB(JMAX(L)-1,L) < ZP(JMAX(L),L)) then
                  IVWALL(L)= 2
               else
                  IVWALL(L)= 1
               endif
            endif
            call SEDTRA(L)
         endif

         ! If IPROFL= 1, Subroutine12 CHANGE computes the change of the bottom elevation, DELZB(j), at node j during the time step DELT determined in this subroutine which also checks whether IEND= 1 and the end of given ITIME is reached. VY(j,L)=total longshore sediment transport rate integrated from TIMEBC(ITIME) to TIMEBC(ITIME+1) used in Subroutine12 CHANGE if IQYDY= 1
         if (IPROFL == 1) then
            call CHANGE(ITIME, L, IEND, 1)
            do 430 J= 1,JMAX(L)
               if (TIME == TIMEBC(ITIME)) then
                  VY(J,L)=0.D0
                  DZX(J,L)=ZB(J,L)
               endif
               ZB(J,L)=ZB(J,L)+DELZB(J,L)
               if (TIME == 0.D0) then
                  VBX(J,L)=0.D0
                  VSX(J,L)=0.D0
                  VBY(J,L)=0.D0
                  VSY(J,L)=0.D0
               endif
               VBX(J,L)=VBX(J,L)+DELT*QBX(J)
               VSX(J,L)=VSX(J,L)+DELT*QSX(J)
               if (IANGLE == 1) then
                  VBY(J,L)=VBY(J,L)+DELT*QBY(J)
                  VSY(J,L)=VSY(J,L)+DELT*QSY(J)
                  VY(J,L)=VY(J,L)+DELT*(QBY(J)+QSY(J))
               endif
430         end do
            if (IVWALL(L) >= 1) then
               if (ZB(JMAX(L),L) < ZP(JMAX(L),L)) then
                  ZB(JMAX(L),L)=ZP(JMAX(L),L)
               endif
            endif
         else
            IEND= 1
            DELT = 1.D0
         endif
         if (IQYDY == 1 .and. IEND == 1) then
            if (L == ILINE) call CHANGE(ITIME, L, IEND, 2)
         endif

         ! If IOVER = 1, store time series of wave overtopping rate and sediment transport rates at landward end node JMAX QO(L)=sum of wave overtopping, overflow and seepage rates
         if (IOVER == 1) then
            if (TIME == TIMEBC(ITIME)) then
               TSQO(L) = 0.D0
               TSQBX(L) = 0.D0
               TSQSX(L) = 0.D0
            endif
            if (IPOND == 1 .and. NOPOND == 0) then
               TSQO(L)=TSQO(L)+DELT*QM
            else
               TSQO(L)=TSQO(L)+DELT*QO(L)
            endif
            TSQBX(L)=TSQBX(L)+DELT*QBX(JMAX(L))
            TSQSX(L)=TSQSX(L)+DELT*QSX(JMAX(L))
         endif

         ! Subroutine OUTPUT stores input and computed results when IEND= 1.
         if (IPROFL == 1) TIME=TIME+DELT
         if (IEND == 1) then
            call OUTPUT(ITIME, L, ITEQO, ICONV)
         endif
         if (IPROFL == 0) TIME=TIMEBC(ITIME+1)

         ! Compute the BSLOPE at time = (TIME+DELT)
         if (IPROFL == 1) then
            do 501 J= 1,JMAX(L)
               if (J == 1) then
                  BSLOPE(1,L)=(ZB(2,L)-ZB(1,L))/DX
               else
                  if (J == JMAX(L)) then
                     BSLOPE(JMAX(L),L)=(ZB(JMAX(L),L)-ZB(JMAX(L)-1,L))/DX
                  else
                     BSLOPE(J,L)=(ZB(J+1,L)-ZB(J-1,L))/DX2
                  endif
               endif
501         end do

            ! Compute new bottom RCREST if IOVER= 1 and IPOND=0
            if (IOVER == 1 .and. IPOND == 0) then
               RCREST(L) = ZB(1,L)
               do 502 J= 2,JMAX(L)
                  DUM = ZB(J,L) - RCREST(L)
                  if (DUM >= 0.D0) then
                     RCREST(L) = ZB(J,L)
                     JCREST(L) = J
                  endif
502            end do
            endif
            
            ! Compute new porous layer thickness HP if IPERM= 1 or if ISEDAV= 1
            if (IPERM == 1 .OR. ISEDAV == 1) then
               do 503 J= 1,JMAX(L)
                  HP(J,L) = ZB(J,L) - ZP(J,L)
                  if (HP(J,L) < 0.D0) HP(J,L)=0.D0
503            end do
            endif
         endif

         ! If IEND = 0, go to 888 for the next landward marching computation
         if (IEND == 0) goto 888

         ! If IEND= 1 and L is less than ILINE, reset time for next cross-shore line
         if (L < ILINE) then
            TIME=TIMEBC(ITIME)
         endif

         ! IEND= 1 and specify the seaward input for the next value of ITIME if ITIME is less than NTIME and L=ILINE.
         if (ITIME < NTIME .and. L == ILINE) then
            ITIME1=ITIME+1
            TP=TPBC(ITIME1)
            HRMS(1)=HRMSBC(ITIME1)
            
            ! NPT=integer used in Subroutine14 SMOOTH
            ! NPE=integer used in Subroutine15 EXTRAPO
            NPT= 1+NINT(HRMS(1)/DX)
            NPE= 1+NINT(HRMS(1)/DX2)
            if (IPROFL == 1 .and. IPERM == 1) then
               NPT=NPT+NINT(HRMS(1)/DXD2)
               if (HRMS(1) < 0.05D0) NPT = NPT + NINT(HRMS(1) / DX)
            endif
            
            ! if (IVWALL(L).EQ.2) NPT=NPT+NINT(HRMS(1)/DX)
            WSETUP(1)=WSETBC(ITIME1)
            ANGLE=WANGBC(ITIME1)
            IANGLE = 1
            if (ANGLE == 0.D0) IANGLE=0
            if (IANGLE == 1) IWCINT=0
            
            ! No wave and current interaction for IANGLE= 1
            WKPO=TWOPI*TWOPI/(GRAV*TP*TP)
            ! where WKPO is the deep water wave number
         endif

998   end do
   ! **************** END OF ILINE COMPUTATION ***************************

999 end do

! **************** END OF TIME MARCHING COMPUTATION ********************

#if defined ARGINOUT 
   nOutSize = JR        ! Is nExpectedRows in CoastalME
#endif

#if defined FILEINOUT || ARGINBOTHOUT   
   do i = 20, 40
      if (i == 39) cycle
      
      ! This is necessary to get CShore to finish writing all output files before returning to the calling program
      flush (i)
      close (i)
   end do
#endif

#if defined EXE
   stop 0
#else
   return
#endif

#if defined EXE
end program CShore
#else
end subroutine CShore
#endif
