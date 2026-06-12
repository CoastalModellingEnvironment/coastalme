#if defined ARGINOUT || ARGINBOTHOUT
!===============================================================================================================================
!
! This is a wrapper for CShore which gets input from arguments which are supplied by a calling program (e.g. CoastalME), and puts output into arguments which are passed back to the calling program
!
!===============================================================================================================================
subroutine CShoreWrapper(In_ILINE, In_IPROFL, In_IPERM, In_IOVER, In_IWCINT, In_IROLL, In_IWIND, In_ITIDE, In_ILAB, In_NWAVE, In_NSURG, In_DX, In_GAMMA, In_TWAVE, In_TPIN, In_HRMSIN, In_WANGIN, In_TSURG, In_SWLIN, In_NBINP, In_XBINP, In_ZBINP, In_FBINP, Out_IError, Out_nOutSize, Out_XYDist, Out_FreeSurfaceStd, Out_WaveSetupSurge, Out_SinWaveAngleRadians, Out_FractionBreakingWaves) bind(c, name = "CShoreWrapper")

   use CShoreShared
   
   ! DFM TODO 070 Make these dynamically allocated. NN increased to 1000 for CoastalME 1.3.3
   integer, parameter :: NN = 1000, NL = 1

   integer(c_int), intent(in) :: In_ILINE       ! Number of cross-shore lines i.e. number of CoastalME profiles
   integer(c_int), intent(in) :: In_IPROFL      ! 0 for fixed bottom profile (assumes ISEDAV = 0), 1 for profile evolution computation (need to input ISEDAV)
                                                ! ISEDAV = 0 for unlimited bottom sediment, ISEDAV = 1 for sediment availability limited by hard bottom
   integer(c_int), intent(in) :: In_IPERM       ! 0 for impermeable bottom, 1 for permeable bottom of stone structure
   integer(c_int), intent(in) :: In_IOVER       ! 0 for no wave overtopping and overflow on crest (assumes IWTRAN = 0), 1 for wave overtopping and overflow (need to input IWTRAN)
                                                ! IWTRAN = 0 for no standing water landward of crest, IWTRAN = 1 for wave transmission into bay or estuary
                                                ! IPOND = 0 for no ponding water landward of JSWL node, IPOND = 1 for ponding water in runnel landward of ridge for IPERM = 0
                                                ! INFILT = 0 for no infiltration landward of sand dune, INFILT = 1 for infiltration if IOVER = 1, IPERM = 0, IPROFL = 1

   integer(c_int), intent(in) :: In_IWCINT      ! 0 for no wave/current interaction, 1 for wave/current interaction in frequency dispersion, momentum and wave action equations
   integer(c_int), intent(in) :: In_IROLL       ! 0 for no roller effects, 1 for roller effects in governing equations
   integer(c_int), intent(in) :: In_IWIND       ! 0 for no wind effects, 1 for wind shear stresses on momentum equations
   integer(c_int), intent(in) :: In_ITIDE       ! 0 for no tidal effect on currents, 1 for longshore and cross-shore tidal currents

   integer(c_int), intent(in) :: In_ILAB        ! 0 for field data set (waves and water levels read separately), 1 for laboratory data set (waves and water levels read together)
   integer(c_int), intent(in) :: In_NWAVE       ! Number of waves at x = 0 starting from time = 0
   integer(c_int), intent(in) :: In_NSURG       ! Number of water levels at x = 0 from time = 0
   
   real(c_double), intent(in) :: In_DX          ! Nodal spacing for input bottom geometry
   real(c_double), intent(in) :: In_GAMMA       ! Breaking wave height-depth ratio (0.5 - 1.0)
   
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_TWAVE  ! Timestep
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_TPIN   ! Wave period
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_HRMSIN ! Hrms
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_WANGIN ! Wave angle
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_TSURG  ! Surge time   
   real(c_double), intent(in), dimension(In_NWAVE+1) :: In_SWLIN  ! Surge level   

   ! The bottom geometry is divided into segments of different inclination and roughness starting from the seaward boundary for cross-shore line L   
   integer(c_int), intent(in) :: In_NBINP        ! Number of input points of bottom elevation ZB(X) for this profile
   
   real(c_double), intent(in), dimension(In_NBINP, In_ILINE) :: In_XBINP    ! Horizontal distance to input bottom point (J) in metres from the seaward boundary
   real(c_double), intent(in), dimension(In_NBINP, In_ILINE) :: In_ZBINP    ! Dimensional vertical coordinate (+ above datum) of input bottom point (J) in metres
   real(c_double), intent(in), dimension(In_NBINP, In_ILINE) :: In_FBINP    ! Bottom friction factor for segment between points (J) and (J+1)
   
   integer(c_int), intent(out) :: Out_IError               ! Error code to be returned to CoastalME
   integer(c_int), intent(out) :: Out_nOutSize             ! Number of valid elements in Out_XYDist, Out_FreeSurfaceStd, Out_SinWaveAngleRadians, Out_FractionBreakingWaves
   
   real(c_double), intent(out), dimension(NN, NL) :: Out_XYDist
   real(c_double), intent(out), dimension(NN, NL) :: Out_FreeSurfaceStd
   real(c_double), intent(out), dimension(NN, NL) :: Out_WaveSetupSurge
   ! real(c_double), intent(out), dimension(NN, NL) :: Out_StormSurge
   real(c_double), intent(out), dimension(NN, NL) :: Out_SinWaveAngleRadians
   real(c_double), intent(out), dimension(NN, NL) :: Out_FractionBreakingWaves  
 
   
   integer :: I, IDUM, J, K, L, NSLAN, NSLAN1, NSURG1, NTIDE, NTIDE1, NWAVE1, NWIND1
   double precision :: ANG, CONVRT, DUM, RATIO, SGM1, SPORO
   
   interface
      subroutine TSINTP(N1, T1, W1, N2, T2, W2)
         integer, intent(in) :: N1, N2
         double precision, intent(in), dimension(:) :: T1, W1, T2
         double precision, intent(out), dimension(:) :: W2
      end subroutine TSINTP
      
      subroutine CShore(NRET)
         integer, intent(inout) :: NRET
      end subroutine CShore
   end interface

   ! Start of executable section
   ! write (*, *) "In CShoreWrapper"
   
   ! Copy the input values to variables in the common module
   ILINE = In_ILINE
   IPROFL = In_IPROFL
   IPERM = In_IPERM
   IOVER = In_IOVER
   IWCINT = In_IWCINT
   IROLL = In_IROLL
   IWIND = In_IWIND
   ITIDE = In_ITIDE
   ILAB = In_ILAB
   NWAVE = In_NWAVE
   NSURG = In_NSURG
   
   DX = In_DX
   GAMMA = In_GAMMA
   
   ! Copy these input arrays to arrays in the common module: doing this also allocates the target arrays
   TWAVE = In_TWAVE
   TPIN = In_TPIN
   HRMSIN = In_HRMSIN
   WANGIN = In_WANGIN
   TSURG = In_TSURG
   SWLIN = In_SWLIN
   
   ! Allocate memory for the cross-shore lines i.e. for the CoastalME profiles 
   call allocate_cross_shore_lines_size_arrays(ILINE)
   
   NBINP(1) = In_NBINP     ! NOTE Will need to change this with multiple profiles
   
   ! We now need to copy these input arrays to arrays in the common module. But can't do a simple array assign, since XBINP, ZBINP and FBINP must be larger than In_XBINP in the first dimension, because element K+1 (where K is the size of the first dimension) gets accessed in various places e.g. BOTTOM(). This is a bodge
   ndim1size = size(In_XBINP, 1)
   ndim2size = size(In_XBINP, 2)

   allocate(XBINP(ndim1size+1, ndim2size), ZBINP(ndim1size+1, ndim2size), FBINP(ndim1size+1, ndim2size))
   do ndim2 = 1, ndim2size
      XBINP(1:ndim1size, ndim2) = In_XBINP(1:ndim1size, ndim2)
      ZBINP(1:ndim1size, ndim2) = In_ZBINP(1:ndim1size, ndim2)
      FBINP(1:ndim1size, ndim2) = In_FBINP(1:ndim1size, ndim2)
   end do
   
!   XBINP = In_XBINP
!   ZBINP = In_ZBINP
!   FBINP = In_FBINP 
      
   ! NOTE: We may wish to try reading in data for all profiles simultaneously. At present, ILINE always == 1   
   !if (ILINE > 2) then
   !   read (11,*) DUM
   !   IQYDY = NINT(DUM)
   ! endif
   
   ! Allocate memory for the cross-shore nodes arrays NOTE Will need to change this with multiple profiles TODO make bounds dynamic?
   call allocate_cross_shore_nodes_size_arrays(NN)
   
   ! Allocate memory for output arrays TODO can this be dynamic?  FIX THIS why is this called more than once per CShore run?
   call allocate_argument_output_arrays(NN, NL)
   
   ! IPROFL == 0 for fixed bottom profile (no input for ISEDAV=0)
   ! IPROFL == 1 for profile evolution computation(input ISEDAV)
   ! ISEDAV == 0 for unlimited bottom sediment
   ! ISEDAV == 1 for sediment availability limited by hard bottom
   
   ! read (11,*) DUM
   ! IPROFL = NINT(DUM)
   
   ISEDAV = 0
   
   if (IPROFL == 1) then
      ! bdj        read (11, 1110) ISEDAV
      read (11,*) DUM 
      ISEDAV = NINT(DUM) 
   endif

   if (IPROFL == 0) then
      ! Bottom of profile is fixed, set IQYDY = 0 to neglect alongshore gradient of longshore sediment transport. If IQYDY = 1, then include alongshore gradient for cross-shore profile evolution
      IQYDY = 0
   endif
   
   ! if (IPROFL == 1) then
   !    Bottom of profile is not fixed, so we will also need a value for ISEDAV
   !    read (11, 1110) ISEDAV
   ! endif
   
   ! INPUT COMPUTATION OPTION: IPERM == 0 for impermeable bottom, IPERM == 1 for permeable bottom of stone structure
   ! read (11,*) DUM 
   ! IPERM = NINT(DUM) 
   
   ! INPUT COMPUTATION OPTIONS
   ! IOVER == 0 for no wave overtopping and overflow(no additional input)
   ! IOVER == 1 for wave overtopping and overflow on crest(input IWTRAN)
   ! IWTRAN == 0 for no standing water landward of crest
   ! IWTRAN == 1 for wave transmission into bay or estuary
   ! IPOND == 0 for no ponding water landward of JSWL node
   ! IPOND == 1 for ponding water in runnel landward of ridge for IPERM=0
   ! INFILT == 0 for no infiltration landward of sand dune
   ! INFILT == 1 for infiltration if IOVER == 1, IPERM == 0, IPROFL == 1
   ! read (11, *) DUM 
   ! IOVER = NINT(DUM) 

   IWTRAN = 0
   IPOND = 0
   INFILT = 0
   
   ! if (IOVER == 1) then
   !   read (11,*) DUM
   !   IWTRAN = NINT(DUM) 
   ! endif
   
   ! if (IOVER == 1 .and. IWTRAN == 0) then
   !   read (11,*) DUM 
   !   IPOND = NINT(DUM) 
   ! endif
   
   ! if (IOVER == 1 .and. IPERM == 0) then
   !    if (IPROFL == 1) then
   !       read (11,*) DUM
   !       INFILT = NINT(DUM) 
   !    endif
   ! endif

   ! INPUT COMPUTATION OPTION: IWCINT == 0 for no wave and current interaction, IWCINT == 1 for wave and current interaction in frequency dispersion, momentum and wave action equations for IANGLE == 0
   ! read (11,*) DUM 
   ! IWCINT = NINT(DUM) 

   ! INPUT COMPUTATION OPTION: IROLL == 0 for no roller, IROLL == 1 for roller effects in governing equations
   ! read (11,*) DUM 
   ! IROLL = NINT(DUM) 
   
   if (IROLL == 1) RBZERO = 0.1D0

   ! INPUT COMPUTATION OPTION: IWIND == 0 for no wind effect, IWIND == 1 for wind shear stresses on momentum equations
   ! read (11,*) DUM 
   ! IWIND = NINT(DUM) 

   ! INPUT COMPUTATION OPTION: ITIDE == 0 for no tidal effect on currents, ITIDE == 1 for longshore and cross-shore tidal currents

   ! read (11,*) DUM 
   ! ITIDE = NINT(DUM)

   ! COMPUTATIONAL INPUT data: DX = nodal spacing for input bottom geometry
   ! read (11,*) DX 

   ! BREAKER RATIO parameter GAMMA == 0.5 to 1.0
   ! read (11,*) GAMMA 
   
   if (IPROFL == 1) then
      ! Bottom of profile is not fixed, so need sediment characteristics
      ! WF     = sediment fall velocity (m/s)
      ! SG     = sediment specific gravity
      ! D50    = median sediment diameter (mm) converted to (m) below
      ! EFFB   = suspension efficiency due to breaking, eB
      ! EFFF   = suspension efficiency due to friction, ef
      ! SLP    = suspended load parameter
      ! SLPOT  = suspended load parameter due to wave overtopping if IOVER = 1
      ! SPORO  = sediment porosity
      ! SHIELD = critical Shields parameter used if D50 is less than CSEDIA
      ! CSTABN = critical stability number (0.6 or 0.7) used if IPERM = 1
      ! CSEDIA = critical sediment diameter to separate sand and stone
      ! TANPHI = tangent (sediment friction angle)
      ! BLP    = bedload parameter
      ! BLD    = BLP/GRAV/(SG-1) used for bedload transport rate
      ! BEDLM  = parameter m for bedload reduction factor BRF for hard bottom used for ISEDAV = 1
      
      ! The following default values are specified here to reduce input error
      SPORO = 0.4D0
      SHIELD = 0.05D0
      ! CSTABN = 0.6D0
      CSTABN = 0.7D0
      ! EFFF = 0.01D0
      ! EFFB = 0.002D0 to 0.01D0
      ! SLP = 0.1D0 to 0.4D0
      ! SLPOT = 0.1D0 to 3.6D0
      ! TANPHI = 0.63D0 for sand
      ! BLP = 0.001D0 to 0.004D0
      
      if (ISEDAV == 1) BEDLM = 1.0D0

      ! read (11, 1150) D50, WF, SG
      
      ! if (IOVER == 0) read (11, 1150) EFFB, EFFF, SLP
      ! if (IOVER == 1) read (11, 1150) EFFB, EFFF, SLP, SLPOT
      
      ! read (11, 1150) TANPHI, BLP      
      if (EFFF < EFFB) then
         write (*,*) 'CShore ERROR: The suspension efficiency parameter due to bottom friction must be >= the suspension efficiency parameter due to wave breaking'
         
         IError = -2
         
         return
      endif
         
      D50 = D50 * 1.D-3
      SGM1 = SG - 1.D0
      SPORO1 = 1.D0 - SPORO
      WFSGM1 = WF * SGM1
      GSGM1 = GRAV * SGM1
      
      if (IPERM == 0) then
         GSD50S = GSGM1 * D50 * SHIELD
         CSEDIA = 2.D0 * D50
      else
         GSD50S = DSQRT(GSGM1 * D50 * CSTABN)
         CSEDIA = 0.5D0 * D50
      endif
      
      BLD = BLP / GSGM1
   endif

   ! Runup wire height RWH (in metres), if IOVER == 1
   ! if (IOVER == 1) read (11, 1130) RWH

   ! Stone or gravel characteristics, if IPERM = 1. SNP = Stone/gravel porosity in porous layer, SDP = Nominal stone/gravel diameter (m)
   ! if (IPERM == 1) read (11, 1150) SNP, SDP

   ! HWDMIN = minimum water depth (m), is used in the wet and dry zone in subroutine WETDRY
   ! D50 = median sediment diameter (m)   
   if (IOVER == 1) then
      HWDMIN = 1.D-3
      if (IPROFL == 1 .and. IPERM == 0) then
         HWDMIN = D50
         if (DX >= 0.05D0) HWDMIN = 1.D-3
      else
         if (IPROFL == 1) HWDMIN = 1.D-4
      endif
   endif

   ! INPUT WAVE AND WATER LEVEL
   ! If NWAVE == NSURG, then NTIME = number of waves and water levels at x = 0 during time = TIMEBC(i) to time = TIMEBC(i+1)
   ! TIMEBC(i) = time in seconds at the beginning of the specified wave and water level
   ! TPBC(i) = spectral peak or wave period in seconds
   ! HRMSBC(i) = root mean square wave height in meters
   ! WSETBC(i) = wave setup in meters
   ! SWLBC(i) = still water level in meters above the datum used for the input bottom profile
   ! WANGBC(i) = incident wave angle in degrees
   ! For IPROFL = 0, use input TIMEBC(i+1) = 1.0, 2.0, ... to identify each combination of waves and still water level
   
   ! Allocate memory for the wave arrays
   call allocate_offshore_wave_size_arrays(NWAVE+1)

   if (ILAB == 1) then
      ! We have laboratory wave data: waves and water levels are read together
      NTIME = NWAVE
      TIMEBC(1) = 0.D0
      
      do 120 I = 1, NTIME
         ! read (11, 1160) TIMEBC(I+1), TPBC(I), HRMSBC(I), WSETBC(I), SWLBC(I), WANGBC(I)
   
         if (WANGBC(I) < -80.D0 .OR. WANGBC(I) > 80.D0) then
            write (*,*) 'CShore ERROR: incident wave angle =', WANGBC(I), ' but incident wave angle must be between -80 and 80 degrees'

            IError = -3
            return
         endif
120   end do

   else
      ! ILAB == 0, so we have field wave data: wave conditions (NWAVE+1) and water level (NSURG+1) at x = 0 vary continously in time starting from time = 0, unlike step changes assumed for laboratory data. Choose the number of step changes to approximate time series
      NTIME = MAX0(NWAVE, NSURG)
      NWAVE1 = NWAVE + 1
      
      do 130 I = 1, NWAVE1
         ! read (11,*) TWAVE(I), TPIN(I), HRMSIN(I), WANGIN(I)
   
         if (NWAVE == NTIME) then
            TIMEBC(I) = TWAVE(I)
         endif
130   end do

      NSURG1 = NSURG + 1
      
      do 131 I = 1, NSURG1
         ! read (11,*) TSURG(I),SWLIN(I)
   
         if (NSURG == NTIME) then
            TIMEBC(I) = TSURG(I)
         endif
131   end do

      if (TWAVE(1) /= 0.D0 .OR. TSURG(1) /= 0.D0) then
         write (*,*) 'CShore ERROR: non-zero start time for offshore wave conditions and water level'

         IError = -4
         return
      endif
      
      if (TWAVE(NWAVE1) /= TSURG(NSURG1)) then
         write (*,*) 'CShore ERROR: different durations for offshore wave conditions and water level'
         IError = -4
         return
      endif

      ! Subroutine TSINTP interpolates input time series at specified time TIMEBC(I) with I = 1, 2, ... (NTIME+1) and generates stepped time series corresponding to input format for the case of NWAVE = NSURG
      call TSINTP(NWAVE, TWAVE, TPIN, NTIME, TIMEBC, TPBC)
      call TSINTP(NWAVE, TWAVE, HRMSIN, NTIME, TIMEBC, HRMSBC)
      call TSINTP(NWAVE, TWAVE, WANGIN, NTIME, TIMEBC, WANGBC)
      call TSINTP(NSURG, TSURG, SWLIN, NTIME, TIMEBC, SWLBC)
      
      ! Wave setup at x = 0 is assumed to be zero
      WSETBC(1:NTIME) = 0.D0
      
      do 132 I = 1, NTIME
         if (WANGBC(I) < -80.D0 .OR. WANGBC(I) > 80.D0) then
            write (*,*) 'CShore ERROR: incident wave angle =', WANGBC(I), ' but incident wave angle must be between -80 and 80 degrees'
            
            IError = -6
            
            return
         endif
132   end do

   endif
   ! End of field data input (for ILAB == 0)
   
   ! Prepare for ITIME == 1 computation. IF IPOND == 1, ponded water level ZW = SWL at time = 0
   TP = TPBC(1)
   HRMS(1) = HRMSBC(1)
   WSETUP(1) = WSETBC(1)
   ANGLE = WANGBC(1)
   
   if (IPOND == 1) ZW = SWLBC(1)

   if (ANGLE == 0.D0) then
      IANGLE = 0
   else
      IANGLE = 1
      IWCINT = 0
   endif
   
   ! Allocate memory for the bottom geometry arrays
   call allocate_bottom_geometry_size_arrays(JMAXMAX, ILINE)

   ! BOTTOM GEOMETRY and POROUS LAYER BOTTOM if IPERM == 1
   ! The bottom geometry is divided into segments of different inclination and roughness starting from seaward boundary for cross-shore line L.
   ! NPINP(L) = number of input points of impermeable hard bottom
   ! ZP(X) along cross-shore line L only if IPERM = 1 or ISEDAV = 1
   ! XPINP(J,L) = horizontal distance of input point J from x = 0
   ! ZPINP(J,L) = dimensional vertical coordinate in meters of porous layer bottom or hard bottom at point (J) with ZPINP(J) equal to or less than ZBINP(J,L) where ZPINP(1,L) = ZBINP(1,L) imposed
   ! IF ISEDAV == 1, an almost vertical impermeable wall can be specified using two points (NPINP-1) and NPINP where
   ! IVWALL(L) = 0 for no vertical wall along cross-shore line L
   ! IVWALL(L) = 1 for vertical wall with sediment in front
   ! IVWALL(L) = 2 for vertical wall exposed to wave action

   do 160 L = 1, ILINE
      ! Point J = 1 has no corresponding friction factor
      ! read (11,*) XBINP(1,L), ZBINP(1,L)
      XBINP(1,L) = 0.D0
      
      do 140 J = 2, NBINP(L)
         ! read (11,*) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L) ! BDJ
   
         ! If IANGLE = 1, the bottom friction factor must be positive
         if (IANGLE == 1) then
            if (FBINP(J-1,L) <= 0.D0) then
               write (*, 2901) FBINP(J-1,L), (J-1), L
               
               IError = -8
               
               return
            endif
         endif
         
         ! Can't have a perfectly horizontal bottom, or we run into numerical difficulty
         if (ZBINP(J,L) == ZBINP(J-1,L)) then
            ZBINP(J-1,L) = ZBINP(J-1,L) - 1.D-4
         end if
140   end do
2901  format(/'Bottom Friction Factor FBINP(J-1,L) =', D11.4, ' for (J-1) =', I4, ' and L =', I3/'For obliquely incident waves, FBINP must be positive')

      DUM = XBINP(NBINP(L), L) / DX
      IDUM = NINT(DUM)
      DUM = DUM - DBLE(IDUM)
      
      if (DUM < 1.D-5) then
         XBINP(NBINP(L),L) = XBINP(NBINP(L),L) + 1.D-4
      end if

      if (IPERM == 1 .OR. ISEDAV == 1) then
         XPINP(1,L) = 0.D0
         ZPINP(1,L) = ZBINP(1,L)
         ! do 150 J = 2, NPINP(L)
         !    read (11,*) XPINP(J,L), ZPINP(J,L)
! 150      end do

         if (XPINP(NPINP(L),L) < XBINP(NBINP(L),L)) XPINP(NPINP(L),L) = XBINP(NBINP(L),L)
      endif
      
      if (L > 1) DYLINE(L-1) = YLINE(L) - YLINE(L-1)
160 end do

   ! WIND SPEED AND DIRECTION IF IWIND == 1
   ! During time = TIMEBC(i) to time = TIMEBC(i+1)
   ! W10(i) = wind speed (m/s) at 10m elevation above mean sea level
   ! WANGLE(i) = wind direction in degrees at 10 m
   ! WINDCD(i) = wind drag coefficient based on Large and Pond (1981)
   ! TWXSTA(i) = cross-shore wind shear stress/specific water weight
   ! TWYSTA(i) = longshore wind shear stress/specific water weight
   ! RATIO = specific water weight/specific air weight

   ! Wind data time series is read in the same way as field data of waves and water level
   if (IWIND == 1) then
      !   read (11,*) DUM      
      ! NWIND = NINT(DUM)
      
      NWIND1 = NWIND + 1
      
      ! do 190 I= 1, NWIND1
      !      read (11,*) TWIND(I), WIND10(I), WINDAN(I)
! 190   end do

      if (TWIND(1) /= 0.D0) then
         write (*,*) 'CShore ERROR: non-zero start time of wind data'

         IError = -9
         
         return
      endif
      
      if (TWIND(NWIND1) /= TIMEBC(NTIME+1)) then
         write (*,*) 'CShore ERROR: end time of wind data is different from the end time of wave and water level data'
         
         IError = -10
         
         return
      endif
      
      call TSINTP(NWIND, TWIND, WIND10, NTIME, TIMEBC, W10)
      call TSINTP(NWIND, TWIND, WINDAN, NTIME, TIMEBC, WANGLE)
      
      RATIO = 837.D0
      CONVRT = 3.14159D0 / 180.D0
      
      do 2005 I = 1, NTIME
         if (W10(I) > 25.D0) write (*,*) 'CShore ERROR: wind speed at 10m =', W10(I), 'but wind speed must be less than 25 m/s for Large and Pond (1981)'

         if (W10(I) < 11.D0) then
            WINDCD(I) = 1.2D-3
         else
            WINDCD(I) = 0.49D-3 + 0.065D-3 * W10(I)
         endif
         
         DUM = (WINDCD(I) / RATIO / GRAV) * W10(I) * W10(I)
         ANG = CONVRT * WANGLE(I)
         TWXSTA(I) = DUM * DCOS(ANG)
         TWYSTA(I) = DUM * DSIN(ANG)
2005   end do

   else
      do 201 I= 1, NTIME
         TWXSTA(I) = 0.D0
         TWYSTA(I) = 0.D0
201   end do
   endif

   ! LANDWARD STILL WATER LEVEL IF IWTRAN == 1. During time = TIMEBC(i) to time = TIMEBC(i+1), SWLAND(i) = still water level in meters above datum landward of emerged structure or dune. If ISWLSL == 0, seaward and landward still water levels are same and SWLAND(i) = SWLBC(i). If ISWLSL == 1, read time series of landward still water level SLANIN(I) at time TSLAND(I) with I= 1, 2, ..., (NSLAN+1)
   if (IWTRAN == 1) then
      ! read (11,*) DUM
      ! ISWLSL = NINT(DUM)
      
      if (ISWLSL == 0) then
         do 300 I = 1, NTIME
            SWLAND(I) = SWLBC(I)
300      end do

      else
         ! read (11,*) DUM
         ! NSLAN = NINT(DUM)
         
         NSLAN1 = NSLAN + 1
         
         !do 301 I = 1, NSLAN1
            ! read (11,*) TSLAND(I), SLANIN(I)
!301      end do

         if (TSLAND(1) /= 0.D0) then
            write (*, *) 'CShore ERROR: non-zero start time of landward SWL'
            
            IError = -11
            
            return
         endif
         
         if (TSLAND(NSLAN1) /= TIMEBC(NTIME+1)) then
            write (*, *) 'CShore ERROR: the end time of landward SWL differs from the end time of other input time series'
            
            IError = -12
            
            return
         endif
         
         call TSINTP(NSLAN, TSLAND, SLANIN, NTIME, TIMEBC, SWLAND)
      endif
   endif

   ! ALONGSHORE WATER LEVEL GRADIENT IF ITIDE == 1. During time == TIMEBC(i) to time == TIMEBC(i+1)
   ! DETADY(i) = alongshore water level gradient for longshore current
   ! DSWLDT(i) = rate of input water level change only for ILAB == 0

   ! Alongshore gradient data time series is read the same way as field surge data
   if (ITIDE == 1) then
      ! read (11,*) DUM
      ! NTIDE = NINT(DUM)
      
      NTIDE1 = NTIDE + 1
      
      ! do 4001 I= 1, NTIDE1
      !  read (11,*) TTIDE(I), DEDYIN(I)
!4001   end do
!1195  format(D11.1,D11.7)

      if (TTIDE(1) /= 0.D0) then
         write (*, *) 'CShore ERROR: non-zero start time of tide data'
         
         IError = -13
         
         return
      endif
      
      if (TTIDE(NTIDE1) /= TIMEBC(NTIME+1)) then
         write (*, *) 'CShore ERROR: end time of tide data differs from the end time of wave and water level data'
         
         IError = -14
         
         return
      endif
      
      call TSINTP(NTIDE, TTIDE, DEDYIN, NTIME, TIMEBC, DETADY)
   endif

   ! If ITIDE == 1 and ILAB == 0, cross-shore water flux associated with continuous input water level change is accounted for in cross-shore current in wet zone
   if (ITIDE == 1 .and. ILAB == 0) then
      do 4100 I = 1, NSURG
         K = I + 1
         DSDTIN(I) = (SWLIN(K) - SWLIN(I)) / (TSURG(K) - TSURG(I))
4100   end do

      DSDTIN(NSURG1) = DSDTIN(NSURG)
      call TSINTP(NSURG, TSURG, DSDTIN, NTIME, TIMEBC, DSWLDT)
   endif

   ! End of former INPUT routine

   ! Run CShore
   call CShore(IError)

   ! Get the results into the argument output variables
   Out_IError = IError

   ! DFM safety check TEST
   if (IError == -1) return

   Out_nOutSize = nOutSize
   
   do I = 1, nOutSize
      do J = 1, ILINE
         Out_XYDist(I, J) = XYDist(I, J)
         Out_FreeSurfaceStd(I, J) = FreeSurfaceStd(I, J)
         Out_WaveSetupSurge(I, J) = WaveSetupSurge(I, J)
         ! Out_StormSurge(I, J) = StormSurge(I, J) 
         Out_SinWaveAngleRadians(I, J) = SinWaveAngleRadians(I, J)
         Out_FractionBreakingWaves(I, J) = FractionBreakingWaves(I, J)
      end do
   end do
   
   return   
end subroutine CShoreWrapper
#endif   

