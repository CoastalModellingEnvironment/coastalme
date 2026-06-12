!===============================================================================================================================
!
! 2011 CShore: Hooyoung ; 2011 April 25
!
!   _____  _____ _                    
!  / ____|/ ____| |                   
! | |    | (___ | |__   ___  _ __ ___ 
! | |     \___ \| '_ \ / _ \| '__/ _ \
! | |____ ____) | | | | (_) | | |  __/
!  \_____|_____/|_| |_|\___/|_|  \___|
!
!
! Nobuhisa Kobayashi:
! Center for Applied Coastal Research, University of Delaware, Newark, Delaware 19716
!
! Cross-shore wave transformation with Brad Johnson in 1998
!
! Cross-shore sediment transport with Yukiko Tega in 2003 and Andres Payo in 2005
!
! Bottom permeability with Lizbeth Meigs in 2004
!
! Roller effects with Haoyu Zhao in 2004
!
! Wave runup and overtopping with Paco de los Santos in 2006 and with Jill Pietropaolo in 2011
!
! Longshore current and sediment transport with Arpit Agarwal in 2005
!
! Longshore bedload transport rate and wind stresses with Andres Payo in 2007
!
! Wave and current interaction and impermeable and permeable wet/dry zone (no sediment and with sediment) with Ali Farhadzadeh in 2008
!
! Sediment transport on hard bottom (limited sediment availability), wave transmission due to wave overtopping, and new input options for storm surge and wave time series as well as for permeable and impermeable bottom profiles with Ali Farhadzadeh in 2009
!
! Calibration, improvement and verification of CShore by Brad Johnson, Mark Gravens and Jens Figlus in 2009
!
! Infiltration landward of dune crest and dip effect above still water shoreline with Jens Figlus in 2010
!
! Onshore ridge migration into ponded runnel, oblique waves on permeable wet/dry zone, and tidal effect on currents with Jens Figlus in 2010
!
! Multiple cross-shore lines for alongshore gradient of longshore sediment transport and its effect on beach profile evolution with Hooyoung Jung in 2011
!
! Converted CShore to a static library, so that it can be called by CoastalME: Dave Favis-Mortlock 2017
!
! Converted the CShore code to (more-or-less) free format Fortran 2003, changed most arrays from static to dynamic, created a wrapper so that CShore is optionally able to get input from and write output to a calling program via arguments: Dave Favis-Mortlock 2018
!
!===============================================================================================================================


!===============================================================================================================================
!
! This module declares variables which are shared between the CShore subroutines
!
!===============================================================================================================================
module CShoreShared
   ! TODO this is crude since every variable in this module is 'global' i.e. is available to every subroutine that uses the module. Really we need "use CShoreShared only <variable list>" in order to be more selective about what is available to individual subroutines
   
   use :: iso_c_binding          ! for C/C++ interopability
   implicit none   
   
   ! Static array sizes used in the Fortran 77 version of CShore
   ! parameter (NN = 10000)      ! Maximum number of cross-shore nodes i.e. input points per profile
   ! parameter (NB = 30000)      ! Maximum number of offshore wave and water level data
   ! parameter (NL = 1)          ! Maximum number of cross-shore lines i.e. number of CoastalME profiles

   ! Wave arrays
   double precision, allocatable, dimension(:) :: TWAVE, TPIN, HRMSIN, WANGIN, TSURG
   
   ! Used when input and output is via arguments to/from a calling program e.g. CoastalME. Note that we are used 2D arrays here, even though the second array dimension is always one at present (since at present CShore only processes one CoastalME coastline-normal profile per CShore run). The second dimension will be used when we process multiple CoastalME profiles per CShore run
   double precision, allocatable, dimension(:, :) :: XYDist, FreeSurfaceStd, SinWaveAngleRadians, FractionBreakingWaves
   double precision, allocatable, dimension(:, :) :: WaveSetupSurge !, StormSurge
   ! Other shared variables
   character(70) :: VER = 'CShore USACE 2011 (CoastalME library version 2018-10-02)'
   integer :: IError = 0, nOutSize = 0
   double precision, allocatable, dimension(:) :: DUMVEC
   double precision, allocatable, dimension(:) :: QTIDE, SMDEDY
   double precision, allocatable, dimension(:) :: SWLIN, TWIND, WIND10, WINDAN, TSLAND, SLANIN, TTIDE, DEDYIN, DSDTIN      
   
   ! /OPTION/  Computation options and time
   integer :: IPROFL, IANGLE, IROLL, IWIND, IPERM, IOVER, IWCINT, ISEDAV, IWTRAN, ILAB, INFILT, IPOND, ITIDE, ILINE, IQYDY
   integer, allocatable, dimension(:) :: IVWALL
   double precision :: TIME   
   
   ! /PERIOD/  Representative period and input wave angle
   double precision :: TP, WKPO, ANGLE
   double precision, allocatable, dimension(:) :: WT
   
   ! /SEAWAV/  Input waves and water levels
   integer :: NWAVE, NSURG, NWIND, NTIME
   double precision, allocatable, dimension(:) :: TIMEBC, TPBC, HRMSBC, WSETBC, SWLBC, WANGBC
   
   ! /PREDIC/  Unknown wave variables predicted by CShore
   double precision, allocatable, dimension(:) :: HRMS, SIGMA, H, WSETUP, SIGSTA
   
   ! /BINPUT/  Input bottom geometry
   integer, allocatable, dimension(:) :: NBINP
   double precision, allocatable, dimension(:) :: XS, YLINE, DYLINE
   double precision, allocatable, dimension(:,:) :: XBINP, ZBINP, FBINP
   
   ! /BPROFL/  Discretized bottom geometry
   integer :: JMAXMAX = 0
   integer, allocatable, dimension(:) :: JMAX, JSWL
   double precision :: DXD2, DXDX, DX2, DX
   double precision, allocatable, dimension(:) :: XB
   double precision, allocatable, dimension(:, :) :: ZB, FB2, SWLDEP, BSLOPE
   
   ! /CONSTA/  Constants
   double precision :: GRAV = 9.81D0, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
   
   ! /LINEAR/  Linear wave values and wave angle quantities
   double precision :: WKP, WKPSIN, FSX, FSY, FE, QWX, QWY
   double precision, allocatable, dimension(:) :: CP, WN, STHETA, CTHETA
   
   ! /FRICTN/  Dimensionless parameters related to bottom friction
   double precision, allocatable, dimension(:) :: GBX, GBY, GF
   
   ! /WBREAK/  Wave breaking quantities and constants
   double precision :: GAMMA, SISMAX
   double precision, allocatable, dimension(:) :: QBREAK, DBSTA, ABREAK

   ! /CRSMOM/  Terms in cross-shore momentum equation
   double precision, allocatable, dimension(:) :: SXXSTA, TBXSTA

   ! /LOGMOM/  Terms in longshore momentum equation
   double precision, allocatable, dimension(:) :: SXYSTA, TBYSTA

   ! /ENERGY/  Terms in energy or wave action equation
   double precision, allocatable, dimension(:) :: EFSTA, DFSTA

   ! /RUNUP/   Parameters for landward computation limit
   integer :: JR
   double precision :: XR, ZR, SSP

   ! /VELOCY/  Mean and standard deviation of horizontal velocities
   double precision, allocatable, dimension(:) :: UMEAN, USTD, USTA, VMEAN, VSTD, VSTA

   ! /SEDINP/  Sediment input parameters
   double precision :: WF, SG, SPORO1, WFSGM1, GSGM1, TANPHI, BSLOP1, BSLOP2, EFFB, EFFF, D50, SHIELD, GSD50S, BLP, SLP, BLD, BEDLM, CSTABN, CSEDIA

   ! /SEDOUT/  Sediment output variables
   double precision, allocatable, dimension(:) :: PS, VS, QSX, QSY, PB, GSLOPE, QBX, QBY, Q

   ! /SEDVOL/  Sediment transport volume per unit width
   double precision, allocatable, dimension(:, :) :: VBX, VSX, VBY, VSY, VY, DZX

   ! /PROCOM/  Beach profile computation variables
   double precision :: DELT
   double precision, allocatable, dimension(:, :) :: DELZB

   ! /ROLLER/  Roller slope,volume flux and related quantities
   double precision :: RBZERO
   double precision, allocatable, dimension(:) :: RBETA, RQ, RX, RY, RE

   ! /POROUS/  Porous flow input and output variables
   integer, allocatable, dimension(:) :: NPINP
   double precision :: WNU, SNP, SDP, ALPHA, BETA1, BETA2, ALSTA, BESTA1, BESTA2
   double precision, allocatable, dimension(:) :: UPMEAN, UPSTD, DPSTA, QP, UPMWD
   double precision, allocatable, dimension(:, :) :: XPINP, ZPINP, ZP, HP

   ! /OVERTF/  Wave overtopping and overflow variables
   integer, allocatable, dimension(:) :: JCREST
   double precision :: RWH, QOTF, SLPOT
   double precision, allocatable, dimension(:) :: RCREST, QO

   ! /WIND/      Wind speed, direction and shear stresses
   double precision, allocatable, dimension(:) :: W10, WANGLE, WINDCD, TWXSTA, TWYSTA

   ! /SWASHP/  Swash parameters for wet and dry zone
   double precision :: AWD, WDN, EWD, CWD, AQWD, BWD, AGWD, AUWD, WPM, ALSTA2, BE2, BE4

   ! /SWASHY/  Computed swash hydrodynamic variables
   integer :: JWD, JDRY
   double precision :: H1
   double precision, allocatable, dimension(:) :: PWET, USWD, HWD, SIGWD, UMEAWD, USTDWD, VMEAWD, VSTDWD, HEWD, UEWD, QEWD

   ! /WATRAN/  Input landward still water level
   integer :: ISWLSL, JSL, JSL1, LANCOM
   double precision, allocatable, dimension(:) :: SWLAND

   ! /COMPAR/  Computational parameters in subroutines
   integer :: NPT, NPE
   double precision :: HWDMIN

   ! /RRPOND/  Variables for ridge and runnel with ponded water
   integer :: JXW, JX2, NOPOND
   double precision :: ZW, QD, QM

   ! /TIDALC/  Tidal input variables for currents
   double precision, allocatable, dimension(:) :: DETADY, DSWLDT

   ! /SERIES/  Time series of wave overtopping and sediment transport rates
   double precision, allocatable, dimension(:) :: TSQO, TSQBX, TSQSX
   
contains


!===============================================================================================================================
!
! Allocates and initialises the arrays which have 'offshore wave' size i.e. maximum number of offshore wave and water level data (these formerly had dimension NB)
!
!===============================================================================================================================
subroutine allocate_offshore_wave_size_arrays(nSize)
   integer :: nSize

#if defined FILEINOUT
   ! If input and output is via arguments, these have already been allocated in CShore_wrapper()
   allocate(TWAVE(nSize), TPIN(nSize), HRMSIN(nSize), WANGIN(nSize), TSURG(nSize))
   TWAVE = 0.0d0
   TPIN = 0.0d0
   HRMSIN = 0.0d0
   WANGIN = 0.0d0
   TSURG = 0.0d0
#endif

   allocate(QTIDE(nSize), SMDEDY(nSize))
   QTIDE = 0.0d0
   SMDEDY = 0.0d0
   
#if defined FILEINOUT
   ! If input and output is via arguments, this has already been allocated in CShore_wrapper()
   allocate(SWLIN(nSize))
   SWLIN = 0.0d0
#endif   
   
   allocate(TWIND(nSize), WIND10(nSize), WINDAN(nSize), TSLAND(nSize), SLANIN(nSize), TTIDE(nSize), DEDYIN(nSize), DSDTIN(nSize))
   TWIND = 0.0d0
   WIND10 = 0.0d0
   WINDAN = 0.0d0
   TSLAND = 0.0d0
   SLANIN = 0.0d0
   TTIDE = 0.0d0
   DEDYIN = 0.0d0
   DSDTIN = 0.0d0

   allocate(TIMEBC(nSize), TPBC(nSize), HRMSBC(nSize), WSETBC(nSize), SWLBC(nSize), WANGBC(nSize))   
   TIMEBC = 0.0d0
   TPBC = 0.0d0
   HRMSBC = 0.0d0
   WSETBC = 0.0d0
   SWLBC = 0.0d0
   WANGBC = 0.0d0
   
   allocate(W10(nSize), WANGLE(nSize), WINDCD(nSize), TWXSTA(nSize), TWYSTA(nSize))
   W10 = 0.0d0
   WANGLE = 0.0d0
   WINDCD = 0.0d0
   TWXSTA = 0.0d0
   TWYSTA = 0.0d0
   
   allocate(SWLAND(nSize))
   SWLAND = 0.0d0
   
   allocate(DETADY(nSize), DSWLDT(nSize))
   DETADY = 0.0d0
   DSWLDT = 0.0d0
   
   return
end subroutine allocate_offshore_wave_size_arrays


!===============================================================================================================================
!
! Allocates and initialises the arrays which have 'cross shore nodes' size i.e. the maximum number of input points per profile (these formerly had dimension NN)
!
!===============================================================================================================================
subroutine allocate_cross_shore_nodes_size_arrays(nSize)
   integer :: nSize

   allocate(DUMVEC(nSize))
   DUMVEC = 0.0d0
   
   allocate(WT(nSize))
   WT = 0.0d0
   
   allocate(HRMS(nSize), SIGMA(nSize), H(nSize), WSETUP(nSize), SIGSTA(nSize))
   HRMS = 0.0d0
   SIGMA = 0.0d0
   H = 0.0d0
   WSETUP = 0.0d0
   SIGSTA = 0.0d0
   
   ! DFM Note that XB must be larger than nSize, because SRSFP() accesses element JR+1 (where JR is the size of the array). This is a bodge
   allocate(XB(nSize+1))
   XB = 0.0d0
   
   allocate(CP(nSize), WN(nSize), STHETA(nSize), CTHETA(nSize))
   CP = 0.0d0
   WN = 0.0d0
   STHETA = 0.0d0
   CTHETA = 0.0d0
   
   allocate(GBX(nSize), GBY(nSize), GF(nSize))
   GBX = 0.0d0
   GBY = 0.0d0
   GF = 0.0d0
   
   allocate(QBREAK(nSize), DBSTA(nSize), ABREAK(nSize))
   QBREAK = 0.0d0
   DBSTA = 0.0d0
   ABREAK = 0.0d0
   
   allocate(SXXSTA(nSize), TBXSTA(nSize))
   SXXSTA = 0.0d0
   TBXSTA = 0.0d0
   
   allocate(SXYSTA(nSize), TBYSTA(nSize))
   SXYSTA = 0.0d0
   TBYSTA = 0.0d0
   
   allocate(EFSTA(nSize), DFSTA(nSize))
   EFSTA = 0.0d0
   DFSTA = 0.0d0
   
   allocate(UMEAN(nSize), USTD(nSize), USTA(nSize), VMEAN(nSize), VSTD(nSize), VSTA(nSize))
   UMEAN = 0.0d0
   USTD = 0.0d0
   USTA = 0.0d0
   VMEAN = 0.0d0
   VSTD = 0.0d0
   VSTA = 0.0d0
   
   allocate(PS(nSize), VS(nSize), QSX(nSize), QSY(nSize), PB(nSize), GSLOPE(nSize), QBX(nSize), QBY(nSize), Q(nSize))
   PS = 0.0d0
   VS = 0.0d0
   QSX = 0.0d0
   QSY = 0.0d0
   PB = 0.0d0
   GSLOPE = 0.0d0
   QBX = 0.0d0
   QBY = 0.0d0
   Q = 0.0d0
   
   allocate(RBETA(nSize), RQ(nSize), RX(nSize), RY(nSize), RE(nSize))
   RBETA= 0.0d0
   RQ = 0.0d0
   RX = 0.0d0
   RY = 0.0d0
   RE = 0.0d0
   
   allocate(UPMEAN(nSize), UPSTD(nSize), DPSTA(nSize), QP(nSize), UPMWD(nSize))
   UPMEAN = 0.0d0
   UPSTD = 0.0d0
   DPSTA = 0.0d0
   QP = 0.0d0
   UPMWD = 0.0d0
   
   allocate(PWET(nSize), USWD(nSize), HWD(nSize), SIGWD(nSize), UMEAWD(nSize), USTDWD(nSize), VMEAWD(nSize), VSTDWD(nSize), HEWD(nSize), UEWD(nSize), QEWD(nSize))
   PWET = 0.0d0
   USWD = 0.0d0
   HWD = 0.0d0
   SIGWD = 0.0d0
   UMEAWD = 0.0d0
   USTDWD = 0.0d0
   VMEAWD = 0.0d0
   VSTDWD = 0.0d0
   HEWD = 0.0d0
   UEWD = 0.0d0
   QEWD = 0.0d0
   
   return
end subroutine allocate_cross_shore_nodes_size_arrays


!===============================================================================================================================
!
! Allocates and initialises the arrays which have 'cross shore lines' size i.e. maximum number of cross-shore lines aka CoastalME profiles (these formerly had dimension NL)
!
!===============================================================================================================================
subroutine allocate_cross_shore_lines_size_arrays(nSize)
   integer :: nSize
   
  allocate(NBINP(nSize))
   NBINP = 0

   allocate(IVWALL(nSize))
   IVWALL = 0
   
   allocate(XS(nSize))
   XS = 0.0d0

   allocate(YLINE(nSize))
   YLINE = 0.0d0

   allocate(DYLINE(nSize))
   DYLINE = 0.0d0
   
   allocate(JMAX(nSize))
   JMAX = 0

   allocate(JSWL(nSize))
   JSWL = 0
   
   allocate(NPINP(nSize))
   NPINP = 0
   
   allocate(JCREST(nSize))
   JCREST = 0
   
   allocate(RCREST(nSize))
   RCREST = 0.0d0

   allocate(QO(nSize))
   QO = 0.0d0
   
   allocate(TSQO(nSize))
   TSQO = 0.0d0

   allocate(TSQBX(nSize))
   TSQBX = 0.0d0

   allocate(TSQSX(nSize))
   TSQSX = 0.0d0
   
   return
end subroutine allocate_cross_shore_lines_size_arrays


!===============================================================================================================================
!
! Allocates and initialises the 2D arrays which have 'bottom geometry' size (these formerly had dimension NN, NL)
!
!===============================================================================================================================
subroutine allocate_bottom_geometry_size_arrays(nSize1, nSize2)
   integer :: nSize1, nSize2
   
#if defined FILEINOUT
   ! If input and output is via arguments, these have already been allocated in CShore_wrapper()
   ! DFM Note that XBINP, ZBINP and FBINP must be larger than nSize1 in the first dimension, because element K+1 (where K is the size of the first dimension) gets accessed in various places e.g. BOTTOM(). Is done for ARGINOUT in CShore_wrapper.f03. This is a bodge
   allocate(XBINP(nSize1+1, nSize2), ZBINP(nSize1+1, nSize2), FBINP(nSize1+1, nSize2))
   XBINP = 0.0d0
   ZBINP = 0.0d0
   FBINP = 0.0d0
#endif

   allocate(VBX(nSize1, nSize2), VSX(nSize1, nSize2), VBY(nSize1, nSize2), VSY(nSize1, nSize2), VY(nSize1, nSize2), DZX(nSize1, nSize2))
   VBX = 0.0d0
   VSX = 0.0d0
   VBY = 0.0d0
   VSY = 0.0d0
   VY = 0.0d0
   DZX = 0.0d0
   
   allocate(DELZB(nSize1, nSize2))
   DELZB = 0.0d0
   
   allocate(XPINP(nSize1, nSize2), ZPINP(nSize1, nSize2), ZP(nSize1, nSize2), HP(nSize1, nSize2))
   XPINP = 0.0d0
   ZPINP = 0.0d0
   ZP = 0.0d0
   HP = 0.0d0
   
   return
end subroutine allocate_bottom_geometry_size_arrays


!===============================================================================================================================
!
! Allocates and initialises the 2D arrays which have 'bottom geometry output' size (these formerly had dimension NN, NL)
!
!===============================================================================================================================
subroutine allocate_bottom_geometry_output_size_arrays(nSize1, nSize2)
   integer :: nSize1, nSize2
   
   ! DFM Note that ZB must be larger than nSize1 in the first dimension, because SRSFP() accesses element JR+1 (where JR is the size of the first dimension). This is a bodge
   allocate(ZB(nSize1+1, nSize2), FB2(nSize1, nSize2), SWLDEP(nSize1, nSize2), BSLOPE(nSize1, nSize2))
   ZB = 0.0d0
   FB2 = 0.0d0
   SWLDEP = 0.0d0
   BSLOPE = 0.0d0
   
   return
end subroutine allocate_bottom_geometry_output_size_arrays


#if defined ARGINOUT || ARGINBOTHOUT
!===============================================================================================================================
!
! Allocates and initialises the 2D arrays which are used for argument output
!
!===============================================================================================================================
subroutine allocate_argument_output_arrays(nSize1, nSize2)
   integer :: nSize1, nSize2

   if (.not. allocated(XYDist)) allocate(XYDist(nSize1, nSize2))
   if (.not. allocated(FreeSurfaceStd)) allocate(FreeSurfaceStd(nSize1, nSize2))
   if (.not. allocated(SinWaveAngleRadians)) allocate(SinWaveAngleRadians(nSize1, nSize2))
   if (.not. allocated(FractionBreakingWaves)) allocate(FractionBreakingWaves(nSize1, nSize2))
   if (.not. allocated(WaveSetupSurge)) allocate(WaveSetupSurge(nSize1, nSize2))
   ! if (.not. allocated(StormSurge)) allocate(StormSurge(nSize1, nSize2))
   
   ! allocate(XYDist(nSize1, nSize2))
   ! allocate(FreeSurfaceStd(nSize1, nSize2))
   ! allocate(SinWaveAngleRadians(nSize1, nSize2))
   ! allocate(FractionBreakingWaves(nSize1, nSize2))
   ! allocate(WaveSetup(nSize1, nSize2))
   ! allocate(StormSurge(nSize1, nSize2))

   return
end subroutine allocate_argument_output_arrays
#endif


!===============================================================================================================================
!
! Deallocates all arrays
!
!===============================================================================================================================
subroutine deallocate_all_arrays()
   deallocate(TWAVE, TPIN, HRMSIN, WANGIN, TSURG)
   deallocate(QTIDE, SMDEDY)
   deallocate(SWLIN, TWIND, WIND10, WINDAN, TSLAND, SLANIN, TTIDE, DEDYIN, DSDTIN)
   deallocate(TIMEBC, TPBC, HRMSBC, WSETBC, SWLBC, WANGBC)
   deallocate(W10, WANGLE, WINDCD, TWXSTA, TWYSTA)
   deallocate(SWLAND)
   deallocate(DETADY, DSWLDT)
   deallocate(DUMVEC)
   deallocate(WT)
   deallocate(HRMS, SIGMA, H, WSETUP, SIGSTA)
   deallocate(XB)
   deallocate(CP, WN, STHETA, CTHETA)
   deallocate(GBX, GBY, GF)
   deallocate(QBREAK, DBSTA, ABREAK)
   deallocate(SXXSTA, TBXSTA)
   deallocate(SXYSTA, TBYSTA)
   deallocate(EFSTA, DFSTA)
   deallocate(UMEAN, USTD, USTA, VMEAN, VSTD, VSTA)
   deallocate(PS, VS, QSX, QSY, PB, GSLOPE, QBX, QBY, Q)
   deallocate(RBETA, RQ, RX, RY, RE)
   deallocate(UPMEAN, UPSTD, DPSTA, QP, UPMWD)
   deallocate(PWET, USWD, HWD, SIGWD, UMEAWD, USTDWD, VMEAWD, VSTDWD, HEWD, UEWD, QEWD)
   deallocate(IVWALL)
   deallocate(NBINP)
   deallocate(XS, YLINE, DYLINE)
   deallocate(JMAX, JSWL)
   deallocate(NPINP)
   deallocate(JCREST)
   deallocate(RCREST, QO)
   deallocate(TSQO, TSQBX, TSQSX)
   deallocate(XBINP, ZBINP, FBINP)   
   deallocate(ZB, FB2, SWLDEP, BSLOPE)
   deallocate(VBX, VSX, VBY, VSY, VY, DZX)
   deallocate(DELZB)
   deallocate(XPINP, ZPINP, ZP, HP)
   
!#if defined ARGINOUT   
!   deallocate(XYDist, FreeSurfaceStd, SinWaveAngleRadians, FractionBreakingWaves, WaveSetup, StormSurge)
!#endif
   
   return
end subroutine deallocate_all_arrays


end module CShoreShared
