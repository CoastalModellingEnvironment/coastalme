!===============================================================================================================================
!
! This subroutine calculates parameters used in other subroutines
!
!===============================================================================================================================
subroutine PARAM
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   ! CONSTANTS and parameter
   ! PI      = 3.14159
   ! TWOPI   = 2.D0 * PI
   ! GRAV    = acceleration due to gravity specified in subr.2 INPUT
   ! SQR2    = Sqrt(2)
   ! SQR8    = Sqrt(8)
   ! SQRG1   = Sqrt(2/PI)
   ! SQRG2   = 2*Sqrt(2/PI)
   ! WKPO    = deep water wave number for the representative period

   PI = 3.14159D0
   TWOPI = 2.D0*PI
   GRAV = 9.81D0
   SQR2 = DSQRT(2.D0)
   SQR8 = DSQRT(8.D0)
   SQRG1= DSQRT(2.D0/PI)
   SQRG2= 2.D0*SQRG1
   WKPO = (TWOPI)**2.D0/(GRAV*TP**2.D0)

   ! POROUS FLOW RESISTANCE parameterS IF IPERM= 1
   ! SNP = stone porosity specified in Subroutine2 INPUT
   ! SDP = nominal stone diameter specified in Subroutine2 INPUT
   ! WNU = kinematic viscosity of water (m*m/s)
   ! WPM = maximum seepage velocity (m/s)
   ! If INFILT= 1, WPM is computed using SNP=SPORO and SDP=D50 of sand
   ! in Subroutine2 INPUT
   if (IPERM == 1 .OR. INFILT == 1) then
      WNU = 1.0D-6
      A = 1000.D0
      B = 5.D0
      if (IPERM == 1) then
         DUMP=SNP
         DUMD=SDP
      endif
      if (INFILT == 1) then
         DUMP= 1.D0-SPORO1
         DUMD=D50
      endif
      C = 1.D0 - DUMP
      ALPHA = A*WNU*C**2.D0/(DUMP*DUMD)**2.D0
      BETA1 = B*C/DUMP**3.D0/DUMD
      BETA2 = 7.5D0*B*C/SQR2/DUMP**2.D0
      ! Need to divide BETA2 by WT(J) in Subroutine9 POFLOW
      ALSTA  = ALPHA/GRAV
      BESTA1 = BETA1/GRAV
      BESTA2 = BETA2/GRAV
      ALSTA2 = ALSTA*ALSTA
      BE2    = 2.D0*BESTA1
      BE4    = 2.D0*BE2
      WPM    = (DSQRT(ALSTA2+BE4)-ALSTA)/BE2
   endif

   ! SWASH parameterS IN WET AND DRY ZONE IF IOVER= 1
   ! AWD = swash velocity parameter
   ! AWD= 2.0 calibrated for structures (IPROFL=0 or IPERM= 1)
   ! AWD= 1.6 calibrated for wave overwash of sand dunes
   ! EWD = duration-based exceedance probability for output
   ! where AWD has not been calibrated extensively and
   ! EWD=0.01-0.02 approximately corresponds to 2-percent exceedance
   ! probability based on individual overtopping events.
   if (IOVER == 1) then
      if (IPROFL == 0 .OR. IPERM == 1) then
         AWD= 2.0D0
      else
         AWD= 1.6D0
      endif
      EWD = 0.015D0
      if (IPERM == 1) EWD=0.01D0
      
      ! The following parameters are constant in Subroutine16 WETDRY
      CWD= 0.75D0*DSQRT(PI)
      AQWD = CWD*AWD
      AGWD = AWD*AWD
      AUWD = 0.5D0*DSQRT(PI)*AWD
      BWD = (2.D0-9.D0*PI/16.D0)*AGWD + 1.D0
   endif

   return
end subroutine PARAM
