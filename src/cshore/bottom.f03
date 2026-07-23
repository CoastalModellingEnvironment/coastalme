!===============================================================================================================================
!
! This subroutine calculates the bottom geometry using input DX between two adjacent nodes along ILINE cross-shore line. Smooth input ZB(J,L) to reduce numerical irregularity
!
!===============================================================================================================================
subroutine BOTTOM
   use CShoreShared
   
   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   interface 
      subroutine SMOOTH(NUM, RAW, F)
         integer, intent(in) :: NUM
         double precision, intent(in), dimension(:) :: RAW
         double precision, intent(out), dimension(:) :: F
      end subroutine SMOOTH
   end interface

   double precision, allocatable, dimension(:) ::  SLOPE, PSLOPE, ZBRAW, ZPRAW

   ! The structure geometry is divided into segments of different inclination and roughness for each cross-shore line L.
   ! NBINP(L)    = number of input bottom points
   ! For segments starting from the seaward boundary:
   ! SLOPE(K)    = slope of segment K(+ upslope, - downslope)
   ! FBINP(K,L)  = bottom friction factor
   ! XBINP(K,L)  = dimensional horizontal distance from seaward boundary to the seaward end of segment K
   ! ZBINP(K,L)  = dimensional vertical coordinate (+ above datum) at the seaward end of segment K
   ! PSLOPE(K)   = slope of porous layer bottom
   ! XPINP(K,L)  = dimensional horizontal distance of porous layer bottom at the seaward end of segment K
   ! ZPINP(K,L)  = dimensional vertical coordinate of porous layer bottom at the seaward end of segment K
   
   do L = 1, ILINE
      DUM = XBINP(NBINP(L), L) / DX
      JMAX(L) = NINT(DUM) + 1
      DUM = DX * DBLE(JMAX(L) - 1) - XBINP(NBINP(L), L)
      
      if (DUM > 0.D0) JMAX(L) = JMAX(L) - 1

      if (JMAX(L) > JMAXMAX) JMAXMAX = JMAX(L)
   end do
   
   ! Allocate memory for local arrays
   allocate(SLOPE(JMAXMAX), PSLOPE(JMAXMAX), ZBRAW(JMAXMAX), ZPRAW(JMAXMAX))
   
   ! Allocate memory for shared arrays
   call allocate_bottom_geometry_output_size_arrays(JMAXMAX, ILINE)

   do 100 L = 1, ILINE
      do 120 K = 1, NBINP(L)-1
         DUM = XBINP(K+1,L) - XBINP(K,L)
         SLOPE(K) = (ZBINP(K+1,L) - ZBINP(K,L)) / DUM
         
!          write (*,*) "K =", K, " SLOPE(K) =", SLOPE(K)
120   end do

      ! No vertical wall at landward end unless IVWALL == 1 or 2
      IVWALL(L) = 0
      if (IPERM == 1 .OR. ISEDAV == 1) then
         do 121 K = 1, NPINP(L)-1
            DUM = XPINP(K+1,L) - XPINP(K,L)
            PSLOPE(K) = (ZPINP(K+1,L) - ZPINP(K,L)) / DUM
121      end do
         if (PSLOPE(NPINP(L)-1) > TANPHI) IVWALL(L) = 1
      endif

      ! INITIAL SHORELINE LOCATION AT DATUM Z == 0
      ! XS(L) == horizontal distance between X == 0 and shoreline
      K = 0
900   continue
      if (K == NBINP(L)) then
         XS(L) = XBINP(NBINP(L),L)
         
         ! write (*,*) "Leaving loop XS(L) =", XS(L)
         goto 901
      endif
      
      K = K + 1      
      ! write (*,*) "K =", K, " L =", L, " NBINP(L) =", NBINP(L), " XS(L) =", XS(L)

      CROSS = ZBINP(K,L) * ZBINP(K+1,L)
      ! write (*,*) "CROSS =", CROSS
      
      IF (CROSS > 0.D0) goto 900

      ! DFM safety check bodge ==============================
      if (SLOPE(K) == 0.D0) SLOPE(K) = 1.0D-9
      ! DFM safety check bodge ==============================

      XS(L) = XBINP(K+1,L) - ZBINP(K+1,L) / SLOPE(K)      
      
901   if (L == 1) then
         DXD2 = DX / 2.D0
         DX2 = 2.D0 * DX
         DXDX = DX * DX

         ! NPT = integer used in SMOOTH()
         ! NPE = integer used in EXTRAPO()
         ! bdj 2014-04-29 changed bottom smoothing to rely on rms time-series (not first value at first time)
         NPT = 1 + NINT(maxval(HRMSBC) / DX)
         NPE = 1 + NINT(maxval(HRMSBC) / DX2)
         !      NPT= 1 + NINT(HRMS(1) / DX)
         !      NPE= 1 + NINT(HRMS(1) / DX2)
         ! bdj 2014-04-29 end
         
         if (IPROFL == 1 .and. IPERM == 1) NPT = NPT + 2 * NINT(SDP / DX)
      endif
      
      ! ... CALCULATE BOTTOM GEOMETRY AT EACH NODE
      ! JMAX(L)     = landward edge node corresponding to maximum node number
      ! XB(J)       = horizontal coordinate of node J where XB(1) = 0
      ! ZB(J,L)     = vertical coordinate of bottom at node J (+ above datum)
      ! BSLOPE(J,L) = bottom slope at node J for cross-shore line L
      ! SLOPE(K)    = tangent of local slope of segment K

      ! INTERPOLATION OF BOTTOM POSITION at XB(J)
      ! RCREST(L)   = crest (highest) elevation above datum Z == 0
      ! JCREST(L)   = nodal location of crest for cross-shore Line L
      ! If IPOND == 1, JCREST(L) = nodal location of ridge crest computed in PONDED()
   
      if (L == 1) JDUM = JMAX(L)
      
      if (JMAX(L) < JDUM) write (*,*) "JMAX(L) < JDUM, goto 130"
      if (JMAX(L) < JDUM) goto 130
      JDUM = JMAX(L)
      
      ! DFM safety check bodge ==============================
      if (JDUM > 1000) JDUM = 1000
      ! DFM safety check bodge ==============================

      do 141 J = 1, JDUM
!      do 141 J = 1, JMAX(L)     ! DFM replaced with line above
         XB(J) = DX * DBLE(J-1)
141   end do

130   continue
      ZBRAW(1) = ZBINP(1,L)
      FB2(1,L) = 0.5D0 * FBINP(1,L)
      RCREST(L) = ZBRAW(1)
      
      ! DFM safety check bodge ==============================
      JDUM = JMAX(L)
      ! DFM safety check bodge ==============================

      if (JDUM > 1000) JDUM = 1000

      do 142 J = 2, JDUM
!      do 142 J = 2, JMAX(L)     ! DFM
         do 143 K = 1, NBINP(L)-1
            if ((XB(J) > XBINP(K,L)) .and. (XB(J) <= XBINP(K+1,L))) then
               ZBRAW(J) = ZBINP(K,L) + (XB(J) - XBINP(K,L)) * SLOPE(K)
               FB2(J,L) = 0.5D0 * FBINP(K,L)
               goto 144
            endif
143      end do

144      DUM = ZBRAW(J) - RCREST(L)
         if (IPROFL == 0 .and. DUM >= 0.D0) then
            RCREST(L) = ZBRAW(J)
            JCREST(L) = J
         endif
         if (IPERM == 1 .OR. ISEDAV == 1) then
            if (J == 2) ZPRAW(1) = ZPINP(1,L)
            
            do 145 K= 1, NPINP(L)-1
               if ((XB(J) > XPINP(K,L)) .and. (XB(J) <= XPINP(K+1,L))) then
                  ZPRAW(J) = ZPINP(K,L) + (XB(J) - XPINP(K,L)) * PSLOPE(K)
                  goto 142
               endif
145         end do
         endif
142   end do

      ! Smooth ZBRAW(J) and ZPRAW(J) J = 1-JMAX(L) using SMOOTH()
      JMAXL = JMAX(L)
      
      call SMOOTH(JMAXL, ZBRAW, SLOPE)
      
      if (IPERM == 1 .or. ISEDAV == 1) call SMOOTH(JMAXL, ZPRAW, PSLOPE)
      
      do 149 J= 1, JMAX(L)
         ZB(J,L) = SLOPE(J)
         if (IPERM == 1 .OR. ISEDAV == 1) ZP(J,L)=PSLOPE(J)
149   end do

      ! Calculate bottom slope and JCREST(if IPROFL= 1) using smoothed ZB(J)
      BSLOPE(1,L) = (ZB(2,L) - ZB(1,L))/DX
      
      JMAXM1 = JMAX(L) - 1
      
      BSLOPE(JMAX(L),L) = (ZB(JMAX(L),L) - ZB(JMAXM1,L))/DX
      
      do 150 J= 2, JMAXM1
         BSLOPE(J,L) = (ZB(J+1,L) - ZB(J-1,L))/DX2
150   end do
      
      if (IPROFL == 1 .and. IPOND == 0) then
         RCREST(L) = ZB(1,L)
         
         do 151 J = 2, JMAX(L)
            DUM = ZB(J,L) - RCREST(L)
            if (DUM >= 0.D0) then
               RCREST(L) = ZB(J,L)
               JCREST(L) = J
            endif
151      end do
      endif

      ! HP(J,L) = vertical thickness of porous or sediment layer
      if (IPERM == 1 .or. ISEDAV == 1) then
         do 210 J = 1, JMAX(L)
            HP(J,L) = ZB(J,L) - ZP(J,L)
            if (HP(J,L) < 0.D0) then
               HP(J,L) = 0.D0
               ZP(J,L) = ZB(J,L)
            endif
210      end do
      endif

100 end do

   deallocate(SLOPE, PSLOPE, ZBRAW, ZPRAW)
   
   return
end subroutine BOTTOM

