!===============================================================================================================================
!
! This subroutine computes ponded water level and zone if IPOND= 1
!
!===============================================================================================================================
subroutine PONDED(L)
   use CShoreShared

   implicit integer (I-N)
   implicit double precision (A-H, O-Z)
   
   integer, intent (in) :: L
   integer :: J, JPEAK, JRUN
   
   ! Compute the following quantities for known bottom profile ZB(J,L) above datum at XB(J) for node J for line L
   !   ZW = ponded water level at time=TIME
   !   JCREST(L) = ridge crest node landward of SWL node JSWL
   !   JXW = seaward end node of ponded water zone
   !   JX2 = landward end node of ponded water zone

   !   For TIME=0.0, ZW=SWLBC(1) as specified in Subroutine2 INPUT
   !   For TIME>0, compute ZW at present time level
   
   if (TIME > 0.D0) then
      if (JX2 > JXW) then
         ZW=ZW+DELT*(QO(L)-QM)/(XB(JX2)-XB(JXW))
      endif
      if (ZW > ZB(JMAX(L),L)) ZW=ZB(JMAX(L),L)
   endif

   ! NOPOND=0 for ponded water in runnel
   ! NOPOND= 1 for submerged ridge and runnel seaward of node JSWL(L) or dry runnel with no ponded water
   JRUN=JMAX(L)
   JPEAK=JMAX(L)
   do 100 J=(JSWL(L)+1),(JMAX(L)-1)
      if (ZB(J-1,L) >= ZB(J,L) .and. ZB(J,L) < ZB(J+1,L)) then
         if (ZB(J,L) < ZB(JRUN,L) .and. ZW > ZB(J,L)) JRUN=J
      endif
      if (ZB(J,L) >= ZB(JPEAK,L)) JPEAK=J
100 end do
   if (JRUN == JMAX(L)) then
      NOPOND= 1
      JCREST(L)=JPEAK
      RCREST(L)=ZB(JCREST(L),L)
      JXW=JSWL(L)
      JX2=JMAX(L)
      ZW=ZB(JSWL(L),L)
      goto 200
   else
      NOPOND=0
   endif
   
   ! For NOPOND= 1, node JCREST(J) is highest bottom elevation and water level ZW is set to be still water level JCREST(L) = node of ridge crest located between nodes JSWL(L) and JRUN if NOPOND=0
   JCREST(L)=JSWL(L)
   
   do 110 J=(JSWL(L)+1),(JRUN-1)
      if (ZB(J-1,L) <= ZB(J,L) .and. ZB(J,L) > ZB(J+1,L)) then
         if (ZB(J,L) > ZB(JCREST(L),L)) JCREST(L)=J
      endif
110 end do

   if (JCREST(L) == JSWL(L)) then
      NOPOND= 1
      JCREST(L)=JPEAK
      RCREST(L)=ZB(JCREST(L),L)
      JXW=JSWL(L)
      JX2=JMAX(L)
      ZW=ZB(JSWL(L),L)
      goto 200
   endif
   
   RCREST(L)=ZB(JCREST(L),L)
   
   ! If ponded water in runnel is full landward of ridge crest, lower ZW to ridge crest elevation
   if (ZW > ZB(JCREST(L),L)) ZW=ZB(JCREST(L),L)

   ! Find nodes JXW and JX2 at water level ZW
   J=JCREST(L)
   
120 if (ZB(J,L) <= ZW) then
      JXW=J
      goto 121
   else
      J=J+1
      if (J == JRUN) then
         JXW=JRUN-1
         goto 121
      endif
      goto 120
   endif
   
121 J=JRUN

125 if (ZB(J,L) > ZW) then
      JX2=J-1
      goto 200
   else
      J=J+1
      if (J == JMAX(L)) then
         JX2=JMAX(L)
         goto 200
      endif
      goto 125
   endif

200 continue

   return      
end subroutine PONDED

