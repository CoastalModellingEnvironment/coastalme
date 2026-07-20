#if defined FILEINOUT || ARGINBOTHOUT
!===============================================================================================================================
!
! This subroutine opens all output files, it is built and called only if BUILDVERSION == FILEINOUT or ARGINBOTHOUT
!
!===============================================================================================================================
subroutine OUTPUT_OPENER
   open (unit = 20, file = 'ODOC', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 21, file = 'OBPROF', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 22, file = 'OSETUP', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 23, file = 'OPARAM', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 24, file = 'OXMOME', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 25, file = 'OYMOME', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 26, file = 'OENERG', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 27, file = 'OXVELO', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 28, file = 'OYVELO', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 29, file = 'OROLLE', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 30, file = 'OBSUSL', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 31, file = 'OPORUS', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 32, file = 'OCROSS', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 33, file = 'OLONGS', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 34, file = 'OSWASH', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 35, file = 'OSWASE', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 36, file = 'OTIMSE', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 37, file = 'OCRVOL', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 38, file = 'OLOVOL', status = 'UNKNOWN', access = 'SEQUENTIAL')
   open (unit = 40, file = 'OMESSG', status = 'UNKNOWN', access = 'SEQUENTIAL')

   return
end subroutine OUTPUT_OPENER
#endif
