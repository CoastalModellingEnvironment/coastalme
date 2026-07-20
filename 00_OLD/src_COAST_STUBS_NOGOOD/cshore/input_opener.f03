#if defined FILEINOUT
!===============================================================================================================================
!
! This subroutine opens the input file, it is built and called only if BUILDVERSION == FILEINOUT
!
!===============================================================================================================================
subroutine INPUT_OPENER
   open (unit = 11, file = 'infile', status = 'OLD', access = 'SEQUENTIAL')
   
   return
end subroutine INPUT_OPENER
#endif

