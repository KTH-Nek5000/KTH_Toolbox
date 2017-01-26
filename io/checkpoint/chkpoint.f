!> @file chkpoint.f
!! @ingroup chkpoint
!! @brief Set of checkpoint routines
!! @details This is a main interface reading/writing runtime parameters
!! and calling proper submodule.
!=======================================================================
!> @brief Read checkpoint parameters
!! @ingroup chkpoint
      subroutine chkpt_param_get()
      implicit none

      include 'SIZE'            !
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'CHKPOINTD'

!     local variables
      integer i_out,ifnd
      real d_out
!-----------------------------------------------------------------------
!     default values
      CHKPTSTEP = 100
      IFCHKPTRST = .FALSE.
!     read the file
      ierr=0
      if (NID.eq.0) then
!     do we restart
        call finiparser_getBool(i_out,'_chkpoint:ifchkptrst',ifnd)
        if (ifnd.eq.1.and.i_out.eq.1) then
           IFCHKPTRST = .TRUE.
        endif
!     checkpoint frequency
        call finiparser_getDbl(d_out,'_chkpoint:chkptstep',ifnd)
        if (ifnd.eq.1) then
           CHKPTSTEP = int(d_out)
        endif
      endif

!     broadcast data
      call bcast(CHKPTSTEP,ISIZE)
      call bcast(IFCHKPTRST,LSIZE)

      return
      end
!=======================================================================
!> @brief Main checkpoint interface
!! @ingroup chkpoint
!! @note This routine should be called in userchk
      subroutine checkpoint_main
      implicit none

      include 'CHKPOINTD'        ! IFCHKPTRST

!     local variables
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
      if(.not.ifcalled) then
        ifcalled=.TRUE.
        call chkpt_init
      endif

      if(ifchkptrst) call chkpt_read

      call chkpt_write

      return
      end
!=======================================================================

