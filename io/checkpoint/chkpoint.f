!> @file chkpoint.f
!! @ingroup chkpoint
!! @brief Set of checkpoint routines
!! @details This is a main interface reading/writing runtime parameters
!! and calling proper submodule.
!=======================================================================
!> @brief Read checkpoint parameters
!! @ingroup chkpoint
!! @param[in]  fid    file unit
!! @note This routine should be called in @ref runprm_in
!! @todo Check iostat value for missing namelist in the file
!! @see @ref readers_writers_page
      subroutine chkpt_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'CHKPOINTD'

!     argument list
      integer fid

!     local variables
      integer ierr

!     namelists
      namelist /CHKPOINT/ CHKPTSTEP, IFCHKPTRST
!-----------------------------------------------------------------------
!     default values
      CHKPTSTEP = 100
      IFCHKPTRST = .FALSE.
!     read the file
      ierr=0
      if (NID.eq.0) then
        rewind(fid)
        read(unit=fid,nml=CHKPOINT,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading CHKPOINT parameters.$')

!     broadcast data
      call bcast(CHKPTSTEP,ISIZE)
      call bcast(IFCHKPTRST,LSIZE)

      return
      end
!=======================================================================
!> @brief Write checkpoint parameters
!! @ingroup chkpoint
!! @param[in]  fid    file unit
!! @note This routine should be called in @ref runprm_out
!! @see @ref readers_writers_page
      subroutine chkpt_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'CHKPOINTD'

!     argument list
      integer fid

!     local variables
      integer ierr

!     namelists
      namelist /CHKPOINT/ CHKPTSTEP, IFCHKPTRST
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=CHKPOINT,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing CHKPOINT parameters.$')

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

