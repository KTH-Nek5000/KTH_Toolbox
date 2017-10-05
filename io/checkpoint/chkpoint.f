!> @file chkpoint.f
!! @ingroup chkpoint
!! @brief Set of checkpoint routines
!! @details This is a main interface reading/writing runtime parameters
!! and calling proper submodule.
!=======================================================================
!> @brief Register checkpointing module
!! @ingroup chkpoint
!! @note This routine should be called in userchk during first step
!!    after call to rprm_dict_get
      subroutine chkpt_register()
      implicit none

      include 'SIZE'
      include 'CHKPOINTD'
      include 'MNTRLP'

!     local variables
      integer lpmid
!-----------------------------------------------------------------------
!     find parent module
      call mntr_mod_is_name_reg(lpmid,'NEK5000')
!     register module
      call mntr_mod_reg(chpt_id,lpmid,chpt_name,'Checkpointing I/O')
      if (lpmid.le.0) then
         call mntr_log(chpt_id,lp_vrb,
     $        'ERROR: parent module ['//'NEK5000'//'] not registered')
      endif

!     register parameters
      call rprm_rp_log_reg(chpt_ifrst_id,chpt_id,'READCHKPT',
     $     'Restat from checkpoint',.false.)

      call rprm_rp_int_reg(chpt_fnum_id,chpt_id,'CHKPFNUMBER',
     $     'Restart file number',1)

      call rprm_rp_int_reg(chpt_step_id,chpt_id,'CHKPINTERVAL',
     $     'Checkpiont saving frequency (number of time steps)',500)

      call rprm_rp_str_reg(chpt_wtime_id,chpt_id,'WALLTIME',
     $     'Simulation wall time','00:00')

!     call submodule registration
      call chkpts_register

      return
      end subroutine
!=======================================================================
!> @brief Initilise checkpointing module
!! @ingroup chkpoint
!! @note This routine should be called in userchk during first step
!!    after call to rprm_dict_get
      subroutine chkpt_init()
      implicit none

      include 'SIZE'
      include 'CHKPOINTD'
      include 'MNTRLP'

!     local variables
      character*132 lkey, c_out
      integer i_out,ifnd, ierr, nhour, nmin
      real d_out
      logical ifsec
!-----------------------------------------------------------------------
!     get runtime parameters
      call rprm_rp_log_get(chpt_ifrst,chpt_ifrst_id)
      call rprm_rp_int_get(chpt_fnum,chpt_fnum_id)
      call rprm_rp_int_get(chpt_step,chpt_step_id)
      call rprm_rp_str_get(chpt_wtimes,chpt_wtime_id)

!     get wall clock
      lkey = trim(adjustl(chpt_wtimes))
!     check string format
      ierr = 0
      if (lkey(3:3).ne.':') ierr = 1
      if (.not.(LGE(lkey(1:1),'0').and.LLE(lkey(1:1),'9'))) ierr = 1
      if (.not.(LGE(lkey(2:2),'0').and.LLE(lkey(2:2),'9'))) ierr = 1
      if (.not.(LGE(lkey(4:4),'0').and.LLE(lkey(4:4),'9'))) ierr = 1
      if (.not.(LGE(lkey(5:5),'0').and.LLE(lkey(5:5),'9'))) ierr = 1

      if (ierr.eq.0) then
         read(lkey(1:2),'(I2)') nhour
         read(lkey(4:5),'(I2)') nmin
         chpt_wtime = 60.0*(nmin +60*nhour)
      else
         call mntr_log(chpt_id,lp_inf,'Wrong wall time format')
      endif

!     call submodule initialisation
      call chkpts_init

      chpt_ifinit=.true.

      return
      end
!=======================================================================
!> @brief Main checkpoint interface
!! @ingroup chkpoint
!! @note This routine should be called in userchk
      subroutine chkpt_main
      implicit none

      include 'SIZE'
      include 'CHKPOINTD'
      include 'MNTRLP'

!-----------------------------------------------------------------------
      if(chpt_ifrst) then
         call chkpts_read
      endif

      call chkpts_write

      return
      end
!=======================================================================

