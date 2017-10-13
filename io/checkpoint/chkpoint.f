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
      include 'RPRMD'
      include 'FRAMELP'
      include 'CHKPOINTD'

!     local variables
      integer lpmid
      integer itest
      real rtest
      logical ltest
      character*20 ctest
!-----------------------------------------------------------------------
!     find parent module
      call mntr_mod_is_name_reg(lpmid,'NEK5000')
!     register module
      call mntr_mod_reg(chpt_id,lpmid,chpt_name,'Checkpointing I/O')
      if (lpmid.le.0) then
         call mntr_log(chpt_id,lp_vrb,
     $        'ERROR: parent module ['//'NEK5000'//'] not registered')
      endif

!     register and set active section
      call rprm_sec_reg(chpt_sec_id,chpt_id,'_'//adjustl(chpt_name),
     $     'Runtime paramere section for checkpoint module')
      call rprm_sec_set_act(.true.,chpt_sec_id)

!     register parameters
      call rprm_rp_reg(chpt_ifrst_id,chpt_sec_id,'READCHKPT',
     $     'Restat from checkpoint',rprm_par_log,0,0.0,.false.,' ')

      call rprm_rp_reg(chpt_fnum_id,chpt_sec_id,'CHKPFNUMBER',
     $     'Restart file number',rprm_par_int,1,0.0,.false.,' ')

      call rprm_rp_reg(chpt_step_id,chpt_sec_id,'CHKPINTERVAL',
     $     'Checkpiont saving frequency (number of time steps)',
     $      rprm_par_int,500,0.0,.false.,' ')

      call rprm_rp_reg(chpt_wtime_id,chpt_sec_id,'WALLTIME',
     $     'Simulation wall time',rprm_par_str,0,0.0,.false.,'00:00')

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
      include 'RPRMD'
      include 'FRAMELP'
      include 'CHKPOINTD'

!     local variables
      integer ierr, nhour, nmin
      integer itmp
      real rtmp
      logical ltmp
      character*20 ctmp

!     functions
      logical chkpts_is_initialised
!-----------------------------------------------------------------------
!     get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_ifrst_id,rprm_par_log)
      chpt_ifrst = ltmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_fnum_id,rprm_par_int)
      chpt_fnum = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_step_id,rprm_par_int)
      chpt_step = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_wtime_id,rprm_par_str)
      chpt_wtimes = ctmp

!     get wall clock
      ctmp = trim(adjustl(chpt_wtimes))
!     check string format
      ierr = 0
      if (ctmp(3:3).ne.':') ierr = 1
      if (.not.(LGE(ctmp(1:1),'0').and.LLE(ctmp(1:1),'9'))) ierr = 1
      if (.not.(LGE(ctmp(2:2),'0').and.LLE(ctmp(2:2),'9'))) ierr = 1
      if (.not.(LGE(ctmp(4:4),'0').and.LLE(ctmp(4:4),'9'))) ierr = 1
      if (.not.(LGE(ctmp(5:5),'0').and.LLE(ctmp(5:5),'9'))) ierr = 1

      if (ierr.eq.0) then
         read(ctmp(1:2),'(I2)') nhour
         read(ctmp(4:5),'(I2)') nmin
         chpt_wtime = 60.0*(nmin +60*nhour)
      else
         call mntr_log(chpt_id,lp_inf,'Wrong wall time format')
      endif

!     call submodule initialisation
      call chkpts_init

!     is everything initialised
      if (chkpts_is_initialised) chpt_ifinit=.true.

      return
      end
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup chkpoint
!! @return chkpt_is_initialised
      logical function chkpt_is_initialised()
      implicit none

      include 'SIZE'
      include 'CHKPOINTD'
!-----------------------------------------------------------------------
      chkpt_is_initialised = chpt_ifinit

      return
      end function
!=======================================================================
!> @brief Main checkpoint interface
!! @ingroup chkpoint
!! @note This routine should be called in userchk
      subroutine chkpt_main
      implicit none

      include 'SIZE'
      include 'CHKPOINTD'
      include 'FRAMELP'

!-----------------------------------------------------------------------
      if(chpt_ifrst) then
         call chkpts_read
      endif

      call chkpts_write

      return
      end
!=======================================================================

