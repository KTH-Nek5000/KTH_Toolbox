!> @file pwit.f
!! @ingroup pwit
!! @brief Set of subroutines to perform power iterations within time stepper
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Register power iteration module
!! @ingroup pwit
!! @note This interface is called by @ref tstpr_register
      subroutine stepper_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'TSTPRD'
      include 'PWITD'

      ! local variables
      integer lpmid, il
      real ltim
      character*2 str

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,pwit_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(pwit_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,tstpr_name)
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//trim(tstpr_name)//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(pwit_id,lpmid,pwit_name,
     $      'Power iterations for time stepper')

      ! register timers
      ! initialisation
      call mntr_tmr_reg(pwit_tmr_ini_id,tstpr_tmr_ini_id,pwit_id,
     $     'PWIT_INI','Power iteration initialisation time',.true.)
      ! submodule operation
      call mntr_tmr_reg(pwit_tmr_evl_id,tstpr_tmr_evl_id,pwit_id,
     $     'PWIT_EVL','Power iteration evolution time',.true.)

      ! register and set active section
      call rprm_sec_reg(pwit_sec_id,pwit_id,'_'//adjustl(pwit_name),
     $     'Runtime paramere section for power iteration module')
      call rprm_sec_set_act(.true.,pwit_sec_id)

      ! register parameters
      call rprm_rp_reg(pwit_l2n_id,pwit_sec_id,'L2N',
     $     'Vector initial norm',rpar_real,0,1.0,.false.,' ')

      ! set initialisation flag
      pwit_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(pwit_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise power iteration module
!! @ingroup pwit
!! @note This interface is called by @ref tstpr_init
      subroutine stepper_init()
      implicit none

      include 'SIZE'
      include 'SOLN'            ! V[XYZ]P, TP
      include 'MASS'            ! BM1
      include 'FRAMELP'
      include 'TSTPRD'
      include 'PWITD'

      ! local variables
      integer itmp, il, set_in
      real rtmp, ltim, lnorm
      logical ltmp
      character*20 ctmp

      ! to get checkpoint runtime parameters
      integer ierr, lmid, lsid, lrpid

      ! functions
      real dnekclock, cnht_glsc2_wt
      logical chkpts_is_initialised
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (pwit_ifinit) then
         call mntr_warn(pwit_id,
     $        'module ['//trim(pwit_name)//'] already initialised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,pwit_l2n_id,rpar_real)
      pwit_l2n = rtmp

      ! check the restart flag
      ! check if checkpointing module was registered and take parameters
      ierr = 0
      call mntr_mod_is_name_reg(lmid,'CHKPT')
      if (lmid.gt.0) then
         call rprm_sec_is_name_reg(lsid,lmid,'_CHKPT')
         if (lsid.gt.0) then
            ! restart flag
            call rprm_rp_is_name_reg(lrpid,lsid,'READCHKPT',rpar_log)
            if (lrpid.gt.0) then
               call rprm_rp_get(itmp,rtmp,ltmp,ctmp,lrpid,rpar_log)
               pwit_ifrst = ltmp
            else
               ierr = 1
               goto 30
            endif
            if (pwit_ifrst) then
               ! checkpoint set number
               call rprm_rp_is_name_reg(lrpid,lsid,'CHKPFNUMBER',
     $              rpar_int)
               if (lrpid.gt.0) then
                  call rprm_rp_get(itmp,rtmp,ltmp,ctmp,lrpid,rpar_int)
                  pwit_fnum = itmp
               else
                  ierr = 1
                  goto 30
               endif
            endif
         else
            ierr = 1
         endif
      else
         ierr = 1
      endif

 30   continue

      ! check for errors
      call mntr_check_abort(pwit_id,ierr,
     $            'Error reading checkpoint parameters')

      ! read checkpoint file
      if (pwit_ifrst) then
         if(.not.chkpts_is_initialised()) call mntr_abort(pwit_id,
     $        'Checkpointing module not initialised')
         set_in = pwit_fnum -1
         call stepper_read(set_in)
      endif

      ! initial growth rate
      pwit_grw = 0.0

      ! normalise vector
      lnorm = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)
      lnorm = sqrt(pwit_l2n/lnorm)
      call cnht_opcmult (VXP,VYP,VZP,TP,lnorm)

      if (tstpr_pr.ne.0) then
         ! normalise pressure ???????????????????????????
         itmp = nx2*ny2*nz2*nelv
         call cmult(PRP,lnorm,itmp)
      endif

      ! make sure the velocity and temperature fields are continuous at
      ! element faces and edges
      call tstpr_dssum

      ! save intial vector
      call cnht_opcopy (pwit_vx,pwit_vy,pwit_vz,pwit_t,VXP,VYP,VZP,TP)

      ! stamp log file
      call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS initialised')
      call mntr_logr(pwit_id,lp_prd,'L2NORM = ',pwit_l2n)

      ! everything is initialised
      pwit_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(pwit_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup pwit
!! @return stepper_is_initialised
      logical function stepper_is_initialised()
      implicit none

      include 'SIZE'
      include 'PWITD'
!-----------------------------------------------------------------------
      stepper_is_initialised = pwit_ifinit

      return
      end function
!=======================================================================
!> @brief Renormalise vector and check convergence.
!! @ingroup pwit
!! @note This interface is defined in @ref tstpr_main
!! @remarks This routine uses global scratch space SCRUZ
      subroutine stepper_vsolve
      implicit none

      include 'SIZE'            ! NIO
      include 'TSTEP'           ! TIME, LASTEP, NSTEPS
      include 'INPUT'           ! IFHEAT
      include 'MASS'            ! BM1
      include 'SOLN'            ! V[XYZ]P, TP
      include 'FRAMELP'
      include 'TSTPRD'
      include 'PWITD'

      ! scratch space
      real  TA1 (LPX1*LPY1*LPZ1*LELV), TA2 (LPX1*LPY1*LPZ1*LELV),
     $     TA3 (LPX1*LPY1*LPZ1*LELV), TAT (LPX1*LPY1*LPZ1*LELT)
      COMMON /SCRUZ/ TA1, TA2, TA3, TAT

      ! local variables
      integer itmp
      real lnorm, grth_old, ltim

      ! functions
      real dnekclock, cnht_glsc2_wt
!-----------------------------------------------------------------------
      ! timing
      ltim=dnekclock()

      ! normalise vector
      lnorm = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)
      lnorm = sqrt(pwit_l2n/lnorm)
      call cnht_opcmult (VXP,VYP,VZP,TP,lnorm)

      if (tstpr_pr.ne.0) then
         ! normalise pressure ???????????????????????????
         itmp = nx2*ny2*nz2*nelv
         call cmult(PRP,lnorm,itmp)
      endif

      ! make sure the velocity and temperature fields are continuous at
      ! element faces and edges
      call tstpr_dssum

      ! compare current and prevoius growth rate
      grth_old = pwit_grw
      pwit_grw = 1.0/lnorm
      grth_old = pwit_grw - grth_old

      ! get L2 norm of the update
      call cnht_opsub3 (TA1,TA2,TA3,TAT,pwit_vx,pwit_vy,pwit_vz,pwit_t,
     $     VXP,VYP,VZP,TP)
      lnorm = cnht_glsc2_wt(TA1,TA2,TA3,TAT,TA1,TA2,TA3,TAT,BM1)
      lnorm = sqrt(lnorm)

      ! log stamp
      call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS: convergence')
      call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
      call mntr_logr(pwit_id,lp_prd,'Growth ',pwit_grw)

      itmp = 0
      if (IFHEAT) itmp = 1

      !write down current field
      call outpost2(VXP,VYP,VZP,PRP,TP,itmp,'PWI')

      ! write down field difference
      call outpost2(TA1,TA2,TA3,PRP,TAT,itmp,'VDF')

      ! check convergence
      if(lnorm.lt.tstpr_tol.and.grth_old.lt.tstpr_tol) then
         call mntr_log(pwit_id,lp_prd,'Reached stopping criteria')
         ! mark the last step
         LASTEP = 1
      else
         ! save current vector and restart stepper
         call cnht_opcopy(pwit_vx,pwit_vy,pwit_vz,pwit_t,VXP,VYP,VZP,TP)
      endif

      ! save checkpoint
      if (LASTEP.eq.1.or.tstpr_cmax.eq.tstpr_vstep) then
         call stepper_write()

         ! mark the last step
         LASTEP = 1
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(pwit_tmr_evl_id,1,ltim)

      if (LASTEP.eq.1) then
         ! final log stamp
         call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS finalised')
         call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
         call mntr_logr(pwit_id,lp_prd,'Growth ',pwit_grw)
      endif

      return
      end subroutine
!=======================================================================
!> @brief Read restart files
!! @ingroup pwit
!! @param[in]  set_in  restart set number
      subroutine stepper_read(set_in)
      implicit none

      include 'SIZE'            ! NIO
      include 'TSTEP'           ! TIME, LASTEP, NSTEPS
      include 'INPUT'           ! IFMVBD, IFREGUO
      include 'FRAMELP'
      include 'CHKPTD'
      include 'CHKPTMSD'
      include 'PWITD'

      ! argument list
      integer set_in

      ! local variables
      integer ifile, step_cnt, fnum
      character*132 fname(chkptms_fmax)
      logical ifreguol
!-----------------------------------------------------------------------
      ! no regular mesh
      ifreguol= IFREGUO
      IFREGUO = .false.

      call mntr_log(pwit_id,lp_inf,'Reading checkpoint snapshot')

      ! initialise I/O data
      call io_init

      ! get set of file names in the snapshot
      ifile = 1
      call chkptms_set_name(fname, fnum, set_in, ifile)

      ! read files
      call chkptms_restart_read(fname, fnum)

      ! put parameters back
      IFREGUO = ifreguol

      return
      end subroutine
!=======================================================================
!> @brief Write restart files
!! @ingroup pwit
      subroutine stepper_write
      implicit none

      include 'SIZE'            ! NIO
      include 'TSTEP'           ! TIME, LASTEP, NSTEPS
      include 'INPUT'           ! IFMVBD, IFREGUO
      include 'FRAMELP'
      include 'CHKPTD'
      include 'CHKPTMSD'
      include 'PWITD'

      ! local variables
      integer ifile, step_cnt, set_out, fnum
      character*132 fname(chkptms_fmax)
      logical ifcoord
      logical ifreguol
!-----------------------------------------------------------------------
      ! no regular mesh
      ifreguol= IFREGUO
      IFREGUO = .false.

      call mntr_log(pwit_id,lp_inf,'Writing checkpoint snapshot')

      ! initialise I/O data
      call io_init

      ! get set of file names in the snapshot
      ifile = 1
      call chkpt_get_fset(step_cnt, set_out)
      call chkptms_set_name(fname, fnum, set_out, ifile)

      ifcoord = .true.
      ! write down files
      call chkptms_restart_write(fname, fnum, ifcoord)

      ! put parameters back
      IFREGUO = ifreguol

      return
      end subroutine
!=======================================================================

