!> @file tstpr.f
!! @ingroup tstpr
!! @brief Set of subroutines to use time steppers for e.g. power
!!    iterations or solution of eigenvalue problem with Arnoldi algorithm
!! @author Adam Peplinski
!! @date Mar 7, 2016
! preprocessing flag for pressure reconstruction
!#define PRS_REC
#undef PRS_REC
!=======================================================================
!> @brief Register time stepper module
!! @ingroup tstpr
!! @note This routine should be called in frame_usr_register
      subroutine tstpr_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'TSTPRD'

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
      call mntr_mod_is_name_reg(lpmid,tstpr_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(tstpr_name)//'] already registered')
         return
      endif

      ! check if conjugated heat transfer module was registered
      call mntr_mod_is_name_reg(lpmid,'CNHT')
      if (lpmid.gt.0)  then
         call mntr_warn(lpmid,
     $        'module ['//'CNHT'//'] already registered')
      else
         call cnht_register()
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(tstpr_id,lpmid,tstpr_name,
     $      'Time stepper')

      ! register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      ! total time
      call mntr_tmr_reg(tstpr_tmr_tot_id,lpmid,tstpr_id,
     $     'TSTPR_TOT','Time stepper total time',.false.)
      lpmid = tstpr_tmr_tot_id
      ! initialisation
      call mntr_tmr_reg(tstpr_tmr_ini_id,lpmid,tstpr_id,
     $     'TSTPR_INI','Time stepper initialisation time',.true.)
      ! submodule operation
      call mntr_tmr_reg(tstpr_tmr_evl_id,lpmid,tstpr_id,
     $     'TSTPR_EVL','Time stepper evolution time',.true.)

      ! register and set active section
      call rprm_sec_reg(tstpr_sec_id,tstpr_id,'_'//adjustl(tstpr_name),
     $     'Runtime paramere section for time stepper module')
      call rprm_sec_set_act(.true.,tstpr_sec_id)

      ! register parameters
      call rprm_rp_reg(tstpr_mode_id,tstpr_sec_id,'MODE',
     $     'Simulation mode',rpar_str,10,0.0,.false.,'DIR')

      call rprm_rp_reg(tstpr_step_id,tstpr_sec_id,'STEPS',
     $     'Length of stepper phase',rpar_int,40,0.0,.false.,' ')

      call rprm_rp_reg(tstpr_cmax_id,tstpr_sec_id,'MAXCYC',
     $     'Max number of stepper cycles',rpar_int,10,0.0,.false.,' ')

      call rprm_rp_reg(tstpr_tol_id,tstpr_sec_id,'TOL',
     $    'Convergence threshold',rpar_real,0,1.0d-6,.false.,' ')

      ! place for submodule registration
      ! register arnoldi or power iterations
      call stepper_register()

      ! set initialisation flag
      tstpr_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(tstpr_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise time stepper module
!! @ingroup tstpr
!! @note This routine should be called in frame_usr_init
      subroutine tstpr_init()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'TSTEP'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'ADJOINT'
      include 'TSTPRD'

      ! local variables
      integer itmp, il
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      ! functions
      real dnekclock, cnht_glsc2_wt
      logical cnht_is_initialised
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (tstpr_ifinit) then
         call mntr_warn(tstpr_id,
     $        'module ['//trim(tstpr_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! intialise conjugated heat transfer
      if (.not.cnht_is_initialised()) call cnht_init

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_mode_id,rpar_str)
      if (trim(ctmp).eq.'DIR') then
        tstpr_mode = 1
      else if (trim(ctmp).eq.'ADJ') then
        tstpr_mode = 2
      else if (trim(ctmp).eq.'OIC') then
        tstpr_mode = 3
      else
        call mntr_abort(tstpr_id,
     $        'wrong simulation mode; possible values: DIR, ADJ, OIC')
      endif

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_step_id,rpar_int)
      tstpr_step = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_cmax_id,rpar_int)
      tstpr_cmax = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_tol_id,rpar_real)
      tstpr_tol = rtmp

      ! check simulation parameters
      if (.not.IFTRAN) call mntr_abort(tstpr_id,
     $   'time stepper requres transient simulation; IFTRAN=.T.')

      if (NSTEPS.eq.0) call mntr_abort(tstpr_id,
     $   'time stepper requres NSTEPS>0')

      if (PARAM(12).ge.0) call mntr_abort(tstpr_id,
     $   'time stepper assumes constant dt')

      if (.not.IFPERT) call mntr_abort(tstpr_id,
     $   'time stepper has to be run in perturbation mode')

      if (IFBASE)  call mntr_abort(tstpr_id,
     $   'time stepper assumes constatnt base flow')

      if (NPERT.ne.1) call mntr_abort(tstpr_id,
     $   'time stepper requires NPERT=1')

      ! initialise cycle counters
      tstpr_istep = 0
      tstpr_vstep = 0

      ! vector length
      tstpr_nv  = NX1*NY1*NZ1*NELV ! velocity single component
      if (IFHEAT) then        !temperature
         tstpr_nt  = NX1*NY1*NZ1*NELT
      else
         tstpr_nt  = 0
      endif
      tstpr_np  = NX2*NY2*NZ2*NELV ! presure

      ! place for submodule initialisation
      ! arnoldi or power iterations
      call stepper_init

      ! zero presure
      call rzero(PRP,tstpr_np)

      ! set initial time
      TIME=0.0

      ! make sure NSTEPS is bigger than the possible number of iterations
      ! in time stepper phase; multiplication by 2 for OIC
      NSTEPS = max(NSTEPS,tstpr_step*tstpr_cmax*2+10)

      IFADJ = .FALSE.
      if (tstpr_mode.eq.2) then
         ! Is it adjoint mode
         IFADJ = .TRUE.
      elseif  (tstpr_mode.eq.3) then
         ! If it is optimal initial condition save initial L2 norm
         tstpr_L2ini = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)

         if (tstpr_L2ini.eq.0.0) call mntr_abort(tstpr_id,
     $   'tstpr_init, tstpr_L2ini = 0')

         call mntr_log(tstpr_id,lp_prd,
     $  'Optimal initial condition; direct phase start')
      endif

      ! set cpfld for conjugated heat transfer
      if (IFHEAT) call cnht_cpfld_set

      ! everything is initialised
      tstpr_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(tstpr_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup tstpr
!! @return tstpr_is_initialised
      logical function tstpr_is_initialised()
      implicit none

      include 'SIZE'
      include 'TSTPRD'
!-----------------------------------------------------------------------
      tstpr_is_initialised = tstpr_ifinit

      return
      end function
!=======================================================================
!> @brief Control time stepper after every nek5000 step and call suitable
!! stepper_vsolve of required submodule
!! @ingroup tstpr
      subroutine tstpr_main()
      implicit none

      include 'SIZE'            ! NIO
      include 'TSTEP'           ! ISTEP, TIME
      include 'INPUT'           ! IFHEAT, IF3D
      include 'MASS'            ! BM1
      include 'SOLN'            ! V[XYZ]P, PRP, TP, VMULT, V?MASK
      include 'ADJOINT'         ! IFADJ
      include 'FRAMELP'
      include 'TSTPRD'

      ! global comunication in nekton
      integer nidd,npp,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      ! local variables
      real grw         ! growth rate
      real ltim        ! timing

      ! functions
      real dnekclock, cnht_glsc2_wt
!-----------------------------------------------------------------------
      if (ISTEP.eq.0) return

      ! step counting
      tstpr_istep = tstpr_istep + 1

      ! stepper phase end
      if (mod(tstpr_istep,tstpr_step).eq.0) then
         ! timing
         ltim = dnekclock()
         ! check for the calculation mode
         if (tstpr_mode.eq.3.and.(.not.IFADJ)) then
            ! optimal initial condition

            call mntr_log(tstpr_id,lp_prd,
     $      'Optimal initial condition; adjoint phase start')

            IFADJ = .TRUE.

            ! iteration count
            tstpr_istep = 0

            ! set time and iteration number
            TIME=0.0
            ISTEP=0

            ! get L2 norm after direct phase
            tstpr_L2dir = cnht_glsc2_wt(VXP,VYP,VZP,TP,
     $         VXP,VYP,VZP,TP,BM1)
            ! normalise vector
            grw = sqrt(tstpr_L2ini/tstpr_L2dir)
            call cnht_opcmult (VXP,VYP,VZP,TP,grw)

#ifdef PRS_REC
            ! normalise pressure ???????????????????????????
            call cmult(PRP,grw,tstpr_np)
#else
            ! zero presure
            call rzero(PRP,tstpr_np)
#endif
            ! set cpfld for conjugated heat transfer
            if (IFHEAT) call cnht_cpfld_set
         else
            !stepper phase counting
            tstpr_istep = 0
            tstpr_vstep = tstpr_vstep +1

            call mntr_logi(tstpr_id,lp_prd,'Finished stepper phase:',
     $           tstpr_vstep)

            if (tstpr_mode.eq.3) then
               ! optimal initial condition
               call mntr_log(tstpr_id,lp_prd,
     $         'Optimal initial condition; rescaling solution')

               ! get L2 norm after direct phase
               tstpr_L2adj = cnht_glsc2_wt(VXP,VYP,VZP,TP,
     $                 VXP,VYP,VZP,TP,BM1)
               ! normalise vector after whole cycle
               grw = sqrt(tstpr_L2dir/tstpr_L2ini)! add direct growth
               call cnht_opcmult (VXP,VYP,VZP,TP,grw)
#ifdef PRS_REC
               ! normalise pressure ???????????????????????????
               call cmult(PRP,grw,tstpr_np)
#endif
            endif

            ! run vector solver (arpack, power iteration)
            call stepper_vsolve

            if (LASTEP.ne.1) then
               ! stepper restart;
               ! set time and iteration number
               TIME=0.0
               ISTEP=0
#ifndef PRS_REC
               ! zero presure ?????????????????????????????????????
               call rzero(PRP,tstpr_np)
#endif
               if (tstpr_mode.eq.3) then
                  ! optimal initial condition
                  call mntr_log(tstpr_id,lp_prd,
     $            'Optimal initial condition; direct phase start')

                  IFADJ = .FALSE.

                  ! get initial L2 norm
                  tstpr_L2ini = cnht_glsc2_wt(VXP,VYP,VZP,TP,
     $                    VXP,VYP,VZP,TP,BM1)
                  ! set cpfld for conjugated heat transfer
                  if (IFHEAT) call cnht_cpfld_set

               endif
            endif

         endif               ! tstpr_mode.eq.3.and.(.not.IFADJ)
         ! timing
         ltim = dnekclock() - ltim
         call mntr_tmr_add(tstpr_tmr_evl_id,1,ltim)
      endif                  ! mod(tstpr_istep,tstpr_step).eq.0

      return
      end subroutine
!=======================================================================
!> @brief Average velocity and temperature at element faces.
!! @ingroup tstpr
      subroutine tstpr_dssum
      implicit none

      include 'SIZE'            ! N[XYZ]1
      include 'INPUT'           ! IFHEAT
      include 'SOLN'            ! V[XYZ]P, TP, [VT]MULT
      include 'TSTEP'           ! IFIELD
      include 'TSTPRD'       ! tstpr_nt

      ! local variables
      integer ifield_tmp
!-----------------------------------------------------------------------
      ! make sure the velocity and temperature fields are continuous at
      ! element faces and edges
      ifield_tmp = IFIELD
      IFIELD = 1
#ifdef AMR
      call amr_oph1_proj(vxp,vyp,vzp,nx1,ny1,nz1,nelv)
#else
      call opdssum(vxp,vyp,vzp)
      call opcolv (vxp,vyp,vzp,vmult)
#endif

      if(IFHEAT) then
         IFIELD = 2
#ifdef AMR
         call h1_proj(tp,nx1,ny1,nz1)
#else
         call dssum(tp,nx1,ny1,nz1)
         call col2 (tp,tmult,tstpr_nt)
#endif
      endif
      IFIELD = ifield_tmp

      return
      end subroutine
!=======================================================================
