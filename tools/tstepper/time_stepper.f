!> @file time_stepper.f
!! @ingroup tstepper
!! @brief Set of subroutines to use time steppers for power iterations or
!!     solution of eigenvalue problem with Arnoldi algorithm
!!
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Read parameters for time stepper and call suitable param_reader
!! of required submodule
!! @ingroup tstepper
!! @param[in]  fid    file unit
!! @note This routine should be called in @ref runprm_in
!! @todo Check iostat value for missing namelist in the file
!! @see @ref readers_writers_page
      subroutine tstpr_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'TIME_STEPPERD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len
      integer itmp(3)

!     namelists
      namelist /TSTEPPER/ TSTMODE, TSTSTEP, TSTCMAX,
     $     TSTTOL, TSTIFUZAWA
!-----------------------------------------------------------------------
!     default values
      TSTMODE = 1
      TSTSTEP = 40
      TSTCMAX = 40
      TSTTOL = 1.0d-6
      TSTIFUZAWA = .TRUE.
!     read the file
      ierr=0
      if (NID.eq.0) then
        rewind(fid)
        read(unit=fid,nml=TSTEPPER,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading TSTEPPER parameters.$')

!     broadcast data
      if (NID.eq.0) then
         itmp(1) = TSTMODE
         itmp(2) = TSTSTEP
         itmp(3) = TSTCMAX
      endif
      len = 3*ISIZE
      call bcast(itmp,len)
      if (NID.ne.0) then
         TSTMODE = itmp(1)
         TSTSTEP = itmp(2)
         TSTCMAX = itmp(3)
      endif

      call bcast(TSTTOL,WDSIZE)
      call bcast(TSTIFUZAWA,LSIZE)

! call submodule runtime parametr reader
      call stepper_param_in(fid)

      return
      end
!=======================================================================
!> @brief Write parameters for time stepper and call suitable param_writer
!! of required submodule
!! @ingroup tstepper
!! @param[in]  fid    file unit
!! @note This routine should be called in @ref runprm_out
!! @see @ref readers_writers_page
      subroutine tstpr_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'TIME_STEPPERD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /TSTEPPER/ TSTMODE, TSTSTEP, TSTCMAX,
     $     TSTTOL, TSTIFUZAWA
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=TSTEPPER,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing TSTEPPER parameters.$')

! call required submodule runtime parameter writer
      call stepper_param_out(fid)

      return
      end
!=======================================================================
!> @brief Main time stepper interface
!! @ingroup tstepper
!! @note This routine should be called in userchk
      subroutine tstepper_main
      implicit none

      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP

!     local variables
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
      if(.not.ifcalled) then
        ifcalled=.TRUE.
        call tstpr_init
      endif

      if (ISTEP.gt.0) call tstpr_solve

      return
      end
!=======================================================================
!> @brief Initialise time stepper and call suitable stepper_vinit
!! of required submodule
!! @ingroup tstepper
      subroutine tstpr_init
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NSTEPS
      include 'INPUT_DEF'
      include 'INPUT'           ! IFPERT, IFBASE, IFTRAN, PARAM, 
                                ! IFHEAT, IFUZAWA
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, PRP
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'TIME_STEPPERD'

!     local variables

!     functions
      real dnekclock, cht_glsc2_wt
!-----------------------------------------------------------------------
!     set parameters
      IFTST = .TRUE.
      if (TSTMODE.eq.0) IFTST = .FALSE.

      if (IFTST) then
!     timing
         TSTIME1=dnekclock()

!     initialise related packages
!     conjuagate heat transfer init
         if (IFHEAT) call cht_init

!     check nek5000 parameters
         if (.not.IFTRAN) then
            if (NIO.eq.0) write(*,*)
     $          'ERROR: time stepper requres transient equations'
            call exitt
         endif

         if (NSTEPS.eq.0) then
            if (NIO.eq.0) write(*,*)
     $           'ERROR: time stepper requires NSTEPS>0'
            call exitt
         endif

         if (PARAM(12).ge.0) then
            if (NIO.eq.0) write(*,*)
     $           'ERROR: time stepper assumes constant dt'
            call exitt
         endif

         if (.not.IFPERT) then
            if (NIO.eq.0) write(*,*)
     $         'ERROR: time stepper has to be run in perturbation mode'
            call exitt
         endif

         if (IFBASE) then
            if (NIO.eq.0) write(*,*)
     $           'ERROR: time stepper assumes constatnt base flow'
            call exitt
         endif

         if (NPERT.ne.1) then
            if (NIO.eq.0) write(*,*)
     $           'ERROR: time stepper requires NPERT=1'
            call exitt
         endif

!     make sure NSTEPS is bigger than the possible number of iterations
!     in time stepper phase
         NSTEPS = max(NSTEPS,TSTSTEP*TSTCMAX+10)

!     initialise cycle counters
         TSTISTEP = 0
         TSTVSTEP = 0

!     for timing
         TSTIMET=0.0
         TSTIMES=0.0
         TSTIMESA=0.0

!     vector length
         NVECAV  = NX1*NY1*NZ1*NELV ! velocity single component
         if(IFHEAT) then        !temperature
            NVECAT  = NX1*NY1*NZ1*NELT
         else
            NVECAT  = 0
         endif
         NVECAP  = NX2*NY2*NZ2*NELV ! presure

!     print info
         if (NIO.eq.0) then
            write(*,*)
            write(*,*) 'TIME STEPPER initialised'
            if (TSTMODE.eq.1) then
               write(*,*) 'DIRECT mode'
            elseif (TSTMODE.eq.2) then
               write(*,*) 'ADJOINT mode'
            elseif (TSTMODE.eq.3) then
               write(*,*) 'Opt. init. cond.'
            endif
            write(*,*)
            write(*,*) 'Parameters:'

            write(*,'(A15,G13.5)') 'TOL = ',TSTTOL
            write(*,'(A15,I13)') 'Nstep = ',TSTSTEP
            write(*,'(A15,I13)') 'Ncmax = ',TSTCMAX
            write(*,*)
         endif

!     place to initialise vector solver (arpack, power iteration ...)
         call stepper_vinit

!     setup simulation parameters
!     set time and iteration number
         TIME=0.0

!     should be the first step of every cycle performed with Uzawa 
!     turned on?
         IFUZAWA = TSTIFUZAWA

!     zero presure
         call rzero(PRP,NVECAP)

         IFADJ = .FALSE.
         if (TSTMODE.eq.2) then
!     Is it adjoint mode
            IFADJ = .TRUE.
         elseif  (TSTMODE.eq.3) then
!     If it is optimal initial condition save initial L2 norm
            TSTL2INI = cht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)

            if (TSTL2INI.eq.0.0) then
               if (NIO.eq.0)
     $              write(*,*) 'ERROR; tstpr_init, TSTL2INI = 0'
               call exitt
            endif
         endif

!     set cpfld for conjugated heat transfer
         if (IFHEAT) call cht_cpfld_set

!     timing
         TSTIME2=dnekclock()
         TSTIMET=TSTIME2-TSTIME1

      endif                     ! IFTST

      if (TSTMODE.eq.3.and.NIO.eq.0) write(*,*)
     $     'TSTEPPER: opt. init. cond. direct phase start'

      return
      end
!=======================================================================
!> @brief Control time stepper after every nek5000 step and call suitable
!! stepper_vsolve of required submodule
!! @ingroup tstepper
      subroutine tstpr_solve
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, TIME
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D, IFUZAWA
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, PRP, TP, VMULT, V?MASK
      include 'ADJOINT_DEF'
      include 'ADJOINT'         ! IFADJ
      include 'TIME_STEPPERD'

!     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL

!     local variables
      real grw                  ! growth rate

!     functions
      real dnekclock, op_glsc2_wt, cht_glsc2_wt
!-----------------------------------------------------------------------
      if (IFTST) then

!     turn off Uzawa after first step
         IFUZAWA = .FALSE.

!     step counting
         TSTISTEP = TSTISTEP + 1

!     stepper phase end
         if (mod(TSTISTEP,TSTSTEP).eq.0) then

!     check for the calculation mode
            if (TSTMODE.eq.3.and.(.not.IFADJ)) then
!     optimal initial condition
!     stamp output
               if (NIO.eq.0) write(*,*)
     $         'TSTEPPER: opt. init. cond. adjoint phase start'

               IFADJ = .TRUE.

!     itaration count
               TSTISTEP = 0

!     should be the first step of every cycle performed with Uzawa 
!     turned on?
               IFUZAWA = TSTIFUZAWA

!     set time and iteration number
               TIME=0.0
               ISTEP=0

!     get L2 norm after direct phase
               TSTL2DIR = cht_glsc2_wt(VXP,VYP,VZP,TP,
     $              VXP,VYP,VZP,TP,BM1)
!     normalise vector
               grw = sqrt(TSTL2INI/TSTL2DIR)
               call cht_opcmult (VXP,VYP,VZP,TP,grw)

!     zero presure
               call rzero(PRP,NVECAP)

!     set cpfld for conjugated heat transfer
               if (IFHEAT) call cht_cpfld_set
            else
!     stepper phase counting
               TSTISTEP = 0
               TSTVSTEP = TSTVSTEP +1

!     timing
               TSTIME1=dnekclock()
!     average stepper phase length
               TSTIMES=TSTIME1-TSTIME2
               TSTIME2=TSTIME1
               TSTIMESA=(TSTIMESA*(TSTVSTEP-1)+TSTIMES)/real(TSTVSTEP)

!     stamp output
               if (NIO.eq.0) write(*,210) TSTVSTEP
 210           FORMAT('TSTEPPER: after',I5,' stepper phase')

               if (TSTMODE.eq.3) then
!     optimal initial condition
!     stamp output
                  if (NIO.eq.0) write(*,*)
     $                 'TSTEPPER: opt. init. cond. soln. rescaling'

!     get L2 norm after direct phase
                  TSTL2ADJ = cht_glsc2_wt(VXP,VYP,VZP,TP,
     $                 VXP,VYP,VZP,TP,BM1)

!     normalise vector after whole cycle
                  grw = sqrt(TSTL2DIR/TSTL2INI)! add direct growth
                  call cht_opcmult (VXP,VYP,VZP,TP,grw)

                  if (NIO.eq.0) then
                     write(*,*) 'Scaling factors:'
                     write(*,'(A15,G13.5)') 'TSTL2INI = ',TSTL2INI
                     write(*,'(A15,2G13.5)') 'TSTL2DIR = ',TSTL2DIR,
     $                    TSTL2DIR/TSTL2INI
                     write(*,'(A15,2G13.5)') 'TSTL2ADJ = ',TSTL2ADJ,
     $                    TSTL2ADJ/TSTL2INI
                  endif
               endif

!     place to run vector solver (arpack, power iteration)
               call stepper_vsolve

               if(LASTEP.eq.1) then
!     finalise stepper
!     timing
                  TSTIME3=dnekclock()
                  TSTIMET=TSTIMET+TSTIME3-TSTIME1

                  call tstpr_end
               else
!     stepper restart;
!     set time and iteration number
                  TIME=0.0
                  ISTEP=0

!     should be the first step of every cycle performed with Uzawa 
!     turned on?
                  IFUZAWA = TSTIFUZAWA

!     zero presure
                  call rzero(PRP,NVECAP)

                  if (TSTMODE.eq.3) then
!     optimal initial condition
!     stamp output
                     if (NIO.eq.0) write(*,*)
     $                'TSTEPPER: opt. init. cond. direct phase start'

                     IFADJ = .FALSE.

!     get initial L2 norm
                     TSTL2INI = cht_glsc2_wt(VXP,VYP,VZP,TP,
     $                    VXP,VYP,VZP,TP,BM1)
                     
!     set cpfld for conjugated heat transfer
                     if (IFHEAT) call cht_cpfld_set

                  endif
               endif

!     timing
               TSTIME3=dnekclock()
               TSTIMET=TSTIMET+TSTIME3-TSTIME1

            endif               ! TSTMODE.eq.3.and.(.not.IFADJ)
         endif                  ! mod(TSTISTEP,TSTSTEP).eq.0
      endif                     ! IFTST

      return
      end
!=======================================================================
!> @brief Finalise time stepper
!! @ingroup tstepper
      subroutine tstpr_end
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TSTEP_DEF'
      include 'TSTEP'           ! LASTEP
      include 'TIME_STEPPERD'
!-----------------------------------------------------------------------
      if (IFTST.and.(LASTEP.eq.1)) then
         if (NIO.eq.0) then
            write(*,*) ''
            write(*,*) 'TSTEPPER: finalize'
            write(*,*) '   Number of stepper cycles   ',TSTVSTEP
            write(*,*) '   Time spent in tstpr_solve  ',TSTIMET
            write(*,*) '   Average stepper phase time ',TSTIMESA
         endif
      endif

      return
      end
!=======================================================================
!> @brief Average velocity and temperature at element faces.
!! @ingroup tstepper
      subroutine tstpr_dssum
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! N[XYZ]1
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, [VT]MULT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! IFIELD
      include 'TIME_STEPPERD'   ! NVECAT

!     local variables
      integer ifield_tmp
!-----------------------------------------------------------------------
!     make sure the velocity and temperature fields are continuous at
!     element faces and edges
      ifield_tmp = IFIELD
      IFIELD = 1
      call opdssum(VXP,VYP,VZP)
      call opcolv (VXP,VYP,VZP,VMULT)
      if(IFHEAT) then
         IFIELD = 2
         call dssum(TP,NX1,NY1,NZ1)
         call col2 (TP,TMULT,NVECAT)
      endif
      IFIELD = ifield_tmp

      return
      end
!=======================================================================
