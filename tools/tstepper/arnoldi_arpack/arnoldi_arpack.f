!> @file arnoldi_arpack.f
!! @ingroup arnoldi_arpack
!! @brief Set of subroutines to solve eigenvalue problem with Arnoldi
!!     algorithm using PARPACK/ARPACK
!!
!! @author Adam Peplinski
!! @date Mar 7, 2016
!     To define ARPACK mode: direct or inverse
!#define ARPACK_DIRECT
#undef ARPACK_DIRECT
!
!     IMPORTANT!! No restart option for serial ARPACK version 
!     (only for PARPACK)
!=======================================================================
!> @brief Read runtime parameters for arnoldi_arpack
!! @ingroup arnoldi_arpack
!! @param[in]  fid    file unit
!! @note This interface is defined in @ref tstpr_param_in
!! @todo Check iostat value for missing namelist in the file
!! @see @ref readers_writers_page
      subroutine stepper_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'ARNOLDI_ARPACKD'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr, len
      integer itmp(4)

!     namelists
      namelist /ARPACK/ ARNKRYLOV, ARNEGV, ARNISTART, ARNISTOP
!-----------------------------------------------------------------------
!     default values
      ARNKRYLOV = 100
      ARNEGV = 10
      ARNISTART = 1
      ARNISTOP = 2
!     read the file
      ierr=0
      if (NID.eq.0) then
        rewind(fid)
        read(unit=fid,nml=ARPACK,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading ARPACK parameters.$')

!     broadcast data
      if (NID.eq.0) then
         itmp(1) = ARNKRYLOV
         itmp(2) = ARNEGV
         itmp(3) = ARNISTART
         itmp(4) = ARNISTOP
      endif
      len = 4*ISIZE
      call bcast(itmp,len)
      if (NID.ne.0) then
         ARNKRYLOV = itmp(1)
         ARNEGV = itmp(2)
         ARNISTART = itmp(3)
         ARNISTOP = itmp(4)
      endif

      return
      end
!=======================================================================
!> @brief Write runtime parameters for arnoldi_arpack
!! @ingroup arnoldi_arpack
!! @param[in]  fid    file unit
!! @note This interface is defined in @ref tstpr_param_out
!! @see @ref readers_writers_page
      subroutine stepper_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'      
      include 'ARNOLDI_ARPACKD'         

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /ARPACK/ ARNKRYLOV, ARNEGV, ARNISTART, ARNISTOP
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=ARPACK,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing ARPACK parameters.$')

      return
      end
!=======================================================================
!> @brief Initialise arnoldi_arpack
!! @ingroup arnoldi_arpack
!! @note This interface is defined in @ref tstpr_init
      subroutine stepper_vinit
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NSTEPS
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN_DEF'
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'CHKPOINTD'        ! IFCHKPTRST
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     ARPACK include file
      INCLUDE 'debug.h'

!     local variables
      integer il

!     functions
      real dnekclock
!-----------------------------------------------------------------------
!     set parameters
!     timing
      ARPTIME1=dnekclock()

!     check nek5000 parameters
!     Simulation with temperature or passive scalars has to be performed
!     in inverse mode due to speciffic inner product.
#ifdef ARPACK_DIRECT
!     standard eigenvalue problem A*x = lambda*x
      if(IFHEAT) then
         if (NIO.eq.0) write(*,*)
     $        'ERROR: IFHEAT requires #undef ARPACK_DIRECT'
         call exitt
      endif
#endif

      if (ARNKRYLOV.gt.LDIMA) then
         if (NIO.eq.0) write(*,*) 'ERROR: ARNKRYLOV too big'
         call exitt
      endif

      if (ARNEGV.ge.(ARNKRYLOV/2)) then
         if (NIO.eq.0) write(*,*)'ERROR: ARNEGV is >ARNKRYLOV/2'
         call exitt
      endif

!     make sure NSTEPS is bigger than the possible number of iteration in arnoldi
      NSTEPS = max(NSTEPS,TSTSTEP*ARNKRYLOV*TSTCMAX+10)

!     for timing
      ARPTIMET=0.0

!     related to restart
      NPARP = 0
      NCARP = 0
      RNMARP= 0.0

!     initialize ARPACK parameters
#ifdef ARPACK_DIRECT
!     standard eigenvalue problem A*x = lambda*x
      BMATARP='I'
#else
!     generalized eigenvalue problem A*x = lambda*B*x
      BMATARP='G'
#endif
!     eigenvalues of largest magnitude
      WHICHARP='LM'

      call izero(IPARP,11)
      call izero(IPNTARP,14)
!     exact shifts with respect to the current Hessenberg matrix
      IPARP(1)=1
!     maximum number of Arnoldi update iterations allowed
      IPARP(3)=TSTCMAX
#ifdef ARPACK_DIRECT
!     A*x = lambda*x
      IPARP(7)=1
#else
!     A*x = lambda*M*x, M symmetric positive definite; BMATARP='G'
      IPARP(7)=2
#endif
!     used size of WORKLA
      NWLARP = (3*ARNKRYLOV+6)*ARNKRYLOV

!     user supplied initial conditions
      INFARP=1

!     get eigenvectors
      RVARP=.true.
!     compute Ritz vectors
      HOWARP='A'
!     select should be specifird for HOWARP='S'

!     no shift
      SIGARP(1) = 0.0
      SIGARP(2) = 0.0

!     vector length

!     single vector length in Krylov space
!     velocity
      NVECAS = NVECAV*NDIM
!     temperature
      if(IFHEAT) then
         NVECAS = NVECAS + NVECAT
      endif
      if (NVECAS.gt.LVAS) then
         if (NIO.eq.0) then
            write(*,*) 'ERROR: NVECAS too big'
            write(*,*) 'NVECAS = ', NVECAS
            write(*,*) 'LVAS   = ', LVAS
         endif
         call exitt
      endif

!     initialize arrays
      call rzero(WORKDA,WDDIMA)
      call rzero(WORKLA,WLDIMA)
      call rzero(WORKEA,WEDIMA)
      call rzero(VBASEA,LVAS*LDIMA)
      call rzero(RESIDA,LVAS)
      call rzero(DRIARP,LDIMA*4)

!     info level from ARPACK
      ndigit = -3
      logfil = 6
      mngets = 0
      mnaitr = 2
      mnapps = 0
      mnaupd = 2
      mnaup2 = 2
      mneupd = 0

!     PLACE FOR RESTART
      if (IFCHKPTRST) then
!     read checkpoint
         call arn_rst_read
      else
!     if no restatrt fill RESIDA with initial conditions
!     V?MASK removes points at the wall and inflow
#ifdef ARPACK_DIRECT
!     A*x = lambda*x
!     velocity
         call col3(RESIDA(1),VXP,V1MASK,NVECAV)
         call col3(RESIDA(1+NVECAV),VYP,V2MASK,NVECAV)
         if (IF3D) call col3(RESIDA(1+2*NVECAV),VZP,V3MASK,NVECAV)
!     no temperature here
#else
!     A*x = lambda*M*x
!     velocity
         call copy(RESIDA(1),VXP,NVECAV)
         call copy(RESIDA(1+NVECAV),VYP,NVECAV)
         if (IF3D) call copy(RESIDA(1+2*NVECAV),VZP,NVECAV)
!     temperature
         if(IFHEAT) call copy(RESIDA(1+NDIM*NVECAV),TP,NVECAT)
#endif

!     initialize rest of variables
!     firts call
         IDOARP=0
      endif

!     ARPACK interface
      call arn_naupd

!     we should start stepper here
      if (IDOARP.ne.-1.and.IDOARP.ne.1) then
         if (NIO.eq.0) then
            write(*,*) 'ARNOLDI_ARPACK: '
            write(*,*) ' Error with arn_naupd, IDOARP = ', IDOARP
         endif
         call exitt
      endif
         
!     print info
      if (NIO.eq.0) then
         write(*,*)
         write(*,*) 'ARPACK initialised'
         write(*,*) 'Parameters:'
         write(*,'(A15,A1)') 'BMAT = ',BMATARP
         write(*,'(A15,A2)') 'WHICH = ',WHICHARP
         write(*,'(A15,G13.5)') 'TOL = ',TSTTOL
         write(*,'(A15,I13)') 'NEV = ',ARNEGV
         write(*,'(A15,I13)') 'NCV = ',ARNKRYLOV
         write(*,'(A15,I13)') 'IPARAM(1) = ',IPARP(1)
         write(*,'(A15,I13)') 'IPARAM(3) = ',IPARP(3)
         write(*,'(A15,I13)') 'IPARAM(7) = ',IPARP(7)
         write(*,'(A15,L)') 'RVEC = ',RVARP
         write(*,'(A15,A1)') 'HOWMNY = ',HOWARP
         if (HOWARP.eq.'S') then
            do il=1,LDIMA
               write(*,100) il,SELARP(il)
            enddo
         endif
      endif

!     timing
      ARPTIME2=dnekclock()
      ARPTIMET=ARPTIME2-ARPTIME1

 100  FORMAT('  SELECT(',I2,') = ',L)

      return
      end
!=======================================================================
!> @brief Create Krylov space, get Ritz values and restart
!!  stepper phase.
!! @ingroup arnoldi_arpack
!! @note This interface is defined in @ref tstpr_solve
      subroutine stepper_vsolve
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'INPUT_DEF'
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN_DEF'
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !

!     functions
      real dnekclock
!-----------------------------------------------------------------------
!     timing
      ARPTIME1=dnekclock()

!     fill work array with velocity
!     V?MASK removes points at the boundary
#ifdef ARPACK_DIRECT
!     A*x = lambda*x
!     velocity
      call col3(WORKDA(IPNTARP(2)),VXP,V1MASK,NVECAV)
      call col3(WORKDA(IPNTARP(2)+NVECAV),VYP,V2MASK,NVECAV)
      if (IF3D)
     $     call col3(WORKDA(IPNTARP(2)+2*NVECAV),VZP,V3MASK,NVECAV)
!     no temperature here
#else
!     velocity
!     A*x = lambda*M*x
      call copy(WORKDA(IPNTARP(2)),VXP,NVECAV)
      call copy(WORKDA(IPNTARP(2)+NVECAV),VYP,NVECAV)
      if (IF3D) call copy(WORKDA(IPNTARP(2)+2*NVECAV),VZP,NVECAV)
!     temperature
      if(IFHEAT) call copy(WORKDA(IPNTARP(2)+NDIM*NVECAV),TP,NVECAT)
#endif

!     ARPACK interface
      call arn_naupd

      if (IDOARP.eq.-2) then
!     checkpoint
         call arn_rst_save
         call arn_end
      elseif (IDOARP.eq.99) then
!     finalise
         call arn_esolve
         call arn_end
      elseif (IDOARP.eq.-1.or.IDOARP.eq.1) then
!     stepper restart, nothing to do
      else
         if (NIO.eq.0) then
            write(*,*) 'ARNOLDI_ARPACK: '
            write(*,*) ' Error with arn_naupd, IDOARP = ', IDOARP
         endif
         call exitt
      endif

!     timing
      ARPTIME2=dnekclock()
      ARPTIMET=ARPTIMET+ARPTIME2-ARPTIME1

      return
      end
!=======================================================================
!> @brief ARPACK postprocessing
!! @ingroup arnoldi_arpack
      subroutine arn_esolve
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NID, LDIMT1
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, DT, LASTEP
      include 'SOLN_DEF'
      include 'SOLN'            ! VX, VY, VZ, VMULT, V?MASK
      include 'INPUT_DEF'
      include 'INPUT'           ! IFXYO,IFPO,IFVO,IFTO,IFPSO,IF3D,IFHEAT
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !

!     local variables
      integer il, iunit, ierror
      real dumm
      logical lifxyo, lifpo, lifvo, lifto, lifpso(LDIMT1)
!     functions
      real dnekclock

!     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
!-----------------------------------------------------------------------
      if (IDOARP.eq.99) then

!     timing
         ARPTIME1=dnekclock()

         if (NIO.eq.0) write(*,405) IPARP(5)
 405     FORMAT('ARNOLDI_ARPACK: ',I4,' eigenvectors converged;
     $postprocessing')

#ifdef MPI
         call pdneupd(NEKCOMM,RVARP,HOWARP,SELARP,DRIARP,DRIARP(1,2),
     $        VBASEA,LVAS,SIGARP(1),SIGARP(2),WORKEA,BMATARP,NVECAS,
     $        WHICHARP,ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,
     $        IPARP,IPNTARP,WORKDA,WORKLA,NWLARP,IERRARP)
#else
         call dneupd(RVARP,HOWARP,SELARP,DRIARP,DRIARP(1,2),
     $        VBASEA,LVAS,SIGARP(1),SIGARP(2),WORKEA,BMATARP,NVECAS,
     $        WHICHARP,ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,
     $        IPARP,IPNTARP,WORKDA,WORKLA,NWLARP,IERRARP)
#endif

!     timing
         ARPTIME2=dnekclock()
         ARPTIMET=ARPTIMET+ARPTIME2-ARPTIME1

         if (IERRARP.eq.0) then

            if (NIO.eq.0) then
               write(*,*) 'ARNOLDI_ARPACK:'
               write(*,*) 'writing eigenvalues and eigenvectors.'
            endif

            ierror=0
!     open file 
            if (NID.eq.0) then
!     find free unit
               call IO_file_freeid(iunit, ierror)
               if (ierror.eq.0) then
                  open (unit=iunit,file='eigenvalues.txt',
     $                 action='write', iostat=ierror)
                  write(unit=iunit,FMT=410,iostat=ierror)
 410              FORMAT(10x,'I',17x,'re(RITZ)',17x,'im(RITZ)',17x,
     $                 'ln|RITZ|',16x,'arg(RITZ)')

               endif
            endif
!     error check
            call err_chk(ierror,'Error opening eigenv. file.$')

!     integration time
            dumm = DT*TSTSTEP
            dumm = 1.0/dumm

!     copy and set output parameters
            lifxyo = IFXYO
            IFXYO = .TRUE.
            lifpo= IFPO
            IFPO = .false.
            lifvo= IFVO
            IFVO = .true.
            lifto= IFTO
            if (IFHEAT) then
               IFTO = .TRUE.
            else
               IFTO = .FALSE.
            endif
            do il=1,LDIMT1
               lifpso(il)= IFPSO(il)
               IFPSO(il) = .false.
            enddo

!     we have to take into account storrage of imaginary and real
!     parts of eigenvectors in arpack.
!     The complex Ritz vector associated with the Ritz value
!     with positive imaginary part is stored in two consecutive
!     columns.  The first column holds the real part of the Ritz
!     vector and the second column holds the imaginary part.  The
!     Ritz vector associated with the Ritz value with negative
!     imaginary part is simply the complex conjugate of the Ritz
!     vector associated with the positive imaginary part.

            ierror=0
            do il=1,IPARP(5)
!     copy eigenvectors to perturbation variables
               call copy(VXP,VBASEA(1,il),NVECAV)
               call copy(VYP,VBASEA(1+NVECAV,il),NVECAV)
               if (IF3D) call copy(VZP,VBASEA(1+2*NVECAV,il),NVECAV)
               if(IFHEAT) then
                  call copy(TP,VBASEA(1+NDIM*NVECAV,il),NVECAT)
                  call outpost2(VXP,VYP,VZP,PRP,TP,1,'egv')
               else
                  call outpost2(VXP,VYP,VZP,PRP,TP,0,'egv')
               endif
!     possible place to test error

!     get growth rate; get eigenvalues of continuous operator
               DRIARP(il,3) = log(sqrt(DRIARP(il,1)**2+
     $              DRIARP(il,2)**2))*dumm
               DRIARP(il,4) = atan2(DRIARP(il,2),DRIARP(il,1))*dumm

               if (NID.eq.0)  write(unit=iunit,FMT=*,iostat=ierror)
     $           il,DRIARP(il,1),DRIARP(il,2),DRIARP(il,3),DRIARP(il,4)
            enddo
!     error check
            call err_chk(ierror,'Error witing to eigenv. file.$')

!     put output variables back
            IFXYO = lifxyo
            IFPO = lifpo
            IFVO = lifvo
            IFTO = lifto
            do il=1,LDIMT1
               IFPSO(il) = lifpso(il)
            enddo

!     close eigenvalue file
            if (NID.eq.0)  close(unit=iunit)

         else                   ! IERRARP
            if (NIO.eq.0) then
               write(*,*) 'ARNOLDI_ARPACK:'
               write(*,*) ' Error with _neupd, info = ', IERRARP
               write(*,*) ' Check the documentation of _naupd.'
            endif
            call exitt
         endif                  ! IERRARP

!     finish run
         LASTEP=1

      endif

      return
      end
!=======================================================================
!> @brief Finalise arnoldi
!! @ingroup arnoldi_arpack
      subroutine arn_end
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID, NDIM, NPERT
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP, NSTEPS, LASTEP
      include 'ARNOLDI_ARPACKD'
!-----------------------------------------------------------------------
      if ((LASTEP.eq.1)) then
         if (NIO.eq.0) then
            write(*,*) ''
            write(*,*) 'ARPACK: finalize'
            write(*,*) '   Time spent in ARPACK      ',ARPTIMET
         endif
      endif

      return
      end
!=======================================================================
!> @brief Interface to pdnaupd
!! @ingroup arnoldi_arpack
      subroutine arn_naupd
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO, NDIM, N[XYZ]1
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN_DEF'
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'MASS_DEF'
      include 'MASS'            ! BM1
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !

!     global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
!-----------------------------------------------------------------------
#ifdef MPI
      call pdnaupd(NEKCOMM,IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $     TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,WORKDA,
     $     WORKLA,NWLARP,INFARP,NPARP,RNMARP,NCARP)
#else
      call dnaupd(IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $     TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,WORKDA,
     $     WORKLA,NWLARP,INFARP)
#endif

!     error check
      if (INFARP.lt.0) then
         if (NIO.eq.0) then
            write(*,*) 'ARNOLDI_ARPACK: '
            write(*,*) ' Error with _naupd, info = ', INFARP
            write(*,*) ' Check the documentation of _naupd.'
         endif
         call exitt
      endif

      if (IDOARP.eq.2) then
         do
!     A*x = lambda*M*x
!     multiply by weights and masks
!     velocity
            call col3(WORKDA(IPNTARP(2)),BM1,V1MASK,NVECAV)
            call col3(WORKDA(IPNTARP(2)+NVECAV),BM1,V2MASK,NVECAV)
            if (IF3D) call col3(WORKDA(IPNTARP(2)+2*NVECAV),
     $           BM1,V3MASK,NVECAV)

!     temperature
            if(IFHEAT) then
               call col3(WORKDA(IPNTARP(2)+NDIM*NVECAV),
     $              BM1,TMASK,NVECAT)

!     coefficients
               call cht_weight_fun (WORKDA(IPNTARP(2)),
     $              WORKDA(IPNTARP(2)+NVECAV),
     $              WORKDA(IPNTARP(2)+2*NVECAV),
     $              WORKDA(IPNTARP(2)+NDIM*NVECAV),1.0)
            endif

            call col2(WORKDA(IPNTARP(2)),WORKDA(IPNTARP(1)),NVECAS)

#ifdef MPI
            call pdnaupd(NEKCOMM,IDOARP,BMATARP,NVECAS,WHICHARP,
     $           ARNEGV,TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,
     $           IPNTARP,WORKDA,WORKLA,NWLARP,INFARP,NPARP,RNMARP,
     $           NCARP)
#else
            call dnaupd(IDOARP,BMATARP,NVECAS,WHICHARP,ARNEGV,
     $           TSTTOL,RESIDA,ARNKRYLOV,VBASEA,LVAS,IPARP,IPNTARP,
     $           WORKDA,WORKLA,NWLARP,INFARP)
#endif

!     error check
            if (INFARP.lt.0) then
               if (NIO.eq.0) then
                  write(*,*) 'ARNOLDI_ARPACK: '
                  write(*,*) ' Error with _naupd, info = ', INFARP
                  write(*,*) ' Check the documentation of _naupd.'
               endif
               call exitt
            endif
            if (IDOARP.ne.2) exit
         enddo
      endif                     ! IDOARP.eq.2

!     restart stepper
      if (IDOARP.eq.-1.or.IDOARP.eq.1) then

         if (NIO.eq.0) write(*,*) 'ARNOLDI_ARPACK: restarting stepper'

!     move renormed data back to nekton
!     velocity
         call copy(VXP,WORKDA(IPNTARP(1)),NVECAV)
         call copy(VYP,WORKDA(IPNTARP(1)+NVECAV),NVECAV)
         if (IF3D) call copy(VZP,WORKDA(IPNTARP(1)+2*NVECAV),NVECAV)
!     temperature
         if(IFHEAT) call copy(TP,WORKDA(IPNTARP(1)+NDIM*NVECAV),NVECAT)

!     make sure the velocity and temperature fields are continuous at
!     element faces and edges
         call tstpr_dssum
      endif                     ! IDOARP.eq.-1.or.IDOARP.eq.1

      return
      end
!=======================================================================

!=======================================================================
