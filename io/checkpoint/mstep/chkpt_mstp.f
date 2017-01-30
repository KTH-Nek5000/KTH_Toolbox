!> @file chkpt_mstp.f
!! @ingroup chkpoint_mstep
!! @brief Set of multi-file checkpoint routines for DNS and  linear
!! @details VERSION FOR DNS AND LINEAR SIMULATIONS
!!    This version supports only proper restart for DNS and for
!!    perturbation mode with single(!!!!!) NPERT=1
!!    perturbation and not advected basefield IFBASE=F. For DNS only
!!    'rs8' files are written. For perturbation mode  files 'rs8'
!!    include perturbation data and the basefield is saved in 'rb8'.
!!
!=======================================================================
!> @brief Initialise multi-file checkpoint routines
!! @ingroup chkpoint_mstep
!! @note This interface is defined in @ref checkpoint_main
      subroutine chkpt_init
      implicit none

      include 'SIZE'            ! NID, NPERT
      include 'TSTEP'           ! IOSTEP, ISTEP
      include 'INPUT'           ! SESSION, IFPERT, IFBASE
      include 'PARALLEL'        ! ISIZE
      include 'CHKPOINTD'       ! CHKPTSTEP, IFCHKPTRST
      include 'CHKPTMSTPD'      ! CHKPTSET_O, CHKPTSET_I, CHKPTNRSF,
                                ! CHKPTNFILE, CHKPTNSET

!     local variables
      integer il, ierr, itmp

      character*132 fname, bname

      character*3 prefix

      integer iunit
!-----------------------------------------------------------------------
!     this should be done only once

!     check checkpoint frequency
      if (CHKPTSTEP.le.0) then
         CHKPTSTEP = NSTEPS-CHKPTNRSF+1
         if (NIO.eq.0) write (*,*) 'Warning; checkpoint_init: ',
     $           ' CHKPTSTEP = 0; resetting to ', CHKPTSTEP
      elseif(mod(NSTEPS,CHKPTSTEP).ne.(CHKPTNRSF-1)) then
         if (NIO.eq.0) write (*,*) 'Warning; checkpoint_init: ',
     $           ' CHKPTSTEP and NSTEPS not optimal'
      endif

!     set negetive value of chkptset_i to mark that restart was not 
!     initialised yet
      CHKPTSET_I = -1

!     check perturbation parameters
      if (IFPERT) then
         if (IFBASE.or.NPERT.gt.1.or.IFMHD) then
            if (NIO.eq.0) write(*,*) 'CHECKPOINT: not supported mode'
            call exitt
         endif
      else                   ! IFPERT
         if (IFMHD) then
            if(NIO.eq.0)
     $           write(*,*) 'CHECKPOINT: not supported mode'
            call exitt
         endif
      endif                  ! IFPERT

!     create chkptrstf name (SESSION.restart)
      bname = trim(adjustl(SESSION))
      fname = trim(bname)
      fname = trim(fname)//'.restart'
      if (len_trim(fname).gt.80) then
         if(NIO.eq.0) write(*,*)
     $      'ERROR; checkpoint: too long restart file name'
         call exitt
      endif
      chkptrstf = trim(fname)

!     create names acording to mfo_open_files
      if (IFCHKPTRST) then

!     create names of restart files
!     get set number from the file SESSION.restart
         ierr=0
         if (NID.eq.0) then
!     find free unit
            call io_file_freeid(iunit, ierr)
!     get file name and open file
            if (ierr.eq.0) then
               open (unit=iunit,file=chkptrstf,status='old',
     $               action='read',iostat=ierr)
               if(ierr.eq.0) read(unit=iunit,fmt=*,iostat=ierr) itmp
               close(unit=iunit,iostat=ierr)
            endif
            if (ierr.eq.0.and.(itmp.lt.0.or.itmp.ge.CHKPTNSET)) ierr=1

         endif
         call err_chk(ierr,'Error reading .restart file.$')

         call bcast(itmp ,ISIZE)
         CHKPTSET_I = itmp

!     create restart file names
!     fill chkptfname array with 'rb8' and 'rs' file names
         if (IFPERT) then
!     assumes only single perturbation
            CHKPTNFILE = 2

!     prefix and name for base flow
            prefix(1:2)='rb'
            call chkpt_fname(CHKPTFNAME(1,1), bname, prefix,
     $           CHKPTSET_I, ierr)

            if (ierr.ne.0) then
               if (NIO.eq.0) write(*,*)
     $            'ERROR; checkpoint: base flow file name error'
               call exitt
            endif

!     for now IFBASE=F only
            do il=2,CHKPTNRSF
               CHKPTFNAME(il,1) = CHKPTFNAME(1,1)
            enddo

!     create prefix and name for perturbation
            prefix(1:2)='rs'
            call chkpt_fname(CHKPTFNAME(1,2), bname, prefix,
     $           CHKPTSET_I, ierr)

            if (ierr.ne.0) then
               if (NIO.eq.0) write(*,*)
     $            'ERROR; checkpoint: perturbation file name error'
               call exitt
            endif

         else                ! DNS
            CHKPTNFILE = 1

!     create prefix and name for DNS
            prefix(1:2)='rs'
            call chkpt_fname(CHKPTFNAME(1,1), bname, prefix,
     $           CHKPTSET_I, ierr)

            if (ierr.ne.0) then
               if (NIO.eq.0) write(*,*)
     $            'ERROR: checkpoint: DNS file name error'
               call exitt
            endif

         endif

      endif                  ! IFCHKPTRST

!     set initial file set number; it is different from chkptset_i 
!     because nek5000 always starts from 1
      CHKPTSET_O = 1

      return
      end
!=======================================================================
!> @brief Generate restart file names
!! @ingroup chkpoint_mstep
!! @param[out] fname  restart file name
!! @param[in]  bname  base name
!! @param[in]  prefix prefix
!! @param[in]  nset   file set number
!! @param[out] ierr   error mark
      subroutine chkpt_fname(fname, bname, prefix, nset, ierr)
      implicit none

      include 'SIZE'            ! NIO
      include 'CHKPTMSTPD'      ! CHKPTNRSF, CHKPTNSET

!     argument list
      character*80 fname(CHKPTNRSF)
      character*132 bname
      character*3 prefix
      integer nset, ierr

!     local variables
      character*132 fnamel   ! local file name
      character*3 prefixl    ! local prefix
      integer itmp

      character*6  str

      character*17 kst
      save         kst
      data         kst / '0123456789abcdefx' /
!-----------------------------------------------------------------------
!     create prefix and name for DNS
      ierr = 0
      prefixl(1:2) = prefix(1:2)
      itmp=min(17,CHKPTNSET*CHKPTNRSF) + 1
      prefixl(3:3)=kst(itmp:itmp)
      call io_mfo_fname(fnamel,bname,prefixl,ierr)
      if (ierr.ne.0) then
         if (NIO.eq.0) write(*,*) 'ERROR: checkpoint: file name error'
         return
      endif

!     is fname too long?
      if (len_trim(fnamel).gt.(80-5)) then
         if (NIO.eq.0)
     $      write(6,*) 'ERROR: checkpoint: too long file name'
         return
      endif

      do itmp=1,CHKPTNRSF
         write(str,'(i5.5)') CHKPTNRSF*nset+itmp
         fname(itmp) = trim(fnamel)//trim(str)//char(0)
      enddo

      return
      end
!=======================================================================
!> @brief Read full file restart set.
!! @ingroup chkpoint_mstep
!! @note This interface is defined in @ref checkpoint_main
      subroutine chkpt_read()
      implicit none

      include 'SIZE'            !
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! IFREGUO
      include 'CHKPOINTD'       ! IFCHKPTRST
      include 'CHKPTMSTPD'      ! CHKPTNRSF, CHKPTNFILE, CHKPTFNAME

!     local variables
      logical lifreguo
      integer ifile, i
      real p67
!-----------------------------------------------------------------------
!     no regular mesh
      lifreguo= IFREGUO
      IFREGUO = .false.

      if (IFCHKPTRST.and.(ISTEP.lt.CHKPTNRSF)) then

         ifile = ISTEP+1  ! istep=0,1,...

         p67 = param(67)
         param(67) = 6.00
         do i=1,CHKPTNFILE
            call chcopy (initc(i),CHKPTFNAME(ifile,i),80)
         enddo
         call bcast  (initc,132*CHKPTNFILE)

         if (IFPERT) then       ! perturbation
!     for IFBASE=F basefield loaded during istep=0 only
!     perturbation from 'rs' files in all steps
            if (IFBASE) then
               call restart_pert(CHKPTNFILE,1)
            elseif (ifile.eq.1) then
               call restart_pert(CHKPTNFILE,1)
            else
               call restart_pert(CHKPTNFILE,2)
            endif
         elseif (IFMHD) then    ! MHD
            call restart(2)
         else                   ! DNS
            call restart(1)
         endif                  ! IFPERT

         param(67)=p67

      endif

!     put parameters back
      IFREGUO = lifreguo

      return
      end
!=======================================================================
!> @brief Write full file restart set
!! @ingroup chkpoint_mstep
!! @note This interface is defined in @ref checkpoint_main.
!! @note This is version of full_restart_save routine. In additional
!! file it saves the file set number to allow for automatic restart.
      subroutine chkpt_write()
      implicit none

      include 'SIZE'            ! NID
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! IFREGUO
      include 'CHKPOINTD'       ! CHKPTSTEP
      include 'CHKPTMSTPD'      ! CHKPTNRSF, CHKPTSET_O, CHKPTNSET

!     local variables
      integer ierr, iotest, iunit
      integer save_size

      logical lifreguo
!-----------------------------------------------------------------------
!     no regular mesh
      lifreguo= IFREGUO
      IFREGUO = .false.

!     save rs.. rb.. files
!     there is some problem with nfld_save; sometimes it would be good
!     to increase its value, but restart_nfld doesn't allow for that
      save_size=8  ! For checkpoint
      call restart_save_pert(CHKPTSTEP,save_size,chkptnrsf)

!     save file set in SESSION.restart
!     we do it after the last file in the set was sucsesfully written
      ierr=0
      if(NID.eq.0) then
         iotest = 0
         if (ISTEP.gt.CHKPTSTEP/2.and.
     $     mod(ISTEP+CHKPTSTEP-iotest,CHKPTSTEP).eq.(CHKPTNRSF-1)) then

            CHKPTSET_O = mod(CHKPTSET_O+1,CHKPTNSET)

!     find free unit
            call io_file_freeid(iunit, ierr)
!     get file name and open file
            if(ierr.eq.0) then
               open (unit=iunit,file=CHKPTRSTF,
     $              action='write',iostat=ierr)
               if(ierr.eq.0)
     $              write(unit=iunit,fmt=*,iostat=ierr) CHKPTSET_O
               close(unit=iunit,iostat=ierr)
            endif
         endif                  ! ISTEP
      endif                     ! NID
      call err_chk(ierr,'Error writing to .restart file.$')

!     put parameters back
      IFREGUO = lifreguo

      return
      end
!=======================================================================
!     VERSION FOR PERTURBATION MODE
!     following two subroutines are modiffications of
!     full_restart_save
!     restart_save
!     from prepost.f.
!     This version supports only proper restart for single(!!!!!)
!     perturbation with not advected basefield IFBASE=F. In this case
!     files rs8 include perturbation and the basefield is saved in rb8.
!     The more general case (npert/=1 and IFBASE/=F) requeres increse
!     of file numbers to be written and e.g. hanges in restart_nfld.
!=======================================================================
!> @brief Save checkpoint variables in a single set of files.
!! @ingroup chkpoint_mstep
!! @param[in] iosave
!! @param[in] save_size
!! @param[in] nfldi
!! @note This is version of restart routine. It
      subroutine restart_save_pert(iosave,save_size,nfldi)
      implicit none

      include 'SIZE'
      include 'RESTART'
      include 'TSTEP'
      include 'INPUT'
      include 'SOLN'

!     argument list
      integer iosave,save_size,nfldi

!     local variables
      character*3 prefix, bprefix

      character*17 kst
      save         kst
      data         kst / '0123456789abcdefx' /

      logical if_full_pres_tmp, ifxyo_tmp

      integer icalld
      save    icalld
      data    icalld  /0/

      integer i2, iosav, iotest, iwdsizo, m1, mfld, mt
      integer nfld, nfld2, npscal1
      real p66
!-----------------------------------------------------------------------
      iosav = iosave

      if (iosav.eq.0) iosav = iostep
      if (iosav.eq.0) return

      iotest = 0
!     if (iosav.eq.iostep) iotest = 1  ! currently spoiled because of
!                                      ! incompatible format of .fld
!                                      ! and multi-file i/o;  the latter
!                                      ! is the only form used for restart

      nfld  = nfldi*2
      nfld2 = nfld/2
      mfld  = min(17,nfld)
!     this is not supported; only second order time integration for MHD?
!     Why only 2 files per set not 3?
!      if (ifmhd) nfld2 = nfld/4

      i2 = iosav/2
      m1 = istep+iosav-iotest
      mt = mod(istep+iosav-iotest,iosav)
      prefix = '   '
      bprefix = '   '

      if (istep.gt.iosav/2  .and.
     $   mod(istep+iosav-iotest,iosav).lt.nfld2) then ! save
         prefix(1:2)='rs'
         prefix(3:3)=kst(mfld+1:mfld+1)
         bprefix(1:2)='rb'
         bprefix(3:3)=kst(mfld+1:mfld+1)

         iwdsizo = wdsizo
         wdsizo  = save_size
         p66 = param(66)
         param(66) = 6          ! force multi-file out
         ifxyo_tmp = IFXYO
         IFXYO = .TRUE.         ! force writing coordinates

         npscal1 = npscal+1
         if (.not.ifheat) npscal1 = 0

         if_full_pres_tmp = if_full_pres
         if (save_size.eq.8) if_full_pres = .true. !Preserve mesh 2 pressure

         if (IFPERT) then

!     save basefiled
!     do this only once
            if (icalld.eq.0) then
               icalld=1
               call outpost2(vx,vy,vz,pr,t,npscal1,bprefix)
            endif

!     save perturbation
            call outpost2(vxp(1,1),vyp(1,1),vzp(1,1),prp(1,1),tp(1,1,1)
     $           ,npscal1,prefix) ! perturbation
         else                   ! DNS
            call outpost2(vx,vy,vz,pr,t,npscal1,prefix) ! basefield
         endif

         wdsizo    = iwdsizo  ! Restore output parameters
         param(66) = p66
         IFXYO = ifxyo_tmp
         if_full_pres = if_full_pres_tmp

      endif

      if_full_pres = .false.
      return
      end
!=======================================================================
!> @brief Load checpoint variables from single set of files.
!! @ingroup chkpoint_mstep
!! @param[in] nfiles     last file number
!! @param[in] nfilstart  first file number
!! @note This is version of restart routine. It excludes .fld format
!! and adds lower loop bound to allow for ifbase=.false.
      subroutine restart_pert(nfiles,nfilstart)
      implicit none

      include 'SIZE'            ! NID
      include 'TSTEP'           ! TIME
      include 'INPUT'           ! PARAM
      include 'RESTART'

!     argument list
      integer nfiles, nfilstart

!     local variables
      logical ifbasefl

      integer ifile,ndumps

      character*132 fname
!     functions
      real glmax
!-----------------------------------------------------------------------
      if(nfiles.lt.1) return
      if(NID.eq.0) write(6,*) 'Reading checkpoint data'

! use new reader (only binary support)
      if (PARAM(67).eq.6.0) then
         do ifile=nfilstart,nfiles
            call sioflag(ndumps,fname,initc(ifile))
            call mfi_pert(fname,ifile)
         enddo
         call setup_convect(3)
         if (nid.ne.0) time=0
         time = glmax(time,1) ! Sync time across processors
         return
      else
         if (NIO.eq.0) write(6,*) 'RESTART_PERT supporst .f only'
         call exitt
      endif

      return
      end
!=======================================================================
!> @brief Load checkpoint variables from single file.
!! @ingroup chkpoint_mstep
!! @param[in] fname  file name
!! @param[in] ifile  field pointer (nonlinear, mhd, perturbation; at the same time file number)
!! @note This is version of mfi including perturbation field.
!! @remark This routine uses global scratch space SCRNS, SCRCG.
      subroutine mfi_pert(fname,ifile)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'RESTART'
      include 'PARALLEL'
      include 'SOLN'
      include 'GEOM'

!     aragument list
      character*132  fname
      integer ifile

!     local variables
      character*132 hdr
      logical if_full_pres_tmp

      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk), pm1(lx1*ly1*lz1,lelv)
      common /scrns/ wk
      common /scrcg/ pm1
      integer e, j, k, ierr

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8
      integer iofldsr
      real tiostart, tio, dnbyte

!     functions
      real dnekclock, glsum
!-----------------------------------------------------------------------
      tiostart=dnekclock()

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      if_full_pres_tmp = if_full_pres
      if (wdsizr.eq.8) if_full_pres = .true. !Preserve mesh 2 pressure

      iofldsr = 0
      if (ifgetxr) then      ! if available
         offs = offs0 + ndim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetx) then
!            if(nio.eq.0) write(6,*) 'Reading mesh'
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.false.)
         else                ! skip the data
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + ndim
      endif

      if (ifgetur) then
         offs = offs0 + iofldsr*stride + ndim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetu) then
!     MHD
            if (ifmhd.and.ifile.eq.2) then
!               if(nio.eq.0) write(6,*) 'Reading B field'
               call mfi_getv(bx,by,bz,wk,lwk,.false.)
!     perturbation mode
!     importatn assumption; there are no MHD perturbation
            elseif (ifpert.and.ifile.ge.2) then
               j=ifile-1  ! pointer to perturbation field
               call mfi_getv(vxp(1,j),vyp(1,j),vzp(1,j),wk,lwk,.false.)
            else
!               if(nio.eq.0) write(6,*) 'Reading velocity field'
               call mfi_getv(vx,vy,vz,wk,lwk,.false.)
            endif
         else
            call mfi_getv(vx,vy,vz,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + ndim
      endif

      if (ifgetpr) then
         offs = offs0 + iofldsr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetp) then
!           if(nid.eq.0) write(6,*) 'Reading pressure field'
            call mfi_gets(pm1,wk,lwk,.false.)
         else
            call mfi_gets(pm1,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + 1
      endif

      if (ifgettr) then
         offs = offs0 + iofldsr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgett) then
!     perturbation mode
            if (ifpert.and.ifile.ge.2) then
               j=ifile-1  ! pointer to perturbation field
               call mfi_gets(tp(1,1,j),wk,lwk,.false.)
            else
!               if(nid.eq.0) write(6,*) 'Reading temperature field'
               call mfi_gets(t,wk,lwk,.false.)
            endif
         else
            call mfi_gets(t,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + 1
      endif

      do k=1,ldimt-1
         if (ifgtpsr(k)) then
            offs = offs0 + iofldsr*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            if (ifgtps(k)) then
!               if(nid.eq.0) write(6,'(A,I2,A)') ' Reading ps',k,' field'
!     perturbation mode
               if (ifpert.and.ifile.ge.2) then
                  j=ifile-1     ! pointer to perturbation field
                  call mfi_gets(tp(1,k+1,j),wk,lwk,.false.)
               else
                  call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.false.)
               endif
            else
               call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.true.)
            endif
            iofldsr = iofldsr + 1
         endif
      enddo
      nbyte = 0
      if(nid.eq.pid0r) nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr

      if (ifgtim) time = timer

      ierr = 0

#ifdef MPIIO
      if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
#else
      if (nid.eq.pid0r) call byte_close(ierr)
#endif
      call err_chk(ierr,'Error closing restart file, in mfi_pert.$')
      tio = dnekclock()-tiostart

      dnbyte = nbyte
      nbyte = glsum(dnbyte,1)
      nbyte = nbyte + iHeaderSize + 4 + isize*nelgr

      if(NIO.eq.0) write(6,7) istep,time,
     &             nbyte/tio/1024/1024/10,
     &             nfiler
    7 format(/,i9,1pe12.4,' done :: Read checkpoint data',/,
     &       30X,'avg data-throughput = ',f7.1,'MBps',/,
     &       30X,'io-nodes = ',i5,/)

!     neither perturbation, nor mhd here
      if (ifaxis) call axis_interp_ic_full_pres(pm1) ! Interpolate to axi mesh
      if (ifgetp) call map_pm1_to_pr_pert(pm1,ifile) ! Interpolate pressure

      if_full_pres = if_full_pres_tmp

      return
      end
!=======================================================================
!> @brief Map loaded variables from velocity to axisymmetric mesh
!! @ingroup chkpoint_mstep
!! @param[in] pm1    pressure loaded form the file
!! @note This is version of axis_interp_ic taking into account fact
!! pressure does not have to be written on velocity mesh.
!! @remark This routine uses global scratch space CTMP0.
      subroutine axis_interp_ic_full_pres(pm1)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'RESTART'
      include 'SOLN'
      include 'GEOM'
      include 'IXYZ'

!     argument list
      real pm1(lx1,ly1,lz1,lelv)

!     scratch space
      real axism1 (lx1,ly1)
      common /ctmp0/ axism1

!     local variables
      integer e
!-----------------------------------------------------------------------
      if (.not.ifaxis) return

      do e=1,nelv
         if (ifrzer(e)) then
           if (ifgetx) then
             call mxm   (xm1(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy  (xm1(1,1,1,e),axism1,nx1*ny1)
             call mxm   (ym1(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy  (ym1(1,1,1,e),axism1,nx1*ny1)
           endif
           if (ifgetu) then
             call mxm    (vx(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy   (vx(1,1,1,e),axism1,nx1*ny1)
             call mxm    (vy(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy   (vy(1,1,1,e),axism1,nx1*ny1)
           endif
           if (ifgetw) then
             call mxm    (vz(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy   (vz(1,1,1,e),axism1,nx1*ny1)
           endif
           if (ifgetp.and.(ifsplit.or.(.not.(if_full_pres)))) then
             call mxm    (pm1(1,1,1,e),nx1,iatlj1,ny1,axism1,ny1)
             call copy   (pm1(1,1,1,e),axism1,nx1*ny1)
           endif
           if (ifgett) then
             call mxm  (t (1,1,1,e,1),nx1,iatlj1,ny1,axism1,ny1)
             call copy (t (1,1,1,e,1),axism1,nx1*ny1)
           endif
           do ips=1,npscal
            is1 = ips + 1
            if (ifgtps(ips)) then
             call mxm (t(1,1,1,e,is1),nx1,iatlj1,ny1,axism1,ny1)
             call copy(t(1,1,1,e,is1),axism1,nx1*ny1)
            endif
           enddo
         endif
      enddo

      return
      end
!=======================================================================
!> @brief Map loaded pressure to the pressure mesh
!! @ingroup chkpoint_mstep
!! @param[in] pm1    pressure loaded form the file
!! @param[in] ifile  field pointer (nonlinear, mhd, perturbation;
!!   at the same time file number)
!! @note This is version of map_pm1_to_pr including perturbation
      subroutine map_pm1_to_pr_pert(pm1,ifile)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'RESTART'
      include 'SOLN'

!     argument list
      real pm1(lx1*ly1*lz1,lelv)
      integer ifile

!     local variables
      logical if_full_pres_tmp

      integer e, nxyz2, j, ie2
!-----------------------------------------------------------------------
      nxyz2 = nx2*ny2*nz2

      if (ifmhd.and.ifile.eq.2) then
         do e=1,nelv
            if (if_full_pres) then
               call copy  (pm(1,1,1,e),pm1(1,e),nxyz2)
            else
               call map12 (pm(1,1,1,e),pm1(1,e),e)
            endif
         enddo
      elseif (ifpert.and.ifile.ge.2) then
         j=ifile-1     ! pointer to perturbation field
         if (ifsplit) then
            call copy (prp(1,j),pm1,nx1*ny1*nz1*nelv)
         else
            do e=1,nelv
               ie2 = (e-1)*nxyz2+1
               if (if_full_pres) then
                  call copy(prp(ie2,j),pm1(1,e),nxyz2)
               else
                  call map12(prp(ie2,j),pm1(1,e),e)
               endif
            enddo
         endif
      elseif (ifsplit) then
         call copy (pr,pm1,nx1*ny1*nz1*nelv)
      else
         do e=1,nelv
            if (if_full_pres) then
               call copy  (pr(1,1,1,e),pm1(1,e),nxyz2)
            else
               call map12 (pr(1,1,1,e),pm1(1,e),e)
            endif
         enddo
      endif
   
      return
      end
!=======================================================================
