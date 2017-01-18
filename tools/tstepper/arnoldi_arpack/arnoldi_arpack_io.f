!> @file arnoldi_arpack_io.f
!! @ingroup arnoldi_arpack
!! @brief Set of checkpoint routines for arnoldi_arpack module.
!!
!=======================================================================
!> @brief Take care of savings for restart
!! @ingroup arnoldi_arpack
      subroutine arn_rst_save
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TSTEP_DEF'
      include 'TSTEP'           ! LASTEP
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
!-----------------------------------------------------------------------
!     save checkpoint for IDOARP=-2
      if (IDOARP.eq.-2) then

         if (NIO.eq.0) then
            write(*,*) 'ARNOLDI_ARPACK: '
            write(*,*) 'IDOARP  = ', IDOARP
            write(*,*) 'Writing checkpoint'
         endif

!     save parameters and WORKLA; independent on processor;
!     serial output
         call arn_write_par('ARP')

!     save big arrays; parallel output
         call mfo_arnv('ARV')

!     this is the last step
         LASTEP=1

      endif

      return
      end
!=======================================================================
!> @brief Read from checkpoints
!! @ingroup arnoldi_arpack
      subroutine arn_rst_read
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NIO
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'         !
!-----------------------------------------------------------------------
      if (NIO.eq.0) then
         write(*,*) 'ARNOLDI_ARPACK: '
         write(*,*) 'Reading checkpoint'
      endif
!     read parameters and WORKLA; independent on processor; serial input
      call arn_read_par('ARP')

!     read big arrays; parallel input
      call mfi_arnv('ARV')

      return
      end
!=======================================================================
!     Subroutines to manipulate files for arnoldi restart
!=======================================================================
!> @brief Write procesor independent data
!! @ingroup arnoldi_arpack
!! @param[in]   prefix    prefix
      subroutine arn_write_par(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character*132 fname

      character*6  str

      integer lwdsizo, ierr

      logical lifreguo
!-----------------------------------------------------------------------
      call nekgsync()

!     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.

      ierr = 0
!     this is done by master node only
      if (NID.eq.0) then
!     create file name
         call IO_mfo_fname(fname,SESSION,prefix,ierr)
         if (ierr.eq.0) then
            write(str,'(i5.5)') ARNISTOP
            fname = trim(fname)//trim(str)

!     open file; only serial
            call IO_mbyte_open_srl(fname,NID,ierr)
         endif

         if (ierr.eq.0) then
            call mfo_arnp

!     close the file; only serial
            call byte_close(ierr)
         endif
      endif

      call err_chk(ierr,'arn_write_par: file error. Abort. $')

!     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo

      return
      end
!=======================================================================
!> @brief Write procesor independent variables
!! @ingroup arnoldi_arpack
      subroutine mfo_arnp
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     local variables
      character*16 hdr
      integer ahdsize
      parameter (ahdsize=16)

      integer ibsw_out, il, itmp(33), ierr

      real*4 test_pattern

      real*4 rtmp4(6), workla4(2*WLDIMA)
      real*8 rtmp8(3), workla8(WLDIMA)
      equivalence (rtmp4,rtmp8)
      equivalence (workla4,workla8)
!-----------------------------------------------------------------------
!     write IDOARP and character varialbes
      call blank(hdr,ahdsize)

      write(hdr,1) IDOARP,BMATARP,WHICHARP,TSTMODE! 14
 1    format('#arp',1x,i2,1x,a1,1x,a2,1x,i1)

      ! if we want to switch the bytes for output
      ! switch it again because the hdr is in ASCII
      call get_bytesw_write(ibsw_out)
      if (ibsw_out.ne.0) call set_bytesw_write(0)

      call byte_write(hdr,ahdsize/4,ierr)

!     write test pattern for byte swap
      test_pattern = 6.54321

      call byte_write(test_pattern,1,ierr)

!     collect and write integer varialbes
       itmp(1) = NVECAS
       itmp(2) = ARNEGV
       itmp(3) = ARNKRYLOV
       itmp(4) = NWLARP
       itmp(5) = INFARP
       itmp(6) = NPARP
       itmp(7) = NCARP
       itmp(8) = TSTSTEP
       do il=1,11
          itmp(8+il) = IPARP(il)
       enddo
       do il=1,14
          itmp(19+il) = IPNTARP(il)
       enddo

       call byte_write(itmp,33,ierr)

!     collect and write real variables
       rtmp8(1) = TSTTOL
       rtmp8(2) = RNMARP
       rtmp8(3) = DT

       call byte_write(rtmp4,6,ierr)

!     write WORKLA
       call copy(workla8,WORKLA,NWLARP)

       call byte_write(workla4,2*NWLARP,ierr)

      return
      end
!=======================================================================
!> @brief Read procesor independent data
!! @ingroup arnoldi_arpack
!! @param[in]   prefix    prefix
      subroutine arn_read_par(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character*132 fname

      character*6  str

      integer lwdsizo, il, ierr

      logical lifreguo
!-----------------------------------------------------------------------
      call nekgsync()

!     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.

      ierr = 0
!     this is done by master node only
      if (NID.eq.0) then
!     create file name
         call IO_mfo_fname(fname,SESSION,prefix,ierr)
         if (ierr.eq.0) then
            write(str,'(i5.5)') ARNISTART
            fname = trim(fname)//trim(str)

!     open file; only serial
            call IO_mbyte_open_srl(fname,NID, ierr)
         endif

         if (ierr.eq.0) then
!     read parameters
            call mfi_arnp

!     close the file; only serial
            call byte_close(ierr)
         endif
      endif

      call err_chk(ierr,'arn_read_par: file error. Abort. $')

      ierr = 0
      if (NID.eq.0) then
!     check and copy parameters
!     is it correct restart
         if (IDOARP0.ne.-2) then
            write(*,*) 'ERROR: arn_read_par, wrong IDOARP0; abort.'
            write(*,*) 'IDOARP0 = ', IDOARP0
            ierr=1
         endif

!     is it the same ARPACK mode
         if (BMATARP0.ne.BMATARP) then
            write(*,*) 'ERROR: arn_read_par, different ARPACK modes.'
            write(*,*) 'BMATARP0 = ', BMATARP0
            write(*,*) 'BMATARP  = ', BMATARP
            ierr=1
         endif

!     do we look for the same eigenvectors
         if (WHICHARP0.ne.WHICHARP) then
            write(*,*) 'ERROR: arn_read_par, different mode selection.'
            write(*,*) 'WHICHARP0 = ', WHICHARP0
            write(*,*) 'WHICHARP  = ', WHICHARP
            ierr=1
         endif

!     is it the same integration mode
         if (TSTMODE0.ne.TSTMODE) then
            write(*,*) 'ERROR: arn_read_par, wrong simulation mode.'
            write(*,*) 'TSTMODE0 = ', TSTMODE0
            write(*,*) 'TSTMODE  = ', TSTMODE
            ierr=1
         endif

!     this should be removed later as it does not allow to change processor number
!     is the length of the vector the same
         if (NVECAS0.ne.NVECAS) then
            write(*,*) 'ERROR: arn_read_par, different vector lenght;
     $abort. Check IFHEAT!'
            write(*,*) 'NVECAS0 = ', NVECAS0
            write(*,*) 'NVECAS  = ', NVECAS
            ierr=1
         endif

!     what is the size of Krylov space
!     would it be possible to change this?; related NPARP, NWLARP,
!     IPNTARP
         if (ARNKRYLOV0.ne.ARNKRYLOV) then
            write(*,*) 'ERROR: arn_read_par, different Krylov space
     $size.'
            write(*,*) 'ARNKRYLOV0 = ', ARNKRYLOV0
            write(*,*) 'ARNKRYLOV  = ', ARNKRYLOV
            ierr=1
         endif

         if (NWLARP0.ne.NWLARP) then
            write(*,*) 'ERROR: arn_read_par, different size of work
     $array'
            write(*,*) 'NWLARP0 = ', NWLARP0
            write(*,*) 'NWLARP  = ', NWLARP
            ierr=1
         endif

!     stopping criterion
         if (TSTTOL0.ne.TSTTOL) then
            write(*,*) 'WARNING: arn_read_par, different stopping
     $criterion'
            write(*,*) 'TSTTOL0 = ', TSTTOL0
            write(*,*) 'TSTTOL  = ', TSTTOL
         endif

!     number of eigenvalues
         if (ARNEGV0.ne.ARNEGV) then
            write(*,*) 'WARNING: arn_read_par, different number of
     $eigenvalues'
            write(*,*) 'ARNEGV0 = ', ARNEGV0
            write(*,*) 'ARNEGV  = ', ARNEGV
         endif

c     stepper phase length
         if (DTARP0.ne.DT) then
            write(*,*) 'WARNING: arn_read_par, different timestep'
            write(*,*) 'DTARP0 = ', DTARP0
            write(*,*) 'DT     = ', DT
         endif

         if (TSTSTEP0.ne.TSTSTEP) then
            write(*,*) 'WARNING: arn_read_par, different number of
     $steps in stepper phase'
            write(*,*) 'TSTSTEP0 = ', TSTSTEP0
            write(*,*) 'TSTSTEP  = ', TSTSTEP
         endif

!     check IPARP
         if (IPARP0(1).ne.IPARP(1)) then
            write(*,*) 'ERROR: arn_read_par, different shift in ARPACK'
            write(*,*) 'IPARP0(1) = ', IPARP0(1)
            write(*,*) 'IPARP(1)  = ', IPARP(1)
            error=1
         endif

         if (IPARP0(3).ne.IPARP(3)) then
            write(*,*) 'WARNING: arn_read_par, different cycle number'
            write(*,*) 'IPARP0(3) = ', IPARP0(3)
            write(*,*) 'IPARP(3)  = ', IPARP(3)
         endif

         if (IPARP0(7).ne.IPARP(7)) then
            write(*,*) 'ERROR: arn_read_par, different ARPACK modes'
            write(*,*) 'IPARP0(7) = ', IPARP0(7)
            write(*,*) 'IPARP(7)  = ', IPARP(7)
            ierr=1
         endif

!     copy rest of parameters
         NPARP = NPARP0
         NCARP = NCARP0
         INFARP= INFARP0
         RNMARP= RNMARP0
         do il=4,11
            IPARP(il) = IPARP0(il)
         enddo
         IPARP(2) = IPARP0(2)
         do il=1,14
            IPNTARP(il) = IPNTARP0(il)
         enddo
      endif                     ! NID

      call err_chk(ierr,'Error restarting arnoldi.$')

!      call nekgsync()

      IDOARP = -2
!     broadcast
      call bcast(NPARP,ISIZE)
      call bcast(NCARP,ISIZE)
      call bcast(INFARP,ISIZE)
      call bcast(IPARP,11*ISIZE)
      call bcast(IPNTARP,14*ISIZE)
      call bcast(RNMARP,WDSIZE)

      call bcast(WORKLA,NWLARP*WDSIZE)

!     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo

      return
      end
!=======================================================================
!> @brief Read procesor independent variables
!! @ingroup arnoldi_arpack
      subroutine mfi_arnp
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'

!     local variables
      character*16 hdr
      character*4 dummy
      integer ahdsize
      parameter (ahdsize=16)

      integer ibsw_out, il, itmp(33), ierr

      real*4 test_pattern

      real*4 rtmp4(6), workla4(2*WLDIMA)
      real*8 rtmp8(3), workla8(WLDIMA)
      equivalence (rtmp4,rtmp8)
      equivalence (workla4,workla8)

      logical if_byte_swap_test, if_byte_sw_loc
!     functions
      integer indx2
!-----------------------------------------------------------------------
!     read IDOARP and character varialbes
      call blank(hdr,ahdsize)

      call byte_read(hdr,ahdsize/4,ierr)

      if (indx2(hdr,132,'#arp',4).eq.1) then
         read(hdr,*) dummy,IDOARP0,BMATARP0,WHICHARP0,TSTMODE0! 14
      else
         write(*,*) 'ERROR: mfi_arnp, wrong header; abort.'
         call exitt
      endif

!     read test pattern for byte swap
      call byte_read(test_pattern,1,ierr)
!     determine endianess
      if_byte_sw_loc = if_byte_swap_test(test_pattern,ierr)

!     read integer varialbes
       call byte_read(itmp,33,ierr)
       if (if_byte_sw) call byte_reverse(itmp,33,ierr)

       NVECAS0 = itmp(1)
       ARNEGV0  = itmp(2)
       ARNKRYLOV0  = itmp(3)
       NWLARP0 = itmp(4)
       INFARP0 = itmp(5)
       NPARP0  = itmp(6)
       NCARP0  = itmp(7)
       TSTSTEP0 = itmp(8)
       do il=1,11
          IPARP0(il) = itmp(8+il)
       enddo
       do il=1,14
          IPNTARP0(il) = itmp(19+il)
       enddo

!     read real variables
       call byte_read(rtmp4,6,ierr)
       if (if_byte_sw) call byte_reverse(rtmp4,6,ierr)

       TSTTOL0 = rtmp8(1)
       RNMARP0 = rtmp8(2)
       DTARP0  = rtmp8(3)

!     read WORKLA
       if (NWLARP0.le.WLDIMA) then
          call byte_read(workla4,2*NWLARP0,ierr)
          if (if_byte_sw) call byte_reverse(workla4,2*NWLARP0,ierr)

          call copy(WORKLA,workla8,NWLARP0)
       else
          write(*,*) 'ARNOLDI_ARPACKD: restart error'
          write(*,*) 'NWLARP0 = ',NWLARP0
          write(*,*) 'WLDIMA  = ',WLDIMA
          call exitt
       endif

      return
      end
!-----------------------------------------------------------------------
!     following subroutines are modiffications of
!
!     mfo_outfld
!     from prepost.f;
!     mfi
!     from ic.f
!=======================================================================
!> @brief Write procesor dependent data (long vectors)
!! @ingroup arnoldi_arpack
!! @param[in]   prefix    prefix
!! @remark This routine uses global scratch space SCRUZ
      subroutine mfo_arnv(prefix)  ! muti-file output
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      include 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzo8
      integer lwdsizo, il, ierr
      integer ioflds, nout
      real dnbyte, tio, tiostart

      character*132  fname

      character*6  str

      logical lifxyo, lifpo, lifvo, lifto, lifreguo, lifpso(LDIMT1)

!     functions
      real dnekclock_sync, glsum

!     scratch space
      real UR1(LXO*LXO*LXO*LELT), UR2(LXO*LXO*LXO*LELT),
     $     UR3(LXO*LXO*LXO*LELT)
      common /SCRUZ/  UR1, UR2, UR3
!-----------------------------------------------------------------------
      tiostart=dnekclock_sync()

!     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.
      lifxyo= IFXYO
      IFXYO = .false.
      lifpo= IFPO
      IFPO = .false.
      lifvo= IFVO
      IFVO = .true.
      lifto= IFTO
      IFTO = .false.
      do il=1,LDIMT1
         lifpso(il)= IFPSO(il)
         IFPSO(il) = .false.
      enddo

      nout = NELT
      NXO  = NX1
      NYO  = NY1
      NZO  = NZ1

!     set offset
      offs0 = iHeaderSize + 4 + isize*nelgt

      ierr = 0
      if (nid.eq.pid0) then         ! open files on i/o nodes
!     create file name
         call IO_mfo_fname(fname,SESSION,prefix,ierr)
         if (ierr.eq.0) then
            write(str,'(i5.5)') ARNISTOP
            fname = trim(fname)//trim(str)
            call mbyte_open(fname,fid0,ierr)
         endif
      endif

      call err_chk(ierr,'Error opening file in mfo_arnv. Abort. $')

      call mfo_write_hdr                     ! create element mapping +
                                             ! write hdr
      nxyzo8  = NXO*NYO*NZO
      strideB = nelB * nxyzo8*WDSIZO
      stride  = nelgt* nxyzo8*WDSIZO

      ioflds = 0
      ! dump all fields based on the t-mesh to avoid different
      ! topologies in the post-processor

!     resid array
      call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,RESIDA)

!     workd array
      do il=0,2
         call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,WORKDA(1+NVECAS*il))
      enddo

!     krylov space
      do il=1,ARNKRYLOV
         call  mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     UR1,UR2,UR3,VBASEA(1,il))
      enddo

      dnbyte = 1.*ioflds*nout*WDSIZO*NXO*NYO*NZO

!     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo
      IFXYO = lifxyo
      IFPO = lifpo
      IFVO = lifvo
      IFTO = lifto
      do il=1,LDIMT1
         IFPSO(i) = lifpso(il)
      enddo

      if (NID.eq.PID0)
#ifdef MPIIO
     &   call byte_close_mpi(ifh_mbyte,ierr)
#else
     &   call byte_close(ierr)
#endif

      tio = dnekclock_sync()-tiostart
      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4. +ISIZE*NELGT
      dnbyte = dnbyte/1024/1024
      if(nid.eq.0) write(6,7)  ISTEP,TIME,dnbyte,dnbyte/tio,
     &             NFILEO
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &       30X,'file size = ',3pG12.2,'MB',/,
     &       30X,'avg data-throughput = ',0pf7.1,'MB/s',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
!=======================================================================
!> @brief Read procesor dependent data (long vectors)
!! @ingroup arnoldi_arpack
!! @param[in]   prefix    prefix
!! @remark This routine uses global scratch space SCRNS, SCRUZ
      subroutine mfi_arnv(prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'

!     argument list
      character*3 prefix

!     local variables
      character*132  fname

      character*6  str

      integer e, il, iofldsr, ierr
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8
      real dnbyte, tio, tiostart

!     functions
      real dnekclock, glsum

!     scratch space
      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk)
      common /scrns/ wk

      real UR1(LX1,LY1,LZ1,LELT), UR2(LX1,LY1,LZ1,LELT),
     $     UR3 (LX1,LY1,LZ1,LELT)
      COMMON /scruz/ UR1, UR2, UR3
!-----------------------------------------------------------------------
      tiostart=dnekclock()

!     create file name
      ierr = 0
      if (nid.eq.pid0r) then         ! open files on i/o nodes
         call IO_mfo_fname(fname,SESSION,prefix,ierr)
         if (ierr.eq.0) then
            write(str,'(i5.5)') ARNISTART
            fname = trim(fname)//trim(str)
         endif
      endif

      call err_chk(ierr,'Error opening file in mfi_arnv. Abort. $')

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping

      offs0   = iHeadersize + 4 + ISIZE*NELGR
      nxyzr8  = NXR*NYR*NZR
      strideB = nelBr* nxyzr8*WDSIZR
      stride  = nelgr* nxyzr8*WDSIZR

!     read arrays
      iofldsr = 0
!     resid array
      call mfi_singlev(iofldsr,offs0,stride,strideB,
     $     UR1,UR2,UR3,RESIDA(1))

!     workd array
      do il=0,2
         call mfi_singlev(iofldsr,offs0,stride,strideB,
     $        UR1,UR2,UR3,WORKDA(1+NVECAS*il))
      enddo

!     krylov space
      do il=1,ARNKRYLOV
         call mfi_singlev(iofldsr,offs0,stride,strideB,
     $        UR1,UR2,UR3,VBASEA(1,il))
      enddo

      nbyte = 0
      if(nid.eq.pid0r) then
         nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr
!     close files
#ifdef MPIIO
         call byte_close_mpi(ifh_mbyte,ierr)
#else
         call byte_close(ierr)
#endif
      endif
      call nekgsync
      tio = dnekclock()-tiostart

      dnbyte = nbyte
      nbyte = glsum(dnbyte,1)
      nbyte = nbyte + iHeaderSize + 4 + isize*nelgr

      if(nid.eq.0) write(6,7) istep,time,
     &             nbyte/tio/1024/1024/10,
     &             nfiler
    7 format(/,i9,1pe12.4,' done :: Read checkpoint data',/,
     &       30X,'avg data-throughput = ',f7.1,'MBps',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
!=======================================================================
!> @brief Write single Krylov vector to the file
!! @ingroup arnoldi_arpack
!! @param[inout] ioflds         Vector counter
!! @param[in]    nout           local number of elements to write
!! @param[in]    offs0          global file offset (header +...)
!! @param[in]    stride         single vector length
!! @param[in]    strideB        space saved for processes with lower nid
!! @param[in]    ur1,ur2,ur3    output arrays
!! @param[in]    vect           Krylov vector
      subroutine mfo_singlev(ioflds,nout,offs0,stride,strideB,
     $     ur1,ur2,ur3,vect)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'

!     argument list
      integer ioflds,nout
      integer*8 offs0,stride,strideB
      real UR1(LXO*LXO*LXO*LELT), UR2(LXO*LXO*LXO*LELT),
     $     UR3(LXO*LXO*LXO*LELT)
      real VECT(LVAS)

!     local variables
      integer*8 offs
!-----------------------------------------------------------------------
      offs = offs0 + ioflds*stride + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)

      call copy(UR1,VECT(1),NVECAV)
      call copy(UR2,VECT(1+NVECAV),NVECAV)
      if (IF3D) call copy(UR3,VECT(1+2*NVECAV),NVECAV)

      call mfo_outv(UR1,UR2,UR3,nout,NXO,NYO,NZO)
      ioflds = ioflds + NDIM

      if (IFHEAT) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         call copy(UR1,VECT(1+NDIM*NVECAV),NVECAT)

         call mfo_outs(UR1,nout,NXO,NYO,NZO)
         ioflds = ioflds + 1
      endif

      return
      end
!=======================================================================
!> @brief Read single Krylov vector from the file
!! @ingroup arnoldi_arpack
!! @param[inout] iofldr         Vector counter
!! @param[in]    offs0          global file offset (header +...)
!! @param[in]    stride         single vector length
!! @param[in]    strideB        space saved for processes with lower nid
!! @param[in]    ur1,ur2,ur3    input arrays
!! @param[out]   vect           Krylov vector
!! @remark This routine uses global scratch space SCRNS
      subroutine mfi_singlev(iofldr,offs0,stride,strideB,
     $     ur1,ur2,ur3,vect)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'RESTART_DEF'
      include 'RESTART'
      include 'TIME_STEPPERD'
      INCLUDE 'ARNOLDI_ARPACKD'
!     argument list
      integer iofldr
      integer*8 offs0,stride,strideB
      real UR1(LX1*LX1*LX1*LELT), UR2(LX1*LX1*LX1*LELT),
     $     UR3(LX1*LX1*LX1*LELT)
      real VECT(LVAS)

      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk)
      common /scrns/ wk

!     local variables
      integer*8 offs
!-----------------------------------------------------------------------
      offs = offs0 + iofldr*stride + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)
      call mfi_getv(UR1,UR2,UR3,wk,lwk,.false.)

      call copy(VECT(1),UR1,NVECAV)
      call copy(VECT(1+NVECAV),UR2,NVECAV)
      if (IF3D) call copy(VECT(1+2*NVECAV),UR3,NVECAV)
      iofldr = iofldr + NDIM

      if (IFHEAT) then
         offs = offs0 + iofldr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         call mfi_gets(UR1,wk,lwk,.false.)
         call copy(VECT(1+NDIM*NVECAV),UR1,NVECAT)
         iofldr = iofldr + 1
      endif

      return
      end
!=======================================================================

!=======================================================================

