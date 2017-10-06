!> @file mntr.f
!! @ingroup monitor
!! @brief Set of monitoring routines for KTH framework
!! @author Adam Peplinski
!! @date Sep 28, 2017
!=======================================================================
!> @brief Initialise monitor by registering Nek5000 and monitor
!! @ingroup monitor
      subroutine mntr_register_mod
      implicit none

      include 'SIZE'
      include 'MNTRD'
      include 'MNTRLP'

!     local variables
      character*2 str
      character*200 lstring
!-----------------------------------------------------------------------
!     first register neck5000
      mntr_nek_id = 1
      mntr_mod_id(mntr_nek_id) = 0
      mntr_mod_name(mntr_nek_id) = mntr_nek_name
      mntr_mod_dscr(mntr_nek_id) = 'CFD solver'
      mntr_mod_num = mntr_mod_num + 1
      mntr_mod_mpos = mntr_mod_mpos + 1

!     next monitor
      mntr_id = 2
      mntr_mod_id(mntr_id) = mntr_nek_id
      mntr_mod_name(mntr_id) = mntr_name
      mntr_mod_dscr(mntr_id) = 'Monitoring module'
      mntr_mod_num = mntr_mod_num + 1
      mntr_mod_mpos = mntr_mod_mpos + 1

!     set log threshold
      mntr_lp_def = lp_inf

!     log changes
      lstring = 'Registered module ['//trim(mntr_mod_name(mntr_nek_id))
      lstring= trim(lstring)//']: '//trim(mntr_mod_dscr(mntr_nek_id))
      call mntr_log(mntr_id,lp_inf,trim(lstring))

      lstring = 'Registered module ['//trim(mntr_mod_name(mntr_id))
      lstring= trim(lstring)//']: '//trim(mntr_mod_dscr(mntr_id))
      call mntr_log(mntr_id,lp_inf,trim(lstring))

      write(str,'(I2)') mntr_lp_def
      call mntr_log(mntr_id,lp_inf,
     $     'Initial log threshold set to: '//trim(str))

      return
      end subroutine
!=======================================================================
!> @brief Register monitor runtime parameters
!! @ingroup monitor
      subroutine mntr_register_par
      implicit none

      include 'SIZE'
      include 'RPRMD'
      include 'MNTRD'
      include 'MNTRLP'

!     local variables
      integer rpid,itmp
      real rtmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
!     register and set active section
      call rprm_sec_reg(mntr_sec_id,mntr_id,'_'//adjustl(mntr_name),
     $     'Runtime paramere section for monitor module')
      call rprm_sec_set_act(.true.,mntr_sec_id)

!     register parameters
      call rprm_rp_reg(mntr_lp_def_id,mntr_sec_id,'LOGLEVEL',
     $     'Logging threshold for toolboxes',rprm_par_int,lp_inf,
     $      0.0,.false.,' ')

      return
      end subroutine
!=======================================================================
!> @brief Initialise monitor module
!! @ingroup monitor
      subroutine mntr_init
      implicit none

      include 'SIZE'
      include 'RPRMD'
      include 'MNTRD'
      include 'MNTRLP'

!     local variables
      integer itmp
      real rtmp
      logical ltmp
      character*20 ctmp
      character*2 str
!-----------------------------------------------------------------------
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,mntr_lp_def_id,rprm_par_int)
      mntr_lp_def = itmp

      write(str,'(I2)') mntr_lp_def
      call mntr_log(mntr_id,lp_inf,
     $     'Reseting log threshold to: '//trim(str))

      return
      end subroutine
!=======================================================================
!> @brief Register new module
!! @ingroup monitor
!! @param[out] mid      current module id
!! @param[in]  pmid     parent module id
!! @param[in]  mname    module name
!! @param[in]  mdscr    module description
      subroutine mntr_mod_reg(mid,pmid,mname,mdscr)
      implicit none

      include 'SIZE'
      include 'PARALLEL'        ! ISIZE
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      integer mid, pmid
      character*(*) mname, mdscr

!     local variables
      character*10  lname
      character*132 ldscr
      integer slen,slena

      integer il, ipos
!-----------------------------------------------------------------------
!     check name length
      slena = len_trim(adjustl(mname))
!     remove trailing blanks
      slen = len_trim(mname) - slena + 1
      if (slena.gt.mntr_lstl_mnm) then
         call mntr_log(mntr_id,lp_deb,
     $        'too long module name; shortenning')
         slena = min(slena,mntr_lstl_mnm)
      endif
      call blank(lname,mntr_lstl_mnm)
      lname= mname(slen:slen+slena- 1)
      call capit(lname,slena)

!     check description length
      slena = len_trim(adjustl(mdscr))
!     remove trailing blanks
      slen = len_trim(mdscr) - slena + 1
      if (slena.ge.mntr_lstl_mds) then
         call mntr_log(mntr_id,lp_deb,
     $        'too long module description; shortenning')
         slena = min(slena,mntr_lstl_mnm)
      endif
      call blank(ldscr,mntr_lstl_mds)
      ldscr= mdscr(slen:slen + slena - 1)

!     find empty space
      ipos = 0

!     to ensure consistency I do it on master and broadcast result
      if (nid.eq.mntr_pid0) then

!     check if module is already registered
         do il=1,mntr_mod_mpos
            if (mntr_mod_id(il).ge.0.and.
     $         mntr_mod_name(il).eq.lname) then
               ipos = -il
               exit
            endif
         enddo

!     find empty spot
         if (ipos.eq.0) then
            do il=1,mntr_id_max
               if (mntr_mod_id(il).eq.-1) then
                  ipos = il
                  exit
               endif
            enddo
         endif
      endif

!     broadcast mid
      call bcast(ipos,isize)

!     error; no free space found
      if (ipos.eq.0) then
         mid = ipos
         call mntr_abort(mntr_id,
     $        'module ['//trim(lname)//'] cannot be registered')
!     module already registered
      elseif (ipos.lt.0) then
         mid = abs(ipos)
         call mntr_abort(mntr_id,
     $    'Module ['//trim(lname)//'] is already registered')
!     new module
      else
         mid = ipos
!        check if parent module is registered
         if (pmid.gt.0) then
            if (mntr_mod_id(pmid).ge.0) then
               mntr_mod_id(ipos) = pmid
            else
               mntr_mod_id(ipos) = 0
               call mntr_log(mntr_id,lp_inf,
     $       "Module's ["//trim(lname)//"] parent not registered.")
            endif
         else
            mntr_mod_id(ipos) = 0
         endif
         mntr_mod_name(ipos)=lname
         mntr_mod_dscr(ipos)=ldscr
         mntr_mod_num = mntr_mod_num + 1
         if (mntr_mod_mpos.lt.ipos) mntr_mod_mpos = ipos
         call mntr_log(mntr_id,lp_inf,
     $       'Registered module ['//trim(lname)//']: '//trim(ldscr))
      endif

      return
      end subroutine
!=======================================================================
!> @brief Check if module name is registered and return its id.
!! @ingroup monitor
!! @param[out] mid      module id
!! @param[in]  mname    module name
      subroutine mntr_mod_is_name_reg(mid,mname)
      implicit none

      include 'SIZE'
      include 'PARALLEL'        ! ISIZE
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      integer mid
      character*(*) mname

!     local variables
      character*10  lname
      character*3 str
      integer slen,slena

      integer il, ipos
!-----------------------------------------------------------------------
!     check name length
      slena = len_trim(adjustl(mname))
!     remove trailing blanks
      slen = len_trim(mname) - slena + 1
      if (slena.gt.mntr_lstl_mnm) then
         call mntr_log(mntr_id,lp_deb,
     $          'too long module name; shortenning')
         slena = min(slena,mntr_lstl_mnm)
      endif
      call blank(lname,mntr_lstl_mnm)
      lname= mname(slen:slen+slena- 1)
      call capit(lname,slena)

!     find module
      ipos = 0

!     to ensure consistency I do it on master and broadcast result
      if (nid.eq.mntr_pid0) then
!     check if module is already registered
         do il=1,mntr_mod_mpos
            if (mntr_mod_id(il).ge.0.and.
     $         mntr_mod_name(il).eq.lname) then
               ipos = il
               exit
            endif
         enddo
      endif

!     broadcast ipos
      call bcast(ipos,isize)

      if (ipos.eq.0) then
         mid = -1
         call mntr_log(mntr_id,lp_inf,
     $        'Module ['//trim(lname)//'] not registered')
      else
         mid = ipos
         write(str,'(I3)') ipos
         call mntr_log(mntr_id,lp_inf,
     $        'Module ['//trim(lname)//'] registered with mid='//str)
      endif

      return
      end subroutine
!=======================================================================
!> @brief Check if module id is registered. This operation is performed locally
!! @ingroup monitor
!! @param[in] mid      module id
      logical function mntr_mod_is_id_reg(mid)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      integer mid
!-----------------------------------------------------------------------
      mntr_mod_is_id_reg = mntr_mod_id(mid).ge.0

      return
      end function
!=======================================================================
!> @brief Get number of registered modules. This operation is performed locally
!! @ingroup monitor
!! @param[out]    nmod     module number
!! @param[out]    mmod     max module id
      subroutine mntr_mod_get_number(nmod,mmod)
      implicit none

      include 'SIZE'
      include 'MNTRD'

!     argument list
      integer nmod, mmod
!-----------------------------------------------------------------------
      nmod = mntr_mod_num
      mmod = mntr_mod_mpos

      return
      end subroutine
!=======================================================================
!> @brief Get module name an parent id for given module id. This operation is performed locally
!! @ingroup monitor
!! @param[out]    pmid     parent module id
!! @param[out]    mname    module name
!! @param[inout]  mid      module id
      subroutine mntr_mod_get_info(mname, pmid,mid)
      implicit none

      include 'SIZE'
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      character*10 mname
      integer mid, pmid

!     local variables
      character*5 str
!-----------------------------------------------------------------------
      if (mntr_mod_id(mid).ge.0) then
         pmid = mntr_mod_id(mid)
         mname = mntr_mod_name(mid)
      else
         mid = -1
         write(str,'(I3)') mid
         call mntr_log(mntr_id,lp_vrb,
     $        'Module id'//trim(str)//' not registered')
      endif

      return
      end subroutine
!=======================================================================
!> @brief Write log message
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
      subroutine mntr_log(mid,priority,logs)
      implicit none

      include 'SIZE'
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      integer mid,priority
      character*(*) logs

!     local variables
      character*200 llogs
      character*5 str
      integer slen, slena
!-----------------------------------------------------------------------
!     check log priority
      if (priority.lt.mntr_lp_def) return

!     done only by master
      if (nid.eq.mntr_pid0) then

!     check description length
         slena = len_trim(adjustl(logs))
!     remove trailing blanks
         slen = len_trim(logs) - slena + 1
         if (slena.ge.mntr_lstl_log) then
            if (mntr_lp_def.le.lp_deb) write(*,*)' ['//mntr_name//'] ',
     $       'too long log string; shortenning'
            slena = min(slena,mntr_lstl_log)
         endif
         call blank(llogs,mntr_lstl_mds)
         llogs= logs(slen:slen + slena - 1)

!     check module id
         if (mntr_mod_id(mid).ge.0) then
!     add module name
            write(*,*) ' ['//trim(mntr_mod_name(mid))//'] '//trim(llogs)
         else
            write(str,'(I3)') mid
            write(*,*) ' ['//trim(mntr_name)//'] ',
     $      ' WARNING: module'//trim(str)//' not registered;'
            write(*,*) 'Log body: '//trim(llogs)
         endif
      endif

      return
      end subroutine
!=======================================================================
!> @brief Write log message from given process
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
!! @param[in] prid      process id
      subroutine mntr_log_local(mid,priority,logs,prid)
      implicit none

      include 'SIZE'
      include 'MNTRD'
      include 'MNTRLP'

!     argument list
      integer mid,priority, prid
      character*(*) logs

!     local variables
      character*200 llogs
      character*5 str
      integer slen, slena
!-----------------------------------------------------------------------
!     check log priority
      if (priority.lt.mntr_lp_def) return

!     done only by given process

!     check description length
      slena = len_trim(adjustl(logs))
!     remove trailing blanks
      slen = len_trim(logs) - slena + 1
      if (slena.ge.mntr_lstl_log) then
         if (mntr_lp_def.le.lp_deb) write(*,*)' ['//mntr_name//'] ',
     $       'too long log string; shortenning'
         slena = min(slena,mntr_lstl_log)
      endif
      call blank(llogs,mntr_lstl_mds)
      llogs= logs(slen:slen + slena - 1)

!     check module id
      if (mntr_mod_id(mid).ge.0) then
!     add module name
       write(*,*) ' ['//trim(mntr_mod_name(mid))//'] nid= ',prid,
     $      ' '//trim(llogs)
      else
         write(str,'(I3)') mid
         write(*,*) ' ['//trim(mntr_name)//'] ',
     $   ' WARNING: module'//trim(str)//' not registered;'
         write(*,*) 'Log body: nid= ',prid,' '//trim(llogs)
      endif

      return
      end subroutine
!=======================================================================
!> @brief Write log message adding single integer
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
!! @param[in] ivar      integer variable
      subroutine mntr_logi(mid,priority,logs,ivar)
      implicit none

!     argument list
      integer mid,priority,ivar
      character*(*) logs

!     local variables
      character*10 str
!-----------------------------------------------------------------------
      write(str,'(I8)') ivar
      call mntr_log(mid,priority,trim(logs)//' '//trim(str))

      return
      end subroutine
!=======================================================================
!> @brief Write log message adding single real
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
!! @param[in] rvar      real variable
      subroutine mntr_logr(mid,priority,logs,rvar)
      implicit none

!     argument list
      integer mid,priority
      character*(*) logs
      real rvar

!     local variables
      character*20 str
!-----------------------------------------------------------------------
      write(str,'(E15.8)') rvar
      call mntr_log(mid,priority,trim(logs)//' '//trim(str))

      return
      end subroutine
!=======================================================================
!> @brief Write log message adding single logical
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
!! @param[in] rvar      logical variable
      subroutine mntr_logl(mid,priority,logs,lvar)
      implicit none

!     argument list
      integer mid,priority
      character*(*) logs
      logical lvar

!     local variables
      character*2 str
!-----------------------------------------------------------------------
      write(str,'(L2)') lvar
      call mntr_log(mid,priority,trim(logs)//' '//trim(str))

      return
      end subroutine
!=======================================================================
!> @brief Write warning message
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] logs      log body
      subroutine mntr_warn(mid,logs)
      implicit none

      include 'MNTRLP'

!     argument list
      integer mid,priority
      character*(*) logs
!-----------------------------------------------------------------------
      call mntr_log(mid,lp_inf,'WARNING: '//logs)
      return
      end subroutine
!=======================================================================
!> @brief Write error message
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] logs      log body
      subroutine mntr_error(mid,logs)
      implicit none

      include 'MNTRLP'

!     argument list
      integer mid
      character*(*) logs
!-----------------------------------------------------------------------
      call mntr_log(mid,lp_err,'ERROR: '//logs)
      return
      end subroutine
!=======================================================================
!> @brief Abort simulation
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] logs      log body
      subroutine mntr_abort(mid,logs)
      implicit none

      include 'MNTRLP'

!     argument list
      integer mid
      character*(*) logs
!-----------------------------------------------------------------------
      call mntr_log(mid,lp_err,'ABORT: '//logs)
      call exitt
      return
      end subroutine
!=======================================================================
!> @brief Abort simulation
!! @ingroup monitor
!! @param[in] mid       module id
!! @param[in] ierr      error flag
!! @param[in] logs      log body
      subroutine mntr_check_abort(mid,ierr,logs)
      implicit none

      include 'MNTRLP'

!     argument list
      integer mid,ierr
      character*(*) logs

!     local variables
      integer itest
      character*5 str
!     functions
      integer iglmax
!-----------------------------------------------------------------------
      itest = iglmax(ierr,1)

      if (itest.gt.0) then
         write(str,'(I3)') itest
         call mntr_log(mid,lp_err,
     $         'ABORT: '//trim(logs)//' ierr='//trim(str))
         call exitt
      endif
      return
      end subroutine
!=======================================================================

