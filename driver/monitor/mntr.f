!> @file mntr.f
!! @ingroup monitor
!! @brief Set of monitoring routines for KTH framework
!! @author Adam Peplinski
!! @date Sep 28, 2017
!=======================================================================
!> @brief Initialise monitor by registering Nek5000 an dmonitor
      subroutine mntr_init
      implicit none

      include 'SIZE'
      include 'MNTRD'
      include 'MNTRLP'

!     local variables
      character*2 str
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
      mntr_lp_def = 0 ! just for now
      write(str,'(I2)') mntr_lp_def
      call mntr_log(mntr_id,lp_ess,'Log threshold set to: '//trim(str))

      return
      end
!=======================================================================
!> @brief Register new module
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

!     error; no space found
      if (ipos.eq.0) then
         mid = ipos
         call mntr_abort(mntr_id,
     $        'module ['//trim(lname)//'] cannot be registered')
!     module already registered
      elseif (ipos.lt.0) then
         mid = abs(ipos)
         call mntr_log(mntr_id,lp_inf,
     $    'Module ['//trim(lname)//'] is already registered')
!     check parent
         if(pmid.ne.mntr_mod_id(mid)) call mntr_warn(mntr_id,
     $      "Module's ["//trim(lname)//"] parent inconsistent")
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
      end
!=======================================================================
!> @brief Check if module is registered and return its id.
!! @param[out] mid      module id
!! @param[in]  mname    module name
      subroutine mntr_mod_is_reg(mid,mname)
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

!     broadcast mid
      call bcast(ipos,isize)

!     error; no space found
      if (ipos.gt.0) then
         mid = ipos
         write(str,'(I3)') ipos
         call mntr_log(mntr_id,lp_vrb,
     $        'Module ['//trim(lname)//'] registered with mid='//str)
!     module already exist
      else
         mid = -1
         call mntr_log(mntr_id,lp_vrb,
     $        'Module ['//trim(lname)//'] not registered')
      endif

      return
      end
!=======================================================================
!> @brief Write log message
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
      end
!=======================================================================
!> @brief Write log message adding single integer
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
      end
!=======================================================================
!> @brief Write log message adding single real
!! @param[in] mid       module id
!! @param[in] priority  log priority
!! @param[in] logs      log body
!! @param[in] rvar      integer variable
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
      end
!=======================================================================
!> @brief Write warning message
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
      end
!=======================================================================
!> @brief Write error message
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
      end
!=======================================================================
!> @brief Abort simulation
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
      end
!=======================================================================
!> @brief Abort simulation
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
      end
!=======================================================================

