!> @file nseb.f
!! @ingroup nseb
!> @brief Adding white noise in a box
!! @author Adam Peplinski
!! @date Feb 2, 2017
!=======================================================================
!> @brief Register noise box module
!! @ingroup nseb
!! @note This routine should be called in frame_usr_register
      subroutine nseb_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'NSEBD'

      ! local variables
      integer lpmid, il, jl
      real ltim
      ! initial values for nseb_rfc
      real rfc(3,ldim)
      character*1 str1, str2

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,nseb_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(nseb_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'Parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(nseb_id,lpmid,nseb_name,
     $                 'Adding white noise in rectangular domain')

      ! register timer
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      call mntr_tmr_reg(nseb_tmr_id,lpmid,nseb_id,
     $      'NSEB_TOT','Noise box total time',.false.)

      ! register and set active section
      call rprm_sec_reg(nseb_sec_id,nseb_id,'_'//adjustl(nseb_name),
     $     'Runtime paramere section for nseb module')
      call rprm_sec_set_act(.true.,nseb_sec_id)

      ! register parameters
      call rprm_rp_reg(nseb_tim_id,nseb_sec_id,'TIME',
     $     'Time to add noise',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(nseb_amp_id,nseb_sec_id,'AMPLITUDE',
     $     'Noise amplitude',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(nseb_bmin_id(1),nseb_sec_id,'BOXMINX',
     $     'Position of lower left box corner; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(nseb_bmin_id(2),nseb_sec_id,'BOXMINY',
     $     'Position of lower left box corner; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(nseb_bmin_id(ndim),nseb_sec_id,
     $     'BOXMINZ','Position of lower left box corner; dimension Z ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(nseb_bmax_id(1),nseb_sec_id,'BOXMAXX',
     $     'Position of upper right box corner; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(nseb_bmax_id(2),nseb_sec_id,'BOXMAXY',
     $     'Position of upper right box corner; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(nseb_bmax_id(ndim),nseb_sec_id,
     $     'BOXMAXZ','Position of upper right box corner; dimension Z ',
     $     rpar_real,0,0.0,.false.,' ')

      ! function coefficients for random number generator
      ! initial values
      rfc(1,1) = 3.0e4
      rfc(2,1) =-1.5e3
      rfc(3,1) = 0.5e5
      rfc(1,2) = 2.3e4
      rfc(2,2) = 2.3e3
      rfc(3,2) =-2.0e5
      if (IF3D) then
         rfc(1,ldim) = 2.e4
         rfc(2,ldim) = 1.e3
         rfc(3,ldim) = 1.e5
      end if
      do il = 1, ldim
         write(str1,'(I1.1)') il
         do jl=1, 3
         write(str2,'(I1.1)') jl
         call rprm_rp_reg(nseb_rfc_id(jl,il),nseb_sec_id,
     $     'FRC'//str2//'_'//str1,
     $     'Function coefficient for random number gnerator ',
     $     rpar_real,0,rfc(jl,il),.false.,' ')
         end do
      end do

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(nseb_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise noise box module
!! @ingroup nseb
!! @note This routine should be called in frame_usr_init
!! @remark This routine uses global scratch space \a SCRUZ
      subroutine nseb_init()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'NSEBD'

      ! local variables
      integer ierr, nhour, nmin, il, jl
      integer itmp
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (nseb_ifinit) then
         call mntr_warn(nseb_id,
     $        'module ['//trim(nseb_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_tim_id,rpar_real)
      nseb_tim = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_amp_id,rpar_real)
      nseb_amp = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmin_id(1),
     $     rpar_real)
      nseb_bmin(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmin_id(2),
     $     rpar_real)
      nseb_bmin(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmin_id(ndim),
     $     rpar_real)
         nseb_bmin(ndim) = rtmp
      endif

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmax_id(1),
     $     rpar_real)
      nseb_bmax(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmax_id(2),
     $     rpar_real)
      nseb_bmax(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_bmax_id(ndim),
     $     rpar_real)
         nseb_bmax(ndim) = rtmp
      endif

      ! function coefficients for random number generator
      do il = 1, ldim
         do jl=1, 3
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,nseb_rfc_id(jl,il),
     $     rpar_real)
         nseb_rfc(jl,il) = rtmp
         end do
      end do

      ! is everything initialised
      nseb_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(nseb_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup nseb
!! @return nseb_is_initialised
      logical function nseb_is_initialised()
      implicit none

      include 'SIZE'
      include 'NSEBD'
!-----------------------------------------------------------------------
      nseb_is_initialised = nseb_ifinit

      return
      end function
!=======================================================================
!> @brief Add noise to velocity field in a box
!! @ingroup nseb
      subroutine nseb_noise_add()
      implicit none

      include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
      include 'TSTEP'           ! TIME, DT
      include 'PARALLEL'        ! LGLEL
      include 'INPUT'           ! IF3D
      include 'SOLN'            ! VX, VY, VZ, VMULT
      include 'GEOM'            ! XM1, YM1, ZM1
      include 'FRAMELP'
      include 'NSEBD'

      ! local variables
      integer iel, ieg, il, jl, kl, nl
      real xl(LDIM)
      logical ifadd
      real fcoeff(3)            !< coefficients for random distribution
      real ltim

      ! functions
      real dnekclock, math_ran_dst
!-----------------------------------------------------------------------
      ! add noise
      if (nseb_amp.gt.0.0) then
         if (nseb_tim.ge.TIME.and.nseb_tim.le.(TIME+DT)) then

            ! timing
            ltim = dnekclock()
            call mntr_log(nseb_id,lp_inf,
     $          "Adding noise to velocity field")

            do iel=1,NELV
               do kl=1,NZ1
                  do jl=1,NY1
                     do il=1,NX1
                        ieg = LGLEL(iel)
                        xl(1) = XM1(il,jl,kl,iel)
                        xl(2) = YM1(il,jl,kl,iel)
                        if (IF3D) xl(NDIM) = ZM1(il,jl,kl,iel)
                        ifadd = .TRUE.
                        do nl=1,NDIM
                           if (xl(nl).lt.nseb_bmin(nl).or.
     $                          xl(nl).gt.nseb_bmax(nl)) then
                              ifadd = .FALSE.
                              exit
                           endif
                        enddo

                        if (ifadd) then
                           VX(il,jl,kl,iel)=VX(il,jl,kl,iel)+nseb_amp*
     $                      math_ran_dst(il,jl,kl,ieg,xl,nseb_rfc(1,1))
                           VY(il,jl,kl,iel)=VY(il,jl,kl,iel)+nseb_amp*
     $                      math_ran_dst(il,jl,kl,ieg,xl,nseb_rfc(1,2))
                           if (IF3D) VZ(il,jl,kl,iel)=VZ(il,jl,kl,iel)+
     $                          nseb_amp*math_ran_dst(il,jl,kl,ieg,xl,
     $                                          nseb_rfc(1,ldim))
                        endif

                     enddo
                  enddo
               enddo
            enddo

            ! face averaging
            call opdssum(VX,VY,VZ)
            call opcolv (VX,VY,VZ,VMULT)

            ! timing
            ltim = dnekclock() - ltim
            call mntr_tmr_add(nseb_tmr_id,1,ltim)

         endif
      endif

      return
      end subroutine
!=======================================================================
