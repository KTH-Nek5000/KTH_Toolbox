!> @file spnb.f
!! @ingroup spnb
!! @brief Sponge/fringe for simple box mesh
!! @author Adam Peplinski
!! @date Feb 1, 2017
!=======================================================================
!> @brief Register sponge box module
!! @ingroup spnb
!! @note This routine should be called in frame_usr_register
      subroutine spnb_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'SPNBD'

      ! local variables
      integer lpmid
      real ltim

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,spnb_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(spnb_name)//'] already registered')
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
      call mntr_mod_reg(spnb_id,lpmid,spnb_name,
     $          'Sponge/fringe for rectangular domain')

      ! register timer
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      call mntr_tmr_reg(spnb_tmr_id,lpmid,spnb_id,
     $     'spnb_INI','Sponge calculation initialisation time',.false.)

      ! register and set active section
      call rprm_sec_reg(spnb_sec_id,spnb_id,'_'//adjustl(spnb_name),
     $     'Runtime paramere section for sponge box module')
      call rprm_sec_set_act(.true.,spnb_sec_id)

      ! register parameters
      call rprm_rp_reg(spnb_str_id,spnb_sec_id,'STRENGTH',
     $     'Sponge strength',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_wl_id(1),spnb_sec_id,'WIDTHLX',
     $     'Sponge left section width; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_wl_id(2),spnb_sec_id,'WIDTHLY',
     $     'Sponge left section width; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(spnb_wl_id(ndim),spnb_sec_id,
     $     'WIDTHLZ','Sponge left section width; dimension Z ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_wr_id(1),spnb_sec_id,'WIDTHRX',
     $     'Sponge right section width; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_wr_id(2),spnb_sec_id,'WIDTHRY',
     $     'Sponge right section width; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(spnb_wr_id(ndim),spnb_sec_id,
     $     'WIDTHRZ','Sponge right section width; dimension Z ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_dl_id(1),spnb_sec_id,'DROPLX',
     $     'Sponge left drop/rise section width; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_dl_id(2),spnb_sec_id,'DROPLY',
     $     'Sponge left drop/rise section width; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(spnb_dl_id(ndim),spnb_sec_id,
     $    'DROPLZ','Sponge left drop/rise section width; dimension Z ',
     $    rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_dr_id(1),spnb_sec_id,'DROPRX',
     $     'Sponge right drop/rise section width; dimension X ',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spnb_dr_id(2),spnb_sec_id,'DROPRY',
     $     'Sponge right drop/rise section width; dimension Y ',
     $     rpar_real,0,0.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(spnb_dr_id(ndim),spnb_sec_id,
     $   'DROPRZ','Sponge right drop/rise section width; dimension Z ',
     $    rpar_real,0,0.0,.false.,' ')

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(spnb_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise sponge box module
!! @ingroup spnb
!! @param[in] lvx, lvy, lvz   velocity field to be stored as reference field
!! @note This routine should be called in frame_usr_init
!! @remark This routine uses global scratch space \a SCRUZ
      subroutine spnb_init(lvx,lvy,lvz)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'FRAMELP'
      include 'SPNBD'

      ! argument list
      real lvx(LX1*LY1*LZ1*LELV),lvy(LX1*LY1*LZ1*LELV),
     $     lvz(LX1*LY1*LZ1*LELV)

      ! local variables
      integer ierr, nhour, nmin
      integer itmp
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      integer ntot, il, jl
      real bmin(LDIM), bmax(LDIM)

      real xxmax, xxmax_c, xxmin, xxmin_c, arg
      real lcoord(LX1*LY1*LZ1*LELV)
      common /SCRUZ/ lcoord

      ! functions
      real dnekclock, glmin, glmax, math_stepf
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (spnb_ifinit) then
         call mntr_warn(spnb_id,
     $        'module ['//trim(spnb_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_str_id,rpar_real)
      spnb_str = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wl_id(1),rpar_real)
      spnb_wl(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wl_id(2),rpar_real)
      spnb_wl(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wl_id(ndim),
     $        rpar_real)
         spnb_wl(ndim) = rtmp
      endif

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wr_id(1),rpar_real)
      spnb_wr(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wr_id(2),rpar_real)
      spnb_wr(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_wr_id(ndim),
     $        rpar_real)
         spnb_wr(ndim) = rtmp
      endif

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dl_id(1),rpar_real)
      spnb_dl(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dl_id(2),rpar_real)
      spnb_dl(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dl_id(ndim),
     $        rpar_real)
         spnb_dl(ndim) = rtmp
      endif

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dr_id(1),rpar_real)
      spnb_dr(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dr_id(2),rpar_real)
      spnb_dr(2) = rtmp

      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spnb_dr_id(ndim),
     $        rpar_real)
         spnb_dr(ndim) = rtmp
      endif

      ! initialise sponge variables

      ! get box size
      ntot = NX1*NY1*NZ1*NELV
      bmin(1) = glmin(XM1,ntot)
      bmax(1) = glmax(XM1,ntot)
      bmin(2) = glmin(YM1,ntot)
      bmax(2) = glmax(YM1,ntot)
      if(IF3D) then
         bmin(NDIM) = glmin(ZM1,ntot)
         bmax(NDIM) = glmax(ZM1,ntot)
      endif

      ! zero spnb_fun
      call rzero(spnb_fun,ntot)

!
!                        spnb_fun
!                        *********
!
!            bmin   xxmin_c         xxmax_c  bmax
!              |     |  xxmin     xxmax |    |
!              |_____|  |          |    |____|
!            ^       \                 /
!            |        \     x->       /
!    spnb_str|         \             /
!            v          \___________/
!              <-------->           <-------->
!                  wl                   wr
!                    <-->           <-->
!                     dl             dr

      if(spnb_str.gt.0.0) then
         call mntr_log(spnb_id,lp_inf,"Sponge turned on")

         ! save reference field
         call copy(spnb_vr(1,1),lvx, ntot)
         call copy(spnb_vr(1,2),lvy, ntot)
         if (IF3D) call copy(spnb_vr(1,NDIM),lvz, ntot)

         ! for every dimension
         do il=1,NDIM

            if (spnb_wl(il).gt.0.0.or.spnb_wr(il).gt.0.0) then
               if (spnb_wl(il).lt.spnb_dl(il).or.
     $              spnb_wr(il).lt.spnb_dr(il)) then
                  call mntr_abort(spnb_id,"Wrong sponge parameters")
               endif

               ! sponge beginning (rise at xmax; right)
               xxmax = bmax(il) - spnb_wr(il)
               ! end (drop at xmin; left)
               xxmin = bmin(il) + spnb_wl(il)
               ! beginnign of constant part (right)
               xxmax_c = xxmax + spnb_dr(il)
               ! beginnign of constant part (left)
               xxmin_c = xxmin - spnb_dl(il)

               ! get spnb_FUN
               if (xxmax.le.xxmin) then
                  call mntr_abort(spnb_id,"Sponge too wide")
               else
                  ! this should be done by pointers, but for now I avoid it
                  if (il.eq.1) then
                     call copy(lcoord,XM1, ntot)
                  elseif (il.eq.2) then
                     call copy(lcoord,YM1, ntot)
                  elseif (il.eq.3) then
                     call copy(lcoord,ZM1, ntot)
                  endif

                  do jl=1,ntot
                     rtmp = lcoord(jl)
                     if(rtmp.le.xxmin_c) then ! constant; xmin
                        rtmp=spnb_str
                     elseif(rtmp.lt.xxmin) then ! fall; xmin
                        arg = (xxmin-rtmp)/(spnb_wl(il)-spnb_dl(il))
                        rtmp = spnb_str*math_stepf(arg)
                     elseif (rtmp.le.xxmax) then ! zero
                        rtmp = 0.0
                     elseif (rtmp.lt.xxmax_c) then ! rise
                        arg = (rtmp-xxmax)/(spnb_wr(il)-spnb_dr(il))
                        rtmp = spnb_str*math_stepf(arg)
                     else    ! constant
                        rtmp = spnb_str
                     endif
                     spnb_fun(jl)=max(spnb_fun(jl),rtmp)
                  enddo

               endif         ! xxmax.le.xxmin

            endif            ! spnb_w(il).gt.0.0
         enddo


      endif

#ifdef DEBUG
      ! for debugging
      ltmp = ifto
      ifto = .TRUE.
      call outpost2(spnb_vr,spnb_vr(1,2),spnb_vr(1,NDIM),spnb_fun,
     $              spnb_fun,1,'spg')
      ifto = ltmp
#endif

      ! is everything initialised
      spnb_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(spnb_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup spnb
!! @return spnb_is_initialised
      logical function spnb_is_initialised()
      implicit none

      include 'SIZE'
      include 'SPNBD'
!-----------------------------------------------------------------------
      spnb_is_initialised = spnb_ifinit

      return
      end function
!=======================================================================
!> @brief Get sponge forcing
!! @ingroup spnb
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine spnb_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      include 'SPNBD'

      ! argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

      ! local variables
      integer iel, ip
!-----------------------------------------------------------------------
      iel=gllel(ieg)
      if (spnb_str.gt.0.0) then
         ip=ix+nx1*(iy-1+ny1*(iz-1+nz1*(iel-1)))

         if (jp.eq.0) then
            ! dns
            ffx = ffx + spnb_fun(ip)*(spnb_vr(ip,1) - vx(ix,iy,iz,iel))
            ffy = ffy + spnb_fun(ip)*(spnb_vr(ip,2) - vy(ix,iy,iz,iel))
            if (if3d) ffz = ffz + spnb_fun(ip)*
     $           (spnb_vr(ip,ndim) - vz(ix,iy,iz,iel))
         else
            ! perturbation
            ffx = ffx - spnb_fun(ip)*vxp(ip,jp)
            ffy = ffy - spnb_fun(ip)*vyp(ip,jp)
            if(if3d) ffz = ffz - spnb_fun(ip)*vzp(ip,jp)
         endif

      endif

      return
      end subroutine
!=======================================================================
