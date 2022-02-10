!> @file trip.f
!! @ingroup trip_line
!! @brief Tripping function supporting conformal and AMR version of nek5000
!! @details The tripping is based on a similar implementation in the
!!   SIMSON code (Chevalier et al. 2007, KTH Mechanics), and is described
!!   in detail in the paper Schlatter & Örlü, JFM 2012, DOI 10.1017/jfm.2012.324.
!! @author Adam Peplinski
!! @date May 03, 2018
!=======================================================================
!> @brief Register tripping module
!! @ingroup trip_line
!! @note This routine should be called in frame_usr_register
      subroutine trip_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'TRIPD'

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
      call mntr_mod_is_name_reg(lpmid,trip_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(trip_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(trip_id,lpmid,trip_name,
     $      'Tripping along the line')

      ! register timer
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      call mntr_tmr_reg(trip_tmr_id,lpmid,trip_id,
     $     'TRIP_TOT','Tripping total time',.false.)

      ! register and set active section
      call rprm_sec_reg(trip_sec_id,trip_id,'_'//adjustl(trip_name),
     $     'Runtime paramere section for tripping module')
      call rprm_sec_set_act(.true.,trip_sec_id)

      ! register parameters
      call rprm_rp_reg(trip_nline_id,trip_sec_id,'NLINE',
     $     'Number of tripping lines',rpar_int,0,0.0,.false.,' ')

      call rprm_rp_reg(trip_tiamp_id,trip_sec_id,'TIAMP',
     $     'Time independent amplitude',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(trip_tdamp_id,trip_sec_id,'TDAMP',
     $     'Time dependent amplitude',rpar_real,0,0.0,.false.,' ')

      do il=1, trip_nline_max
         write(str,'(I2.2)') il

         call rprm_rp_reg(trip_spos_id(1,il),trip_sec_id,'SPOSX'//str,
     $     'Starting point X',rpar_real,0,0.0,.false.,' ')
         
         call rprm_rp_reg(trip_spos_id(2,il),trip_sec_id,'SPOSY'//str,
     $     'Starting point Y',rpar_real,0,0.0,.false.,' ')

         if (IF3D) then
            call rprm_rp_reg(trip_spos_id(ldim,il),trip_sec_id,
     $           'SPOSZ'//str,'Starting point Z',
     $           rpar_real,0,0.0,.false.,' ')
         endif
        
         call rprm_rp_reg(trip_epos_id(1,il),trip_sec_id,'EPOSX'//str,
     $     'Ending point X',rpar_real,0,0.0,.false.,' ')
         
         call rprm_rp_reg(trip_epos_id(2,il),trip_sec_id,'EPOSY'//str,
     $     'Ending point Y',rpar_real,0,0.0,.false.,' ')

         if (IF3D) then
            call rprm_rp_reg(trip_epos_id(ldim,il),trip_sec_id,
     $           'EPOSZ'//str,'Ending point Z',
     $           rpar_real,0,0.0,.false.,' ')
         endif

         call rprm_rp_reg(trip_smth_id(1,il),trip_sec_id,'SMTHX'//str,
     $     'Smoothing length X',rpar_real,0,0.0,.false.,' ')
         
         call rprm_rp_reg(trip_smth_id(2,il),trip_sec_id,'SMTHY'//str,
     $     'Smoothing length Y',rpar_real,0,0.0,.false.,' ')

         if (IF3D) then
            call rprm_rp_reg(trip_smth_id(ldim,il),trip_sec_id,
     $           'SMTHZ'//str,'Smoothing length Z',
     $           rpar_real,0,0.0,.false.,' ')
         endif

         call rprm_rp_reg(trip_lext_id(il),trip_sec_id,'LEXT'//str,
     $        'Line extension',rpar_log,0,0.0,.false.,' ')
      
         call rprm_rp_reg(trip_rota_id(il),trip_sec_id,'ROTA'//str,
     $        'Rotation angle',rpar_real,0,0.0,.false.,' ')
         call rprm_rp_reg(trip_nmode_id(il),trip_sec_id,'NMODE'//str,
     $     'Number of Fourier modes',rpar_int,0,0.0,.false.,' ')
         call rprm_rp_reg(trip_tdt_id(il),trip_sec_id,'TDT'//str,
     $     'Time step for tripping',rpar_real,0,0.0,.false.,' ')
      enddo

      ! set initialisation flag
      trip_ifinit=.false.
      
      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(trip_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise tripping module
!! @ingroup trip_line
!! @note This routine should be called in frame_usr_init
      subroutine trip_init()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'FRAMELP'
      include 'TRIPD'

      ! local variables
      integer itmp
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      integer il, jl, kl
      real rtmpv(ldim)

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (trip_ifinit) then
         call mntr_warn(trip_id,
     $        'module ['//trim(trip_name)//'] already initiaised.')
         return
      endif
      
      ! timing
      ltim = dnekclock()

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_nline_id,rpar_int)
      trip_nline = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_tiamp_id,rpar_real)
      trip_tiamp = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_tdamp_id,rpar_real)
      trip_tdamp = rtmp
      do il=1,trip_nline
         do jl=1,LDIM
            call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_spos_id(jl,il),
     $           rpar_real)
            trip_spos(jl,il) = rtmp
            call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_epos_id(jl,il),
     $           rpar_real)
            trip_epos(jl,il) = rtmp
            call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_smth_id(jl,il),
     $           rpar_real)
            trip_smth(jl,il) = abs(rtmp)
         enddo
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_lext_id(il),
     $        rpar_log)
         trip_lext(il) = ltmp
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_rota_id(il),
     $        rpar_real)
         trip_rota(il) = rtmp
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_nmode_id(il),
     $        rpar_int)
         trip_nmode(il) = itmp
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,trip_tdt_id(il),
     $        rpar_real)
         trip_tdt(il) = rtmp
      enddo

      ! check simulation dimension
      if (.not.IF3D) call mntr_abort(trip_id,
     $        '2D simulation is not supported.')

      ! get line versor, inverse line lengths and scaled smoothing lengths
      do il=1,trip_nline
         call mntr_logi(trip_id,lp_inf,'Line info; line nr: ',il)
         trip_ilngt(il) = 0.0

         do jl=1,LDIM
            ! the last (third) versor is parallel to the line
            trip_vrs(jl,ldim,il) = trip_epos(jl,il) - trip_spos(jl,il)
            trip_ilngt(il) = trip_ilngt(il) + trip_vrs(jl,ldim,il)**2
         enddo
         trip_ilngt(il) = sqrt(trip_ilngt(il))
         call mntr_logr(trip_id,lp_inf,'Line length: ',trip_ilngt(il))
         if (trip_ilngt(il).gt.0.0) then
            trip_ilngt(il) = 1.0/trip_ilngt(il)
            do jl=1,LDIM
               trip_vrs(jl,ldim,il) = trip_vrs(jl,ldim,il)*
     $            trip_ilngt(il)
            enddo
         else
            call mntr_abort(trip_id,
     $        'Line with zero lenght is not supported.')
         endif
         ! the rest of versors given by cross product starting with the
         ! second versor
         call rzero(rtmpv,ldim)
         ! the first versor guess depends on the last versor coordinates
         if (trip_vrs(1,ldim,il).lt.0.95) then
            rtmpv(1) = 1.0
         else
            rtmpv(2) = 1.0
         endif
         ! the second versor
         call cross(trip_vrs(1,2,il),trip_vrs(1,ldim,il),rtmpv)
         ! the first versor
         call cross(trip_vrs(1,1,il),trip_vrs(1,2,il),
     $              trip_vrs(1,ldim,il))
         ! correct versor length
         do jl = 1, ldim
            rtmp = 0.0
            do kl = 1, ldim
               rtmp = rtmp + trip_vrs(kl,jl,il)*trip_vrs(kl,jl,il)
            end do
            if (rtmp.gt.0.0) then
               rtmp = 1.0/sqrt(rtmp)
               do kl = 1, ldim
                  trip_vrs(kl,jl,il) = trip_vrs(kl,jl,il)*rtmp
               end do
            else
               call mntr_abort(trip_id,
     $               'Line versor with zero lenght.')
            end if
         end do
         ! rotate the first and second versors along the third versor
         do jl = 1, 2
            call math_rot3da(rtmpv,trip_vrs(1,jl,il),trip_vrs(1,ldim,il)
     $       ,trip_rota(il))
            do kl = 1, ldim
               trip_vrs(kl,jl,il) = rtmpv(kl)
            end do
         end do

         ! stump the log
         do jl = 1, ldim
            call mntr_logi(trip_id,lp_inf,'Line versor: ',jl)
            call mntr_logrv(trip_id,lp_inf,'Coordinates:',
     $           trip_vrs(1,jl,il),ldim)
         end do

         ! rescale smoothing lengths
         do jl=1,LDIM
           trip_smth(jl,il) = trip_smth(jl,il)*trip_ilngt(il)
         enddo
         ! get inverse smoothing lenght
         do jl=1,LDIM
            if (trip_smth(jl,il).gt.0.0) then
               trip_ismth(jl,il) = 1.0/trip_smth(jl,il)
            else
               trip_ismth(jl,il) = 1.0
            endif
         enddo
      enddo

      ! get 1D projection and array mapping
      call trip_1dprj

      ! initialise random generator seed and number of time intervals
      do il=1,trip_nline
         trip_seed(il) = -32*il
         trip_ntdt(il) = 1 - trip_nset_max
         trip_ntdt_old(il) = trip_ntdt(il)
      enddo
      
      ! generate random phases (time independent and time dependent)
      call trip_rphs_get

      ! get forcing
      call trip_frcs_get(.true.)
      
      ! everything is initialised
      trip_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(trip_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup trip_line
!! @return trip_is_initialised
      logical function trip_is_initialised()
      implicit none

      include 'SIZE'
      include 'TRIPD'
!-----------------------------------------------------------------------
      trip_is_initialised = trip_ifinit

      return
      end function
!=======================================================================
!> @brief Update tripping
!! @ingroup trip_line
      subroutine trip_update()
      implicit none

      include 'SIZE'
      include 'TRIPD'

      ! local variables
      real ltim
      
      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()      

      ! update random phases (time independent and time dependent)
      call trip_rphs_get

      ! update forcing
      call trip_frcs_get(.false.)

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(trip_tmr_id,1,ltim)

      return
      end subroutine      
!=======================================================================
!> @brief Compute tripping forcing
!! @ingroup trip_line
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine trip_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'TRIPD'

      ! argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

      ! local variables
      integer ipos,iel,il
      real ffn
!-----------------------------------------------------------------------
      iel=GLLEL(ieg)

      do il= 1, trip_nline
         ffn = trip_fsmth(ix,iy,iz,iel,il)
         if (ffn.gt.0.0) then
            ipos = trip_map(ix,iy,iz,iel,il)
            ffn = trip_ftrp(ipos,il)*ffn

            ! I assume forcing direction is given by the second versor
            ffx = ffx + ffn*trip_vrs(1,2,il)
            ffy = ffy + ffn*trip_vrs(2,2,il)
            ffz = ffz + ffn*trip_vrs(ldim,2,il)
         endif
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Reset tripping
!! @ingroup trip_line
      subroutine trip_reset()
      implicit none

      include 'SIZE'
      include 'TRIPD'

      ! local variables
      real ltim
      
      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()      

      ! get 1D projection and array mapping
      call trip_1dprj
      
      ! update forcing
      call trip_frcs_get(.true.)

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(trip_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Get 1D projection, array mapping and forcing smoothing
!! @ingroup trip_line
!! @details This routine supports straight lines given by their starting
!!    and ending points. Additional flagg allows to introuduce forcing
!!    periodicity or contain it between starting and ending points + smooting
!!    lenght in z
!! @remark This routine uses global scratch space \a CTMP0 and \a CTMP1
      subroutine trip_1dprj()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TRIPD'

      ! global memory access
      real lcoord(LX1*LY1*LZ1*LELT)
      common /CTMP0/ lcoord
      integer lmap(LX1*LY1*LZ1*LELT), lmap_el(4,LX1*LY1*LZ1*LELT)
      common /CTMP1/ lmap, lmap_el

      ! local variables
      integer nptot, itmp
      integer il, jl, kl, ll, ml, nl
      real rota, rtmp, epsl
      parameter (epsl = 1.0e-10)
      real rtmpv(ldim), rtmpc(ldim)
      ! functions
      real dot
!-----------------------------------------------------------------------
      nptot = NX1*NY1*NZ1*NELV
      
      ! for each line
      do il=1,trip_nline
         ! reset mapping array
         call ifill(trip_map(1,1,1,1,il),-1,nptot)
         ! initialise number of points per line
         trip_npoint(il) = 0
         ! initialize smoothing factor
         call rzero(trip_fsmth(1,1,1,1,il),nptot)
         ! initialize projected point position
         call rzero(trip_prj(1,il),nptot)

         ! Projection onto the line
         ! count points on the line
         itmp = 0
         do jl = 1, nelv
            do kl = 1, lz1
               do ll = 1, ly1
                  do ml = 1, lx1
                     ! get point position relative to the line start
                     rtmpv(1) = xm1(ml,ll,kl,jl)-trip_spos(1,il)
                     rtmpv(2) = ym1(ml,ll,kl,jl)-trip_spos(2,il)
                     rtmpv(ldim) = zm1(ml,ll,kl,jl)-trip_spos(ldim,il)
                     ! get point coordinates in the local line system
                     do nl = 1, ldim
                        rtmpc(nl) = dot(rtmpv,trip_vrs(1,nl,il),ldim)
                        rtmpc(nl) = rtmpc(nl)*trip_ilngt(il)
                     end do
                     ! distance from the line
                     ! 2D
                     rtmp = (rtmpc(1)*trip_ismth(1,il))**2+
     $                      (rtmpc(2)*trip_ismth(2,il))**2
                     ! do we extend a line beyond its ends
                     if (.not.trip_lext(il)) then
                        if (rtmpc(ldim).lt.0.0) then
                           rtmp = rtmp +
     $                            (rtmpc(ldim)*trip_ismth(ldim,il))**2
                        elseif (rtmpc(ldim).gt.1.0) then
                           rtmp = rtmp +
     $                      ((rtmpc(ldim)-1.0)*trip_ismth(ldim,il))**2
                        end if
                     end if

                     ! get smoothing profile
                     ! Gauss; cannot be used with lines not extended beyond their ending points
                     !trip_fsmth(itmp,jtmp,ktmp,eltmp,il) = exp(-4.0*rtmp)
                     ! limited support
                     if (rtmp.lt.1.0) then
                        trip_fsmth(ml,ll,kl,jl,il) =
     $                       exp(-rtmp)*(1-rtmp)**2
                        ! add the point to the list
                        itmp = itmp + 1
                        ! save data
                        ! coordinate along the line in line length unit
                        lcoord(itmp) = rtmpc(ldim)
                        ! point position
                        lmap_el(1,itmp) = jl
                        lmap_el(2,itmp) = kl
                        lmap_el(3,itmp) = ll
                        lmap_el(4,itmp) = ml
                     else
                        trip_fsmth(ml,ll,kl,jl,il) = 0.0
                     endif
                  end do
               end do
            end do
         end do

         if (itmp.ge.1) then
            ! point sorting acording to the last coordinate
            call sort(lcoord,lmap,itmp)

            ! identify unique points
            trip_npoint(il) = 1
            trip_prj(trip_npoint(il),il) = lcoord(1)
            ! generate mapping
            trip_map(lmap_el(4,1),lmap_el(3,1),
     $          lmap_el(2,1),lmap_el(1,1),il) = trip_npoint(il)
            do jl = 2, itmp
               ! compare positions along the line
               if((lcoord(jl)-trip_prj(trip_npoint(il),il)).gt.
     $              max(epsl,abs(epsl*lcoord(jl)))) then
                  trip_npoint(il) = trip_npoint(il) + 1
                  trip_prj(trip_npoint(il),il) = lcoord(jl)
               endif
               ! generate mapping
               trip_map(lmap_el(4,jl),lmap_el(3,jl),
     $          lmap_el(2,jl),lmap_el(1,jl),il) = trip_npoint(il)
            end do
         end if
      enddo

      return
      end subroutine      
!=======================================================================
!> @brief Generate set of random phases
!! @ingroup trip_line
      subroutine trip_rphs_get
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'TRIPD'
      
      ! local variables
      integer il, jl, kl
      integer itmp
      real trip_ran2

#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif
!-----------------------------------------------------------------------
      ! time independent part
      if (trip_tiamp.gt.0.0.and..not.trip_ifinit) then
         do il = 1, trip_nline
            do jl=1, trip_nmode(il)
               trip_rphs(jl,1,il) = 2.0*pi*trip_ran2(il)
            enddo
         enddo
      endif

      ! time dependent part
      do il = 1, trip_nline
         itmp = int(time/trip_tdt(il))
         call bcast(itmp,ISIZE) ! just for safety
         do kl= trip_ntdt(il)+1, itmp
            do jl= trip_nset_max,3,-1
               call copy(trip_rphs(1,jl,il),trip_rphs(1,jl-1,il),
     $              trip_nmode(il))
            enddo
            do jl=1, trip_nmode(il)
               trip_rphs(jl,2,il) = 2.0*pi*trip_ran2(il)
            enddo
         enddo
         ! update time interval
         trip_ntdt_old(il) = trip_ntdt(il)
         trip_ntdt(il) = itmp
      enddo

#ifdef DEBUG
      ! for testing
      ! to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='trp_rps.txt'//str1//'i'//str2)

      do il=1,trip_nmode(1)
         write(iunit,*) il,trip_rphs(il,1:4,1)
      enddo

      close(iunit)
#endif

      return
      end subroutine
!=======================================================================
!> @brief A simple portable random number generator
!! @ingroup trip_line
!! @details  Requires 32-bit integer arithmetic. Taken from Numerical
!!   Recipes, William Press et al. Gives correlation free random
!!   numbers but does not have a very large dynamic range, i.e only
!!   generates 714025 different numbers. Set seed negative for
!!   initialization
!! @param[in]   il      line number
!! @return      ran
      real function trip_ran2(il)
      implicit none

      include 'SIZE'
      include 'TRIPD'
      
      ! argument list
      integer il

      ! local variables
      integer iff(trip_nline_max), iy(trip_nline_max)
      integer ir(97,trip_nline_max)
      integer m,ia,ic,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /trip_nline_max*0/
!-----------------------------------------------------------------------
      ! initialise
      if (trip_seed(il).lt.0.or.iff(il).eq.0) then
         iff(il)=1
         trip_seed(il)=mod(ic-trip_seed(il),m)
         do j=1,97
            trip_seed(il)=mod(ia*trip_seed(il)+ic,m)
            ir(j,il)=trip_seed(il)
         end do
         trip_seed(il)=mod(ia*trip_seed(il)+ic,m)
         iy(il)=trip_seed(il)
      end if
      
      ! generate random number
      j=1+(97*iy(il))/m
      iy(il)=ir(j,il)
      trip_ran2=iy(il)*rm
      trip_seed(il)=mod(ia*trip_seed(il)+ic,m)
      ir(j,il)=trip_seed(il)

      end function
!=======================================================================
!> @brief Generate forcing along 1D line
!! @ingroup trip_line
!! @param[in] ifreset    reset flag
      subroutine trip_frcs_get(ifreset)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'TRIPD'

      ! argument list
      logical ifreset

#ifdef TRIP_PR_RST
      ! variables necessary to reset pressure projection for P_n-P_n-2
      integer nprv(2)
      common /orthbi/ nprv

      ! variables necessary to reset velocity projection for P_n-P_n-2
      include 'VPROJ'
#endif      
      ! local variables
      integer il, jl, kl, ll
      integer istart
      real theta0, theta
      logical ifntdt_dif

#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif
!-----------------------------------------------------------------------
      ! reset all
      if (ifreset) then
         if (trip_tiamp.gt.0.0) then
            istart = 1
         else
            istart = 2
         endif
         do il= 1, trip_nline
            do jl = istart, trip_nset_max
               call rzero(trip_frcs(1,jl,il),trip_npoint(il))
               do kl= 1, trip_npoint(il)
                  theta0 = 2*pi*trip_prj(kl,il)
                  do ll= 1, trip_nmode(il)
                     theta = theta0*ll
                     trip_frcs(kl,jl,il) = trip_frcs(kl,jl,il) +
     $                    sin(theta+trip_rphs(ll,jl,il))
                  enddo
               enddo
            enddo
         enddo
         ! rescale time independent part
         if (trip_tiamp.gt.0.0) then
            do il= 1, trip_nline
               call cmult(trip_frcs(1,1,il),trip_tiamp,trip_npoint(il))
            enddo
         endif
      else
         ! reset only time dependent part if needed
         ifntdt_dif = .FALSE.
         do il= 1, trip_nline
            if (trip_ntdt(il).ne.trip_ntdt_old(il)) then
               ifntdt_dif = .TRUE.
               do jl= trip_nset_max,3,-1
                  call copy(trip_frcs(1,jl,il),trip_frcs(1,jl-1,il),
     $                 trip_npoint(il))
               enddo
               call rzero(trip_frcs(1,2,il),trip_npoint(il))
               do jl= 1, trip_npoint(il)
                  theta0 = 2*pi*trip_prj(jl,il)
                  do kl= 1, trip_nmode(il)
                     theta = theta0*kl
                     trip_frcs(jl,2,il) = trip_frcs(jl,2,il) +
     $                    sin(theta+trip_rphs(kl,2,il))
                  enddo
               enddo
            endif
         enddo
         if (ifntdt_dif) then
#ifdef TRIP_PR_RST
            ! reset projection space
            ! pressure
            if (int(PARAM(95)).gt.0) then
               PARAM(95) = ISTEP
               nprv(1) = 0      ! veloctiy field only
            endif
            ! velocity
            if (int(PARAM(94)).gt.0) then
               PARAM(94) = ISTEP!+2
               ivproj(2,1) = 0
               ivproj(2,2) = 0
               if (IF3D) ivproj(2,3) = 0
            endif
#endif
         endif
      endif
      
      ! get tripping for current time step
      if (trip_tiamp.gt.0.0) then
         do il= 1, trip_nline
           call copy(trip_ftrp(1,il),trip_frcs(1,1,il),trip_npoint(il))
         enddo
      else
         do il= 1, trip_nline
            call rzero(trip_ftrp(1,il),trip_npoint(il))
         enddo
      endif
      ! interpolation in time
      do il = 1, trip_nline
         theta0= time/trip_tdt(il)-real(trip_ntdt(il))
         if (theta0.gt.0.0) then
            theta0=theta0*theta0*(3.0-2.0*theta0)
            !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
            do jl= 1, trip_npoint(il)
               trip_ftrp(jl,il) = trip_ftrp(jl,il) +
     $              trip_tdamp*((1.0-theta0)*trip_frcs(jl,3,il) +
     $              theta0*trip_frcs(jl,2,il))
            enddo
         else
            theta0=theta0+1.0
            theta0=theta0*theta0*(3.0-2.0*theta0)
            !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
            do jl= 1, trip_npoint(il)
               trip_ftrp(jl,il) = trip_ftrp(jl,il) +
     $              trip_tdamp*((1.0-theta0)*trip_frcs(jl,4,il) +
     $              theta0*trip_frcs(jl,3,il))
            enddo
         endif
      enddo

#ifdef DEBUG
      ! for testing
      ! to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='trp_fcr.txt'//str1//'i'//str2)

      do il=1,trip_npoint(1)
         write(iunit,*) il,trip_prj(il,1),trip_ftrp(il,1),
     $        trip_frcs(il,1:4,1)
      enddo

      close(iunit)
#endif
      
      return
      end subroutine
!=======================================================================
