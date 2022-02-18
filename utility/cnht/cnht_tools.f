!> @file cnht_tools.f
!! @ingroup cnht
!! @brief Set of utilities related to conjugated heat transfer to build
!!    single scalar product for velocity nad temperature.
!! @author Clio Saglietti, Adam Peplinski
!! @date Mar 4, 2019
!=======================================================================
!> @brief Register conjugated heat transfer tools module
!! @ingroup cnht
!! @note This routine should be called in frame_usr_register
      subroutine cnht_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'CNHTD'

      ! local variables
      integer lpmid, il
      character*2 str
!-----------------------------------------------------------------------
      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,cnht_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(cnht_name)//'] already registered')
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
      call mntr_mod_reg(cnht_id,lpmid,cnht_name,
     $      'Conjugated heat transfer tools')

      ! register and set active section
      call rprm_sec_reg(cnht_sec_id,cnht_id,'_'//adjustl(cnht_name),
     $     'Runtime paramere section for conj. heat trans. tool module')
      call rprm_sec_set_act(.true.,cnht_sec_id)

      ! register parameters
      call rprm_rp_reg(cnht_sc_id,cnht_sec_id,'SCLN',
     $     'Norm scaling factor',rpar_real,0,3.36558,.false.,' ')

      call rprm_rp_reg(cnht_sv_id,cnht_sec_id,'SCLV',
     $     'Velocity scaling factor (Pareto curve)',
     $     rpar_real,0,0.5,.false.,' ')

      call rprm_rp_reg(cnht_st_id,cnht_sec_id,'SCLT',
     $     'Temperature scaling factor (Pareto curve)',
     $     rpar_real,0,0.5,.false.,' ')

      call rprm_rp_reg(cnht_gx_id,cnht_sec_id,'GRX',
     $     'X component of gravitational field',
     $     rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(cnht_gy_id,cnht_sec_id,'GRY',
     $     'Y component of gravitational field',
     $     rpar_real,0,1.0,.false.,' ')

      if (IF3D) call rprm_rp_reg(cnht_gz_id,cnht_sec_id,'GRZ',
     $     'Z component of gravitational field',
     $     rpar_real,0,0.0,.false.,' ')

      ! set initialisation flag
      cnht_ifinit=.false.

      return
      end subroutine
!=======================================================================
!> @brief Initilise conjugated heat transfer tools  module
!! @ingroup cnht
!! @note This routine should be called in frame_usr_init
      subroutine cnht_init()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'
      include 'SOLN'
      include 'ADJOINT'
      include 'CNHTD'

      ! local variables
      integer itmp, il
      real rtmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (cnht_ifinit) then
         call mntr_warn(cnht_id,
     $        'module ['//trim(cnht_name)//'] already initiaised.')
         return
      endif

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_sc_id,rpar_real)
      cnht_sc = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_sv_id,rpar_real)
      cnht_sv = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_st_id,rpar_real)
      cnht_st = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_gx_id,rpar_real)
      cnht_gx = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_gy_id,rpar_real)
      cnht_gz = rtmp
      if (IF3D) then
         call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cnht_gz_id,rpar_real)
         cnht_gz = rtmp
      endif

      ! Rayleight and Prandtl numbers
      cnht_Ra = abs(PARAM(2))
      cnht_Ra = abs(PARAM(1))
      BETA_B = cnht_Ra

      ! gravity
      G_ADJ(1) = cnht_gx
      G_ADJ(2) = cnht_gy
      G_ADJ(3) = cnht_gz

      ! temperature gradient
      call gradm1(DTDX,DTDY,DTDZ,T)

      ! everything is initialised
      cnht_ifinit=.true.

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup cnht
!! @return cnht_is_initialised
      logical function cnht_is_initialised()
      implicit none

      include 'SIZE'
      include 'CNHTD'
!-----------------------------------------------------------------------
      cnht_is_initialised = cnht_ifinit

      return
      end function
!=======================================================================
!> @brief Calcualte forcing ralted to conjugated heat transfer
!! @ingroup cnht
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine cnht_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! IF3D, IFHEAT, CPFLD
      include 'PARALLEL'        ! GLLEL
      include 'TSTEP'           ! IFIELD
      include 'SOLN'            ! JP, T, TP
      include 'ADJOINT'         ! IFADJ, G_ADJ, DTD[XYZ]
      include 'CNHTD'           ! CHGR[XYZ]

!     argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

!     local variables
      integer iel, ip
      real rtmp
!-----------------------------------------------------------------------
      if (IFHEAT) then
         iel=GLLEL(ieg)
         if (JP.eq.0) then
            rtmp = T(ix,iy,iz,iel,IFIELD)/CPFLD(1,2)
            ffx = ffx + cnht_gx*rtmp
            ffy = ffy + cnht_gy*rtmp
            if (IF3D) ffz = ffz + cnht_gz*rtmp
         else
            ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1)))
            if (.not.IFADJ) then
               rtmp = TP(ip,IFIELD,JP)/CPFLD(1,2)
               ffx = ffx + G_ADJ(1)*rtmp
               ffy = ffy + G_ADJ(2)*rtmp
               if (IF3D) ffz = ffz + G_ADJ(3)*rtmp
            else
               ffx = ffx - DTDX(ip)*TP(ip,IFIELD,JP)
               ffy = ffy - DTDY(ip)*TP(ip,IFIELD,JP)
               if (IF3D) ffz = ffz - DTDZ(ip)*TP(ip,IFIELD,JP)
            end if
         end if
      endif

      return
      end subroutine
!=======================================================================
!> @brief Set cpfld coefficient for given type of simulation
!! @ingroup cnht
      subroutine cnht_cpfld_set()
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! CPFLD, PARAM
      include 'ADJOINT'         ! IFADJ
      include 'CNHTD'          ! cnht_Ra, cnht_Ra
!-----------------------------------------------------------------------
      if (IFHEAT) then
         if (IFADJ) then
            CPFLD(1,1)=cnht_Ra/sqrt(cnht_Ra)
            CPFLD(1,2)=1.0

            CPFLD(2,1)=1.0/sqrt(cnht_Ra)
            CPFLD(2,2)=1.0
         else
            CPFLD(1,1)=1.0/sqrt(cnht_Ra)
            CPFLD(1,2)=1.0/cnht_Ra

            CPFLD(2,1)=1.0/sqrt(cnht_Ra)
            CPFLD(2,2)=1.0
         endif
      else
         if (PARAM(2).lt.0.0) then
            CPFLD(1,1) = -1.0/PARAM(2)
         else
            CPFLD(1,1) = PARAM(2)
         endif

         if (PARAM(1).lt.0.0) then
            CPFLD(1,2) = -1.0/PARAM(1)
         else
            CPFLD(1,2) = PARAM(1)
         endif
      endif

      return
      end subroutine
!=======================================================================
!> @brief Zero velocity and temperature vectors
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
      subroutine cnht_oprzero (a1,a2,a3,a4)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1)

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call rzero(a1,ntotv)
         call rzero(a2,ntotv)
         if(IF3D) call rzero(a3,ntotv)
      endif
      if (IFHEAT) call rzero(a4,ntott)
      return
      end subroutine
!=======================================================================
!> @brief Copy vectors A=B (velocity and temperature)
!! @ingroup cnht
!! @param[out] a1, a2, a3    vlocity field 1
!! @param[out] a4            temperature field 1
!! @param[in]  b1, b2, b3    vlocity field 2
!! @param[in]  b4            temperature field 2
      subroutine cnht_opcopy (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call copy(a1,b1,ntotv)
         call copy(a2,b2,ntotv)
         if(IF3D) call copy(a3,b3,ntotv)
      endif
      if (IFHEAT) call copy(a4,b4,ntott)

      return
      end subroutine
!=======================================================================
!> @brief Add velocity and temperature vectors A = A+B
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
      subroutine cnht_opadd2 (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call add2(a1,b1,ntotv)
         call add2(a2,b2,ntotv)
         if(IF3D) call add2(a3,b3,ntotv)
      endif
      if (IFHEAT) call add2(a4,b4,ntott)

      return
      end subroutine
!=======================================================================
!> @brief Subtract vectors A = A-B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
      subroutine cnht_opsub2 (a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call sub2(a1,b1,ntotv)
         call sub2(a2,b2,ntotv)
         if(IF3D) call sub2(a3,b3,ntotv)
      endif
      if (IFHEAT) call sub2(a4,b4,ntott)

      return
      end subroutine
!=======================================================================
!> @brief Subtract vectors A = B-C (velocity and temperature)
!! @ingroup cnht
!! @param[out] a1, a2, a3    vlocity field 1
!! @param[out] a4            temperature field 1
!! @param[in]  b1, b2, b3    vlocity field 2
!! @param[in]  b4            temperature field 2
!! @param[in]  c1, c2, c3    vlocity field 3
!! @param[in]  c4            temperature field 3
      subroutine cnht_opsub3 (a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real c1(1),c2(1),c3(1),c4(1)

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call sub3(a1,b1,c1,ntotv)
         call sub3(a2,b2,c2,ntotv)
         if(IF3D) call sub3(a3,b3,c3,ntotv)
      endif
      if (IFHEAT) call sub3(a4,b4,c4,ntott)
      return
      end subroutine
!=======================================================================
!> @brief Multiply vector by constant A = c*A (single coeff. for velocity
!!    and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity fields
!! @param[inout] a4            temperature field
!! @param[in]    const         coefficient
      subroutine cnht_opcmult (a1,a2,a3,a4,const)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1)
      real const

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call cmult(a1,const,ntotv)
         call cmult(a2,const,ntotv)
         if(IF3D) call cmult(a3,const,ntotv)
      endif
      if (IFHEAT) call cmult(a4,const,ntott)
      return
      end subroutine
!=======================================================================
!> @brief Multiply vector by constant A = c*A with separate const. for
!!    velocity and temperature
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity fields
!! @param[inout] a4            temperature field
!! @param[in]    const1        velocity coefficient
!! @param[in]    const2        temperature coefficient
      subroutine cnht_opcmult2c (a1,a2,a3,a4,const1, const2)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1)
      real const1, const2

!     local variables
      integer ntotv, ntott
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         call cmult(a1,const1,ntotv)
         call cmult(a2,const1,ntotv)
         if(IF3D) call cmult(a3,const1,ntotv)
      endif
      if (IFHEAT) call cmult(a4,const2,ntott)
      return
      end subroutine
!=======================================================================
!> @brief  Vector summation with scaling A = A+c*B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
!! @param[in]    coeff         scaling coefficient
      subroutine cnht_opadd2cm (a1,a2,a3,a4,b1,b2,b3,b4,coeff)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real coeff

!     local variables
      integer ntotv, ntott
      integer il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         if (IF3D) then
            do il=1,ntotv
               a1(il) = a1(il) + b1(il)*coeff
               a2(il) = a2(il) + b2(il)*coeff
               a3(il) = a3(il) + b3(il)*coeff
            enddo
         else
            do il=1,ntotv
               a1(il) = a1(il) + b1(il)*coeff
               a2(il) = a2(il) + b2(il)*coeff
            enddo
         endif
      endif
      if (IFHEAT) then
         do il=1,ntott
            a4(il) = a4(il) + b4(il)*coeff
         enddo
      endif
      return
      end subroutine
!=======================================================================
!> @brief  Vector subtraction with scaling A = A-c*B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
!! @param[in]    coeff         scaling coefficient
      subroutine cnht_opsub2cm (a1,a2,a3,a4,b1,b2,b3,b4,coeff)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
      real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
      real coeff

!     local variables
      integer ntotv, ntott
      integer il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      if (IFFLOW) then
         if (IF3D) then
            do il=1,ntotv
               a1(il) = a1(il) - b1(il)*coeff
               a2(il) = a2(il) - b2(il)*coeff
               a3(il) = a3(il) - b3(il)*coeff
            enddo
         else
            do il=1,ntotv
               a1(il) = a1(il) - b1(il)*coeff
               a2(il) = a2(il) - b2(il)*coeff
            enddo
         endif
      endif
      if (IFHEAT) then
         do il=1,ntott
            a4(il) = a4(il) - b4(il)*coeff
         enddo
      endif
      return
      end subroutine
!=======================================================================
!> @brief Weigth velocity and temperature fields
!! @ingroup cnht
!! @param[inout] lvx, lvy, lvz    vlocity fields
!! @param[inout] lt               temperature field
!! @param[in]    coeff            scaling coefficient
      subroutine cnht_weight_fun (lvx,lvy,lvz,lt,coeff)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'MASS'            ! VOLVM1, VOLTM1
      include 'CNHTD'           ! cnht_sc,cnht_sv,cnht_st

      ! argument list
      real lvx(1),lvy(1),lvz(1),lt(1)
      real coeff

      ! local variables
      real f1, f2
!-----------------------------------------------------------------------
      f1=cnht_sv/VOLVM1/coeff
      f2=cnht_st*cnht_sc/VOLTM1/coeff

      !rescale
      call cnht_opcmult2c (lvx,lvy,lvz,lt,f1,f2)

      return
      end subroutine
!=======================================================================
!> @brief Global inner product of velocity and temperature fields
!! @ingroup cnht
!! @param[in] b1, b2, b3    vlocity field 1
!! @param[in] b4            temperature field 1
!! @param[in] x1, x2, x3    vlocity field 2
!! @param[in] x4            temperature field 2
!! @param[in] wt            mass matrix
!! @return cnht_glsc2_wt
      real function cnht_glsc2_wt (b1,b2,b3,b4,x1,x2,x3,x4,wt)
      implicit none

      include 'SIZE'            ! N[XYZ]1, NEL[VT]
      include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
      include 'MASS'            ! VOLVM1, VOLTM1
      include 'CNHTD'           ! cnht_sc,cnht_sv,cnht_st

!     argument list
      real b1(1),b2(1),b3(1),b4(1),x1(1),x2(1),x3(1),x4(1),wt(1)

!     local variables
      integer ntotv, ntott
      real sum, f1, f2
      integer il
!     functions
      real glsum
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT

      ! scaling factor velocity vs temperature
      ! veorsion for newton
      !  f1 = coeff_v
      !  f2 = coeff_T
      ! version for oic
      f1=cnht_sv/VOLVM1
      f2=cnht_st*cnht_sc/VOLTM1

      sum = 0.
      if (IFFLOW) then          !if vel
         if (IFHEAT) then       !if temp & vel
            if (IF3D) then
               do il=1,ntotv
                  sum = sum + wt(il)*(f1*(b1(il)*x1(il)+b2(il)*x2(il)
     &                 +b3(il)*x3(il))+f2*b4(il)*x4(il))
               end do
            else
               do il=1,ntotv
                  sum =sum + wt(il)*(f1*(b1(il)*x1(il)+b2(il)*x2(il))
     &                 +f2*b4(il)*x4(il))
               end do
            end if

            ! for conjugate heat transfer
            if (ntott.gt.ntotv) then
               do il=ntotv+1,ntott
                  sum = sum + wt(il)*f2*b4(il)*x4(il)
               end do
            end if
        else                   !just vel
           if (IF3D) then
              do il=1,ntotv
                 sum = sum + wt(il)*f1*(b1(il)*x1(il)+
     $                b2(il)*x2(il)+b3(il)*x3(il))
              end do
           else
              do il=1,ntotv
                 sum = sum + wt(il)*f1*(b1(il)*x1(il)+b2(il)*x2(il))
              end do
           end if
        end if
      else                      !just temp
         if (IFHEAT) then
            do il=1,ntott
               sum = sum + wt(il)*(f2*b4(il)*x4(il))
            end do
         end if
      end if
      
      cnht_glsc2_wt = glsum(sum,1)
      
      return
      end function
!=======================================================================
