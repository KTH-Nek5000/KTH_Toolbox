!> @file spongeb.f
!! @ingroup sponge_box
!! @brief Sponge/fringe for simple box mesh
!! @author Adam Peplinski
!! @date Feb 1, 2017
!=======================================================================
!> @brief Set sponbge_box parameters
!! @ingroup sponge_box
      subroutine spng_box_param_get()
      implicit none

      include 'SIZE'
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'SPONGEBD'

!     local variables
!     to operate with dictionary
      integer i_out,ifnd
      real d_out
!     for broadcasting
      integer llen, il, ip, lsent
      parameter (lsent = (1+3*LDIM))
      real rtmp(lsent)
!-----------------------------------------------------------------------
!     default values
      SPNG_STR = 1.0
      SPNG_W(1) = 5.0
      SPNG_WL(1) = 2.0
      SPNG_WR(1) = 2.0
      do il=2,NDIM
         SPNG_W(il) = 0.0
         SPNG_WL(il) = 0.0
         SPNG_WR(il) = 0.0
      enddo

!     check dictionary
      if (NID.eq.0) then
!     sponge strength
        call finiparser_getDbl(d_out,'_spongeb:strength',ifnd)
        if (ifnd.eq.1) then
           SPNG_STR = d_out
        endif

!     sponge total width
!     X
        call finiparser_getDbl(d_out,'_spongeb:widthx',ifnd)
        if (ifnd.eq.1) then
           SPNG_W(1) = d_out
        endif
!     Y
        call finiparser_getDbl(d_out,'_spongeb:widthy',ifnd)
        if (ifnd.eq.1) then
           SPNG_W(2) = d_out
        endif
!     Z
        if (NDIM.eq.3) then
           call finiparser_getDbl(d_out,'_spongeb:widthz',ifnd)
           if (ifnd.eq.1) then
              SPNG_W(NDIM) = d_out
           endif
        endif

!     sponge left section width
!     X
        call finiparser_getDbl(d_out,'_spongeb:widthlx',ifnd)
        if (ifnd.eq.1) then
           SPNG_WL(1) = d_out
        endif
!     Y
        call finiparser_getDbl(d_out,'_spongeb:widthly',ifnd)
        if (ifnd.eq.1) then
           SPNG_WL(2) = d_out
        endif
!     Z
        if (NDIM.eq.3) then
           call finiparser_getDbl(d_out,'_spongeb:widthlz',ifnd)
           if (ifnd.eq.1) then
              SPNG_WL(NDIM) = d_out
           endif
        endif

!     sponge right section width
!     X
        call finiparser_getDbl(d_out,'_spongeb:widthrx',ifnd)
        if (ifnd.eq.1) then
           SPNG_WR(1) = d_out
        endif
!     Y
        call finiparser_getDbl(d_out,'_spongeb:widthry',ifnd)
        if (ifnd.eq.1) then
           SPNG_WR(2) = d_out
        endif
!     Z
        if (NDIM.eq.3) then
           call finiparser_getDbl(d_out,'_spongeb:widthrz',ifnd)
           if (ifnd.eq.1) then
              SPNG_WR(NDIM) = d_out
           endif
        endif
      endif

!     broadcast data
      if (NID.eq.0) then
         ip = 1
         rtmp(ip) = SPNG_STR
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_W(il)
         enddo
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_WL(il)
         enddo
         do il=1,NDIM
            ip = ip +1
            rtmp(ip) = SPNG_WR(il)
         enddo
      endif
      llen = lsent*WDSIZE
      call bcast(rtmp,llen)
      if (NID.ne.0) then
         ip = 1
         SPNG_STR = rtmp(ip)
         do il=1,NDIM
            ip = ip +1
            SPNG_W(il) = rtmp(ip)
         enddo
         do il=1,NDIM
            ip = ip +1
            SPNG_WL(il) = rtmp(ip)
         enddo
         do il=1,NDIM
            ip = ip +1
            SPNG_WR(il) = rtmp(ip)
         enddo
      endif

      return
      end
!=======================================================================
!> @brief Init sponge variables in the box domain
!! @ingroup sponge_box
!! @param[in] lvx, lvy, lvz   initial
      subroutine spng_box_init(lvx,lvy,lvz)
      implicit none

      include 'SIZE'
      include 'GEOM'            ! [XYZ]M1
      include 'INPUT'           ! IF3D
      include 'SPONGEBD'        !

!     argument list
!     reference velocity field
      real lvx(LX1*LY1*LZ1*LELV),lvy(LX1*LY1*LZ1*LELV),
     $     lvz(LX1*LY1*LZ1*LELV)

!     local variables
      integer ntot, il, jl
      real bmin(LDIM), bmax(LDIM)

      real rtmp, xxmax, xxmax_c, xxmin, xxmin_c, arg
      real lcoord(LX1*LY1*LZ1*LELV)

      logical ltmp

!     functions
      real glmin, glmax, mth_stepf
!-----------------------------------------------------------------------
!     get box size
      ntot = NX1*NY1*NZ1*NELV
      bmin(1) = glmin(XM1,ntot)
      bmax(1) = glmax(XM1,ntot)
      bmin(2) = glmin(YM1,ntot)
      bmax(2) = glmax(YM1,ntot)
      if(IF3D) then
         bmin(NDIM) = glmin(ZM1,ntot)
         bmax(NDIM) = glmax(ZM1,ntot)
      endif

!     zero SPNG_FUN
      call rzero(SPNG_FUN,ntot)
            
      if(SPNG_STR.gt.0.0) then

!     stamp the file
         if (NIO.eq.0) then
            write(6,*) 'SPONGE TURNED ON'
            write(6,*) 'sponge strenght = ' , SPNG_STR
            write(6,*) 'sponge width = ', SPNG_W
            write(6,*) 'sponge drop width = ', SPNG_WL
            write(6,*) 'sponge rise width = ', SPNG_WR
         endif

!     save reference field
         call copy(SPNG_VR(1,1),lvx, ntot)
         call copy(SPNG_VR(1,2),lvy, ntot)
         if (IF3D) call copy(SPNG_VR(1,NDIM),lvz, ntot)

!     for every dimension
         do il=1,NDIM

            if (SPNG_W(il).gt.0.0) then
               if (SPNG_W(il).lt.(SPNG_WL(il)+SPNG_WR(il)).or.
     $              SPNG_WL(il)+SPNG_WR(il).eq.0.0) then
                  if (NIO.eq.0) then
                     write(6,*) 'ERROR; wrong sponge parameters'
                  endif
                  call exitt
               endif
     
               rtmp = SPNG_W(il)/(SPNG_WL(il)+SPNG_WR(il))
!     sponge beginning (rise at xmax; right)
               xxmax = bmax(il) - SPNG_WR(il)*rtmp
!     end (drop at xmin; left)
               xxmin = bmin(il) + SPNG_WL(il)*rtmp
!     beginnign of constant part (right)
               xxmax_c = xxmax + SPNG_WR(il)
!     beginnign of constant part (left)
               xxmin_c = xxmin - SPNG_WL(il)

!     get SPNG_FUN
               if (xxmax.le.xxmin) then
                  if (NIO.eq.0) write(6,*) 'ERROR; sponge to wide'
                  call exitt
               else
!     this should be done by pointers, but for now I avoid it
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
                        rtmp=SPNG_STR
                     elseif(rtmp.lt.xxmin) then ! fall; xmin
                        arg = (xxmin-rtmp)/SPNG_WL(il)
                        rtmp = mth_stepf(arg)
                     elseif (rtmp.le.xxmax) then ! zero
                        rtmp = 0.0
                     elseif (rtmp.lt.xxmax_c) then ! rise
                        arg = (rtmp-xxmax)/SPNG_WR(il)
                        rtmp = mth_stepf(arg)
                     else    ! constant
                        rtmp = SPNG_STR
                     endif
                     SPNG_FUN(jl)=max(SPNG_FUN(jl),rtmp)
                  enddo

               endif         ! xxmax.le.xxmin

            endif            ! SPNG_W(il).gt.0.0
         enddo

      endif                  ! SPNG_STR.gt.0.0

#ifdef DEBUG
!     for debugging
      ltmp = ifto
      ifto = .TRUE.
      call outpost2(SPNG_VR,SPNG_VR(1,2),SPNG_VR(1,NDIM),SPNG_FUN,
     $              SPNG_FUN,1,'spg')
      ifto = ltmp
#endif

      return
      end

!=======================================================================
!> @brief Calcualte sponge forcing
!! @ingroup sponge_box
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine spng_box_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      include 'SPONGEBD'

!     argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

!     local variables
      integer e, ip
      real rtmp
!-----------------------------------------------------------------------
      if (SPNG_STR.gt.0.0) then
         e=GLLEL(ieg)
         ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(e-1)))

         if (JP.eq.0) then
!     dns
            ffx = ffx + SPNG_FUN(ip)*(SPNG_VR(ip,1) - VX(ix,iy,iz,e))
            ffy = ffy + SPNG_FUN(ip)*(SPNG_VR(ip,2) - VY(ix,iy,iz,e))
            if (IF3D) ffz = ffz + SPNG_FUN(ip)*
     $           (SPNG_VR(ip,NDIM) - VZ(ix,iy,iz,e))
         else
!     perturbation
            ffx = ffx - SPNG_FUN(ip)*VXP(ip,JP)
            ffy = ffy - SPNG_FUN(ip)*VYP(ip,JP)
            if(IF3D) ffz = ffz - SPNG_FUN(ip)*VZP(ip,JP)
         endif

      endif

      return
      end
!=======================================================================
