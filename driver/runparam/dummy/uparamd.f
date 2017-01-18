!> @file uparamd.f
!! @ingroup upardummy
!! @brief Dummy routines for user runtime parameter reading/writing.
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Dummy replacement for user runtime parameter reader.
!! @ingroup upardummy
!! @param[in]  iunit    file unit
!! @note This interface is defined in @ref runprm_in
!! @see @ref readers_writers_page
      subroutine user_param_in(iunit)
      implicit none
!     argument list
      integer iunit
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for user runtime parameter writer.
!! @ingroup upardummy
!! @param[in]  iunit    file unit
!! @note This interface is defined in @ref runprm_out
!! @see @ref readers_writers_page
      subroutine user_param_out(iunit)
      implicit none
!     argument list
      integer iunit
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
