!> @file chkptdm.f
!! @ingroup chkptdm
!! @brief Dummy routines for checkpointing.
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Dummy replacement for checkpoint registration.
!! @ingroup chkptdm
      subroutine chkpts_register
      implicit none
!-----------------------------------------------------------------------
      return
      end subroutine
!=======================================================================
!> @brief Dummy replacement for checkpoint initialisation.
!! @ingroup chkptdm
      subroutine chkpts_init
      implicit none
!-----------------------------------------------------------------------
      return
      end subroutine
!=======================================================================
!> @brief Dummy replacement for check of module initialisation
!! @ingroup chkptdm
!! @return chkpts_is_initialised
      logical function chkpts_is_initialised()
      implicit none
!-----------------------------------------------------------------------
      chkpts_is_initialised = .true.

      return
      end function
!=======================================================================
!> @brief Dummy replacement for checkpoint reader.
!! @ingroup chkptdm
      subroutine chkpts_read
      implicit none
!-----------------------------------------------------------------------
      return
      end subroutine
!=======================================================================
!> @brief Dummy replacement for checkpoint writer.
!! @ingroup chkptdm
      subroutine chkpts_write
      implicit none
!-----------------------------------------------------------------------
      return
      end subroutine
!=======================================================================
