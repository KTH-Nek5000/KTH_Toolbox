!> @file chkptd.f
!! @ingroup chkptdummy
!! @brief Dummy routines for checkpointing.
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Dummy replacement for checkpoint registration.
!! @ingroup chkptdummy
      subroutine chkpts_register
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for checkpoint initialisation.
!! @ingroup chkptdummy
      subroutine chkpts_init
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for checkpoint reader.
!! @ingroup chkptdummy
      subroutine chkpts_read
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for checkpoint writer.
!! @ingroup chkptdummy
      subroutine chkpts_write
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
