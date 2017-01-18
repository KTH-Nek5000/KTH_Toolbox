!> @file chkptd.f
!! @ingroup chkptdummy
!! @brief Dummy routines for checkpointing.
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Dummy replacement for checkpoint initialisation.
!! @ingroup chkptdummy
      subroutine chkpt_init
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for checkpoint reader.
!! @ingroup chkptdummy
      subroutine chkpt_read
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
!> @brief Dummy replacement for checkpoint writer.
!! @ingroup chkptdummy
      subroutine chkpt_write
      implicit none
!-----------------------------------------------------------------------
      return
      end
!=======================================================================
