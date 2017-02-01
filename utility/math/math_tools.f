!> @file math_tools.f
!! @ingroup math
!! @brief Set of math related tools for KTH modules
!! @author Adam Peplinski
!! @date Jan 31, 2017
!=======================================================================
!> @brief Step function
!! @ingroup math
!! @details Continuous step function:
!!  \f{eqnarray*}{
!!    stepf(x) = \left\{ \begin{array}{ll}
!!  0 &\mbox{ if $x \leq x_{min}$} \\
!!  \left(1+e^{\left((x-1)^{-1} + x^{-1}\right)}\right)^{-1} &\mbox{ if $x \leq x_{max}$} \\
!!  1 &\mbox{ if $x >  x_{max}$}
!!       \end{array} \right.
!!  \f}
!!  with \f$ x_{min} = 0.02\f$ and \f$ x_{max}=0.98\f$
!! @param[in] x       function argument
!! @return mth_stepf
      real function mth_stepf(x)
      implicit none

!     argument list
      real x

!     local variables
      real xdmin, xdmax
      parameter (xdmin = 0.001, xdmax = 0.999)
!-----------------------------------------------------------------------
!     get function vale
      if (x.le.xdmin) then
         mth_stepf = 0.0
      else
         if (x.le.xdmax) then
            mth_stepf = 1./( 1. + exp(1./(x - 1.) + 1./x) )
         else
            mth_stepf = 1.
         end if
      end if

      end function mth_stepf
!=======================================================================
