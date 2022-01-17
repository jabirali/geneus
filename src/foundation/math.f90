!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This file provides a common interface to a large library of mathematical
!> functions and subroutines. See the documentation of individual interfaces
!> below for more information about the contents of the mathematical library.

module math_m
  use :: basic_m
  use :: calculus_m
  use :: matrix_m
  private

  !------------------------------------------------------------
  ! Public interface
  !------------------------------------------------------------

  ! Floating-point precision
  public :: sp, dp, qp, wp

  ! Numerical constants
  public :: inf, eps, pi

  ! Basic routines
  public :: unitvector, nonzero, re, im, cx, arg

  ! Matrix routines
  public :: trace, diag, inverse, identity, commutator, anticommutator

  ! Calculus routines
  public :: mean, differentiate, integrate, interpolate, linspace

  !------------------------------------------------------------
  ! Interfaces for working with matrices
  !------------------------------------------------------------

  interface trace
    !! Public interface for functions that calculate a matrix trace.
    module procedure matrix_trace     ! Trace of a matrix
  end interface

  interface inverse
    !! Public interface for functions that calculate a matrix inverse.
    module procedure matrix_inverse_re, &  ! Inverse of a real matrix
                     matrix_inverse_cx     ! Inverse of a complex matrix
  end interface
  
  interface diag
    !! Public interface for functions that deal with matrix diagonals.
    module procedure matrix_diag, &   ! Construct a diagonal matrix
                     vector_diag      ! Extract a matrix diagonal
  end interface

  !------------------------------------------------------------
  ! Interfaces for calculus routines
  !------------------------------------------------------------

  interface mean
    !! Public interface for routines that calculate the mean value.
    module procedure mean_array_re, &               ! Mean of a real array
                     mean_array_cx                  ! Mean of a complex array
  end interface

  interface differentiate
    !! Public interface for various differentiation routines.
    module procedure differentiate_array_re, &      ! Derivative of a real array
                     differentiate_array_cx         ! Derivative of a complex array
  end interface

  interface integrate
    !! Public interface for various integration routines.
    module procedure integrate_array_re, &          ! Trapezoid integration of a real array
                     integrate_array_cx, &          ! Trapezoid integration of a complex array
                     integrate_range_re, &          ! Interpolated integration of a real array
                     integrate_range_cx             ! Interpolated integration of a complex array
  end interface

  interface interpolate
    !! Public interface for various interpolation routines.
    module procedure interpolate_point_re, &        ! Interpolate a scalar function from a real array to a point
                     interpolate_point_cx, &        ! Interpolate a scalar function from a complex array to a point
                     interpolate_array_re, &        ! Interpolate a scalar function from a real array to an array
                     interpolate_array_cx, &        ! Interpolate a scalar function from a complex array to an array
                     interpolate_point_matrix_re    ! Interpolate a matrix function from a real array array to a point
  end interface 

  interface linspace
    !! Public interface for routines that initialize arrays.
    module procedure linspace_array_re              ! Linearly spaced real numbers
  end interface

end module
