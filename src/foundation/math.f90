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
    !!  Public interface for functions that calculate a matrix trace.
        module procedure matrix_trace
    end interface

    interface inverse
    !!  Public interface for functions that calculate a matrix inverse.
        module procedure matrix_inverse_re, matrix_inverse_cx
    end interface

    interface diag
    !!  Public interface for functions that deal with matrix diagonals.
        module procedure matrix_diag, vector_diag
    end interface

    !------------------------------------------------------------
    ! Interfaces for calculus routines
    !------------------------------------------------------------

    interface mean
    !!  Public interface for routines that calculate the mean value.
        module procedure mean_array_re, mean_array_cx
    end interface

    interface differentiate
    !!  Public interface for various differentiation routines.
        module procedure differentiate_array_re, differentiate_array_cx
    end interface

    interface integrate
    !!  Public interface for various integration routines.
        module procedure integrate_array_re, integrate_array_cx, &
                         integrate_range_re, integrate_range_cx
    end interface

    interface interpolate
    !!  Public interface for various interpolation routines.
        module procedure interpolate_point_re, interpolate_point_cx, &
                         interpolate_array_re, interpolate_array_cx, &
                         interpolate_point_matrix_re
    end interface

    interface linspace
    !!  Public interface for routines that initialize arrays.
        module procedure linspace_array_re
    end interface

end module
