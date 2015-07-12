! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers using Fortran 2008 syntax. To declare the floating point precision of a variable, use real(sp),
! real(dp), or real(qp) as the type of the variable. For more information about this topic, see http://fortranwiki.org.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

module module_precision
  use, intrinsic :: iso_fortran_env
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
end module 
