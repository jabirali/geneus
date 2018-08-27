!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This file defines functions that perform some common matrix operations.

module evaluate_m
  use :: math_m
  use :: stdio_m
  private

  ! Declare which routines to export
  public :: evaluate

  ! Declare public interfaces
  interface evaluate
    !! Public interface for various routines that evaluate mathematical expressions.
    module procedure evaluate_scalar_value,  evaluate_scalar_field, &
                     evaluate_vector_value,  evaluate_vector_field, &
                     evaluate_logical_value, evaluate_integer_value
  end interface
contains
  impure subroutine evaluate_logical_value(expression, value)
    !! This subroutine takes a scalar logical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_logical_value('F', output)
    !!     call evaluate_logical_value('T', output)
    !!

    character(*), intent(in)  :: expression !! Either the character 'T' or 'F'
    logical,      intent(out) :: value      !! Result of parsing the expression

    select case(expression)
      case ('T', 't')
        value = .true.
      case ('F', 'f')
        value = .false.
      case default
        call error('Invalid logical expression: "' // trim(expression) // '"')
    end select
  end subroutine

  impure subroutine evaluate_integer_value(expression, value)
    !! This subroutine takes a scalar integer expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_integer_value('10', output)
    !!

    character(*), intent(in)  :: expression !! String containing an integer
    integer,      intent(out) :: value      !! Result of parsing the expression
    integer                   :: iostat     !! Error status from the parsing

    read(expression,*,iostat=iostat) value

    if (iostat /= 0) then
      call error('Invalid integer expression: "' // trim(expression) // '"')
    end if
  end subroutine

  impure subroutine evaluate_scalar_value(expression, value)
    !! This subroutine takes a scalar mathematical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('0',                    output)
    !!     call evaluate_scalar_value('sin(0.3*pi)*exp(-pi)', output)
    !!
    use :: fparser

    character(*), intent(in)  :: expression  !! Scalar-valued mathematical expression
    real(wp),     intent(out) :: value       !! Result of parsing the expression

    ! Make sure the expression is non-empty
    if (scan(expression, '0123456789pi') <= 0) then
      call error('Invalid scalar expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(1)
    call parsef(1, expression, ['pi'])

    ! Evaluate the parsed function
    value = evalf(1, [pi])
  end subroutine

  impure subroutine evaluate_scalar_field(expression, domain, value)
    !! This subroutine takes a scalar mathematical function of some variable 'z' as input,  along 
    !! with an array with discrete values for that variable 'z'.  It parses the provided function,
    !! evaluates it at each 'z'-value in the array, and then returns the discretized scalar field.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('0',                    input(1:n), output(1:n))
    !!     call evaluate_scalar_value('sin(pi*z)*exp(-pi*z)', input(1:n), output(1:n))
    !!
    use :: fparser

    character(*),           intent(in)  :: expression   !! Scalar-valued function of position 'z'
    real(wp), dimension(:), intent(in)  :: domain       !! Domain of the independent variable 'z'
    real(wp), dimension(:), allocatable :: value        !! Result of evaluating the field at each point of the domain
    integer                  :: n

    ! Make sure the expression is non-empty
    if (scan(expression, '0123456789piz') <= 0) then
      call error('Invalid scalar expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(1)
    call parsef(1, expression, ['pi', 'z '])

    ! Allocate memory for the output
    allocate(value(size(domain)))

    ! Evaluate the parsed function
    do n=1,size(domain)
      value(n) = evalf(1, [pi, domain(n)])
    end do
  end subroutine

  impure subroutine evaluate_vector_value(expression, value)
    !! This subroutine takes a vector mathematical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('[0,0,0]',                     output(1:3))
    !!     call evaluate_scalar_value('[sin(0.3*pi),0,cos(0,3*pi)]', output(1:3))
    !!
    use :: fparser

    character(*),           intent(in)  :: expression  !! Vector-valued mathematical expression
    real(wp), dimension(3), intent(out) :: value       !! Result of parsing the expression
    integer,  dimension(4)              :: sep

    ! Find the vector delimiters
    sep(1) = scan(expression, '[', back=.false.)
    sep(2) = scan(expression, ',', back=.false.)
    sep(3) = scan(expression, ',', back=.true. )
    sep(4) = scan(expression, ']', back=.true. )

    ! Make sure the expressions are non-empty
    if (sep(1) <= 0 .or. any(sep(2:4)-sep(1:3) <= 1)) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if
    if (scan(expression(sep(1)+1:sep(2)-1), '0123456789pi') <= 0 .or. &
        scan(expression(sep(2)+1:sep(3)-1), '0123456789pi') <= 0 .or. &
        scan(expression(sep(3)+1:sep(4)-1), '0123456789pi') <= 0) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(3)
    call parsef(1, expression(sep(1)+1:sep(2)-1), ['pi'])
    call parsef(2, expression(sep(2)+1:sep(3)-1), ['pi'])
    call parsef(3, expression(sep(3)+1:sep(4)-1), ['pi'])

    ! Evaluate the parsed function
    value(1) = evalf(1, [pi])
    value(2) = evalf(2, [pi])
    value(3) = evalf(3, [pi])
  end subroutine

  impure subroutine evaluate_vector_field(expression, domain, value)
    !! This subroutine takes a vector mathematical function of some variable 'z' as input,  along 
    !! with an array with discrete values for that variable 'z'.  It parses the provided function,
    !! evaluates it at each 'z'-value in the array, and then returns the discretized scalar field.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('[0,0,0]',                     input(1:n), output(1:3,1:n))
    !!     call evaluate_scalar_value('[sin(pi*z/2),0,cos(pi*z/2)]', input(1:n), output(1:3,1:n))
    !!
    use :: fparser

    character(*),             intent(in)  :: expression   !! Vector-valued function of position 'z'
    real(wp), dimension(:),   intent(in)  :: domain       !! Domain of the independent variable 'z'
    real(wp), dimension(:,:), allocatable :: value        !! Result of evaluating the field at each point of the domain
    integer,  dimension(4)                :: sep
    integer                               :: n

    ! Allocate memory for the output
    allocate(value(3,size(domain)))

    ! Find the vector delimiters
    sep(1) = scan(expression, '[', back=.false.)
    sep(2) = scan(expression, ',', back=.false.)
    sep(3) = scan(expression, ',', back=.true. )
    sep(4) = scan(expression, ']', back=.true. )

    ! Make sure the expressions are non-empty
    if (sep(1) <= 0 .or. any(sep(2:4)-sep(1:3) <= 1)) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if
    if (scan(expression(sep(1)+1:sep(2)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(2)+1:sep(3)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(3)+1:sep(4)-1), '0123456789piz') <= 0) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(3)
    call parsef(1, expression(sep(1)+1:sep(2)-1), ['pi', 'z '])
    call parsef(2, expression(sep(2)+1:sep(3)-1), ['pi', 'z '])
    call parsef(3, expression(sep(3)+1:sep(4)-1), ['pi', 'z '])

    ! Evaluate the parsed function
    do n=1,size(domain)
      value(1,n) = evalf(1, [pi, domain(n)])
      value(2,n) = evalf(2, [pi, domain(n)])
      value(3,n) = evalf(3, [pi, domain(n)])
    end do
  end subroutine
end module
