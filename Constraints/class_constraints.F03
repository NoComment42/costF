

module class_constraints
implicit none

type Constraint
        Character(len=120), private :: constraint_label
        Character(len=120), private :: constraint_type
        Character(len=120), private :: function_type
        procedure(gen_func), pointer, private :: func
        integer, private :: nParams
        real, private :: t ! target
        real, allocatable, dimension(:), private :: params
        Character(len=10), allocatable, dimension(:), private :: pLabels
    
    contains
        procedure, public :: getConstraintLabel
        procedure, public :: getConstraintType
        procedure, public :: getFunctionType
        procedure, public :: getNParams
        procedure, public :: getParameters
        procedure, public :: printConstraint
        
        procedure, public :: eval
        procedure, private :: eval_equal

        procedure, private :: init_equal

        procedure, private :: init_gen_exp
        procedure, private :: generalized_exp
        procedure, private :: init_inv_poly
        procedure, private :: inverse_poly
end type Constraint


interface Constraint
    procedure constructor
end interface Constraint

! general function interface to represent all the possible f functions
! basically defines what input/output should look like
abstract interface 
    function gen_func(this,x) result(v)
    import Constraint
    implicit none
        class(Constraint), intent(in) :: this
        real, intent(in) :: x
        real :: v
    end function gen_func
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(Constraint) function constructor (c_type,f_type,params,label)
    implicit none
        Character(len=*), intent(in) :: c_type, f_type
        Character(len=*), optional :: label
        real, dimension(:) :: params
        Character(len=120) :: tlabel
  
        if ( present(label) ) then
            tlabel = label
        else
            tlabel = "Not Labeled"
        end if

        select case(trim(c_type))
            case("equal")
                call constructor%init_equal(c_type,f_type,params,tlabel)
            case DEFAULT
                write(*,*) "Not a recognized constraint"
        end select
    end function constructor


    subroutine init_equal (this,c_type,f_type,params,label)
    implicit none
        class(Constraint) :: this
        Character(len=*), intent(in) :: c_type, f_type, label
        real, dimension(:) :: params

        this%constraint_label = label
        this%constraint_type = c_type
        this%function_type = f_type
        
        select case(trim(f_type))
            case("generalized_exp")
                this%func => generalized_exp
                call this%init_gen_exp(params)
            case("inverse_poly")
                this%func => inverse_poly
                call this%init_inv_poly(params)
            case DEFAULT
                this%func => NULL()
                write(*,*) "Not recognized function type, function pointer &
                           &is still set to NULL"
        end select

    end subroutine init_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function getConstraintLabel (this) result(answer)
    implicit none
        class(Constraint), intent(in) :: this
        Character(len=120) :: answer
        answer = this%constraint_label
    end function getConstraintLabel

    function getConstraintType (this) result(answer)
    implicit none
        class(Constraint), intent(in) :: this
        Character(len=120) :: answer 
        answer = this%constraint_type
    end function getConstraintType

    function getFunctionType (this) result(answer)
    implicit none
        class(Constraint), intent(in) :: this
        Character(len=120) answer
        answer = this%function_type
    end function getFunctionType


    function getNParams (this) result(np)
    implicit none
        class(Constraint), intent(in) :: this
        integer :: np
        np = this%nParams
    end function getNParams

    subroutine getParameters (this,p,l)
    implicit none
        class(Constraint), intent(in) :: this
        real, allocatable, dimension(:) :: p
        Character(len=10), allocatable, dimension(:) :: l

        if( allocated(p) ) deallocate(p)
        allocate( p(size(this%params)) )
        p = this%params
        
        if( allocated(l) ) deallocate(l)
        allocate( l(size(this%pLabels)) )
        l = this%pLabels

    end subroutine getParameters


    function printConstraint (this,label,npts,tin) result(err)
    implicit none
        class(Constraint), intent(in) :: this
        Character(len=*), intent(in) :: label
        integer, intent(in) :: npts
        real, optional :: tin
        integer :: err

        integer :: i
        real :: x, x1, x2, dx
        real :: tol
        integer,parameter :: maxC = 25
        real, parameter :: step = 1e-1
        real, parameter :: delx = 1e-6
        
        tol = 1e-3
        if( present(tin) ) tol = tin

        !write(*,*)
        !write(*,*)
        !write(*,*)
        call findBounds()

        !write(*,*) 'Bounds :',x1,x2,abs(derivF(x2))

        dx = (x2-x1)/float(npts-1)
        open(33,file=trim(label),status='unknown',action='write')
        
        x = x1-dx
        do i=1,npts
            x = x + dx
            write(33,'(3(ES15.7))') x, this%eval(x),derivF(x)
        end do

        close(33)

        err = 0
    contains
        
        subroutine findBounds()
        implicit none
            logical :: run    
            integer :: c
            real :: xm, fm, fp
        
            fp = 0.0
            x1 = this%t
            do while( abs(fp) < 5.0*tol )
                x1 = x1+step
                fp = derivF(x1)
            end do
            x2 = x1+step
    
            c = 0
            run = .TRUE.
            do while( run .and. (c .le. maxC) )
                fp = derivF(x2)
                if ( abs(fp) < tol ) then
                    xm = 0.5*(x1+x2)
                    fm = derivF(fm)
                    if ( abs(fm) < tol ) then
                        x2 = xm
                    else
                        run = .FALSE.
                    end if
                else
                    x1 = x2
                    x2 = x1 + step
                end if
                c = c+1
            end do
            x1 = 2.0*this%t - x2

        end subroutine findBounds
        
        function derivF (lX) result(ldf)
        implicit none
            real, intent(in) :: lX
            real :: ldf
            ldf = ( this%eval(lX+delx) - this%eval(lX) )/delx
        end function derivF
    end function printConstraint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! f Functions and their initializers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generalized exponential function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_gen_exp(this,p) 
    implicit none
        class(Constraint) :: this
        real, dimension(:), intent(in) :: p

        real :: expo,v,s
 
        if( size(p) .ne. 3 ) then
            write(*,*) "inverse polynomial requires 3 parameters"
            return
        end if

        this%nParams = 2
        allocate( this%params(this%nParams) )
        allocate( this%pLabels(this%nParams) )

        this%pLabels(1) = "b"
        this%pLabels(2) = "n"
       
        v = p(1)
        s = p(2)
        this%t = p(3)
        expo = -2.0*v*s/log(2.0)
        this%params(1) = expo ! b
        this%params(2) = log(2.0)/( v**expo ) ! n

    end subroutine init_gen_exp

    function generalized_exp (this,x) result(v)
    implicit none
        class(Constraint), intent(in) :: this
        real, intent(in) :: x
        real :: v 
        real :: n,b

        b = this%params(1)
        n = this%params(2)
        v = exp(-n*x**b)
    end function generalized_exp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse polynomial function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine init_inv_poly(this,p)
    implicit none
        class(Constraint) :: this
        real, dimension(:), intent(in) :: p

        real :: expo,v,s

        if( size(p) .ne. 3 ) then
            write(*,*) "inverse polynomial requires 3 parameters"
            return
        end if
        
        this%nParams = 2
        allocate( this%params(this%nParams) )
        allocate( this%pLabels(this%nParams) )

        this%pLabels(1) = "b"
        this%pLabels(2) = "n"

        v = p(1)
        s = p(2)
        this%t = p(3)
        expo = 4.0*v*s
        this%params(1) = -expo ! b
        this%params(2) = v**expo ! n
    end subroutine init_inv_poly

    function inverse_poly (this,x) result(v)
    implicit none
        class(Constraint), intent(in) :: this
        real, intent(in) :: x
        real :: v 

        real :: n,b

        b = this%params(1)
        n = this%params(2)
        v = 1.0/( n*x**b + 1 )
    end function inverse_poly


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the constraint value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function eval(this,x) result(v)
    implicit none
        class(Constraint), intent(in) :: this
        real, intent(in) :: x
        real :: v
        real :: g

        g = this%func(abs(x-this%t))   
        select case( trim(this%constraint_type) )
            case("equal")
                v = this%eval_equal(g)
            case DEFAULT
                v = 0.0
        end select

    end function eval

    function eval_equal(this,x) result(v)
    implicit none
        class(Constraint), intent(in) :: this
        real, intent(in) :: x
        real :: v
        v = 1.0 - x
    end function eval_equal


end module class_constraints
