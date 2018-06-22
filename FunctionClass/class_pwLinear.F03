module class_pwLinear
use class_pwFunc
implicit none
private

public :: pwLinear

type, extends(pwFunc) :: pwLinear
        real, private, allocatable :: a(:,:)
    contains
        procedure, private :: initLinear1
        procedure, private :: initLinear2
        procedure, public :: evalF
        procedure, public ::  integrate
        procedure, private :: findEle
end type pwLinear

interface pwLinear
    procedure  constructor 
end interface pwLinear

contains

    type(pwLinear) function constructor(pts,aTemp,ypts)
    implicit none
        real, intent(in) :: pts(:)
        real, optional, intent(in) :: aTemp(:,:)
        real, optional, intent(in) :: ypts(:)
        if( present(aTemp) ) then
            call constructor%initLinear1(pts,aTemp)
        else if( present(ypts) ) then
            call constructor%initLinear2(pts,ypts)
        else
            write(*,*) "ERROR: failed to provide enough info to initialize &
                        &function pwLinear."
        end if
    end function constructor

    subroutine initLinear1(this, pts, aTemp)
    implicit none
        class(pwLinear) :: this
        real, intent(in) :: pts(:)
        real, intent(in) :: aTemp(:,:)

        if( size(pts)-1 /= size(aTemp,1) ) then
            write(*,*) 'ERROR: In Linear constructor number of knots does &
                         &not match the number of elements'
            return
        end if

        if( size(aTemp,2) .ne. 2 ) then
            write(*,*) 'ERROR: In Linear constructor number of paramters is &
                        & not equal to 2, which is required for linear elements'
        end if

        this%numPts = size(pts)

        allocate(this%pts(this%numPts))
        this%pts = pts

        allocate(this%a(this%numPts-1,2))
        this%a = aTemp

    end subroutine

    subroutine initLinear2(this, pts, ypts)
    implicit none
        class(pwLinear) :: this
        real, intent(in) :: pts(:), ypts(:)

        integer :: e

        this%numPts = size(pts)
        allocate( this%pts(this%numPts) )
        this%pts = pts

        allocate( this%a(this%numPts-1,2) )

        ! must calculate the linear elemets
        do e=1,this%numPts-1
            this%a(e,1) = (ypts(e+1) - ypts(e))/(pts(e+1) - pts(e))
            this%a(e,2) = ypts(e) - this%a(e,1)*pts(e)
        end do


    end subroutine initLinear2


    real function evalF(this, x)
    implicit none
        class(pwLinear), intent(in) :: this
        real, intent(in) :: x
        integer :: j

        j = findEle(this,x)

        if( j /= -1 ) then
            evalF = this%a(j,1)*x + this%a(j,2)
            !write(*,*) x, evalF, this%a(j,1), this%a(j,2)
        else
            write(*,*) 'ERROR: x not in range. '
            evalF = 0
        end if

    end function 

    real function integrate(this, tx1, tx2)
    implicit none
        class(pwLinear), intent(in):: this
        real, intent(in) :: tx1, tx2
        integer :: e1, e2, e
        integer :: p
        integer :: step
        real :: x1, x2

        e1 = findEle(this,tx1)
        e2 = findEle(this,tx2)

        step = sign(1,e2-e1)

        p = 1 ! need a default value and I think this is safe
        x1 = tx1
        if( e1 == e2 ) then
            x2 = tx2 
        else if( step < 0 ) then
            p = e1
            x2 = this%pts(p)
        else if( step > 0 ) then
            p = e1+1
            x2 = this%pts(p)
        end if
        
!        write(*,*) 'x1, x2, step: ',x1,x2,step

        integrate = 0.0
        do e = e1, e2, step
            integrate = integrate + integral()
!            write(*,*) 'sum ', integrate, integral()
            p = p + step
            x1 = x2
            if( e+step == e2 ) then
                x2 = tx2
            else if( e+step > e2 ) then
                continue
            else
                x2 = this%pts(p)
            end if
        end do
    
    contains
        real function integral()
        implicit none
            integral = 0.5*this%a(e,1)*(x2**2 - x1**2) &
                         + this%a(e,2)*(x2-x1)
        end function integral
    end function


    integer function findEle(this,x)
    implicit none
        class(pwLinear) :: this
        real, intent(in) :: x
        integer :: i
        real :: t1, t2

        findEle = -1
        do i=1,this%numPts-1
            t1 = this%pts(i) - x
            t2 = this%pts(i+1) - x
            if( sign(1.0,t1) /= sign(1.0,t2) .or. &
                t1 == 0.0 .or. t2 == 0.0 ) then
                findEle = i
            end if
        end do
    end function

end module class_pwLinear
