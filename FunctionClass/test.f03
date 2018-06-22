program test
use class_pwFunc
use class_pwLinear
implicit none

    type(pwLinear), target :: func1, func2

    real :: ans, x1, x2
    real, allocatable :: pts(:), ypts(:)
    real :: a(2,2)

    allocate(pts(3))

    pts(1) = 1.0
    pts(2) = 4.0
    pts(3) = 6.0

    a(1,1) = -1.0
    a(1,2) =  8.0

    a(2,1) =  2.0
    a(2,2) = -4.0

    func1 = pwLinear(pts,a)

    pts(1) = 1.0
    pts(2) = 3.0
    pts(3) = 6.0

    allocate(ypts(3))
    ypts(1) = 2.0
    ypts(2) = 6.0
    ypts(3) = 1.0

    func2 = pwLinear(pts,ypts=ypts)

    deallocate(pts)
    deallocate(ypts)    

    x1 = 2.1
    x2 = 4.5

    ans = func1%integrate(x1, x2)
    write(*,*) 'ans ', ans

    ans = func2%integrate(x1, x2)
    write(*,*) 'ans ', ans

    ans = func1%integrate(x1, x2) - func2%integrate(x1,x2)
    write(*,*) 'ans ', ans

end program test
