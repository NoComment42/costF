module bisection
use class_pwFunc
use class_pwLinear
    
    integer, private :: eLevel = 3
    real, private :: tol = 1e-5

    integer, private :: numR
    real, allocatable, private :: roots(:)
    logical, private :: OK = .FALSE.
    class(pwFunc), pointer, private :: func1, func2

    private :: bisect, add

contains

    subroutine init(f1,f2)
    implicit none
        class(pwFunc), target, intent(in) :: f1,f2
        eLevel = 3
        numR = 0
        if(allocated(roots)) deallocate(roots)
        allocate(roots(1))
        roots = 0.0

        func1 => f1
        func2 => f2
        OK = .TRUE.
    end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    real function f(x) 
    implicit none
        real, intent(in) :: x

        if(OK) then
            f = func1%evalF(x) - func2%evalF(x)
        else
            write(*,*) 'ERROR: need to call init before any other &
                        &procedure'
            f = 0
        end if

    end function f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    subroutine startBisection(x1,x2,nr)
    implicit none
        real, intent(in) :: x1, x2
        integer, intent(out) :: nr

        real :: y1, y2
        
        if(OK) then
            !write(*,*) 'All set, start bisection'
            y1 = f(x1)
            y2 = f(x2)
            call bisect(x1,y1,x2,y2,0)
        else
            write(*,*) 'ERROR: need to call init before any other &
                        &procedure'
        end if
        nr = numR
    end subroutine startBisection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    recursive subroutine bisect(x1,y1,x2,y2,err)
    implicit none
        real, intent(in) :: x1, x2
        real, intent(in) :: y1, y2
        integer, intent(in) :: err

        integer :: lerr
        real :: mpt
        real :: ym

        mpt = 0.5*(x1+x2)
        ym = f(mpt)

        !write(*,*) mpt, abs(ym)
        if( abs(ym) < tol ) then
            call add(mpt,ym)
            return
        end if

        if( ym*y1 < 0 ) then
            call bisect(x1,y1,mpt,ym,0)
        else
            lerr = err + 1
            if( lerr < eLevel ) then
                call bisect(x1,y1,mpt,ym,lerr)
            end if
        end if

        if( ym*y2 < 0 ) then
            call bisect(mpt,ym,x2,y2,0)
        else
            lerr = err + 1
            if( lerr < eLevel ) then
                call bisect(mpt,ym,x2,y2,lerr)
            end if
        end if

    end subroutine bisect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
    subroutine add(x,y)
    implicit none

        real, intent(in) :: x, y

        integer :: i, j
        real, allocatable :: temp(:)


        numR = numR+1

        if( numR > size(roots,1) ) then
            allocate(temp(numR))
            temp(1:numR-1) = roots
            deallocate(roots)
            allocate(roots(numR+10))
            roots = 0.0
            roots(1:numR) = temp
            deallocate(temp)
        end if

        j = numR
        do i=1,numR-1
            if( roots(i) > x ) then
                j = i
                exit
            end if
        end do


        if( j .ne. numR ) then
            allocate(temp(numR))
            temp(j:numR) = roots(j:numR)
            roots(j) = x
            roots(j+1:numR+1) = temp(j:numR)
            deallocate(temp)
        else
            roots(j) = x
        end if


    end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    integer function getNumRoots() result(r)
    implicit none
        if(OK) then
            r = numR
        else
            write(*,*) "ERROR: must call init first before calling &
                        &other procdeures"
            r = 0
        end if
    end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine getRoots(r) 
    implicit none
        real :: r(:)
        
        if(OK) then
            if( size(r,1) == numR ) then
                r = roots(1:numR)
            else
                write(*,*) 'ERROR: r must be of length numR'
            end if
        else
            write(*,*) "ERROR: must call init first before calling &
                        &other procdeures"
        end if

    end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine getELevel(nl)
    implicit none
        integer :: nl
        nl = eLevel
    end subroutine getELevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine setELevel(nL)
    implicit none
        integer :: nL
        eLevel = nL
    end subroutine


end module bisection

