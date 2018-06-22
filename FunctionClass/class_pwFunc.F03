! Main purpose of this class is to set up the basic procedures and 
! variables that all piecewise function type objects should posess.

module class_pwFunc
implicit none

type, abstract :: pwFunc
    integer :: numPts
    real, allocatable :: pts(:)
contains
    procedure(evaluate), deferred :: evalF
    procedure(integrateF), deferred :: integrate
end type 

interface
    real function evaluate(this, x)
        import pwFunc
        class(pwFunc), intent(in) :: this
        real, intent(in) :: x
    end function evaluate

    real function integrateF(this,tx1,tx2)
        import pwFunc
        class(pwFunc), intent(in) :: this
        real, intent(in) :: tx1, tx2
    end function integrateF
end interface


end module class_pwFunc
