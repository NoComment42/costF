program test
use class_constraints
implicit none

    type(Constraint) :: c1, c2
    integer :: i
    integer :: nParams
    real, dimension(3) :: inp
    real, allocatable, dimension(:) :: params
    Character(len=120) :: fname
    Character(len=10), allocatable, dimension(:) :: pLabels
    inp(1) = 0.5
    inp(2) = -2.0
    inp(3) = 1.0
    c1 = Constraint("equal","generalized_exp",inp)
    write(*,*) "constraint label : ", trim(c1%getConstraintLabel())
    write(*,*) "constraint type  : ", trim(c1%getConstraintType())
    write(*,*) "function used    : ", trim(c1%getFunctionType())
    write(*,*)

    nParams = c1%getNParams()
    allocate( params(nParams) )
    call c1%getParameters(params,pLabels)
    write(*,*) "number of parameters for the function : ", nParams
    write(*,*) "Parameters are"
    do i=1,nParams
        write(*,*) " ", pLabels(i)," : ", params(i)
    end do

    write(*,*)
    write(*,*) "eval : ", c1%eval(1.4)
    write(*,*)
    write(*,*)

    c2 = Constraint("equal","inverse_poly",inp,"Second Constraint")
    write(*,*) "constraint label : ", trim(c2%getConstraintLabel())
    write(*,*) "constraint type  : ", trim(c2%getConstraintType())
    write(*,*) "function used    : ", trim(c2%getFunctionType())
    write(*,*)

    nParams = c2%getNParams()
    if( allocated(params) ) deallocate(params)
    allocate( params(nParams) )
    call c2%getParameters(params,pLabels)
    write(*,*) "number of parameters for the function : ", nParams
    write(*,*) "Parameters are"
    do i=1,nParams
        write(*,*) " ", pLabels(i)," : ", params(i)
    end do

    write(*,*)
    write(*,*) "eval : ", c2%eval(1.4)

    fname = "c1.dat"
    i = c1%printConstraint(fname,500)
    fname = "c2.dat"
    i = c2%printConstraint(fname,500)

end program test
