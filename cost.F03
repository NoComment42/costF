program test
use bisection
use class_constraints
implicit none

interface
    subroutine readFile(fname,xpts,ypts,xcol,ycol)
    implicit none
        character(len=120),intent(in) :: fname
        real, allocatable, intent(inout) :: xpts(:),ypts(:)
        integer, intent(inout) :: xcol, ycol
    end subroutine readFile

    subroutine printCurve( filename, tag, f, resolution )
    use class_pwFunc
    implicit none
        character(len=120), intent(in) :: filename, tag
        class(pwFunc), target, intent(in) :: f
        integer , intent(in) :: resolution
        class(pwFunc), pointer :: func
        integer :: i,j
        real :: x1, step
    endsubroutine
    

    subroutine difference( func1, func2, diff, error  )
    use bisection
    implicit none
        type(pwLinear) :: func1, func2
        real :: diff
        integer :: error
        integer :: numRoots
        integer :: i
        real :: sta, stp
        real, allocatable :: roots(:)
        real :: x1, x2, xm
    end subroutine difference


    subroutine setToZero(xpts,ypts,func1) 
    use class_pwLinear
    implicit none
        type(pwLinear), intent(in) :: func1
        real, allocatable, intent(inout) :: xpts(:),ypts(:)
    end subroutine

    
    subroutine calcDiff(func1,func2,diff)
    use bisection
    implicit none
        type(pwLinear), intent(in) :: func1, func2
        real, intent(inout) :: diff
        integer :: nl, error
        logical :: run = .TRUE.
    end subroutine 
end interface

!   parameters
    integer, parameter :: resolution = 4
    integer, parameter :: nTerms = 4
    real, parameter :: M = 0.1
    real, parameter :: gamma = 1.4


!   piecewise functions
    type(pwLinear) :: cpfunc, cpexp
    type(pwLinear) :: cffunc, cfexp
    type(pwLinear) :: zerof

!   constraints
    type(Constraint), allocatable, dimension(:) :: constraints
    integer :: err

!   general counters
    integer :: i

!   for file reading
    integer :: xcol, ycol
    character(len=120) :: fname, tag
    character(len=120) :: pathCpExp, pathCfExp, pathCpSim, pathCfSim
 
!   roots
    integer :: numRoots
    real :: sta, stp
    real, allocatable :: pts(:), ypts(:), expR(:), simR(:)

    real :: cost = 0.0
    real, dimension(nTerms) :: convert, diff



!   set conversion factors
    convert = 1.0

!   set diff to zero
    diff = 0.0

!   open output file
    open(22,file='diffCalc.out',status='unknown',action='write')

!   set pathnames to experiment and simulation results
#ifndef EXPDIR
#define EXPDIR "."
#endif
    write(22,*) 'Using EXPDIR: ' // EXPDIR
    pathCpExp = EXPDIR // '/noflow_cp.dat'
    pathCfExp = EXPDIR // '/noflow_cf.dat'

#ifndef SIMDIR
#define SIMDIR "."
#endif
    write(22,*) 'Using SIMDIR: ' // SIMDIR
    pathCpSim = SIMDIR // '/boundary.dat'
    pathCfSim = SIMDIR // '/boundary.dat'

#ifndef OUTDIR
#define OUTDIR "."
#endif
    write(22,*) 'Using OUTDIR: '//OUTDIR

!   Cp from experiment
    write(22,*) 'Read experiment Cp'
    xcol = 1
    ycol = 2
    call readFile(pathCpExp,pts,ypts,xcol,ycol)

    cpexp = pwLinear(pts,ypts=ypts)

    fname = OUTDIR // '/cp_exp.curve'
    tag = '# Cp Exp'
    call printCurve( fname, tag, cpexp, resolution )

!   Cf from experiment
    write(22,*) 'Read experiment Cf'
    xcol = 1
    ycol = 2
    call readFile(pathCfExp,pts,ypts,xcol,ycol)

    cfexp = pwLinear(pts,ypts=ypts)

    fname = OUTDIR // '/cf_exp.curve'
    tag = '# Cf Exp'
    call printCurve( fname, tag, cfexp, resolution )


!   Cp from simulation
    write(22,*) 'Read simulation Cp'
    xcol = 4
    ycol = 8
    call readFile(pathCpSim,pts,ypts,xcol,ycol)

    ! This one needs to be converted from P/Pinf to Cp
    ypts(:) = (2.0/gamma)*(ypts(:) - 1.0)/(M**2)

    cpfunc = pwLinear(pts,ypts=ypts)

    fname = OUTDIR // '/cp_sim.curve'
    tag = '# Cp Simulation'
    call printCurve( fname, tag, cpfunc, resolution )

!   Cf from experiment
    write(22,*) 'Read simulation Cf'
    xcol = 4
    ycol = 10
    call readFile(pathCfSim,pts,ypts,xcol,ycol)

    cffunc = pwLinear(pts,ypts=ypts)

    fname = OUTDIR // '/cf_sim.curve'
    tag = '# Cf Simulation'
    call printCurve( fname, tag, cffunc, resolution )

    write(22,*) "set a zero'd function"

!   Line of zeros
    call setToZero(pts,ypts,cffunc)
    zerof = pwLinear(pts,ypts=ypts)



!   start compiling cost function result

    write(22,*)
    write(22,*)
    write(22,*) 'Get Cp difference'
    call calcDiff(cpfunc,cpexp,diff(1))

    write(22,*)
    write(22,*) 'Get Cf difference'
    call calcDiff(cffunc,cfexp,diff(2))



    !di = resolution*size(pts,1)
    !step = (stp-sta)/float(di)
    !
    !open(12,file=OUTDIR//'/functions.dat',status='unknown',action='write')
    !open(13, file=OUTDIR//'/diff.curve', status='unknown',action='write')
    !write(13,*) '# Diff'
    !do i=1,di
    !    xm = sta + float(i)*step
    !    y1 = cpfunc%evalF(xm)
    !    y2 = cpexp%evalF(xm)
    !    yd = abs(f(xm))
    !    write(12,'(4(es15.7,X))') xm, y1, y2, yd
    !    write(13,*) xm, yd
    !end do
    !close(13)
    !close(12)
 
 
    write(22,*)
    write(22,*) 'Find separation/reattachment points from experiment'

    call init(cfexp,zerof)
    call setELevel(5)
    sta = cfexp%pts(1)
    stp = cfexp%pts(cfexp%numPts)
    call startBisection(sta,stp,numRoots)
    allocate(expR(numRoots))
    call getRoots(expR)
    write(22,*) 'number of points:', numRoots
    write(22,*) expR

    write(22,*)
    write(22,*) 'Find separation/reattachment points from simulation'

    call init(cffunc,zerof)
    call setELevel(5) ! set this somewhat high
    sta = cffunc%pts(1)
    stp = cffunc%pts(cffunc%numPts)
    call startBisection(sta,stp,numRoots)
    allocate(simR(numRoots))
    call getRoots(simR)
    write(22,*) 'number of points:', numRoots
    write(22,*) simR


    write(22,*)
    write(22,*)
    write(22,*) 'Setup constraint equations'

!   set constraints up
    allocate( constraints(nTerms) )
    constraints(1) = Constraint("equal","generalized_exp",[.5,-0.75,0.])
    constraints(2) = Constraint("equal","generalized_exp",[0.01,-3.,0.])
    constraints(3) = Constraint("equal","generalized_exp",[0.5,-1.,expR(1)])
    constraints(4) = Constraint("equal","generalized_exp",[0.5,-1.,expR(2)])
    err = constraints(1)%printConstraint(OUTDIR//"/cp_constraint.dat",100)
    err = constraints(2)%printConstraint(OUTDIR//"/cf_constraint.dat",500)
    err = constraints(3)%printConstraint(OUTDIR//"/sep_constraint.dat",100)
    err = constraints(4)%printConstraint(OUTDIR//"/attach_constraint.dat",100)

    numRoots = min(size(expR),size(simR))
    do i=1,numRoots
        diff(i+2) = abs(expR(i)-simR(i))
    end do

    write(22,*) 
    write(22,*)
    cost = 0.0
    do i=1,nTerms
        write(22,'("Diff: ",2(ES15.7,X))') diff(i), convert(i)*diff(i)
        cost = cost + diff(i)
    end do

    write(22,*)
    write(22,'("Cost: ",ES15.7)') cost

    diff(3) = simR(1)
    diff(4) = simR(2)

    write(22,*) 
    write(22,*) 
    write(22,*) 
    do i=1,nTerms
        write(22,*) constraints(i)%eval(diff(i))
    end do

    close(22)

end program test

subroutine readFile(fname,xpts,ypts,xcol,ycol)
implicit none

    character(len=120),intent(in) :: fname
    real, allocatable, intent(inout) :: xpts(:),ypts(:)
    integer, intent(inout) :: xcol, ycol

    integer :: i, err, npts, header
    real :: x
    real, allocatable :: tempx(:,:), tempy(:,:)
!    logical :: endHeader = .FALSE.
    
    open(12,file=trim(fname),status='old',action='read')
    if( allocated(xpts)) deallocate(xpts)
    if( allocated(ypts)) deallocate(ypts)
    
    header = 0
    npts = 0 
    do 
        read(12,'(f10.7)',iostat=err) x
        if( err == 0 ) then
            npts = npts+1
!            endHeader = .TRUE.
        else if( err < 0 ) then
            exit
        else if( err > 0 ) then
!            if(endHeader) exit
            header = header + 1
        end if
    end do

    !write(*,*) 'npts ', npts

    allocate(xpts(npts))
    allocate(ypts(npts))
    allocate(tempx(npts,xcol))
    allocate(tempy(npts,ycol-xcol))
    rewind(12)

    do i=1,header
        read(12,*)
    end do
    
    do i=1,npts
        read(12,*,iostat=err) tempx(i,:), tempy(i,:)
        if( err .ne. 0 ) then
            write(*,*) "ERROR: Didn't finish reading file "//trim(fname), err
            exit
        end if
    end do
    close(12)
    xpts(:) = tempx(:,xcol)
    ypts(:) = tempy(:,ycol-xcol)

end subroutine readFile

subroutine printCurve( filename, tag, f, resolution )
use class_pwFunc
implicit none

    character(len=120), intent(in) :: filename, tag
    class(pwFunc), target, intent(in) :: f
    integer , intent(in) :: resolution

    class(pwFunc), pointer :: func

    integer :: i,j
    real :: x1, step

    func => f

    open(16,file=trim(filename),status='unknown',action='write')
    write(16,*) trim(tag)
    do i=1,func%numpts-1
        step = (func%pts(i+1)-func%pts(i))/float(resolution-1)
        do j=0,resolution-1
            x1 = func%pts(i) + float(j)*step
            write(16,*) x1,func%evalF(x1)
        end do
    end do
    close(16)
end subroutine printCurve

subroutine difference( func1, func2, diff, error  )
use bisection
implicit none

    type(pwLinear) :: func1, func2
    real :: diff
    integer :: error


    integer :: numRoots
    integer :: i
    real :: sta, stp
    real, allocatable :: roots(:)
    real :: x1, x2, xm
    real :: td

    diff = 0.0
    error = 0

    sta = max(func1%pts(1),func2%pts(1))
    stp = min(func1%pts(func1%numPts),func2%pts(func2%numPts))
    write(22,*) sta, stp
    call startBisection(sta,stp,numRoots)
    write(22,*) 'number of Roots ', numRoots
 
    if( numRoots .ne. 0 ) then
        allocate(roots(numRoots))
        call getRoots(roots)
        write(22,*) roots
 
        x1 = sta
        x2 = roots(1)
        diff = 0.0
        do i=1,numRoots+1
!            write(22,*) 'x1 x2', x1, x2
            xm = 0.5*(x1+x2)
            if( func1%evalF(xm) > func2%evalF(xm) ) then
                td = func1%integrate(x1,x2) - func2%integrate(x1,x2)
            else
                td = func2%integrate(x1,x2) - func1%integrate(x1,x2)
            end if

            if( td .lt. 0.0 ) then
                error = 1
                diff = td
                return
            end if
            diff = diff + td

!            write(22,*) td, diff   
            x1 = x2
            if( i+1 > numRoots ) then
                x2 = stp
            else
                x2 = roots(i+1)
            end if
        end do
 
    else
        x1 = sta
        x2 = stp
        diff = 0.0
 
        xm = 0.5*(x1+x2)
        if( func1%evalF(xm) > func2%evalF(xm) ) then
            td = func1%integrate(x1,x2) - func2%integrate(x1,x2)
        else
            td = func2%integrate(x1,x2) - func1%integrate(x1,x2)
        end if

        if( td .lt. 0.0 ) then
            error = 1
            diff = td
            return
        end if

        diff = diff + td
!        write(22,*) td, diff

    end if
 

end subroutine difference


subroutine setToZero( xpts, ypts, func1 ) 
use class_pwLinear
implicit none
    type(pwLinear), intent(in) :: func1
    real, allocatable, intent(inout) :: xpts(:),ypts(:)

    if(allocated(xpts)) deallocate(xpts)
    if(allocated(ypts)) deallocate(ypts)
    
    allocate(xpts(2))
    allocate(ypts(2))
    
    xpts(1) = func1%pts(1)
    xpts(2) = func1%pts(func1%numPts)
    ypts = 0.0

end subroutine setToZero


subroutine calcDiff( func1, func2, diff )
use bisection
implicit none
    
    type(pwLinear), intent(in) :: func1, func2
    real, intent(inout) :: diff

    integer :: nl, error
    logical :: run

    call init(func1,func2)
    call getELevel(nl)

    ! need a loop to test for the number of 
    ! levels the bisection algorithm needs, if a difference
    ! calculation comes back negative, there was a missed root.
    run = .TRUE.
    do  while( run .eqv. .TRUE. )
        call difference(func1,func2,diff,error)
        write(22,'("Diff: ",ES15.7," Error: ",I1)') diff, error
        if( error .eq. 0 ) then
            run = .FALSE.
        else
            nl = nl + 2
            write(22,'("Changed failure levels to ",I3)') nl
            call setELevel(nl)
        end if
    end do


end subroutine calcDiff

