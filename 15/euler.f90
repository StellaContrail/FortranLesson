module expressions
    implicit none
contains
    double precision function f(x,y)
        double precision,intent(in) :: x, y
        f = cos(x) * y
    end function
    subroutine next(x, y,dx)
        double precision,intent(in) :: dx, x
        double precision,intent(inout) :: y
        y = y + f(x, y) * dx
    end subroutine
end module
program euler
    use expressions
    implicit none
    double precision x, y, delta
    integer :: i, n = 100
    integer,parameter :: fn = 10
    open(fn, file='output_euler.txt')
    
    delta = 10d0 / dble(n)
    ! Initialize Condition
    x = 0d0
    y = 1d0
    write (fn, *) x, y
    do i = 1, n
        x = delta * dble(i)
        call next(x, y, delta)
        write (fn, *) x, y
    end do
end program euler