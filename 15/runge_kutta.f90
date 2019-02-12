module expressions
    implicit none
contains
    double precision function f(x,y)
        double precision,intent(in) :: x, y
        f = cos(x) * y
    end function
    subroutine next(y,dx,k1,k2,k3,k4)
        double precision,intent(in) :: dx, k1, k2, k3, k4
        double precision y
        y = y + dx * (k1 + 2d0*k2 + 2d0*k3 + k4) / 6d0 
    end subroutine
end module
program runge_kutta
    use expressions
    implicit none
    double precision x, y, delta, k1, k2, k3, k4
    integer :: i, n = 1000
    integer,parameter :: fn = 10
    open(fn, file='output_rk.txt')
    
    delta = 10d0 / dble(n)
    ! Initialize Condition
    x = 0d0
    y = 1d0
    write (fn, *) x, y 
    do i = 1, n
        x = delta * dble(i)
        k1 = f(x, y)
        k2 = f(x+0.5d0*delta, y+0.5d0*delta*k1)
        k3 = f(x+0.5d0*delta, y+0.5d0*delta*k2)
        k4 = f(x+delta, y+delta*k3)
        call next(y, delta, k1, k2, k3, k4)
        write (fn, *) x, y
    end do
end program runge_kutta