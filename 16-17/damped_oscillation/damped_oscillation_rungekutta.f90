module initialize_values
    implicit none
    double precision,parameter :: observeTime = 15.0d0
    double precision,parameter :: startHeight = 1.0d0
    double precision,parameter :: k = 1.0d0
    double precision,parameter :: m = 1.0d0
    ! b = 2\sqrt(km)が境目で過減衰・臨界減衰・減衰振動が指定できる
    double precision,parameter :: b = 1.0d0
end module

module expressions
    use initialize_values
    implicit none
contains
    double precision function f(x, v)
        double precision, intent(in) :: x, v
        f = -(k/m) * x - (b/m) * v
    end function
    !Runge-Kutta method
    subroutine rk4(t, x, v, dt)
        double precision,intent(inout) :: t, x, v
        double precision,intent(in) :: dt
        double precision kx(4), kv(4)
        kx(1) = v
        kv(1) = f(x, v)

        kx(2) = v + 0.5d0 * dt * kv(1)
        kv(2) = f(x + 0.5d0 * dt * kx(1), v + 0.5d0 * dt * kv(1))

        kx(3) = v + 0.5d0 * dt * kv(2)
        kv(3) = f(x + 0.5d0 * dt * kx(2), v + 0.5d0 * dt * kv(2))

        kx(4) = v + dt * kv(3)
        kv(4) = f(x + dt * kx(3), v + dt * kv(3))

        x = x + (kx(1)+2d0*kx(2)+2d0*kx(3)+kx(4)) * dt / 6.0d0
        v = v + (kv(1)+2d0*kv(2)+2d0*kv(3)+kv(4)) * dt / 6.0d0

        t = t + dt
    end subroutine
end module

program damped_oscillation_rungekutta
    use initialize_values
    use expressions
    implicit none
    double precision x, v, t, dt
    integer :: n = 1000, i
    integer,parameter :: fo = 10
    open(fo, file='output_rk.txt')
    t = 0.0d0
    x = startHeight
    v = -2.0d0
    dt = observeTime / dble(n)
    do i = 0, n
        call rk4(t, x, v, dt)
        write(fo, *) t, x, v
    end do
    close(fo)
end program damped_oscillation_rungekutta