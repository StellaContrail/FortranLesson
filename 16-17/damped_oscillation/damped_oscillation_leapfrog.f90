module initialize_values
    implicit none
    double precision,parameter :: observeTime = 15.0d0
    double precision,parameter :: startHeight = 1.0d0
    double precision,parameter :: k = 1.0d0
    double precision,parameter :: m = 1.0d0
    double precision,parameter :: b = 2.0d0
end module

module expressions
    use initialize_values
    implicit none
contains
    subroutine euler(t, x, v, dt)
        double precision t, x, v, dt
        x = v * dt + x
        v = v - (k/m) * x * dt - (b/m) * v * dt
        t = t + dt
    end subroutine
end module

program damped_oscillation_leapfrog
    use initialize_values
    use expressions
    implicit none
    double precision x, v, t, dt
    integer :: n = 1000, i
    integer,parameter :: fo = 10
    dt = observeTime / dble(n)
    open(fo, file='output_lf.txt')
    t = 0.0d0
    x = startHeight
    v = -2.0d0
    v = v - ((k/m)*x+(b/m)*v)*dt*0.5d0
    do i = 1, n-1
        call euler(t, x, v, dt)
        write(fo, *) t, x, v
    end do
    v = v - ((k/m)*x+(b/m)*v)*dt*0.5d0
    call euler(t, x, v, dt)
    write(fo, *) t, x, v

    close(fo)
end program damped_oscillation_leapfrog