module initialize_values
    double precision,parameter :: observeTime = 6.2831853071795862
    double precision,parameter :: startXPosition = 1.0d0
    double precision,parameter :: startYPosition = 0.0d0
    double precision,parameter :: startXVelocity = 0.0d0
    double precision,parameter :: startYVelocity = 1.0d0
    double precision,parameter :: GM = 1.0d0
    double precision,parameter :: m  = 1.0d0
end module initialize_values

module expressions
    use initialize_values
    implicit none
contains
    double precision function f(r, dt)
        double precision,intent(in) :: r, dt
        f = (GM * m * dt) / r**3
    end function
    subroutine lp(x, y, v_x, v_y, dt)
        double precision,intent(inout) :: x, y, v_x, v_y
        double precision,intent(in) :: dt
        double precision r, a
        r = sqrt(x**2 + y**2)
        a = f(r, dt)
        v_x = v_x - x * a * 0.5d0
        v_y = v_y - y * a * 0.5d0
    end subroutine
    subroutine euler(t, x, y, v_x, v_y, dt)
        double precision,intent(inout) :: t, x, y, v_x, v_y
        double precision,intent(in) :: dt
        double precision r, a
        r = sqrt(x**2 + y**2)
        a = f(r, dt)

        x = x + v_x * dt
        y = y + v_y * dt
    
        v_x = v_x - x * a
        v_y = v_y - y * a
    
        t = t + dt
    end subroutine
end module expressions

program gravitational_field_leapfrog
    use initialize_values
    use expressions
    implicit none
    double precision x, y, v_x, v_y, t, dt
    integer :: n = 2**15, i
    integer,parameter :: fo = 10
    dt = observeTime / dble(n)
    open(fo, file='output_leapfrog.txt')
    t = 0.0d0
    x = startXPosition
    y = startYPosition
    v_x = startXVelocity
    v_y = startYVelocity

    call lp(x, y, v_x, v_y, dt)
    do i = 1, n-1
        call euler(t, x, y, v_x, v_y, dt)
        write(fo, *) x, y, t
    end do
    call euler(t, x, y, v_x, v_y, dt)
    write(fo, *) x, y, t

    close(fo)
end program gravitational_field_leapfrog