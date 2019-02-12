!Values to initialze
module initialize_values
    double precision,parameter :: observeTime = 6.2831853071795862 !Time to observe the orbit
    double precision,parameter :: startXPosition = 1.0d0
    double precision,parameter :: startYPosition = 0.0d0
    double precision,parameter :: startXVelocity = 0.0d0
    double precision,parameter :: startYVelocity = 1.0d0
    double precision,parameter :: GM = 1.0d0 !G=Newtonian Constant of Gravitation, M=Mass of the center Planet
    double precision,parameter :: m = 1.0d0 ! m = mass of the circulating planet
end module initialize_values

module expressions
    use initialize_values
    implicit none
contains
    double precision function f(r, dt)
        double precision,intent(in) :: r, dt
        f = (GM * m * dt) / r**3
    end function
    subroutine euler(t, x, y, v_x, v_y, dt)
        double precision,intent(inout) :: t, x, y, v_x, v_y
        double precision,intent(in) :: dt
        double precision x0, y0, r, a
        x0 = x
        y0 = y
        r = sqrt(x0**2 + y0**2)
        a = f(r, dt)

        x = x0 + v_x * dt
        y = y0 + v_y * dt
    
        v_x = v_x - x0 * a
        v_y = v_y - y0 * a
    
        t = t + dt
    end subroutine
end module expressions

program gravitational_field_euler
    use initialize_values
    use expressions
    implicit none
    double precision x, y, v_x, v_y, t, dt
    integer :: n = 2**15, i
    integer,parameter :: fo = 10
    dt = observeTime / dble(n)
    open(fo, file='output_euler.txt')
    t = 0.0d0
    x = startXPosition
    y = startYPosition
    v_x = startXVelocity
    v_y = startYVelocity

    do i = 0, n
        call euler(t, x, y, v_x, v_y, dt)
        write(fo, *) x, y, t
    end do
    close(fo)
end program gravitational_field_euler