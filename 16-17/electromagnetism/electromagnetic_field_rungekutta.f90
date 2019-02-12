module initialize_values
    double precision,parameter :: observeTime = 6.5
    double precision,parameter :: startXPosition = 10.0d0
    double precision,parameter :: startYPosition = 0.0d0
    double precision,parameter :: startXVelocity = -2.5d0
    double precision,parameter :: startYVelocity = 0.0d0
    double precision,parameter :: k = 1.0d0 !8.9876d9
    double precision,parameter :: QS = 3.0d0 ! Squared Electrical Charge
    double precision,parameter :: m = 1.0d0
end module initialize_values

module expressions
    use initialize_values
    implicit none
contains
    !Total Energy
    double precision function energy(v_x, v_y, r)
        double precision,intent(in) :: v_x, v_y, r
        energy = 0.5d0*m*(v_x**2+v_y**2)- k*QS / r
    end function
    !Linear momentum that charged particle has
    function linear_momentum(v_x, v_y)
        double precision,intent(in) :: v_x, v_y
        double precision linear_momentum(2)
        linear_momentum(1) = m * v_x
        linear_momentum(2) = m * v_y
    end function

    !x, y => old position
    !c => variable of direction in which to solve, x or y
    double precision function f(x, y, c)
        double precision,intent(in) :: x, y, c
        double precision r
        r = sqrt(x**2 + y**2)
        f = k*QS/ (m*r**3)
        f = f * c
    end function
    !Runge-Kutta method
    subroutine rk4(t, x, y, v_x, v_y, dt)
        double precision,intent(inout) :: t, x, y, v_x, v_y
        double precision,intent(in) :: dt
        double precision kx(4), kvx(4), ky(4), kvy(4)

        kx(1) = v_x
        ky(1) = v_y
        kvx(1) = f(x, y, x)
        kvy(1) = f(x, y, y)

        kx(2) = v_x + 0.5d0 * dt * kvx(1)
        ky(2) = v_y + 0.5d0 * dt * kvy(1)
        kvx(2) = f(x + 0.5d0 * dt * kx(1), y + 0.5d0 * dt * ky(1), x + 0.5d0 * dt * kx(1))
        kvy(2) = f(x + 0.5d0 * dt * kx(1), y + 0.5d0 * dt * ky(1), y + 0.5d0 * dt * ky(1))

        kx(3) = v_x + 0.5d0 * dt * kvx(2)
        ky(3) = v_y + 0.5d0 * dt * kvy(2)
        kvx(3) = f(x + 0.5d0 * dt * kx(2), y + 0.5d0 * dt * ky(2), x + 0.5d0 * dt * kx(2))
        kvy(3) = f(x + 0.5d0 * dt * kx(2), y + 0.5d0 * dt * ky(2), y + 0.5d0 * dt * ky(2))

        kx(4) = v_x + dt * kvx(3)
        ky(4) = v_y + dt * kvy(3)
        kvx(4) = f(x + dt * kx(3), y + dt * ky(3), x + dt * kx(3))
        kvy(4) = f(x + dt * kx(3), y + dt * ky(3), y + dt * ky(3))

        x = x + (kx(1)+2*kx(2)+2*kx(3)+kx(4)) * dt / 6.0d0
        y = y + (ky(1)+2*ky(2)+2*ky(3)+ky(4)) * dt / 6.0d0
        v_x = v_x + (kvx(1)+2*kvx(2)+2*kvx(3)+kvx(4)) * dt / 6.0d0
        v_y = v_y + (kvy(1)+2*kvy(2)+2*kvy(3)+kvy(4)) * dt / 6.0d0

        t = t + dt
    end subroutine
end module expressions

program electromagnetic_field_rungekutta
    use initialize_values
    use expressions
    implicit none
    double precision x, y, v_x, v_y, t, dt
    integer :: n = 1000, i, j
    integer,parameter :: fo = 10

    write (*, '(A, I0, /)') "Precision (Number of calculation loop) is now set to ", n
    write (*, *) "Observe Time: ", observeTime, " seconds"
    write (*, '(A, F0.5, A, F0.5, A)') "Start Position(m): (", startXPosition, ",", startYPosition, ")"
    write (*, '(A, F0.5, A, F0.5, A)') "Start Velocity(m/s): (", startXVelocity, ",", startYVelocity, ")"
    write (*, '(A, F0.5)') "Q^2(C^2): ", QS
    write (*, '(A, F0.5, /)') "m(kg): ", m

    open(fo, file='output.txt')
    dt = observeTime / dble(n)
    do j = -20, 20
        t = 0.0d0
        x = startXPosition
        y = startYPosition + dble(j) * 0.1d0
        v_x = startXVelocity
        v_y = startYVelocity
        write (* , *) "Start position (", x, y, ")"
        do i = 0, n
            call rk4(t, x, y, v_x, v_y, dt)
            write(fo, *) x, y, t
        end do
        write (fo, *)
    end do

    write (*, *) "Output File >> ", "output.txt"
    close(fo)
end program electromagnetic_field_rungekutta