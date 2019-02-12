program mc_integration
    implicit none
    double precision :: x, y = 0.0d0, ans
    integer :: n = 2**20, i

    call random_seed
    do i = 1, n
        call random_number(x)
        y = y + x*(1.0d0-x)
    end do
    ans = 6.0d0*y/dble(n)
    write (*, *) "Result :", ans
    write (*, *) "Ideal  :", 1.0d0
    write (*, *) "Error  :", 1.0d0 - ans
end program