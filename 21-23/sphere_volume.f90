module constant
    implicit none
    double precision,parameter :: pi = acos(-1.0d0)
end module

program sphere_volume
    use constant
    implicit none
    double precision volume_mc, volume_true
    integer :: n = 0
    write (*, '(a)', advance='no') "Dimension of the sphere :"
    read (*, *) n
    volume_mc = solve_volume(n, 2**25)
    volume_true = solve_true_volume(dble(n))
    write (*, *) "Solved by Monte Carlo :", volume_mc
    write (*, *) "Solved by algebra     :", volume_true
    write (*, *) "Error :", abs(volume_mc - volume_true)
contains
    double precision function solve_volume(n, count)
        integer,intent(in) :: n, count
        double precision,allocatable :: x(:)
        double precision :: sum = 0.0d0
        integer i, j
        integer :: hit = 0
        allocate(x(n))

        do j = 1, count
            call random_number(x)
            do i = 1, n
                sum = sum + x(i)**2
            end do
            if (sum < 1.0d0) then
                hit = hit + 1
            end if
            sum = 0.0d0
        end do
        deallocate(x)
        solve_volume = 2**n * (hit/dble(count))
    end function
    double precision function solve_true_volume(n)
        double precision,intent(in) :: n
        solve_true_volume = pi**(n/2.0d0)/gamma_func(0.5d0*n+1)
    end function
    double precision function gamma_func(n)
        double precision,intent(in) :: n
        integer i
        gamma_func = 1.0d0
        if (mod(n, 0.5d0) > 0.0d0) stop "Further calculation needed in gamma_func"
        if (n - int(n) < 1.0d-6) then
            do i = int(n)-1, 2, -1
                gamma_func = gamma_func * i
            end do
        else
            do i = int(n), 1, -1
                gamma_func = gamma_func * (i-0.5d0)
            end do
            gamma_func = gamma_func * sqrt(pi)
        end if
    end function
end program