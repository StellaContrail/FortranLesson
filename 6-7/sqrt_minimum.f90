program sqrt_minimum
    implicit none
    integer :: i, n = 2**15 !n is precision
    double precision :: a, x, y, dx, ymin = 1.0d10, xx
    integer,parameter :: fo = 10
    open(fo, file='output_sqrt_minimum_sqr_y.txt')

    write(*,'(a)', advance='no') "Sqrt of :"
    read(*,*) a
    if (a < 0) stop "Number must be greater than zero" 
    dx = a / dble(n)
    do i = 0, n
        x = dx * dble(i)
        y = (x**2 - a)**2
        write (fo, *) x, y
        if( y < ymin ) then
            ymin = y
            xx = x
        end if
    end do
    write (*, *) "Ans :", xx
    close(fo)

    write (*, *)
    write (*, *) "All answers with various precision(2^3 ~ 2^31)"
    call lists_of_sqrt(a)
contains
    subroutine lists_of_sqrt(a)
        double precision,intent(in) :: a
        integer i, n, j
        double precision :: x, y, dx, ymin = 1.0d10, xx
        do j = 3, 31
            ymin = 1.0d10
            n = 2**j
            dx = a / dble(n)
            do i = 0, n
                x = dx * dble(i)
                y = x**2 - a
                y = y**2
                if( y < ymin ) then
                    ymin = y
                    xx = x
                end if
            end do
            write(*,'(A,I2,A,F0.14)') "[Precision: 2^", j, "] ", xx
        end do
    end subroutine 
end program sqrt_minimum