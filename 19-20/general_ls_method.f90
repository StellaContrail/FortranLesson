module extension
    implicit none
    private :: gauss_jordan
contains
    ! i th term of f(x)
    double precision function f(i, x)
        integer,intent(in) :: i
        double precision,intent(in) :: x
        f = x**i
    end function
    ! ∂f(x)/∂a_i = g(i, x)
    double precision function g(i, x)
        integer,intent(in) :: i
        double precision,intent(in) :: x
        g = x**i
    end function
    subroutine invert_matrix(a)
        double precision a(0:, 0:)
        double precision, allocatable :: b(:, :), a0(:, :)
        integer i, len_max, len_min
        len_min = lbound(a,1)
        len_max = ubound(a,1)
        allocate( a0(len_min:len_max,len_min:len_max), b(len_min:len_max, len_min:len_max) )
        a0 = a
        b = 0.0d0
        do i = len_min, len_max
            b(i, i) = 1.0d0
        end do
        do i = len_min, len_max
            a = a0
            call gauss_jordan(a, b(:,i), len_min, len_max)
        end do
        a(:, :) = b(:, :)
        deallocate(b)
    end subroutine

    subroutine gauss_jordan(a, b, len_min, len_max)
        integer, intent(in) :: len_min, len_max
        double precision,intent(inout) :: a(len_min:len_max,len_min:len_max), b(len_min:len_max)
        integer i, j, k
        double precision ar
        do k = len_min, len_max
            if(a(k, k) == 0.0d0) stop 'pivot = 0'
            ar = 1.0d0 / a(k, k)
            a(k, k) = 1.0d0
            do j = k + 1, len_max
                a(k, j) = ar * a(k, j)
            end do
            b(k) = ar * b(k)
            do i = len_min, len_max
                if(i /= k) then
                    do j = k + 1, len_max
                        a(i, j) = a(i, j) - a(i, k) * a(k, j)
                    end do
                    b(i) = b(i) - a(i,k) * b(k)
                    a(i, k) = 0.0d0
                end if
            end do
        end do
    end subroutine gauss_jordan
end module

program least_squares_method
    use extension
    implicit none
    integer,parameter :: fi = 10
    integer i, k, j
    integer,parameter :: n = 20, m = 2
    double precision x(1:n), y(1:n), yerr(1:n), d(0:m, 0:m), b(0:m), a(0:m), aerr(0:m)
    open(fi, file='data.txt')
    read (fi, *)

    x(:) = 0.0d0
    y(:) = 0.0d0
    yerr(:) = 0.0d0
    d(:, :) = 0.0d0
    b(:) = 0.0d0
    a(:) = 0.0d0
    aerr(:) = 0.0d0

    do i = 1, n
        read (fi, *) x(i), y(i), yerr(i)
    end do

    do j = 0, m
        do k = 0, m
            do i = 1, n
                d(k, j) = d(k, j) + f(j,x(i))*g(k,x(i)) / yerr(i)**2
            end do
        end do
    end do

    do k = 0, m
        do i = 1, n
            b(k) = b(k) + y(i)*g(k,x(i)) / yerr(i)**2
        end do
    end do
    
    call invert_matrix(d)
    
    do j = 0, m
        do k = 0, m
            a(j) = a(j) + d(j, k)*b(k)
        end do
    end do

    do k = 0, m
        aerr(k) = sqrt(d(k, k))
    end do

    do i = 0, m
        write (*, '(a, i0, a, f0.10)') "a(", i, ") = ", a(i)
    end do
    do i = 0, m
        write (*, '(a, i0, a, f0.10)') "aerr(", i, ") = ", aerr(i)
    end do
    close(fi)
end program least_squares_method