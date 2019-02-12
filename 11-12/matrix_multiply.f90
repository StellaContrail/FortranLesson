!Exercise #11
program matrix_multiply
    implicit none
    integer n, i, j
    double precision,allocatable :: a(:, :), b(:, :), c(:, :)
    double precision sum
    integer matrix_size(2)

    call load2d(10, '3d_matrix_data.txt', a, b)

    matrix_size = shape(a)
    n = matrix_size(1)
    allocate(c(n, n))
    
    write (*, '(A,I0,A,I0,A)') "INPUT [", n, "x", n, "]"
    write (*, *) "A = "
    do i = 1, n
        write (*, *) a(i, :)
    end do
    write (*, *) "B = "
    do i = 1, n
        write (*, *) b(i, :)
    end do
    write (*, *) "RESULT"
    c = matmul(a, b)
    write (*, *) "A x B = "
    do i = 1, n
        write (*, *) c(i, :)
    end do

    do j = 1, n
        do i = 1, n
            sum = sum + a(i, j) * b(i, j)
        end do
    end do
    write (*, *) "A * B = "
    write (*, *) sum
    deallocate(a, b, c)
contains
    subroutine load2d(number, name, a, b)
        integer,intent(in) :: number
        character(*) name
        double precision,allocatable :: a(:,:), b(:,:)
        integer matrix_size, i
        open(number, file=name)
        read (number, *) matrix_size
        allocate(a(matrix_size,matrix_size), b(matrix_size,matrix_size))
        read (number, *)
        do i = 1, matrix_size
            read (number, *) a(i, :)
        end do
        read (number, *)
        do i = 1, matrix_size
            read (number, *) b(i, :)
        end do
    end subroutine
end program matrix_multiply