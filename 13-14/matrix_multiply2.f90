!Exercise #14
program matrix_multiply
    implicit none
    integer n, i, j
    double precision,allocatable :: a(:, :), b(:, :), c(:, :)
    double precision inner_product
    integer matrix_size(2)

    call load2d(10, '3d_matrix_data.txt', a, b)

    matrix_size = shape(a)
    n = matrix_size(1)
    allocate(c(n, n))
    c = matmul(a, b)
    do j = 1, n
        do i = 1, n
            inner_product = inner_product + a(i, j) * b(i, j)
        end do
    end do

    call save2d(6, '', a, b, c, inner_product)
    call save2d(10, 'output_matrix_multiply.txt', a, b, c, inner_product)

    deallocate(a, b, c)
contains
    subroutine load2d(number, name, a, b)
        integer,intent(in) :: number
        character(*) name
        double precision,allocatable :: a(:,:), b(:,:)
        integer matrix_size, i, j
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
        close(number)
    end subroutine
    subroutine save2d(number, name, a, b, c, inner_product)
        integer,intent(in) :: number
        character(*) name
        double precision,intent(in) :: a(:,:), b(:,:), c(:,:)
        double precision inner_product
        integer i

        if (name /= "") open(number, file=name)

        write (number, '(A,I0,A,I0,A)') "INPUT [", n, "x", n, "]"
        write (number, *) "A = "
        do i = 1, n
            write (number, *) a(i, :)
        end do
        write (number, *) "B = "
        do i = 1, n
            write (number, *) b(i, :)
        end do
        write (number, *) "RESULT"
        write (number, *) "A x B = "
        do i = 1, n
            write (number, *) c(i, :)
        end do
        write (number, *) "A * B = "
        write (number, *) inner_product

        if (name /= "") close(number)
    end subroutine
end program matrix_multiply