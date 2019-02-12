!Exercise #8
program inner_product
    implicit none
    integer i
    double precision,allocatable :: u(:), v(:)
    double precision :: result = 0.0d0

    !Assign numbers into each element of arrays
    call load1d(10, '1d_matrix_data.txt', u, v)
    write (*, *) "Matrix Data Fetched from ", "1d_matrix_data.txt"
    write (*, *) "Array u :", u
    write (*, *) "Array v :", v
    do i = 1, size(u)
        result = result + u(i) * v(i)
    end do
    write (*, *) "<u,v> = ", result
    deallocate(u, v)
    
contains
    subroutine load1d(number, name, a, b)
        integer,intent(in) :: number
        character(*) name
        double precision,allocatable :: a(:), b(:)
        integer matrix_size
        open(number, file=name)
        read (number, *) matrix_size
        allocate(a(matrix_size), b(matrix_size))
        read (number, *)
        read (number, *) a
        read (number, *)
        read (number, *) b
        close(number)
    end subroutine
end program inner_product