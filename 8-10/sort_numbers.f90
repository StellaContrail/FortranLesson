program sort_numbers
    implicit none
    integer i
    double precision,allocatable :: a(:)

    call loadList(10, 'num_list.txt', a)
    call sort(a)
    write (*, *) "RESULT"
    do i = 1, size(a)
        write (*, '(A,I5.5,A,f0.6)') "[", i, "] = ", a(i)
    end do
    call savelist(10, 'num_list_output.txt', a)
    write (*, *) "Output >> num_list_output.txt"
    deallocate(a)
contains
    subroutine loadlist(number, name, a)
        integer,intent(in) :: number
        character(*) name
        double precision,allocatable :: a(:)
        integer :: lines = 0, io = 0, i
        open(number, file=name)
        do
            read (number, *, iostat=io)
            if (io /= 0) exit
            lines = lines + 1
        end do
        allocate(a(lines))
        rewind(number)
        do i = 1, lines
            read (number, *) a(i)
        end do
        close(number)
    end subroutine
    subroutine savelist(number, name, a)
        integer,intent(in) :: number
        character(*) name
        double precision a(:)
        integer i
        open(number, file=name)
        do i = 1, size(a)
            write (number, *) a(i)
        end do
        close(number)
    end subroutine
    !選択ソート
    subroutine sort(a)
        double precision a(:)
        double precision min_value
        integer i, n, j, min_point
        n = size(a)
        do i = 1, n - 1
            min_point = i
            min_value = a(i)
            do j = i + 1, n
                if (a(j) < min_value) then
                    min_value = a(j)
                    min_point = j
                end if
            end do
            a(min_point) = a(i)
            a(i) = min_value
        end do
    end subroutine
end program sort_numbers