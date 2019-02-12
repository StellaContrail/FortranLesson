program sum_separately
    implicit none
    integer :: i, n, sum = 0
    
    write (*, '(a)', advance='no') "Count up from zero to :"
    read (*, *) n

    do i = 0, n, 2
        sum = sum + i
    end do
    do i = 1, n, 2
        sum = sum + i
    end do

    write (*, *) sum
end program sum_separately
