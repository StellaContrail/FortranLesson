!Exercise #10
program analytics
    implicit none
    integer,parameter :: fn = 10
    integer :: io = 0, lines = 0, i
    double precision :: average = 0.0d0, variance = 0.0d0, sqerr=0.0d0
    double precision,allocatable :: a(:)

    open(fn, file='statistical_data.txt')
    !ファイルの行数を調べる
    do
        read(fn, *, iostat=io)
        if (io /= 0) exit
        lines = lines + 1
    end do
    rewind(fn)
    allocate(a(lines))
    do i = 1, lines
        read (fn, *) a(i)
    end do
    
    !平均値を求める
    average = sum(a) / dble(lines)
    write (*, *) "Average = ", average
    !分散を求める
    do i = 1, lines
        sqerr = sqerr + (a(i) - average)**2
    end do
    variance = sqerr / dble(lines)
    write (*, *) "Variance = ", variance
    
    deallocate(a)
    close(fn)
end program analytics