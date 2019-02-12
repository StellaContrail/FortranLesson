program gaussrnd
    implicit none
    integer,parameter :: fo = 10, m = 1000, n = 500 !m-># of RandomNumbers, n->Samples
    double precision :: rnd(2), x(2), r, th, rand(2*m) = 0.0d0, interval, min, sum_x(2) = 0d0
    double precision,parameter :: pi2=8.0d0*atan(1.0d0), a=1.0d0
    integer,allocatable :: seed(:)
    integer seedsize, i
    open(fo, file='gauss_distribution.txt')

    ! 毎回違う乱数を生成
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do i = 1, seedsize
        call system_clock(count=seed(i))
    end do
    call random_seed(put=seed(:))

    do i = 1, m
        call random_number(rnd)
        r=sqrt(-log(rnd(1))/a) 
        th=pi2*rnd(2)
        x(1)=r*cos(th)
        x(2)=r*sin(th)
        sum_x(:) = sum_x(:) + x(:)**2d0
        rand(2*i - 1) = x(1)
        rand(2*i) = x(2)
    end do

    if (any(rand < 0)) then
        min = minval(rand, rand < 0)
    else
        min = minval(rand)
    end if
    !階級差を求める
    interval = (maxval(rand) - min) / dble(n)

    do i = 1, n
        !階級に入る度数を調べて出力
        write (fo, *) interval*(2d0*i-1d0)/2d0+min, count(interval*(i-1d0) <= (rand - min) .and. (rand - min) <= interval*i)
    end do

    write (*, *) "Variance (x_1 : x_2)"
    write (*, *) sum_x(:) / dble(m)
    write (*, *) "Error from theoretical value"
    write (*, *) 0.5d0 - sum_x(:) / dble(m)
    deallocate(seed)
    close(fo)
end program gaussrnd