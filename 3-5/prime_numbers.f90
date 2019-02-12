program prime_numbers
    implicit none
    integer :: j = 0, n
    write (*, '(a)', advance='no') "Trial Numbers (0 for infinity):"
    read (*, *) n
    if (n == 0) then
        do
            j = j + 1
            if (isprime(j) .eqv. .true.) write (*, *) j
        end do
    else
        do j = 1, n
            if (isprime(j) .eqv. .true.) write (*, *) j
        end do
    end if
    
contains
    logical function isprime(n)
        integer,intent(in) :: n
        integer i
        isprime = .false.

        if (n <= 1) return
        if (mod(n,2) == 0) return
        !合成数xはp≦√xを満たす素因子pをもつ
        !sqrt(real(n)) always returns sqrt of n
        !since sqrt of n is not always integer, say, for instance, we define for statement from 3 to sqrt of 11 (3.3166)
        !then this for statement will increase i from 3 to 3
        do i = 3, int(sqrt(real(n))), 2
            if (mod(n,i) == 0) return
        enddo
        isprime = .true.
    end function isprime
end program prime_numbers