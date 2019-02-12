program sort_variables
    implicit none
    double precision :: a=3, b=5, c=1
    integer i
    write (*, *) "Original:"
    write (*, *) "a = ", a, "b = ", b, "c = ", c
    call sort(a, b, c)
    write (*, *) "Result"
    write (*, *) "a = ", a, "b = ", b, "c = ", c
contains
    subroutine swap(a, b)
        double precision :: a, b
        double precision :: temp
        temp = a
        a = b
        b = temp
    end subroutine
    subroutine sort(a, b, c)
        double precision,intent(in) :: a, b, c
        if (a < b) then
            !If a < b is true, there're two possibilities that a < b < c, a < c < b, or c < a < b exist. 
            !We cope with second and third variables later, so all we have to do is just judge the two leading variables and swap them (or do not).
            if (c < a) call swap(a, c)
            !Then we can understand that second and third variables are always greater than the a
        else
            !If b <= a is true, there're two possibilities that b < a < c, b < c < a, or c < b < a exist.
            !We again cope with second and third vars later, and we just leave them and judge the two leading vars and swap them or do not.
            if (b < c) then
                call swap(a, b)
            else
                call swap(a, c)
            end if
            !Then we can understand that "Second" and "Third" vars are always greater than the a
        end if
    !In the cases of former if statement as well as the latter one, the only remaining procedure is to judge second and third variables.
        if (c < b) call swap(b, c)
    end subroutine
end program sort_variables