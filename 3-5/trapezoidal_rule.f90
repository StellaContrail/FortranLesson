program trapezoidal_rule
  implicit none
  double precision :: dx, x, y, s = 0.0d0
  double precision,parameter :: pi = acos(-1.0d0)
  integer :: i, n = 0
  write (*, '(a)', advance='no') "Precision (higher is better) :"
  read (*, *) n
  if (n < 1) stop "precision n must be greater than zero"
  dx = pi / dble(n)
  do i = 0, n
     x = dx * dble(i)
     y = sin(x) ! 被積分関数の定義
     if (i == 0 .or. i == n) then
        s = s + 0.5d0 * y
     else
        s = s + y
     end if
  end do
  s = s * dx
  write (*, *) "S =", s
end program trapezoidal_rule