!ニュートン法を用いてx=(a)^(1/k)の近似解を求める
program newton_method
  implicit none
  double precision a, k, x1, x2
  double precision,parameter :: threshold = 1.0d-6
  integer i, max
  data i,max/1,100/
  write (*, *) "SOLVE(a >= 0 .and. k >= 2) : x = a^(1/k)"
  write (*, '(a)', advance='no') "Enter a :"
  read (*, *) a
  if (a < 0) stop "a must be positive"
  write (*, '(a)', advance='no') "Enter k :"
  read (*, *) k
  if (k < 2) stop "k must be more than 2"
  
  !f(x1)*f''(x1)>0となるようにxを選べば収束速度が速いと言われている
  !この場合はx1^k>aが成り立つように選べばよい
  x1 = a
  if (a == 0) then
     x2 = 0
  else
     do i = 1, max
        x2 = x1 - (x1**k - a)/(k * x1**(k-1))
        if (abs(x1 - x2) < threshold) exit
        x1 = x2
     end do
  end if
  write (*, *) "Approximate answer is", x2
  write (*, *) "Trial times", i
end program newton_method