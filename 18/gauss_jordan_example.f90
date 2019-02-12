module subprog
    implicit none
contains
    !Input -> a0, b, n
    !Output -> x
    !a0A=b ~= aA=x
    subroutine gauss_jordan(a0, x, b, n)
        integer, intent(in) :: n
        real(8), intent(in) :: a0(n,n), b(n)
        real(8),intent(out) :: x(n)
        integer i, j, k
        real(8) ar, a(n,n)
        
        !aA=x : A is variables matrix
        a(:,:) = a0(:,:)
        x(:) = b(:)
        
        !k行目を対象にする
        do k = 1, n
            !係数行列の対角要素が0であればアプリケーションを止める
            if(a(k,k) == 0.0d0) stop 'pivot = 0'
            !係数を1にするための倍数を調べる（割り算は掛け算よりもコストがかかるため先に計算をしておく）
            ar = 1.0d0 / a(k, k)
            !k行目のk番目の係数を1にする（アルゴリズムの高速化のため）
            a(k,k) = 1.0d0
            !係数行列のk行目をa(k, k)で割る
            do j = k+1, n
                a(k,j) = ar*a(k,j)
            end do
            !右辺も同様に。
            x(k) = ar*x(k)

            !k行以外の全ての行に対しての操作
            do i = 1,n
                if(i /= k) then
                    !k列目以前は計算したのでk+1から始める（k列目ではk行目以外は0になるはずなので飛ばす）
                    do j = k+1,n
                        !k列目の係数を0になるようにしたため、他の要素についても同じ数を掛けてやらなければいけない。
                        !i行目のk列目の係数を0になるようにするにはa(i,k)を掛ければよかった。
                        !そのため、i行目のj列目の係数に対してはa(k,j)*a(i,k)を掛けたものを引けば良い。
                        a(i,j) = a(i,j)-a(i,k)*a(k,j) 
                    end do
                    !右辺も同様にする。
                    x(i) = x(i)-a(i,k)*x(k)
                    !k列目ではk行目以外は0になるはずなので代入のみする。
                    a(i, k) = 0.0d0
                end if
            end do
        end do
    end subroutine gauss_jordan

    !行列aと行列bにランダムな数字を代入する
    subroutine set_random_mat(a, b, n)
        integer, intent(in) :: n
        real(8), intent(out) :: a(n, n), b(n)
        call random_seed !日時より乱数を計算し、それによりシードを決める
        call random_number(a)
        call random_number(b)
    end subroutine set_random_mat
end module subprog

program main
    use subprog
    implicit none
    real(8), allocatable :: a(:,:), b(:), x(:), r(:)
    integer n, i
    !次元を設定
    n = 3
    allocate (a(n,n), b(n), x(n), r(n))

    call set_random_mat(a,b,n)
    call gauss_jordan(a,x,b,n)

    !組み込み関数でAxを求め、当初の右辺と比べる
    r(:) = b(:)-matmul(a,x)

    write(*,*) 'A: '
    do i=1,n
        write(*,*) a(i,1:n)
    end do
    write(*,*) 'b: ', b(1:n)
    write(*,*) 'x: ', x(1:n)
    write(*,*) 'Ax: ', matmul(a,x)
    write(*,*) 'b-Ax: ', r(1:n)
    deallocate (a,b,x,r) 
end program main