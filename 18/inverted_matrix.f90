module subprog
    implicit none
contains
    !Input -> a0, b, n
    !Output -> x
    !a0A=b ~= aA=x
    subroutine gauss_jordan(a, b, n)
        integer, intent(in) :: n
        double precision,intent(inout) :: a(n,n), b(n)
        integer i, j, k
        double precision ar

        !k行目を対象にする
        do k = 1, n
            !係数行列の対角要素が0であればアプリケーションを止める
            if(a(k, k) == 0.0d0) stop 'pivot = 0'
            !係数を1にするための倍数を調べる（割り算は掛け算よりもコストがかかるため先に計算をしておく）
            ar = 1.0d0 / a(k, k)
            !k行目のk番目の係数を1にする（アルゴリズムの高速化のため）
            a(k, k) = 1.0d0
            !係数行列のk行目をa(k, k)で割る
            do j = k + 1, n
                a(k, j) = ar * a(k, j)
            end do
            !右辺も同様に。
            b(k) = ar*b(k)

            !k行以外の全ての行に対しての操作
            do i = 1, n
                if(i /= k) then
                    !k列目以前は計算したのでk+1から始める（k列目ではk行目以外は0になるはずなので飛ばす）
                    do j = k + 1, n
                        !k列目の係数を0になるようにしたため、他の要素についても同じ数を掛けてやらなければいけない。
                        !i行目のk列目の係数を0になるようにするにはa(i,k)を掛ければよかった。
                        !そのため、i行目のj列目の係数に対してはa(k,j)*a(i,k)を掛けたものを引けば良い。
                        a(i, j) = a(i, j) - a(i, k) * a(k, j) 
                    end do
                    !右辺も同様にする。
                    b(i) = b(i) - a(i,k) * b(k)
                    !k列目ではk行目以外は0になるはずなので代入のみする。
                    a(i, k) = 0.0d0
                end if
            end do
        end do
    end subroutine gauss_jordan

    !行列aと行列bにランダムな数字を代入する
    subroutine set_random_mat(a, b, n)
        integer, intent(in) :: n
        double precision, intent(out) :: a(n, n), b(n,n)
        integer i, j
        double precision rand
        call random_seed !日時より乱数を計算し、それによりシードを決める
        do j = 1, n
            do i = 1, n
                call random_number(rand)
                a(i,j) = int(rand * 10)
            end do
        end do
        !call random_number(a)
        b(:,:) = 0
        do i = 1, n
            b(i, i) = 1
        end do
    end subroutine set_random_mat
end module subprog

program main
    use subprog
    implicit none
    double precision, allocatable :: a(:,:), a0(:,:), b(:,:), b0(:,:), r(:,:)
    integer n, i
    !次元を設定
    n = 3
    allocate (a(n,n), a0(n,n), b(n,n), b0(n,n), r(n,n))
    call set_random_mat(a, b, n)
    a0(:,:) = a(:,:)
    b0(:,:) = b(:,:)

    do i = 1, n
        a(:,:) = a0(:,:)
        call gauss_jordan(a, b(:,i), n)
    end do
    do i = 1, n
        !組み込み関数でAxを求め、当初の右辺と比べる
        r(:,i) = b(:,i) - matmul(a, b0(:,i))
    end do

    write(*,*) 'A: '
    do i=1,n
        write(*,*) a0(i, 1:n)
    end do
    write(*,*) 'A^-1: '
    do i=1,n
        write(*,*) b(i, 1:n)
    end do
    deallocate (a, a0, b, b0, r) 
end program main