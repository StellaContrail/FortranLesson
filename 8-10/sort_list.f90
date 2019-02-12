!ファイル内のリストをアイテムごとIndexに従って入れ替える。
!Indexは小数点でも動く
!それぞれのアイテムの最大バイト数は256 bytes
!ファイルの中のFORMATは以下の通り
!一行目：2 (Index以外のElementの数)
!二行目：
!三～N行目：[INDEX] [DATA1] [DATA2]
program sort_list
    implicit none
    type record
        double precision a
        character(len=256),allocatable :: body(:)
    end type record
    type(record),allocatable :: records(:)
    integer :: i, io = 0, lines = 0, elements_num
    integer,parameter :: fi = 10, fo = 11
    open(fi, file='list.txt')
    open(fo, file='list_output.txt')

    !要素数を読み取る
    read(fi, *) elements_num
    read(fi, *)
    !ファイルの行数を調べる
    do
        read(fi, *, iostat=io)
        if (io /= 0) exit
        lines = lines + 1
    end do
    rewind(fi)
    read(fi, *)
    read(fi, *)
    write (fo, '(I0)'), elements_num
    write (fo, *)

    allocate(records(lines))
    do i = 1, lines
        allocate(records(i)%body(elements_num))
    end do
    write (*, *) "Data Fetched from -> ", "list.txt"
    write (*, '(A,I0,A)') "SOURCE (including ", elements_num, " elements)"
    do i = 1, lines
        read(fi, *) records(i)%a, records(i)%body
        write (*, *) records(i)%a, records(i)%body
    end do
    call sort(records)

    write (*, *) "RESULT"
    do i = 1, lines
        write (*, *) records(i)%a, records(i)%body
        write (fo, *) records(i)%a, records(i)%body
    end do
    write (*, *) "Result output >> ", "list_output.txt"

    close(fi)
    close(fo)
    deallocate(records)
contains
    !選択ソート
    subroutine sort(records)
        type(record) records(:)
        double precision min_value
        integer i, n, j, min_point
        character(len=256),allocatable :: min_body(:)
        
        allocate(min_body(size(records(1)%body)))
        n = size(records)

        do i = 1, n - 1
            min_point = i
            min_value = records(i)%a
            do j = i + 1, n
                if (records(j)%a < min_value) then
                    min_value = records(j)%a
                    min_body = records(j)%body
                    min_point = j
                end if
            end do
            records(min_point)%a = records(i)%a
            records(i)%a = min_value
            records(min_point)%body = records(i)%body
            records(i)%body = min_body
        end do
        deallocate(min_body)
    end subroutine
end program sort_list