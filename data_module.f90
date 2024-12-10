module data_module
    implicit none
contains
    subroutine Turbulent_Intensit_reading(filename, TI, nrows, ncols)
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=8), allocatable, intent(out) :: TI(:,:)
        integer, intent(in) :: nrows, ncols
        integer :: i, ios, unit_id
        real(kind=8), allocatable :: temp(:)
        character(len=256) :: line

        unit_id = 20

        ! 打開檔案
        open(newunit=unit_id, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "Error: Unable to open file!"
            stop
        endif

        ! 跳過表頭
        read(unit_id, "(A)", iostat=ios) line
        if (ios /= 0) then
            print *, "Error: Unable to read header!"
            close(unit_id)
            stop
        endif

        ! 分配陣列
        allocate(temp(ncols))
        allocate(TI(nrows, ncols))

        ! 讀取數據
        do i = 1, nrows
            read(unit_id, *, iostat=ios) temp
            if (ios /= 0) then
                print *, "Error reading data at row:", i
                print *, "Line content:", line
                exit
            endif
            TI(i, :) = temp
        end do

        ! 關閉檔案
        close(unit_id)
    end subroutine Turbulent_Intensit_reading
end module data_module
