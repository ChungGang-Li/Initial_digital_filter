program Box

implicit none

integer :: n_cube, n_cellx, n_celly, n_cellz, n_cube_new
integer :: i,j,k,l,q,step,icube,ii,jj,kk,icount
real :: rdum,dx,dy,dz,x0,y0,z0,SML
real :: box_xp,box_xm, box_yp, box_ym, box_zp, box_zm
character(len=10) :: id_str

real, allocatable :: qcel(:,:,:,:,:)
real, allocatable :: qcel_rms(:,:,:,:,:)
real, allocatable :: x(:,:,:,:)
real, allocatable :: y(:,:,:,:)
real, allocatable :: z(:,:,:,:)
integer, allocatable :: mapping(:)
logical, allocatable :: check(:)

real, allocatable :: vertical(:,:)

integer :: p3dq(5) = (/1,2,3,4,5/)
 
	!Define the position of the box
	box_xm = 0.0
	box_xp = 0.0
	box_ym = 0.0
	box_yp = 0.0
	box_zm = 0.0
	box_zp = 0.0
  
  SML = 1.0e-8

	write(*,*) "Start reading mesh data..."
	open(unit=7, form='unformatted', file='mesh.g' ,convert='little_endian')	   
	read(7) n_cube
	read(7) (n_cellx, n_celly, n_cellz, i = 1, n_cube)
	allocate(x(n_cellx, n_celly, n_cellz, n_cube))
	allocate(y(n_cellx, n_celly, n_cellz, n_cube))
	allocate(z(n_cellx, n_celly, n_cellz, n_cube))
			
	x=0.0
	y=0.0
	z=0.0		
	do i = 1, n_cube
		read(7) &		
           ((( x(j,k,l,i),& 
			j=1,n_cellx),& 
			k=1,n_celly),& 
			l=1,n_cellz),&
            ((( y(j,k,l,i),&
			j=1,n_cellx),& 
			k=1,n_celly),& 
			l=1,n_cellz),&
            ((( z(j,k,l,i),& 
			j=1,n_cellx),& 
			k=1,n_celly),& 
			l=1,n_cellz)
	end do	
		
	close(7)
	write(*,*) "End reading mesh data"
	
	n_cube_new = 0
	allocate(mapping(n_cube))
	allocate(check(n_cube))
	
	check = .true. 
	
	write(*,*) "Start mapping..."
	write(*,*) 'Original cube :',n_cube
		
    
    
  
	
	write(*,*) "Start reading Turbulent_Intensity ..."
 
  open ( unit=7, form='unformatted', file='Turbulent_Intensity.q', convert='little_endian')
  read(7) n_cube
  read(7) (n_cellx, n_celly, n_cellz, i = 1, n_cube)
  
  allocate(qcel(5, n_cellx, n_celly, n_cellz, n_cube))
  
  qcel = 0.0
  
  do i = 1, n_cube
    read(7) rdum, rdum, rdum, rdum
    read(7) &
        ((((qcel(p3dq(q), j, k, l, i), &
        j = 1, n_cellx), &
        k = 1, n_celly), &
        l = 1, n_cellz), &
        q = 1, 5)
  end do	
  
  close(7)
  
  
  write(*,*) "====Reading ENDING===="
  
  
  
  write(*,*) "Identifying vertical points (Y_axis) ..."
  
	! If one of the cells of the cube is inside the box then extract
	do i = 1, n_cube
  
    dx = x(2,1,1,i)-x(1,1,1,i)
    dz = z(1,1,2,i)-z(1,1,1,i)
    
    j = 1
    k = 1
    l = 1

    if(check(i)) then
      if (x(j,k,l,i)+(n_cellx-1)*dx    .GE. SML .AND.&
          x(j,k,l,i)                   .LE. SML .AND.&
          z(j,k,l,i)+(n_cellz-1)*dz    .GE. SML .AND.&
          z(j,k,l,i)                   .LE. SML) then
      
        n_cube_new = n_cube_new + 1
        mapping(n_cube_new) = i
        check(i) = .false.
        
      end if
    end if	
    
	end do

  
  
  
  
  icount = 0
  allocate(vertical(n_cellx*n_cube_new,7))
  vertical = 0.0
  
  !! ==== vertical(Y_axis,Uavg,Vavg,Wavg,Urms,Vrms,Wrms) ==== !!
  
  do icube = 1, n_cube_new
  
    i = mapping(icube)
    
    x0 = x(1,1,1,i)
    y0 = y(1,1,1,i)
    z0 = z(1,1,1,i)
    
		do k = 1, n_celly	
    
       icount = icount+1
    
       vertical(icount,1) = y(1,k,1,i)
       
		end do
    
	end do
  
  
  call quicksort(vertical,1,n_cellx*n_cube_new)
  
  write(*,*) "Averaging in streamwise and spanwise directions ..."
        
  do icube = 1, n_cellx*n_cube_new
    
    y0 = vertical(icube,1)
    
    icount = 0
    do i = 1, n_cube
    
        dy = y(1,2,1,i)-y(1,1,1,i)
        
      ! if (y(1,1,1,i)+(n_celly-1)*dy    .GE. y0 .AND.&
          ! y(1,1,1,i)                   .LE. y0) then
      
        do j = 1, n_cellx
          do k = 1, n_celly
            do l = 1, n_cellz	
            
                if ( abs(y(j,k,l,i)-y0) .LE. SML ) then
                
                  icount = icount + 1
                  
				  vertical(icube,2) = vertical(icube,2)+qcel(1, j, k, l, i)
                  vertical(icube,3) = vertical(icube,3)+qcel(2, j, k, l, i)
				  vertical(icube,4) = vertical(icube,4)+qcel(3, j, k, l, i)
				  vertical(icube,5) = vertical(icube,5)+qcel(4, j, k, l, i)
				  vertical(icube,6) = vertical(icube,6)+qcel(5, j, k, l, i)
                  
                end if
                
            end do
          end do
        end do
        
      ! endif
      
    end do
  
    vertical(icube,2) = vertical(icube,2)/(icount*1.0)
    vertical(icube,3) = vertical(icube,3)/(icount*1.0)
    vertical(icube,4) = vertical(icube,4)/(icount*1.0)
    vertical(icube,5) = vertical(icube,5)/(icount*1.0)
    vertical(icube,6) = vertical(icube,6)/(icount*1.0)
  
  enddo
  
  
  write(*,*) "Output to the file"
  
  open(40, FILE='Statistics.dat', FORM='FORMATTED', status='replace', action='write')
  
	  write(40,*) " Y_axis, Uave, U'U', V'V', W'W', U'V' "
  
      do k = 1, n_cellx*n_cube_new	
      
       IF(K > 1 .and. vertical(k,1) .eq. vertical(k-1,1)) CYCLE
        
       write(40,*) vertical(k,1),vertical(k,2),vertical(k,3),vertical(k,4),&
                vertical(k,5),vertical(k,6)
       
      end do
      
  close(40)
  
  
  
    
end program box


recursive subroutine quicksort(a, first, last)

  implicit none
  real*4  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
  
end subroutine quicksort

