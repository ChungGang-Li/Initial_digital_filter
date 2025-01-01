
module subs
 implicit none
contains

subroutine prim(qprim, qcel, n_cube, n_cellx, n_celly, n_cellz)

real, intent(out) :: qprim(:,:,:,:,:)
real, intent(in) :: qcel(:,:,:,:,:)
real :: rho, u, v, w, E, T, temp
integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz
integer :: i, j, k, l

temp = (1.4-1.0)/287.0

do i = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
				rho = qcel(1,j, k, l, i)
				u = qcel(2,j, k, l, i)/qcel(1,j, k, l, i)
				v = qcel(3,j, k, l, i)/qcel(1,j, k, l, i)
				w = qcel(4,j, k, l, i)/qcel(1,j, k, l, i)
				E = qcel(5,j, k, l, i)/qcel(1,j, k, l, i)
				
				T = (E - 0.5*(u**2.0 + v**2.0 + w**2.0))*temp
				
				qprim(1, j, k, l, i) = rho
				qprim(2, j, k, l, i) = u
				qprim(3, j, k, l, i) = v
				qprim(4, j, k, l, i) = w
				qprim(5, j, k, l, i) = T
				
			end do
		end do	
	end do
end do


end subroutine prim

subroutine add(qsum, qprim, n_cube, n_cellx, n_celly, n_cellz)

real, intent(inout) :: qsum(:,:,:,:,:)
real, intent(in) :: qprim(:,:,:,:,:)
integer :: icube, j, k, l, q
integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz

do icube = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
				do q = 1, 5
					qsum(q, j, k, l, icube) = qsum(q, j, k, l, icube) + &
					qprim(q, j, k, l, icube)
				end do
			end do
		end do
	end do
end do	
end subroutine add	

subroutine avg( field, scalar, n_cube, n_cellx, n_celly, n_cellz)

!real, intent(out) :: This_field
real, intent(inout) :: field(:,:,:,:,:)
real, intent(in) :: scalar
integer :: icube, j, k, l, q
integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz

do icube = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
				do q = 1, 5
					field(q, j, k, l, icube) = field(q, j, k, l, icube)* scalar
					
				end do
				
			end do
		end do
	end do
end do	
end subroutine avg	

subroutine fluctuate(this_field, field, avg_field, n_cube, n_cellx, n_celly, n_cellz)
real, intent(out) :: this_field(:,:,:,:,:)
real, intent(in) :: field(:,:,:,:,:), avg_field(:,:,:,:,:)
integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz
integer :: icube, j, k, l, q
do icube = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
				do q = 1, 5
					this_field(q, j, k, l, icube) = field(q, j, k, l, icube) - &
					avg_field(q, j, k, l, icube)
					
				end do
			end do
		end do
	end do
end do	
end subroutine fluctuate


subroutine Turbulent_Intensity(this_field, avg_field, field, n_cube, n_cellx, n_celly, n_cellz)
	real, intent(out) :: this_field(:,:,:,:,:)
	real, intent(in) :: field(:,:,:,:,:), avg_field(:,:,:,:,:)
	integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz
	integer :: icube, j, k, l, q
	do icube = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
					! 1: u_avg; 2: u'u'; 3: v'v'; 4: w'w'; 4: u'v'
					this_field(1, j, k, l, icube) = avg_field(2, j, k, l, icube)
					
					this_field(2, j, k, l, icube) = field(2, j, k, l, icube) * &
					field(2, j, k, l, icube)
					this_field(3, j, k, l, icube) = field(3, j, k, l, icube) * &
					field(3, j, k, l, icube)
					this_field(4, j, k, l, icube) = field(4, j, k, l, icube) * &
					field(4, j, k, l, icube)
					this_field(5, j, k, l, icube) = field(2, j, k, l, icube) * &
					field(3, j, k, l, icube)
					
				end do
			end do
		end do
	end do
	
end subroutine Turbulent_Intensity


subroutine EHF(this_field, field, n_cube, n_cellx, n_celly, n_cellz)
	real, intent(out) :: this_field(:,:,:,:,:)
	real, intent(in) :: field(:,:,:,:,:)
	integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz
	integer :: icube, j, k, l, q
	do icube = 1, n_cube
	do l = 1, n_cellz
		do k = 1, n_celly
			do j = 1, n_cellx
					! 2: u'T'; 3: v'T'; 4: w'T'
					this_field(2, j, k, l, icube) = field(2, j, k, l, icube) * &
					field(5, j, k, l, icube)
					this_field(3, j, k, l, icube) = field(3, j, k, l, icube) * &
					field(5, j, k, l, icube)
					this_field(4, j, k, l, icube) = field(4, j, k, l, icube) * &
					field(5, j, k, l, icube)
					
				end do
			end do
		end do
	end do
	
	
end subroutine EHF

subroutine SS(this_field, field, avg_field, x, y, z, n_cube, n_cellx, n_celly, n_cellz)
	
	real, intent(out) :: this_field(:,:,:,:,:)
	real, intent(in) :: field(:,:,:,:,:), avg_field(:,:,:,:,:) !field: fluctuating values
	real, intent(in) :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
	integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz
	integer :: i,j,k,l,ii,jj,kk,ll
	real :: inv_dx, inv_dy, inv_dz
	real :: dupdx, dupdy, dupdz, dvpdx, dvpdy, dvpdz, dwpdx, dwpdy, dwpdz
	real :: temp, x1,x2, y1,y2, z1,z2
	
		do ii = 1, n_cube
		x1 = x(1,1,1,ii)
		x2 = x(2,1,1,ii)
		
		y1 = y(1,1,1,ii)
		y2 = y(1,2,1,ii)
		
		z1 = z(1,1,1,ii)
		z2 = z(1,1,2,ii)
		
		! inv_dx = 1d0 / abs(x(1,1,1,ii)-x(2,1,1,ii))
		! inv_dy = 1d0 / abs(y(1,1,1,ii)-y(1,2,1,ii))
		! inv_dz = 1d0 / abs(z(1,1,1,ii)-z(1,1,2,ii))
		
		inv_dx = 1d0 / abs(x1-x2)
		inv_dx = 1d0 / abs(y1-y2)
		inv_dx = 1d0 / abs(z1-z2)
		
		do ll = 1, n_cellz
			do kk = 1, n_celly
				do jj = 1, n_cellx			
					if(jj .eq. 1) then
						dupdx = (field(2, jj + 1, kk, ll, ii) - field(2, jj, kk, ll, ii)) * inv_dx
						dvpdx = (field(3, jj + 1, kk, ll, ii) - field(3, jj, kk, ll, ii)) * inv_dx
						dwpdx = (field(4, jj + 1, kk, ll, ii) - field(4, jj, kk, ll, ii)) * inv_dx
					elseif(jj .eq. n_cellx) then
						dupdx = (field(2, jj, kk, ll, ii) - field(2, jj - 1, kk, ll, ii)) * inv_dx
						dvpdx = (field(3, jj, kk, ll, ii) - field(3, jj - 1, kk, ll, ii)) * inv_dx
						dwpdx = (field(4, jj, kk, ll, ii) - field(4, jj - 1, kk, ll, ii)) * inv_dx
					else
						dupdx = (field(2, jj + 1, kk, ll, ii) - field(2, jj - 1, kk, ll, ii)) * 0.5d0 * inv_dx
						dvpdx = (field(3, jj + 1, kk, ll, ii) - field(3, jj - 1, kk, ll, ii)) * 0.5d0 * inv_dx
						dwpdx = (field(4, jj + 1, kk, ll, ii) - field(4, jj - 1, kk, ll, ii)) * 0.5d0 * inv_dx
					end if

					if(kk .eq. 1) then
						dupdy = (field(2, jj, kk + 1, ll, ii) - field(2, jj, kk, ll, ii) ) * inv_dy
						dvpdy = (field(3, jj, kk + 1, ll, ii) - field(3, jj, kk, ll, ii) ) * inv_dy
						dwpdy = (field(4, jj, kk + 1, ll, ii) - field(4, jj, kk, ll, ii) ) * inv_dy
					elseif(kk .eq. n_celly) then
						dupdy = (field(2, jj, kk, ll, ii)- field(2, jj, kk - 1, ll, ii) ) * inv_dy
						dvpdy = (field(3, jj, kk, ll, ii)- field(3, jj, kk - 1, ll, ii) ) * inv_dy
						dwpdy = (field(4, jj, kk, ll, ii)- field(4, jj, kk - 1, ll, ii) ) * inv_dy
					else
						dupdy = (field(2, jj, kk + 1, ll, ii) - field(2, jj, kk - 1, ll, ii)) * 0.5d0 * inv_dy
						dvpdy = (field(3, jj, kk + 1, ll, ii) - field(3, jj, kk - 1, ll, ii)) * 0.5d0 * inv_dy
						dwpdy = (field(4, jj, kk + 1, ll, ii) - field(4, jj, kk - 1, ll, ii)) * 0.5d0 * inv_dy
					end if

					if(ll .eq. 1) then
						dupdz = (field(2, jj, kk, ll + 1, ii) - field(2, jj, kk, ll, ii)) * inv_dz
						dvpdz = (field(3, jj, kk, ll + 1, ii) - field(3, jj, kk, ll, ii)) * inv_dz
						dwpdz = (field(4, jj, kk, ll + 1, ii) - field(4, jj, kk, ll, ii)) * inv_dz
					elseif(ll .eq. n_cellz) then
						dupdz = (field(2, jj, kk, ll, ii) - field(2, jj, kk, ll - 1, ii)) * inv_dz
						dvpdz = (field(3, jj, kk, ll, ii) - field(3, jj, kk, ll - 1, ii)) * inv_dz
						dwpdz = (field(4, jj, kk, ll, ii) - field(4, jj, kk, ll - 1, ii)) * inv_dz
					else
						dupdz = (field(2, jj, kk, ll + 1, ii) - field(2, jj, kk, ll - 1, ii)) * 0.5d0 * inv_dz
						dvpdz = (field(3, jj, kk, ll + 1, ii) - field(3, jj, kk, ll - 1, ii)) * 0.5d0 * inv_dz
						dwpdz = (field(4, jj, kk, ll + 1, ii) - field(4, jj, kk, ll - 1, ii)) * 0.5d0 * inv_dz
					end if
						
						
										
						this_field(1, jj, kk, ll, ii) = 2.0d0* dupdx**2 + dvpdx**2 + dwpdx**2 + dupdy**2 &
						+ 2.0d0*dvpdy**2 + dwpdy**2 + dupdz**2 + dvpdz**2 + 2.0d0*dwpdz**2 &
						+ 2.0d0*(dupdy*dvpdx + dupdz*dwpdx + dvpdz*dwpdy)
										
					end do
				end do					
			end do		
		end do	


	
	! if (allocated(up)) then
        ! deallocate(up)
    ! end if
	
	! if (allocated(vp)) then
        ! deallocate(vp)
    ! end if
	
	! if (allocated(wp)) then
        ! deallocate(wp)
    ! end if
	
	! deallocate(up)
	! deallocate(vp)
	! deallocate(wp)
	

end subroutine SS

subroutine Kol_epsilon_EHF(this_field, avg_field, mu, p0, Pr, rho0, R, Cv, n_cube, n_cellx, n_celly, n_cellz)

real, intent(inout) :: this_field(:,:,:,:,:)
real, intent(in) ::  avg_field(:,:,:,:,:) !,field(:,:,:,:,:)
real, intent(in) :: mu, p0, Pr, rho0, R, Cv
integer, intent(in) :: n_cube, n_cellx, n_celly, n_cellz

integer :: icube, j, k, l, n_q
real ::  T, temp, Cp, NuT 

temp = p0/(rho0*R) 
Cp = k*Cv

    do icube = 1, n_cube
       do l = 1, n_cellz
          do k = 1, n_celly
             do j = 1, n_cellx
                
				T = avg_field(5, j, k, l, icube)

				NuT = (mu*(T / temp)**1.5 * (temp + 110.4)/(T + 110.4)) &
				/avg_field(1, j, k, l, icube)
				
				this_field(5, j, k, l, icube) = NuT * this_field(1, j, k, l, icube)  
				! Dissipation rate of Turbulence Kinetic Energy(Epsilon) is stored in (5, ...)
				
				this_field(1, j, k, l, icube) = sqrt(NuT)/((this_field(1, j, k, l, icube))**0.25)  
				! Kolmogorov length scale (Eta) is stored in (1, ...)
				
				this_field(2, j, k, l, icube) = this_field(2, j, k, l, icube) *Cp &
				* avg_field(1, j, k, l, icube)
				
				this_field(3, j, k, l, icube) = this_field(2, j, k, l, icube) *Cp &
				* avg_field(1, j, k, l, icube)
				
				this_field(4, j, k, l, icube) = this_field(2, j, k, l, icube) *Cp &
				* avg_field(1, j, k, l, icube)
                        
             end do
          end do
       end do
    end do

end subroutine Kol_epsilon_EHF

end module subs

program Kol

use subs

implicit none

integer :: n_cube, n_cellx, n_celly, n_cellz
integer :: mn_cube, mn_cellx, mn_celly, mn_cellz
integer :: i,j,k,l,h,q,step1,step2,icube,ii,jj,kk,ll,hh, icount, nstep
real :: rdum
real :: scalar
character(len=10) :: id_str1, id_str2
real, allocatable :: qcel(:,:,:,:,:)
real, allocatable :: qprim(:,:,:,:,:)
!real, allocatable :: qsum(:,:,:,:,:)
real, allocatable :: qavg(:,:,:,:,:)
real, allocatable :: qfc(:,:,:,:,:)
real, allocatable :: qKol_EHF(:,:,:,:,:)
real, allocatable :: q_out(:,:,:,:,:)

integer :: p3dq(5) = (/1,2,3,4,5/)

real, allocatable :: x(:,:,:,:)  ! grid nodes
real, allocatable :: y(:,:,:,:)  ! grid nodes
real, allocatable :: z(:,:,:,:)  ! grid nodes

real :: mu, p0, Pr, rho0, R, Cv


!!!! ================ USER INPUT ================ !!!!


integer :: step_start = 201!start step
integer :: step_end = 400 !end step
integer :: N_output_gap = 1 !interval of each step 

	icount = 0
	
	mu = 0.0000185d0 
	p0 = 101300.0d0
	Pr = 1.4d0
	rho0 = 1.1842d0
	R = 287.0d0
	Cv = 717.5d0
	

!!!! ================ USER INPUT ================ !!!!
  

	write(*,*) "Start reading mesh data..."
	open ( unit=1, form='unformatted', file='mesh.g', convert='little_endian') 
	read(1) mn_cube
	read(1) (mn_cellx, mn_celly, mn_cellz, i = 1, mn_cube)
	!allocate(qsum(5, mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(qavg(5, mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(qfc(5, mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(qKol_EHF(5, mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(q_out(5, mn_cellx, mn_celly, mn_cellz, mn_cube))
	
	allocate(x(mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(y(mn_cellx, mn_celly, mn_cellz, mn_cube))
	allocate(z(mn_cellx, mn_celly, mn_cellz, mn_cube))
	
	
	!qsum = 0.0
	qavg = 0.0
	qfc = 0.0
	qKol_EHF = 0.0
	q_out = 0.0
	
	x=0.0
	y=0.0
	z=0.0
	
	do h = 1, mn_cube
		read(1) &		
           ((( x(j,k,l,h), j=1,mn_cellx), k=1,mn_celly), l=1,mn_cellz),&
           ((( y(j,k,l,h), j=1,mn_cellx), k=1,mn_celly), l=1,mn_cellz),&
           ((( z(j,k,l,h), j=1,mn_cellx), k=1,mn_celly), l=1,mn_cellz)
	end do
	
	close(1)

	write(*,*) 'Starting getting Average results'

	do step1 = step_start,step_end,N_output_gap
		write(*,*) 'step = ',step1
		write(id_str1,'(i10.10)') step1	! an integer of width 10 with zeros at the left
		open ( unit=2, form='unformatted', file='field_'//id_str1//'.q', convert='little_endian') 
		read(2) n_cube
		read(2) (n_cellx, n_celly, n_cellz, i = 1, n_cube)  
		
		allocate(qcel(5, n_cellx, n_celly, n_cellz, n_cube))
		allocate(qprim(5, n_cellx, n_celly, n_cellz, n_cube))

		qcel = 0.0
		qprim = 0.0
		
		do i = 1, n_cube
			read(2) rdum, rdum, rdum, rdum
			read(2) ((((qcel(q, j, k, l, i), &
			j = 1, n_cellx), k = 1, n_celly), l = 1, n_cellz), &
			q = 1, 5)
		end do	
		close(2)
		
		call prim(qprim, qcel, n_cube, n_cellx, n_celly, n_cellz)
		call add(qavg, qprim, n_cube, n_cellx, n_celly, n_cellz)
		
		icount = icount + 1
		
		if (allocated(qcel)) then
			deallocate(qcel)
		end if
		if (allocated(qprim)) then
			deallocate(qprim)
		end if
	end do
		
	write(*,*) 'the average number: ' , icount
	scalar = 1.0/icount
	
	call avg(qavg, scalar, n_cube, n_cellx, n_celly, n_cellz)
	
	open(10, FILE='qavg.q', FORM='unformatted', &
        status='replace', action='write', convert='little_endian')
	write(10) n_cube
	write(10) (n_cellx, n_celly, n_cellz, i = 1, n_cube)
	do i = 1, n_cube
		write(10) real(1d0), real(1d0), real(100d0), real(0.001)
		write(10) &	
         ((((qavg(p3dq(q), j, k, l, i), &
         j = 1, n_cellx), &
         k = 1, n_celly), &
         l = 1, n_cellz), &
         q = 1, 5)
	end do	
	write(*,*)" Output end "
	close(10)
	
	write(*,*) 'Average results OVER'
	
	!----------------------------------------------------------------!
	
	icount = 0
	write(*,*) 'Starting getting Fluctuating results'
	
	
	do step2 = step_start,step_end,N_output_gap
		
		write(*,*) 'step = ',step2
		write(id_str2,'(i10.10)') step2	! an integer of width 10 with zeros at the left
		open ( unit=3, form='unformatted', file='field_'//id_str2//'.q', convert='little_endian') 
		read(3) n_cube
		read(3) (n_cellx, n_celly, n_cellz, i = 1, n_cube)  
		
		allocate(qcel(5, n_cellx, n_celly, n_cellz, n_cube))
		allocate(qprim(5, n_cellx, n_celly, n_cellz, n_cube))

		qcel = 0.0
		qprim = 0.0
		
		do i = 1, n_cube
			read(3) rdum, rdum, rdum, rdum
			read(3) ((((qcel(q, j, k, l, i), &
			j = 1, n_cellx), k = 1, n_celly), l = 1, n_cellz), &
			q = 1, 5)

		end do	
		close(3)
		
		call prim(qprim, qcel, n_cube, n_cellx, n_celly, n_cellz)
		call fluctuate(qfc, qprim, qavg, n_cube, n_cellx, n_celly, n_cellz)
		call Turbulent_Intensity(qKol_EHF, qavg, qfc, n_cube, n_cellx, n_celly, n_cellz)
		call add(q_out, qKol_EHF, n_cube, n_cellx, n_celly, n_cellz)
		
		icount = icount + 1
		
		if (allocated(qcel)) then
			deallocate(qcel)
		end if
		if (allocated(qprim)) then
			deallocate(qprim)
		end if
		
	end do
	
	write(*,*) 'the average number of Turbulent_Intensity: ' , icount
	scalar = 1.0/icount
	call avg(q_out, scalar, n_cube, n_cellx, n_celly, n_cellz)
	write(*,*) 'Average over'
	
	!----------------------------------------------------------------!
	
	open(10, FILE='Turbulent_Intensity.q', FORM='unformatted', &
        status='replace', action='write', convert='little_endian')
	write(10) n_cube
	write(10) (n_cellx, n_celly, n_cellz, i = 1, n_cube)
	do i = 1, n_cube
		write(10) real(1d0), real(1d0), real(100d0), real(0.001)
		write(10) &	
         ((((q_out(p3dq(q), j, k, l, i), &
         j = 1, n_cellx), &
         k = 1, n_celly), &
         l = 1, n_cellz), &
         q = 1, 5)
	end do	
	write(*,*)" Output end "
	close(10)
	
	
	
	
	
	

	! do step2 = step_start,step_end,N_output_gap
		
		! write(*,*) 'step = ',step2
		! write(id_str2,'(i10.10)') step2	! an integer of width 10 with zeros at the left
		! open ( unit=3, form='unformatted', file='field_'//id_str2//'.q', convert='little_endian') 
		! read(3) n_cube
		! read(3) (n_cellx, n_celly, n_cellz, i = 1, n_cube)  
		
		! allocate(qcel(5, n_cellx, n_celly, n_cellz, n_cube))
		! allocate(qprim(5, n_cellx, n_celly, n_cellz, n_cube))

		! qcel = 0.0
		! qprim = 0.0
		
		! do i = 1, n_cube
			! read(3) rdum, rdum, rdum, rdum
			! read(3) ((((qcel(q, j, k, l, i), &
			! j = 1, n_cellx), k = 1, n_celly), l = 1, n_cellz), &
			! q = 1, 5)

		! end do	
		! close(3)
		
		! call prim(qprim, qcel, n_cube, n_cellx, n_celly, n_cellz)
		! call fluctuate(qfc, qprim, qavg, n_cube, n_cellx, n_celly, n_cellz)
		! call EHF(qKol_EHF, qfc, n_cube, n_cellx, n_celly, n_cellz)
		! call SS(qKol_EHF, qfc, qavg, x, y, z, n_cube, n_cellx, n_celly, n_cellz)
		! call add(q_out, qKol_EHF, n_cube, n_cellx, n_celly, n_cellz)
		
		! icount = icount + 1
		
		! if (allocated(qcel)) then
			! deallocate(qcel)
		! end if
		! if (allocated(qprim)) then
			! deallocate(qprim)
		! end if
		
	! end do
	
	! write(*,*) 'the average number of Kol_EHF: ' , icount
	! scalar = 1.0/icount
	! call avg(q_out, scalar, n_cube, n_cellx, n_celly, n_cellz)
	! write(*,*) 'Average over'
	
	! call Kol_epsilon_EHF(q_out, qavg, mu, p0, Pr, rho0, R, Cv, n_cube, n_cellx, n_celly, n_cellz)
	
	! write(*,*) "Kol EHF Calculation over, outputing Results" 
	
	! !----------------------------------------------------------------!
	
	! open(10, FILE='Kol_EHF.q', FORM='unformatted', &
        ! status='replace', action='write', convert='little_endian')
	! write(10) n_cube
	! write(10) (n_cellx, n_celly, n_cellz, i = 1, n_cube)
	! do i = 1, n_cube
		! write(10) real(1d0), real(1d0), real(100d0), real(0.001)
		! write(10) &	
         ! ((((q_out(p3dq(q), j, k, l, i), &
         ! j = 1, n_cellx), &
         ! k = 1, n_celly), &
         ! l = 1, n_cellz), &
         ! q = 1, 5)
	! end do	
	! write(*,*)" Output end "
	! close(10)
	
	
end program Kol
