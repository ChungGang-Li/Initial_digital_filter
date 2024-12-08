program cns_example
  use bcm
  !use engine

  implicit none

  type(mesh_t) :: mesh
  type(flowfield_t) :: field
  type(field_bc_t) :: bc                       
  type(field_eval_t) :: eval
  type(icns_t) :: icns
  type(cns_t) :: cns
  type(param_t) :: param
  type(geom_t) :: geom
  type(cad_t), allocatable, target :: cad(:)
  type(openshell_t), allocatable, target :: opsh(:)
  type(lagfield_t), allocatable, target :: lag(:)
  type(moving_param_t), allocatable, target :: mparam(:)
  type(moving_t) :: moving
  character(len=FNAME_LEN) :: fname,probe_param
  real(kind=dp) :: Pfix(2)
  real(kind=dp) :: vel1(3),vel2(3),vel3(3),vel4(3),vel5(3),&
                   vel6(3),vel7(3),vel8(3),vel9(3),vel10(3)
 
  procedure(bc_apply) :: CNS_absorbing_BC
  procedure(bc_apply), pointer :: bcp
  
  !type(lkin_procedure_t), allocatable, target :: lkin(:)
  
  call parameters_set_default(param)

  call bcm_init
  
  ! probe_param = 'probe_param.dat'
  ! call parse_parameters(trim(probe_param),param)

  param%dt = 0.001
  param%CFL = 1000d0
 
  
  param%dt_dim = param%dt
  param%max_iter = 10
  param%lag_inner_iter = 1
  param%Simulation_type = 0
  
  ! param%NCx = -1
  ! param%gravity_c = 0.5993
   param%solver_wall_func = 4
   ! param%solver_turb = 9
   ! param%ADA_freq = 10
  
  
  param%output_average = .true.
  param%output_rms = .true.
  param%T_average = 1000   

  param%monitor_freq = 10
  param%sample_freq = 1
  param%clcd_freq = 10

  param%rho0 = 1.1842d0
  param%p0 = 101300d0
  param%pinf = 101300d0
  param%uinf = 1.776d0
  param%eps = 1.0d0
  
  param%e = 0.0001d0
 ! param%e = 1.0d0
  param%k = 1.4d0
  param%R = 287d0
  param%dtau = param%dt
  param%Pr = 0.72d0
  param%mu = 0.0000185d0
  param%cv = 717.5d0  
  param%T_end =36

  param%jstop = '72:00:00'
  param%deadend = .true.
  
  
  param%lag_IBM_type = 2
  param%lag_thermal_bc_type = 1
  param%output_clear = .false.
  

  param%rel_tol = -0.0001d0
  param%abs_tol = -0.0001d0
  param%intf_diagonal = .false.
  
  
  param%restart_geom = .false.
  param%split_cad_input = .false.
  
  
  param%binary = .false.
  param%plot3d = .true.

  param%rank_dir = .false.
  param%output_lag_stl = .false.
  param%output_lag_pcd = .false.
  param%output_uns = .false.
  param%output_map = .false.
  param%outpath = './output'
  param%output_P3Dcnt = .true.
  
  
  ! param%rho_limit = 0.5
  ! param%rholimit = 10.5
  ! param%u_limit = -10d0
  ! param%ulimit = 10d0
  ! param%v_limit = -10d0
  ! param%vlimit = 10d0
  ! param%w_limit = -10d0
  ! param%wlimit = 10d0
  ! param%Prs_limit = 0.5*param%p0
  ! param%Prslimit = 22.5*param%p0
  

  ! Load mesh, and allocate a 4 volume wide halo region
  fname = '32cell.bin'
  call mesh_init(mesh, param, fname, 4)

  ! Create an compressible flow field for the given mesh
  call flow_field_init(field, mesh, FIELD_CNS)

  ! Initialize boundary conditions
  
  call flow_field_add_scratch(field, 8)
  
  field%scratch(6)%q(1, :, :, :, :) = param%CFL
  field%scratch(6)%q(2, :, :, :, :) = 1d0                                         
  
  
  call lag_parameter_set_default(param)
  


  ! Initialize boundary conditions
  bcp => CNS_absorbing_BC
  call field_bc_init(bc, mesh)
  
  
  call field_bc_set(bc, field%type, XM, bcp)
  call field_bc_set(bc, field%type, XP, bcp)
  call field_bc_set(bc, field%type, YP, bcp)
  call field_bc_set(bc, field%type, YM, bcp)
  call field_bc_set(bc, field%type, ZP, bcp)
  call field_bc_set(bc, field%type, ZM, bcp)
  

  ! Set initial conditions
  ! field = (/param%rho0, &
	        ! param%rho0*param%uinf, &
			   ! 0.0d0, &
			   ! 0.0d0, &
			   ! param%pinf/(param%K-1d0)+0.5d0*param%rho0*param%uinf*param%uinf /)
			   
  ! Set initial conditions
  ! field = (/param%rho0, &
	        ! param%rho0*param%uinf, &
			   ! 0.0d0, &
			   ! 0.0d0, &
			   ! param%pinf/(param%K-1d0)+0.5d0*param%rho0*param%uinf*param%uinf /)
			   
  call own_user_flow_field_init(field,param)
  
         

   ! call restart_init(param, P3DCNT, 'field_cnt_0000004000.q')
 
   ! call restart_init(param, P3DNOD, 'field_0000004000.q')

  ! Set geometry
  ! link geometry pointer
  geom%cad => cad
  geom%opsh => opsh
  geom%lag => lag
  
  
  param%lag_COMM_TYPE = 7
  
  
  fname = "Channel.stl"
  call geom_init(geom, mesh, param, fname, cad, lag(:))
  geom%lag%T_fix = 298.0592 
                            
  call cns_init(cns, field, geom, bc, eval, param, lag(:))

  call cns_set_flux(cns, ROE_LM)
  call cns_set_viscous(cns, CE2ND)
  call cns_set_tstep(cns, LUSGS_SLTS_MIX)
  
  
  ! Initialize flow field evaluations
 call field_eval_init(eval, param, './clcd_cns.txt', VOLUME)
!  call field_eval_init(eval, param, './clcd_iblank.txt', IBLANK)
!  call field_eval_init(eval, param, './clcd_openshell.txt', OPNSHL)
!  call field_eval_init(eval, param, './clcd_lagrange.txt', LAGRAN)
!  call field_eval_init(eval, param, './clcd_lvolume.txt', LVOLUME)
!  call field_eval_init(eval, param, './clcd_osurf.txt', OSURF)
  !!!!
  
  
  call cns_solve(cns, param)
  
  call cns_free(cns)
  
  call bcm_finalize
 
 
end program cns_example





!! ------------------------------------------------------------------------- !!
!! ---------------------- User defined function (UDF) ---------------------- !!


subroutine CNS_absorbing_BC(bc, field, param, t, pos, nil0, nil1, nil2)

	use boundary
	use flow_field
	use parameters

	implicit none

    
	type(boundary_t), intent(in) :: bc    
	type(flowfield_t), intent(inout) :: field
	type(param_t), intent(in) :: param
	integer, intent(in), optional :: pos
	real(kind=dp), intent(in) :: t
	integer, intent(in), optional :: nil1, nil2
	integer :: ierr
	real(kind=dp), intent(in), optional :: nil0
	integer :: i, j, k ,l, m, q, icube
	real(kind=dp) :: rho,u,v,w,p,vv,c
	real(kind=dp) :: dx, dy, dz
	real(kind=dp) :: U_total, V_total, Volume
	real(kind=dp) :: gU_total, gV_total, gVolume
	real(kind=dp) :: f, fave, u_tau, Uave, heigh
	
	Uave = param%uinf*param%rho0
	
	U_total = 0d0
	V_total = 0d0
	
	!$omp single   
	
	do icube = 1, field%bcm%n_cube
		
      dx = field%bcm%dxcube(icube)
      dy = field%bcm%dycube(icube)
      dz = field%bcm%dzcube(icube)
    
		do l = n_band + 1, n_band + field%bcm%n_cellz
			do k = n_band + 1, n_band + field%bcm%n_celly
				do j = n_band + 1, n_band + field%bcm%n_cellx
				
					u = field%qcnt(2, j, k, l, icube)
								
					U_total = U_total+u*dx*dy*dz
					V_total = V_total+dx*dy*dz
    
				enddo
			enddo
		enddo
		
	enddo
	
	call MPI_AllReduce(U_total, gU_total, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
	   
	call MPI_AllReduce(V_total, gV_total, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
	
	U_total = gU_total/gV_total
	
	f = (Uave-U_total)/param%dt
	
	do icube = 1, field%bcm%n_cube
		
		do l = n_band + 1, n_band + field%bcm%n_cellz
			do k = n_band + 1, n_band + field%bcm%n_celly
				do j = n_band + 1, n_band + field%bcm%n_cellx
				
					u = field%qcnt(2, j, k, l, icube)/field%qcnt(1, j, k, l, icube)
								
					field%scratch(1)%q(1, j, k, l, icube) = 0.0d0 
					field%scratch(1)%q(2, j, k, l, icube) = f
					field%scratch(1)%q(3, j, k, l, icube) = 0.0d0 
					field%scratch(1)%q(4, j, k, l, icube) = 0.0d0 
					field%scratch(1)%q(5, j, k, l, icube) = f*u
					
				enddo
			enddo
		enddo
		
	enddo
	
	!$omp end single
	
end subroutine CNS_absorbing_BC



subroutine own_user_flow_field_init(field, param)

  use flow_field
  use field_intf			
  use parameters
  
  implicit none

  type(flowfield_t), intent(inout) ::field
  type(param_t), intent(in) :: param
  integer i, j, k, l, ii, jj, kk, ll, icube
  
  real(kind=dp) :: u,v,w,dx,dy,dz
  real(kind=dp) :: Re_tau, heigh, U_tau, tau_w, y_plus, Upro, yy, intensity, pii, temp
  
  integer M, N, My, Mz, nx, ny, nz, nxx, nyy, nzz, step
  integer :: offset_j, offset_k, offset_l									 
  
  real(8), allocatable :: U1_x(:,:,:,:)
  real(8), allocatable :: U2_x(:,:,:,:)
  real(8), allocatable :: U1_y(:,:,:,:)
  real(8), allocatable :: U2_y(:,:,:,:)
  real(8), allocatable :: U1_Z(:,:,:,:)
  real(8), allocatable :: U2_z(:,:,:,:)
  real(8), allocatable :: randomx(:,:,:,:),randomy(:,:,:,:),randomz(:,:,:,:)
  real(8), allocatable :: Bi(:,:), Bj(:,:), Bk(:,:)
  real(8), allocatable :: filter_kernel(:,:,:,:)
  real(8), allocatable :: u_x(:,:,:)
  real(8), allocatable :: v_x(:,:,:)
  real(8), allocatable :: w_x(:,:,:)
  
  
  nx = 2;		 
  ny = 2;
  nz = 2;
  
  intensity = 20;
  
  M = field%bcm%n_cellx
  My = field%bcm%n_celly
  Mz = field%bcm%n_cellz
  
  nyy = 2 * ny
  nzz = 2 * nz
  nxx = 2 * nx
  
  pii = 3.14159
  
  allocate(U1_x(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(U2_x(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(U1_y(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(U2_y(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(U1_Z(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(U2_z(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(randomx(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(randomy(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  allocate(randomz(2*nxx+M, 2*nyy+My, 2*nzz+Mz, field%bcm%n_cube))
  
  allocate(Bi(-nxx:nxx, field%bcm%n_cube), Bj(-nyy:nyy, field%bcm%n_cube), Bk(-nzz:nzz, field%bcm%n_cube))
  allocate(filter_kernel(-nxx:nxx, -nyy:nyy, -nzz:nzz, field%bcm%n_cube))
  
  allocate(u_x(My, Mz, field%bcm%n_cube), v_x(My, Mz, field%bcm%n_cube), w_x(My, Mz, field%bcm%n_cube))
  
  
  call random_number(U1_x)
  call random_number(U2_x)
  call random_number(U1_y)
  call random_number(U2_y)
  call random_number(U1_z)
  call random_number(U2_z)
  
  U1_x = max(U1_x, 1.0d-10)
  U2_x = max(U2_x, 1.0d-10)
  U1_y = max(U1_y, 1.0d-10)
  U2_y = max(U2_y, 1.0d-10)
  U1_z = max(U1_z, 1.0d-10)
  U2_z = max(U2_z, 1.0d-10)

  
  randomx = intensity * sqrt(-2.0d0 * log(U1_x)) * cos(2.0d0 * pii * U2_x)
  randomy = intensity * sqrt(-2.0d0 * log(U1_y)) * cos(2.0d0 * pii * U2_y)
  randomz = intensity * sqrt(-2.0d0 * log(U1_z)) * cos(2.0d0 * pii * U2_z)
  
  deallocate(U1_x, U2_x, U1_y, U2_y, U1_z, U2_z)
  
  
  !$omp do
  
   do i = 1, field%bcm%n_cube
      do l = n_band + 1, n_band + field%bcm%n_cellz
        do k = n_band + 1, n_band + field%bcm%n_celly
           do j = n_band + 1, n_band + field%bcm%n_cellx

			    jj = j-n_band+nx
				kk = k-n_band+ny
				ll = l-n_band+nz
		   
				field%scratch(2)%q(1, j, k, l, i) = randomx(jj, kk, ll, i)
			    field%scratch(2)%q(2, j, k, l, i) = randomy(jj, kk, ll, i)
			    field%scratch(2)%q(3, j, k, l, i) = randomz(jj, kk, ll, i)
            end do
         end do
      end do
   end do
	!$omp end do
   
   
  !$omp single
  call field_intf_q(field,1,3,field%scratch(2)%qp)
  call field_intf_finalize(field,1,3,field%scratch(2)%qp)
  
   if(param%intf_diagonal) then

       call field_intf_q_edge(field,1,3,field%scratch(2)%qp)
       call field_intf_finalize_edge(field,1,3,field%scratch(2)%qp)
       call field_intf_q_corner(field,1,3,field%scratch(2)%qp)
       call field_intf_finalize_corner(field,1,3,field%scratch(2)%qp)

    end if

	!$omp end single			 


   offset_j = (2*nxx+M  - field%bcm%m_cellx) / 2
   offset_k = (2*nyy+My - field%bcm%m_celly) / 2
   offset_l = (2*nzz+Mz - field%bcm%m_cellz) / 2
   
  !$omp do
   do i = 1, field%bcm%n_cube
      do l = 1, 2*nzz+Mz
         do k = 1, 2*nyy+My
            do j = 1, 2*nxx+M
   
               ! 對應索引
               if (j > offset_j .and. j <= field%bcm%m_cellx + offset_j .and. &
                   k > offset_k .and. k <= field%bcm%m_celly + offset_k .and. &
                   l > offset_l .and. l <= field%bcm%m_cellz + offset_l) then
   
                  randomx(j, k, l, i) = field%scratch(2)%q(1, j - offset_j, k - offset_k, l - offset_l, i)
                  randomy(j, k, l, i) = field%scratch(2)%q(2, j - offset_j, k - offset_k, l - offset_l, i)
                  randomz(j, k, l, i) = field%scratch(2)%q(3, j - offset_j, k - offset_k, l - offset_l, i)
               else
                  randomx(j, k, l, i) = 0.0d0
                  randomy(j, k, l, i) = 0.0d0
                  randomz(j, k, l, i) = 0.0d0
               endif
            end do
         end do
      end do
   end do
	!$omp end do
   
  
  !$omp do	  
  do icube = 1, field%bcm%n_cube
	  
     do i = -nxx, nxx
		Bi(i, icube) = exp(-pii * i**2 / (2.0d0 * nxx**2))
     end do
     
	 do j = -nyy, nyy
		Bj(j, icube) = exp(-pii * j**2 / (2.0d0 * nyy**2))
     end do
     
	 do k = -nzz, nzz
		Bk(k, icube) = exp(-pii * k**2 / (2.0d0 * nzz**2))
     end do
	 
	 temp = 0d0
	 
	 temp = temp + sum(Bi(:, icube)**2) * sum(Bj(:, icube)**2) * sum(Bk(:, icube)**2)
	 temp = max(sqrt(temp),1e-10)
	 
	 do i = -nxx, nxx
       Bi(i, icube) = Bi(i, icube) / temp
     end do
     do i = -nyy, nyy
       Bj(i, icube) = Bj(i, icube) / temp
     end do
     do i = -nzz, nzz
       Bk(i, icube) = Bk(i, icube) / temp
     end do
	 
	   do i = -nxx, nxx
		  do j = -nyy, nyy
			 do k = -nzz, nzz
			 
				filter_kernel(i, j, k, icube) = Bi(i, icube) * Bj(j, icube) * Bk(k, icube)
				
			 end do
		  end do
	   end do
	   
   end do
   !$omp end do
   
   
   
	!$omp do	  
   do icube = 1, field%bcm%n_cube
    
		do step = 1, M
		
		  do k = 1, My
			do l = 1, Mz
			
				 u_x(k, l, icube) = sum(filter_kernel(:, :, :, icube) * randomx(step:step+2*nxx, k:k+2*nyy, l:l+2*nzz, icube))
				 v_x(k, l, icube) = sum(filter_kernel(:, :, :, icube) * randomy(step:step+2*nxx, k:k+2*nyy, l:l+2*nzz, icube))
				 w_x(k, l, icube) = sum(filter_kernel(:, :, :, icube) * randomz(step:step+2*nxx, k:k+2*nyy, l:l+2*nzz, icube))
				 
			end do
		  end do


		  do k = n_band + 1, n_band + field%bcm%n_celly
           do l = n_band + 1, n_band + field%bcm%n_cellz
		   
			  field%qcnt(2, step+n_band, k, l, icube) = u_x(k-n_band, l-n_band, icube)
			  field%qcnt(3, step+n_band, k, l, icube) = v_x(k-n_band, l-n_band, icube)
			  field%qcnt(4, step+n_band, k, l, icube) = w_x(k-n_band, l-n_band, icube)
			  
			end do
		  end do
		  
		  
		end do
		
		
   end do
   
   !$omp end do

  
  
  
  
  Re_tau = 180d0
  heigh = 0.05
  U_tau = Re_tau*param%mu/param%rho0/(0.5*heigh)
  tau_w = U_tau*U_tau*param%rho0

  !$omp do
    do i = 1, field%bcm%n_cube
    
      dx = field%bcm%dxcube(i)
      dy = field%bcm%dycube(i)
      dz = field%bcm%dzcube(i)
    
      do l = 1, field%bcm%m_cellz
         do k = 1, field%bcm%m_celly
            do j = 1, field%bcm%m_cellx
            
                  yy = 0.5*heigh-abs(field%bcm%ycube(i)+(k-4-0.5)*dy)
                  
                  y_plus = param%rho0*U_tau*yy/param%mu
                  
                  if (y_plus < 10d0) then
                  
                    Upro = tau_w/param%mu*yy
                  
                  else
                  
                    Upro = (1d0/0.4*log(y_plus)+5.5)*U_tau
                    
                  endif
                  
                  
                  u = field%qcnt(2, j, k, l, i)/param%rho0+Upro
				  v = field%qcnt(3, j, k, l, i)/param%rho0
                  w = field%qcnt(4, j, k, l, i)/param%rho0
            
                  field%qcnt(1, j, k, l, i) = param%rho0
                  field%qcnt(2, j, k, l, i) = param%rho0*u
                  field%qcnt(3, j, k, l, i) = param%rho0*v
                  field%qcnt(4, j, k, l, i) = param%rho0*w
                  field%qcnt(5, j, k, l, i) = param%P0 / (param%k - 1d0) + 0.5*param%rho0*(u*u+v*v+w*w)
				  
            end do
         end do
      end do
   end do
  !$omp end do
  
  deallocate(randomx, randomy, randomz)
  deallocate(Bi, Bj, Bk)
  deallocate(filter_kernel)
  deallocate(u_x, v_x, w_x)
  
  
end subroutine own_user_flow_field_init



subroutine Inlet_gen(Upro)

  implicit none

  real(8), intent(inout) :: Upro
  integer icube, j, k, l
  real(8) :: X
  
  real(8) :: dx,dy,dz

  ! call random_seed()
  call random_number(X)
 
  Upro = Upro + Upro*(X-0.5)*0.44
  
end subroutine Inlet_gen	 




!! ---------------------- User defined function (UDF) ---------------------- !!
!! ------------------------------------------------------------------------- !!
