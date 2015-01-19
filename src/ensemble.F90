module ensemble_variables
 use vertical
 use global
 use writeout_variables
 
 implicit none
 
 type send_ens
  real(kind=8), allocatable :: z(:)
 end type send_ens
 
 type synced_ens
  type(send_ens) :: send
  integer :: send_req
  type(send_ens), allocatable :: recv(:)
  integer, allocatable :: recv_req(:)
  type(var) :: z_send
  type(var) :: z
  type(out_var) :: out
  integer :: tag
 end type
 
 
 type unsynced_ens
  type(send_ens) :: send
  integer :: send_req
  type(send_ens), allocatable :: recv(:)
  integer, allocatable :: recv_req(:)
  type(var) :: z_send
  type(var) :: z
  type(out_var) :: out
  integer :: tag
 end type
  
 type unsynced_ensfin
  type(var) :: z
  type(out_var) :: out
 end type
 
 type synced_ensfin
  type(var) :: z
  type(out_var) :: out
 end type
 
 integer, save :: ens_max_tag=0
  
 type(unsynced_ens) :: a_h_h(nz), a_h_u(nz), a_h_v(nz), a_h_z(nz)
 type(synced_ens) :: b_u_u(nz), b_v_v(nz)
 type(unsynced_ens) :: b_uu_h(nz), b_vv_h(nz), b_uv_z(nz), b_uv_h(nz)
 type(unsynced_ens) :: a_dd_h(0:nz), a_dpx_u(0:nz), a_dpy_v(0:nz)
 type(synced_ens) :: a_d_h(0:nz)
 type(unsynced_ens) :: a_p_h(0:nz), a_px_u(0:nz), a_py_v(0:nz)
 type(unsynced_ens) :: b_lapu_u(nz), b_lapv_v(nz)
 type(unsynced_ens) :: b_lap2u_u(nz), b_lap2v_v(nz)
 type(unsynced_ensfin) :: bp_uu_h(nz), bp_vv_h(nz), bp_uv_z(nz), bp_uv_h(nz)
 type(synced_ensfin) :: ap_dd_h(0:nz)
 type(unsynced_ensfin) :: ap_dpx_u(0:nz), ap_dpy_v(0:nz)
 type(unsynced_ensfin) :: tend_u1(nz), tend_u2(nz), tend_u3(nz), tend_u4(nz), tend_u5(nz), &
      tend_u6(nz), tend_u7(nz), tend_u8(nz), tend_u9(nz), tend_u10(nz)
 type(unsynced_ensfin) :: tend_v1(nz), tend_v2(nz), tend_v3(nz), tend_v4(nz), tend_v5(nz),  &
      tend_v6(nz), tend_v7(nz), tend_v8(nz), tend_v9(nz), tend_v10(nz)
 type(synced_ensfin) :: aq_huu_h(nz), aq_hvv_h(nz), aq_dd_h(0:nz)
 type(synced_ensfin) :: h_bp_uu_h(nz), h_bp_vv_h(nz)
 type(unsynced_ensfin) :: aq_huv_z(nz), h_bp_uv_z(nz)

end module

module ensemble

 implicit none
 
 contains
 
 subroutine init_synced_ens(dat,grd,vec)
  use ensemble_variables
  use parallel
  use writeout
  use grid
  
  implicit none
  
  type(synced_ens), intent(out) :: dat
  integer, intent(in) :: grd(:,:)
  logical, intent(in) :: vec
  integer :: i
  
  dat%tag=20000+ens_max_tag
  ens_max_tag=ens_max_tag+1
  
  if (ubound(grd,1)-mx == 1) then
   if (ubound(grd,2)-my == 1) then
    call create_zparam(dat%z_send)
   else
    call create_uparam(dat%z_send)
   end if
  else
   if (ubound(grd,2)-my == 1) then
    call create_vparam(dat%z_send)
   else
    call create_hparam(dat%z_send)
   end if
  end if
  
  allocate(dat%send%z(count(dat%z_send%z%core)))
  dat%send_req=mpi_request_null
  
  if (ens_name == proc_master) then
  
   call create_var(dat%z,ubound(grd,1)-mx,ubound(grd,2)-my,vec)
   call init_writeouts(dat%out,grd)
  
   allocate(dat%recv(0:ens_num-1))
   do i=0,ens_num-1
    allocate(dat%recv(i)%z(count(dat%z%z%core)))
   end do
   allocate(dat%recv_req(0:ens_num-1))
   dat%recv_req=mpi_request_null
  end if

 end subroutine
 
 subroutine init_unsynced_ens(dat,grd)
  use ensemble_variables
  use parallel
  use writeout
  use grid
  
  implicit none
  
  type(unsynced_ens), intent(out) :: dat
  integer, intent(in) :: grd(:,:)
  integer :: i
  
  dat%tag=20000+ens_max_tag
  ens_max_tag=ens_max_tag+1
  
  if (ubound(grd,1)-mx == 1) then
   if (ubound(grd,2)-my == 1) then
    call create_zparam(dat%z_send)
   else
    call create_uparam(dat%z_send)
   end if
  else
   if (ubound(grd,2)-my == 1) then
    call create_vparam(dat%z_send)
   else
    call create_hparam(dat%z_send)
   end if
  end if
  
  allocate(dat%send%z(count(dat%z_send%z%core)))
  dat%send_req=mpi_request_null
  
  if (ens_name == proc_master) then
  
   if (ubound(grd,1)-mx == 1) then
    if (ubound(grd,2)-my == 1) then
     call create_zparam(dat%z)
    else
     call create_uparam(dat%z)
    end if
   else
    if (ubound(grd,2)-my == 1) then
     call create_vparam(dat%z)
    else
     call create_hparam(dat%z)
    end if
   end if
   call init_writeouts(dat%out,grd)
   
   allocate(dat%recv(0:ens_num-1))
   do i=0,ens_num-1
    allocate(dat%recv(i)%z(count(dat%z%z%core)))
   end do
   allocate(dat%recv_req(0:ens_num-1))
   dat%recv_req=mpi_request_null
  end if

 end subroutine
 
 
 subroutine init_ensemble()
  use parallel
  use ensemble_variables
  use writeout
  use grid
  
  implicit none
  
  integer :: i, j
  
  
  
  do i=0,nz
   call init_unsynced_ens(a_dd_h(i),hgrid_out)
   call init_synced_ens(a_d_h(i),hgrid_out,.false.)
   call init_unsynced_ens(a_dpx_u(i),ugrid_out)
   call init_unsynced_ens(a_dpy_v(i),vgrid_out)
   call init_unsynced_ens(a_p_h(i),hgrid_out)
   call init_unsynced_ens(a_px_u(i),ugrid_out)
   call init_unsynced_ens(a_py_v(i),vgrid_out)
   
   call create_var(ap_dd_h(i)%z,0,0,.false.)
   call create_var(aq_dd_h(i)%z,0,0,.false.)
   call init_writeouts(ap_dd_h(i)%out,hgrid_out)
   call init_writeouts(aq_dd_h(i)%out,hgrid_out)

   call create_uparam(ap_dpx_u(i)%z)
   call create_vparam(ap_dpy_v(i)%z)
   call init_writeouts(ap_dpx_u(i)%out,ugrid_out)
   call init_writeouts(ap_dpy_v(i)%out,vgrid_out)

  end do
  do i=1,nz
  
   
   call init_unsynced_ens(a_h_h(i),hgrid_out)
   call init_unsynced_ens(a_h_u(i),ugrid_out)
   call init_unsynced_ens(a_h_v(i),vgrid_out)
   call init_unsynced_ens(a_h_z(i),zgrid_out)
   call init_synced_ens(b_u_u(i),ugrid_out,.true.)
   call init_synced_ens(b_v_v(i),vgrid_out,.true.)
   call init_unsynced_ens(b_uu_h(i),hgrid_out)
   call init_unsynced_ens(b_vv_h(i),hgrid_out)
   call init_unsynced_ens(b_uv_h(i),hgrid_out)
   call init_unsynced_ens(b_uv_z(i),zgrid_out)
   call init_unsynced_ens(b_lapu_u(i),ugrid_out)
   call init_unsynced_ens(b_lapv_v(i),vgrid_out)
   call init_unsynced_ens(b_lap2u_u(i),ugrid_out)
   call init_unsynced_ens(b_lap2v_v(i),vgrid_out)
   
   call create_hparam(bp_uu_h(i)%z)
   call create_hparam(bp_vv_h(i)%z)
   call create_zparam(bp_uv_z(i)%z)
   call create_hparam(bp_uv_h(i)%z)
   call create_var(aq_huu_h(i)%z,0,0,.false.)
   call create_var(aq_hvv_h(i)%z,0,0,.false.)
   call create_var(h_bp_uu_h(i)%z,0,0,.false.)
   call create_var(h_bp_vv_h(i)%z,0,0,.false.)
   call create_zparam(aq_huv_z(i)%z)
   call create_zparam(h_bp_uv_z(i)%z)
   
   call init_writeouts(bp_uu_h(i)%out,hgrid_out)
   call init_writeouts(bp_vv_h(i)%out,hgrid_out)
   call init_writeouts(bp_uv_z(i)%out,zgrid_out)
   call init_writeouts(bp_uv_h(i)%out,hgrid_out)
   
   call init_writeouts(aq_huu_h(i)%out,hgrid_out)
   call init_writeouts(aq_hvv_h(i)%out,hgrid_out)
   call init_writeouts(h_bp_uu_h(i)%out,hgrid_out)
   call init_writeouts(h_bp_vv_h(i)%out,hgrid_out)
   call init_writeouts(aq_huv_z(i)%out,zgrid_out)
   call init_writeouts(h_bp_uv_z(i)%out,zgrid_out)
   
   call create_uparam(tend_u1(i)%z)
   call init_writeouts(tend_u1(i)%out,ugrid_out)
   call create_uparam(tend_u2(i)%z)
   call init_writeouts(tend_u2(i)%out,ugrid_out)
   call create_uparam(tend_u3(i)%z)
   call init_writeouts(tend_u3(i)%out,ugrid_out)
   call create_uparam(tend_u4(i)%z)
   call init_writeouts(tend_u4(i)%out,ugrid_out)
   call create_uparam(tend_u5(i)%z)
   call init_writeouts(tend_u5(i)%out,ugrid_out)
   call create_uparam(tend_u6(i)%z)
   call init_writeouts(tend_u6(i)%out,ugrid_out)
   call create_uparam(tend_u7(i)%z)
   call init_writeouts(tend_u7(i)%out,ugrid_out)
   call create_uparam(tend_u8(i)%z)
   call init_writeouts(tend_u8(i)%out,ugrid_out)
   call create_uparam(tend_u9(i)%z)
   call init_writeouts(tend_u9(i)%out,ugrid_out)
   call create_uparam(tend_u10(i)%z)
   call init_writeouts(tend_u10(i)%out,ugrid_out)
   call create_vparam(tend_v1(i)%z)
   call init_writeouts(tend_v1(i)%out,vgrid_out)
   call create_vparam(tend_v2(i)%z)
   call init_writeouts(tend_v2(i)%out,vgrid_out)
   call create_vparam(tend_v3(i)%z)
   call init_writeouts(tend_v3(i)%out,vgrid_out)
   call create_vparam(tend_v4(i)%z)
   call init_writeouts(tend_v4(i)%out,vgrid_out)
   call create_vparam(tend_v5(i)%z)
   call init_writeouts(tend_v5(i)%out,vgrid_out)
   call create_vparam(tend_v6(i)%z)
   call init_writeouts(tend_v6(i)%out,vgrid_out)
   call create_vparam(tend_v7(i)%z)
   call init_writeouts(tend_v7(i)%out,vgrid_out)
   call create_vparam(tend_v8(i)%z)
   call init_writeouts(tend_v8(i)%out,vgrid_out)
   call create_vparam(tend_v9(i)%z)
   call init_writeouts(tend_v9(i)%out,vgrid_out)
   call create_vparam(tend_v10(i)%z)
   call init_writeouts(tend_v10(i)%out,vgrid_out)
  end do
  
 end subroutine
 
 subroutine send_unsynced_ens(dat)
  use ensemble_variables
  use parallel
  
  implicit none
  
  type(unsynced_ens), intent(inout) :: dat
  integer :: i,j,k
  
  
  k=0
  call mpi_wait(dat%send_req,mpi_status_ignore,stat)
  do j=1,dat%z_send%ny
   do i=1,dat%z_send%nx
    if (dat%z_send%z(i,j)%core) then
     k=k+1
     dat%send%z(k)=dat%z_send%z(i,j)%z
    end if
   end do
  end do
  
  if (ens_name == proc_master) then
   do k=0,ens_num-1
    call mpi_wait(dat%recv_req(k),mpi_status_ignore,stat)
    call mpi_irecv(dat%recv(k)%z(1),ubound(dat%recv(k)%z),mpi_double_precision,  &
       k*ens_images+proc_name-ens_master,dat%tag,  &
       mpi_comm_world,dat%recv_req(k),stat)
   end do
  end if
  call mpi_isend(dat%send%z(1),ubound(dat%send%z),mpi_double_precision,  &
       proc_name-ens_master,dat%tag,  &
       mpi_comm_world,dat%send_req,stat)
       
 end subroutine
 
 subroutine recv_unsynced_ens(dat)
  use ensemble_variables
  use parallel
  
  implicit none
  
  type(unsynced_ens), intent(inout) :: dat
  integer :: i,j,k,l
  
  
  if (ens_name == proc_master) then
   dat%z%z%z=0.0d0
   do l=0,ens_num-1
    k=0
    call mpi_wait(dat%recv_req(l),mpi_status_ignore,stat)
    do j=1,dat%z%ny
     do i=1,dat%z%nx
      if (dat%z%z(i,j)%core) then
       k=k+1
       dat%z%z(i,j)%z=dat%z%z(i,j)%z+dat%recv(l)%z(k)
      end if
     end do
    end do
   end do
   dat%z%z%z=dat%z%z%z/dble(ens_num)
  end if
  
 end subroutine
 
 
 subroutine send_synced_ens(dat)
  use ensemble_variables
  use parallel
  
  implicit none
  
  type(synced_ens), intent(inout) :: dat
  integer :: i,j,k
  
  
  k=0
  call mpi_wait(dat%send_req,mpi_status_ignore,stat)
  do j=1,dat%z_send%ny
   do i=1,dat%z_send%nx
    if (dat%z_send%z(i,j)%core) then
     k=k+1
     dat%send%z(k)=dat%z_send%z(i,j)%z
    end if
   end do
  end do
  
  if (ens_name == proc_master) then
   do k=0,ens_num-1
    call mpi_wait(dat%recv_req(k),mpi_status_ignore,stat)
    call mpi_irecv(dat%recv(k)%z(1),ubound(dat%recv(k)%z),mpi_double_precision,  &
       k*ens_images+proc_name-ens_master,dat%tag,  &
       mpi_comm_world,dat%recv_req(k),stat)
   end do
  end if
  call mpi_isend(dat%send%z(1),ubound(dat%send%z),mpi_double_precision,  &
       proc_name-ens_master,dat%tag,  &
       mpi_comm_world,dat%send_req,stat)
       
 end subroutine
 
 subroutine recv_synced_ens(dat)
  use ensemble_variables
  use parallel
  use sync
  
  implicit none
  
  type(synced_ens), intent(inout) :: dat
  integer :: i,j,k,l
  
  
  if (ens_name == proc_master) then
   dat%z%z%z=0.0d0
   do l=0,ens_num-1
    k=0
    call mpi_wait(dat%recv_req(l),mpi_status_ignore,stat)
    do j=1,dat%z%ny
     do i=1,dat%z%nx
      if (dat%z%z(i,j)%core) then
       k=k+1
       dat%z%z(i,j)%z=dat%z%z(i,j)%z+dat%recv(l)%z(k)
      end if
     end do
    end do
   end do
   dat%z%z%z=dat%z%z%z/dble(ens_num)
   call start_sync(dat%z)
  end if
  
  
 end subroutine
  
 
 subroutine integerate_ens()
  use ensemble_variables
  use variables
  use operate
  use sync
  use solver_variables
  use density
  
  implicit none
  
  integer :: i, j
  
  if (rigid_lid) then
   true_d(nz)%z%z=-pres%z%z/gp(nz)
   true_p(nz)%z%z=-pres%z%z
  else
   true_d(nz)%z%z=depth(nz)%z%z
   true_p(nz)%z%z=depth(nz)%z%z*gp(nz)
  end if
  call start_sync(true_d(nz))
  call start_sync(true_p(nz))
  do i=nz-1,0,-1
   true_d(i)%z%z=depth(i)%z%z
   true_p(i)%z%z=0
   do j=nz,i+1,-1
    true_p(i)%z%z=true_p(i)%z%z+(true_d(j)%z%z-true_d(i)%z%z)*gp(j)
   end do
   call start_sync(true_d(i))
   call start_sync(true_p(i))
  end do
  
  
  do i=1,nz
   a_h_h(i)%z_send%z%z=h(i)%z%z
   call send_unsynced_ens(a_h_h(i))
   a_h_u(i)%z_send%z%z=h_u(i)%z%z
   call send_unsynced_ens(a_h_u(i))
   a_h_v(i)%z_send%z%z=h_v(i)%z%z
   call send_unsynced_ens(a_h_v(i))
   a_h_z(i)%z_send%z%z=h_z(i)%z%z
   call send_unsynced_ens(a_h_z(i))
   b_u_u(i)%z_send%z%z=h_u(i)%z%z*u(i)%z%z
   call send_synced_ens(b_u_u(i))
   b_v_v(i)%z_send%z%z=h_v(i)%z%z*v(i)%z%z
   call send_synced_ens(b_v_v(i))
   b_uu_h(i)%z_send%z%z=h(i)%z%z*Ax(u(i)%z%z**2)
   call send_unsynced_ens(b_uu_h(i))
   b_vv_h(i)%z_send%z%z=h(i)%z%z*Ay(v(i)%z%z**2)
   call send_unsynced_ens(b_vv_h(i))
   b_lapu_u(i)%z_send%z%z=h_u(i)%z%z*avaru%z%z*lapu(i)%z%z
   call send_unsynced_ens(b_lapu_u(i))
   b_lapv_v(i)%z_send%z%z=h_v(i)%z%z*avarv%z%z*lapv(i)%z%z
   call send_unsynced_ens(b_lapv_v(i))
   b_lap2u_u(i)%z_send%z%z=h_u(i)%z%z*bvaru%z%z*(sGxx(lapu(i))+Gy(sGy(lapu(i)))) 
   call send_unsynced_ens(b_lap2u_u(i))
   b_lap2v_v(i)%z_send%z%z=h_v(i)%z%z*bvarv%z%z*(sGyy(lapv(i))+Gx(sGx(lapv(i)))) 
   call send_unsynced_ens(b_lap2v_v(i))
   b_uv_z(i)%z_send%z%z=h_z(i)%z%z*sAy(u(i))*sAx(v(i))
   call send_unsynced_ens(b_uv_z(i))
   b_uv_h(i)%z_send%z%z=h(i)%z%z*Ax(u(i)%z%z)*Ay(v(i)%z%z)
   call send_unsynced_ens(b_uv_h(i))
  end do
  
   
  do i=0,nz
   a_dd_h(i)%z_send%z%z=true_d(i)%z%z**2
   call send_unsynced_ens(a_dd_h(i))
   a_dpx_u(i)%z_send%z%z=sAx(true_d(i))*sGx(true_p(i))
   call send_unsynced_ens(a_dpx_u(i))
   a_dpy_v(i)%z_send%z%z=sAy(true_d(i))*sGy(true_p(i))
   call send_unsynced_ens(a_dpy_v(i))
   a_d_h(i)%z_send%z%z=true_d(i)%z%z
   call send_synced_ens(a_d_h(i))
   a_p_h(i)%z_send%z%z=true_p(i)%z%z
   call send_unsynced_ens(a_p_h(i))
   a_px_u(i)%z_send%z%z=sGx(true_p(i))
   call send_unsynced_ens(a_px_u(i))
   a_py_v(i)%z_send%z%z=sGy(true_p(i))
   call send_unsynced_ens(a_py_v(i))
  end do
  
 end subroutine
  
  
  
 subroutine init_ens_write()
  use ensemble_variables
  use variables
  use operate
  use sync
  use density
  use writeout
  use parallel
  
  implicit none
  
  integer :: i
 
  if (ens_name == proc_master) then
 
   do i=1,nz
    call recv_unsynced_ens(a_h_h(i))
    call init_var_write(a_h_h(i)%z%z,a_h_h(i)%out)
    
    call recv_unsynced_ens(a_h_u(i))
    call init_var_write(a_h_u(i)%z%z,a_h_u(i)%out)
    
    call recv_unsynced_ens(a_h_v(i))
    call init_var_write(a_h_v(i)%z%z,a_h_v(i)%out)
    
    call recv_unsynced_ens(a_h_z(i))
    call init_var_write(a_h_z(i)%z%z,a_h_z(i)%out)
    
    call recv_synced_ens(b_u_u(i))
    b_u_u(i)%z%z%z=b_u_u(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call start_sync(b_u_u(i)%z)
    aq_huu_h(i)%z%z%z=a_h_h(i)%z%z%z*Ax(b_u_u(i)%z%z%z**2)
    call start_sync(aq_huu_h(i)%z)
    call init_var_write(b_u_u(i)%z%z,b_u_u(i)%out)
    
    call recv_synced_ens(b_v_v(i))
    b_v_v(i)%z%z%z=b_v_v(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call start_sync(b_v_v(i)%z)
    aq_hvv_h(i)%z%z%z=a_h_h(i)%z%z%z*Ay(b_v_v(i)%z%z%z**2)
    call start_sync(aq_hvv_h(i)%z)
    call init_var_write(b_v_v(i)%z%z,b_v_v(i)%out)
    
    call recv_unsynced_ens(b_uu_h(i))
    b_uu_h(i)%z%z%z=b_uu_h(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_h(i)%z%z%z),(a_h_h(i)%z%z%z == 0.0d0))
    call init_var_write(b_uu_h(i)%z%z,b_uu_h(i)%out)
    
    call recv_unsynced_ens(b_vv_h(i))
    b_vv_h(i)%z%z%z=b_vv_h(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_h(i)%z%z%z),(a_h_h(i)%z%z%z == 0.0d0))
    call init_var_write(b_vv_h(i)%z%z,b_vv_h(i)%out)
    
    call recv_unsynced_ens(b_uv_z(i))
    b_uv_z(i)%z%z%z=b_uv_z(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_z(i)%z%z%z),(a_h_z(i)%z%z%z == 0.0d0))
    call init_var_write(b_uv_z(i)%z%z,b_uv_z(i)%out)
    
    call recv_unsynced_ens(b_uv_h(i))
    b_uv_h(i)%z%z%z=b_uv_h(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_h(i)%z%z%z),(a_h_h(i)%z%z%z == 0.0d0))
    call init_var_write(b_uv_h(i)%z%z,b_uv_h(i)%out)
    
    call recv_unsynced_ens(b_lapu_u(i))
    b_lapu_u(i)%z%z%z=b_lapu_u(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(b_lapu_u(i)%z%z,b_lapu_u(i)%out)
    
    call recv_unsynced_ens(b_lapv_v(i))
    b_lapv_v(i)%z%z%z=b_lapv_v(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(b_lapv_v(i)%z%z,b_lapv_v(i)%out)
    
    call recv_unsynced_ens(b_lap2u_u(i))
    b_lap2u_u(i)%z%z%z=b_lap2u_u(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(b_lap2u_u(i)%z%z,b_lap2u_u(i)%out)
    
    call recv_unsynced_ens(b_lap2v_v(i))
    b_lap2v_v(i)%z%z%z=b_lap2v_v(i)%z%z%z*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(b_lap2v_v(i)%z%z,b_lap2v_v(i)%out)
    
   end do
   do i=1,nz
    
    bp_uu_h(i)%z%z%z=b_uu_h(i)%z%z%z-Ax(b_u_u(i)%z%z%z**2)
    h_bp_uu_h(i)%z%z%z=a_h_h(i)%z%z%z*bp_uu_h(i)%z%z%z
    call start_sync(h_bp_uu_h(i)%z)
    call init_var_write(bp_uu_h(i)%z%z,bp_uu_h(i)%out)
    
    bp_vv_h(i)%z%z%z=b_vv_h(i)%z%z%z-Ay(b_v_v(i)%z%z%z**2)
    h_bp_vv_h(i)%z%z%z=a_h_h(i)%z%z%z*bp_vv_h(i)%z%z%z
    call start_sync(h_bp_vv_h(i)%z)
    call init_var_write(bp_vv_h(i)%z%z,bp_vv_h(i)%out)
    
    aq_huv_z(i)%z%z%z=a_h_z(i)%z%z%z*sAy(b_u_u(i)%z)*sAx(b_v_v(i)%z)
    bp_uv_z(i)%z%z%z=b_uv_z(i)%z%z%z-sAy(b_u_u(i)%z)*sAx(b_v_v(i)%z)
    h_bp_uv_z(i)%z%z%z=a_h_z(i)%z%z%z*bp_uv_z(i)%z%z%z
    call init_var_write(bp_uv_z(i)%z%z,bp_uv_z(i)%out)
    
    bp_uv_h(i)%z%z%z=b_uv_h(i)%z%z%z-Ax(b_u_u(i)%z%z%z)*Ay(b_v_v(i)%z%z%z)
    call init_var_write(bp_uv_h(i)%z%z,bp_uv_h(i)%out)
    
   end do
   
   
   do i=0,nz
    
    call recv_synced_ens(a_d_h(i))
    a_d_h(i)%z%z%z=a_d_h(i)%z%z%z
    call start_sync(a_d_h(i)%z)
    aq_dd_h(i)%z%z%z=a_d_h(i)%z%z%z**2
    call start_sync(aq_dd_h(i)%z)
    call init_var_write(a_d_h(i)%z%z,a_d_h(i)%out)
    
    call recv_unsynced_ens(a_p_h(i))
    call init_var_write(a_p_h(i)%z%z,a_p_h(i)%out)
    
    call recv_unsynced_ens(a_dd_h(i))
    call init_var_write(a_dd_h(i)%z%z,a_dd_h(i)%out)
    
    call recv_unsynced_ens(a_dpx_u(i))
    call init_var_write(a_dpx_u(i)%z%z,a_dpx_u(i)%out)
    
    call recv_unsynced_ens(a_dpy_v(i))
    call init_var_write(a_dpy_v(i)%z%z,a_dpy_v(i)%out)
    
    call recv_unsynced_ens(a_px_u(i))
    call init_var_write(a_px_u(i)%z%z,a_px_u(i)%out)
    
    call recv_unsynced_ens(a_py_v(i))
    call init_var_write(a_py_v(i)%z%z,a_py_v(i)%out)
    
   end do
   
   do i=0,nz
   
   
    ap_dd_h(i)%z%z%z=a_dd_h(i)%z%z%z-a_d_h(i)%z%z%z**2
    call start_sync(ap_dd_h(i)%z)
    call init_var_write(ap_dd_h(i)%z%z,ap_dd_h(i)%out)
        
    ap_dpx_u(i)%z%z%z=a_dpx_u(i)%z%z%z-sAx(a_d_h(i)%z)*a_px_u(i)%z%z%z
    call init_var_write(ap_dpx_u(i)%z%z,ap_dpx_u(i)%out)
    
    ap_dpy_v(i)%z%z%z=a_dpy_v(i)%z%z%z-sAy(a_d_h(i)%z)*a_py_v(i)%z%z%z
    call init_var_write(ap_dpy_v(i)%z%z,ap_dpy_v(i)%out)
    
   end do
   
   do i=1,nz
    tend_u1(i)%z%z%z=-sGx(aq_huu_h(i)%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u1(i)%z%z,tend_u1(i)%out)
       
    tend_u2(i)%z%z%z=-Gy(aq_huv_z(i)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u2(i)%z%z,tend_u2(i)%out)
       
    tend_u3(i)%z%z%z=Ay(fz%z%z)*sAxy(b_v_v(i)%z)
    call init_var_write(tend_u3(i)%z%z,tend_u3(i)%out)
    
    tend_u4(i)%z%z%z=-0.5d0*sum(gp(i:nz))*(sGx(aq_dd_h(i)%z)-sGx(aq_dd_h(i-1)%z))*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u4(i)%z%z,tend_u4(i)%out)
       
    tend_u5(i)%z%z%z=-(sAx(a_d_h(i)%z)*a_px_u(i)%z%z%z-  &
       sAx(a_d_h(i-1)%z)*a_px_u(i-1)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u5(i)%z%z,tend_u5(i)%out)
       
    if (i == nz) then
     tend_u6(i)%z%z%z=utau%z%z*merge(0.0d0,1.0d0/a_h_u(i)%z%z%z,(a_h_u(i)%z%z%z == 0.0d0))
    else
     tend_u6(i)%z%z%z=0.0d0
    end if
    call init_var_write(tend_u6(i)%z%z,tend_u6(i)%out)
    
    tend_u7(i)%z%z%z=-sGx(h_bp_uu_h(i)%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u7(i)%z%z,tend_u7(i)%out)
       
    tend_u8(i)%z%z%z=-Gy(h_bp_uv_z(i)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u8(i)%z%z,tend_u8(i)%out)
       
    tend_u9(i)%z%z%z=-0.5d0*sum(gp(i:nz))*(sGx(ap_dd_h(i)%z)-sGx(ap_dd_h(i-1)%z))*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u9(i)%z%z,tend_u9(i)%out)
       
    tend_u10(i)%z%z%z=-(ap_dpx_u(i)%z%z%z-ap_dpx_u(i-1)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_u(i)%z%z%z),(a_h_u(i)%z%z%z == 0.0d0))
    call init_var_write(tend_u10(i)%z%z,tend_u10(i)%out)
    
    
       
    tend_v1(i)%z%z%z=-Gx(aq_huv_z(i)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v1(i)%z%z,tend_v1(i)%out)
    
    tend_v2(i)%z%z%z=-sGy(aq_hvv_h(i)%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v2(i)%z%z,tend_v2(i)%out)
       
    tend_v3(i)%z%z%z=-Ax(fz%z%z)*sAyx(b_u_u(i)%z)
    call init_var_write(tend_v3(i)%z%z,tend_v3(i)%out)
    
    tend_v4(i)%z%z%z=-0.5d0*sum(gp(i:nz))*(sGy(aq_dd_h(i)%z)-sGy(aq_dd_h(i-1)%z))  &
           *merge(0.0d0,1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v4(i)%z%z,tend_v4(i)%out)
       
    tend_v5(i)%z%z%z=-(sAy(a_d_h(i)%z)*a_py_v(i)%z%z%z-&
          sAy(a_d_h(i-1)%z)*a_py_v(i-1)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v5(i)%z%z,tend_v5(i)%out)
       
    if (i == nz) then
     tend_v6(i)%z%z%z=vtau%z%z*merge(0.0d0,1.0d0/a_h_v(i)%z%z%z,(a_h_v(i)%z%z%z == 0.0d0))
    else
     tend_v6(i)%z%z%z=0.0d0
    end if
    call init_var_write(tend_v6(i)%z%z,tend_v6(i)%out)
    
    tend_v7(i)%z%z%z=-sGy(h_bp_vv_h(i)%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v7(i)%z%z,tend_v7(i)%out)
       
    tend_v8(i)%z%z%z=-Gx(h_bp_uv_z(i)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v8(i)%z%z,tend_v8(i)%out)
       
    tend_v9(i)%z%z%z=-0.5d0*sum(gp(i:nz))*(sGy(ap_dd_h(i)%z)-sGy(ap_dd_h(i-1)%z))*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v9(i)%z%z,tend_v9(i)%out)
       
    tend_v10(i)%z%z%z=-(ap_dpy_v(i)%z%z%z-ap_dpy_v(i-1)%z%z%z)*merge(0.0d0,    &
       1.0d0/(a_h_v(i)%z%z%z),(a_h_v(i)%z%z%z == 0.0d0))
    call init_var_write(tend_v10(i)%z%z,tend_v10(i)%out)
    
   end do
   
  end if
   
   
  end subroutine
 
  
 
 subroutine write_ensemble(n)
  use ensemble_variables
  use variables
  use time
  use writeout
  use sync
  use operate
  use parallel
  use density
  
  implicit none
  
  integer, intent(in) :: n
  character(32) :: filename
  integer :: i, num
  logical :: dir_e
  
  if (ens_name == proc_master) then
   num = (n-1)/wstep-1
   if (proc_name == ens_master) then
    write (filename, "(a12)") './data/ens/.'
    inquire( file=filename, exist=dir_e )
    do while (.not.(dir_e))   
     write (filename, "(a17)") 'mkdir -p data/ens'
     call system(filename,stat)
     write (filename, "(a12)") './data/ens/.'
     inquire( file=filename, exist=dir_e )
    end do
    write (filename, "(a11,i4.4,a2)") './data/ens/', num, '/.'
    inquire( file=filename, exist=dir_e )
    do while (.not.(dir_e))   
     write (filename, "(a18,i4.4)") 'mkdir -p data/ens/', num
     call system(filename)
     write (filename, "(a11,i4.4,a2)") './data/ens/', num, '/.'
     inquire( file=filename, exist=dir_e )
    end do
   
    do i=1,nz
    
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_h_h_', i, '.dat'
     call end_var_write(a_h_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_h_u_', i, '.dat'
     call end_var_write(a_h_u(i)%out,filename,0)
   
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_h_v_', i, '.dat'
     call end_var_write(a_h_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_h_z_', i, '.dat'
     call end_var_write(a_h_z(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/b_u_u_', i, '.dat'
     call end_var_write(b_u_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/b_v_v_', i, '.dat'
     call end_var_write(b_v_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/b_uu_h_', i, '.dat'
     call end_var_write(b_uu_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/b_vv_h_', i, '.dat'
     call end_var_write(b_vv_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/b_uv_z_', i, '.dat'
     call end_var_write(b_uv_z(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/b_uv_h_', i, '.dat'
     call end_var_write(b_uv_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/bp_uu_h_', i, '.dat'
     call end_var_write(bp_uu_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/bp_vv_h_', i, '.dat'
     call end_var_write(bp_vv_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/bp_uv_z_', i, '.dat'
     call end_var_write(bp_uv_z(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/bp_uv_h_', i, '.dat'
     call end_var_write(bp_uv_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/b_lapu_u_', i, '.dat'
     call end_var_write(b_lapu_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/b_lapv_v_', i, '.dat'
     call end_var_write(b_lapv_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a11,i1.1,a4)") 'data/ens/', num, '/b_lap2u_u_', i, '.dat'
     call end_var_write(b_lap2u_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a11,i1.1,a4)") 'data/ens/', num, '/b_lap2v_v_', i, '.dat'
     call end_var_write(b_lap2v_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u1_', i, '.dat'
     call end_var_write(tend_u1(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u2_', i, '.dat'
     call end_var_write(tend_u2(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u3_', i, '.dat'
     call end_var_write(tend_u3(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u4_', i, '.dat'
     call end_var_write(tend_u4(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u5_', i, '.dat'
     call end_var_write(tend_u5(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u6_', i, '.dat'
     call end_var_write(tend_u6(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u7_', i, '.dat'
     call end_var_write(tend_u7(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u8_', i, '.dat'
     call end_var_write(tend_u8(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_u9_', i, '.dat'
     call end_var_write(tend_u9(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/tend_u10_', i, '.dat'
     call end_var_write(tend_u10(i)%out,filename,0)
     
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v1_', i, '.dat'
     call end_var_write(tend_v1(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v2_', i, '.dat'
     call end_var_write(tend_v2(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v3_', i, '.dat'
     call end_var_write(tend_v3(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v4_', i, '.dat'
     call end_var_write(tend_v4(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v5_', i, '.dat'
     call end_var_write(tend_v5(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v6_', i, '.dat'
     call end_var_write(tend_v6(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v7_', i, '.dat'
     call end_var_write(tend_v7(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v8_', i, '.dat'
     call end_var_write(tend_v8(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/tend_v9_', i, '.dat'
     call end_var_write(tend_v9(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/tend_v10_', i, '.dat'
     call end_var_write(tend_v10(i)%out,filename,0)
    end do
       
    
    do i=0,nz
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/a_dd_h_', i, '.dat'
     call end_var_write(a_dd_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/a_dpx_u_', i, '.dat'
     call end_var_write(a_dpx_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/a_dpy_v_', i, '.dat'
     call end_var_write(a_dpy_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_d_h_', i, '.dat'
     call end_var_write(a_d_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a7,i1.1,a4)") 'data/ens/', num, '/a_p_h_', i, '.dat'
     call end_var_write(a_p_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/a_px_u_', i, '.dat'
     call end_var_write(a_px_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a8,i1.1,a4)") 'data/ens/', num, '/a_py_v_', i, '.dat'
     call end_var_write(a_py_v(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a9,i1.1,a4)") 'data/ens/', num, '/ap_dd_h_', i, '.dat'
     call end_var_write(ap_dd_h(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/ap_dpx_u_', i, '.dat'
     call end_var_write(ap_dpx_u(i)%out,filename,0)
     
     write (filename, "(a9,i4.4,a10,i1.1,a4)") 'data/ens/', num, '/ap_dpy_v_', i, '.dat'
     call end_var_write(ap_dpy_v(i)%out,filename,0)
     
    end do
        
   end if
   
  end if
  
  
 end subroutine
 
end module
