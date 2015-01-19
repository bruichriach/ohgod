

MODULE writeout_variables
 use writeout_grid
 use params
 use global

 IMPLICIT NONE
 
 type send_dat
  real(kind=8), allocatable :: dat(:)
 end type send_dat
  
 
 type out_var
  integer :: tag
  integer :: req
  integer :: nx, ny
  real(kind=8), pointer :: dat(:,:)
  type(point), pointer :: p(:,:)
  integer, allocatable :: req_master(:)
  real(kind=8), allocatable :: send(:)
  type(send_dat), allocatable :: recv(:)
  real(kind=8), allocatable :: out(:,:)
  integer, pointer :: grid(:,:)
 end type out_var
  
 
 
 type(out_var), save :: h_out(nz), u_out(nz), v_out(nz), zeta_out(nz)
 type(out_var), save :: q_h_out(nz), q_z_out(nz)
 type(out_var), save :: ke_out(nz), ape_out(nz), depth_out(0:nz), mont_out(nz)
 type(out_var), save :: pres_out
 type(out_var), save :: div_u_out
 type(out_var), save :: utau_out, vtau_out
 type(out_var), save :: tendh_out(nz), tendu_out(nz), tendv_out(nz)
 type(out_var), save :: smagu_out(nz), smagv_out(nz)

end module
 
MODULE writeout

 implicit none
 
 CONTAINS
 
 subroutine write_config(dat,grid,varname)
  use writeout_grid
  use writeout_variables
  use variables
  use parallel
  
  implicit none
 
  type(var), intent(in) :: dat
  integer, intent(in) :: grid(:,:)
  character(*), intent(in) :: varname
  type(out_var) :: out
  character(32) :: filename
  
  
   call mpi_barrier(mpi_comm_world,stat)
   call init_writeouts(dat%z,dat%p,out,grid)
   call init_var_write(out)
   write (filename, "(a14,i4.4,a4)") 'global/'//trim(varname)//'_', ens_name, '.dat'
   call end_var_write(out,adjustl(filename),0)
   call mpi_barrier(mpi_comm_world,stat)
  
 end subroutine 
 
 subroutine read_var(dat,layer,varname,folder,x,y)
  use writeout_grid
  use writeout_variables
  use parallel
  use variables
  use grid
  
  implicit none
  
  type(var), intent(inout) :: dat
  character(*), intent(in) :: varname, folder
  type(out_var) :: out
  character(32) :: filename, format, foldername, fullname
  integer, intent(in) :: layer,x,y
  integer :: i,j,k,l
  logical :: file_exist
  
   
  
  
  write (foldername, "(a28,i4.4)") trim(folder)//'/', ens_name
  write (filename, "(a27,i1.1,a4)") trim(varname)//'_',layer,'.dat'
  write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
  inquire( file=fullname, exist=file_exist)
  
  if (file_exist) then
  
  call init_readins(out,1,x,y)
   
  
  

  if (proc_name == ens_master) then
   open(unit=10,form='unformatted',file=adjustl(fullname),status='unknown')
   write (format,"(a2,i3,a8)") '( ', mx+x, 'd24.16 )'
   do j=1,my+y
    read(10) (out%out(i,j),i=1,mx+x)
   end do
   close(10)
    
  
  
   do k=0,ens_images-1 
    l=0
    do j=1,my+y
     do i=1,mx+x
      if ((init_grid(i,j) == ens_master+k).or.   &
         (init_grid(i-x,j) == ens_master+k).or.   &
         (init_grid(i,j-y) == ens_master+k).or.   &
         (init_grid(i-x,j-y) == ens_master+k)) then
       l=l+1
       out%recv(k)%dat(l)=out%out(i,j)
      end if
     end do
    end do
    
    call mpi_isend(out%recv(k)%dat(1),ubound(out%recv(k)%dat),mpi_double_precision,  &
       ens_master+k,out%tag,  &
       mpi_comm_world,out%req_master(k),stat)
   end do
  end if
  call mpi_irecv(out%send(1),ubound(out%send),mpi_double_precision,  &
       ens_master,out%tag,  &
       mpi_comm_world,out%req,stat)
       
  
  k=0
  call mpi_wait(out%req,mpi_status_ignore,stat)
  do j=1,ubound(dat%z,2)
   do i=1,ubound(dat%z,1)
    if (dat%p(i,j)%core) then
      k=k+1
      dat%z(i,j)=out%send(k)
    end if
   end do
  end do
  
  if (proc_name == ens_master) then
   call mpi_waitall(ens_images, out%req_master,mpi_statuses_ignore,stat)
   print *, 'read in complete for', filename
  end if
  
  
  end if
  

  
 end subroutine
 
 subroutine write_main(folder)
  use writeout_grid
  use writeout_variables
  use parallel
  use variables
  use solver_variables
 
  implicit none
 
  character(*), intent(in) :: folder
  character(32) :: filename
  integer :: k
  logical :: dir_e
  
  if (proc_name == ens_master) then
   write (filename, "(a1,i4.4)") '/', ens_name
   inquire( file=folder//filename, exist=dir_e )
   do while (.not.(dir_e))   
    call system('mkdir -p '//folder//filename,stat)
    inquire( file=folder//filename, exist=dir_e )
   end do
  end if
  
  do k=1,nz
   call init_var_write(h_out(k))
   call init_var_write(u_out(k))
   call init_var_write(v_out(k))
  end do
  call init_var_write(pres_out)
  do k=1,nz
   write (filename, "(a1,i4.4,a3,i1.1,a4)") '/', ens_name, '/h_', k, '.dat'
   call end_var_write(h_out(k),folder//filename,0)
   write (filename, "(a1,i4.4,a3,i1.1,a4)") '/', ens_name, '/u_', k, '.dat'
   call end_var_write(u_out(k),folder//filename,0)
   write (filename, "(a1,i4.4,a3,i1.1,a4)") '/', ens_name, '/v_', k, '.dat'
   call end_var_write(v_out(k),folder//filename,0)
  end do
   write (filename, "(a1,i4.4,a6,i1.1,a4)") '/', ens_name, '/pres_', 0, '.dat'
  call end_var_write(pres_out,folder//filename,0)
  
 end subroutine 
  
 subroutine init_writeouts(dat,p,out,grid)
  use writeout_grid
  use writeout_variables
  use parallel
  use global
  
  implicit none
  
 
  real(kind=8), target, intent(in) :: dat(:,:)
  type(point), target, intent(in) :: p(:,:)
  integer, target, intent(in) :: grid(:,:)
  type(out_var), intent(inout) :: out
  integer :: k
  
  max_write_tag=max_write_tag+1
  
  out%nx = ubound(dat,1)
  out%ny = ubound(dat,2)
  
  out%dat => dat
  out%p => p
  
  out%grid => grid
  allocate(out%out(ubound(grid,1),ubound(grid,2)))
  out%out(1:ubound(grid,1),1:ubound(grid,2))=0
  allocate(out%send(count((grid == proc_name))))
  out%req=mpi_request_null
  if (proc_name == ens_master) then
   allocate(out%req_master(0:ens_images-1))
   out%req_master(0:ens_images-1)=mpi_request_null
   allocate(out%recv(0:ens_images-1))
   do k=0,ens_images-1
    allocate(out%recv(k)%dat(count((grid == ens_master+k))))
   end do
  end if
  out%tag=10000+max_write_tag
   
  
  
 end subroutine
 
 
  
 subroutine init_readins(out,tag,i,j)
  use writeout_grid
  use writeout_variables
  use parallel
  use grid
  
  implicit none
 
  integer, intent(in) :: tag
  type(out_var), intent(inout) :: out
  integer :: k
  integer, intent(in) :: i,j
  
  
  allocate(out%out(mx+i,my+j))
  out%out=0
  allocate(out%send(count((init_grid(1:mx+i,1:my+j) == proc_name).or.   &
                          (init_grid(1-i:mx,1:my+j) == proc_name).or.   &
                          (init_grid(1:mx+i,1-j:my) == proc_name).or.   &
                          (init_grid(1-i:mx,1-j:my) == proc_name))))
  out%req=mpi_request_null
  if (proc_name == ens_master) then
   allocate(out%req_master(0:ens_images-1))
   out%req_master(0:ens_images-1)=mpi_request_null
   allocate(out%recv(0:ens_images-1))
   do k=0,ens_images-1
    allocate(out%recv(k)%dat(count((init_grid(1:mx+i,1:my+j) == ens_master + k).or.   &
                          (init_grid(1-i:mx,1:my+j) == ens_master + k).or.   &
                          (init_grid(1:mx+i,1-j:my) == ens_master + k).or.   &
                          (init_grid(1-i:mx,1-j:my) == ens_master + k))))
   end do
  end if
  out%tag=10000+tag
   
  
  
 end subroutine
 
 
 subroutine init_var_write(out)
  use writeout_grid
  use writeout_variables
  use variables
  use parallel

  IMPLICIT NONE

  type(out_var), intent(inout) :: out
  integer :: i,j,k
  

  k=0
  call mpi_wait(out%req,mpi_status_ignore,stat)
  do j=1,out%ny
   do i=1,out%nx
    if (out%p(i,j)%core) then
     if (out%grid(out%p(i,j)%i,out%p(i,j)%j) ==  proc_name) then
      k=k+1
      out%send(k)=out%dat(i,j)
     end if
    end if
   end do
  end do
  
  if (proc_name == ens_master) then
   do k=0,ens_images-1
    call mpi_irecv(out%recv(k)%dat(1),ubound(out%recv(k)%dat),mpi_double_precision,  &
       ens_master+k,out%tag,  &
       mpi_comm_world,out%req_master(k),stat)
   end do
  end if
  call mpi_isend(out%send(1),ubound(out%send),mpi_double_precision,  &
       ens_master,out%tag,  &
       mpi_comm_world,out%req,stat)
  
 end subroutine
 
 subroutine end_var_write(out,filename,fmt)
  use writeout_grid
  use writeout_variables
  use parallel
 
  implicit none
  
  character(*), intent(in) :: filename
  type(out_var), intent(inout) :: out
  integer :: i,j,k,l
  integer, intent(in) :: fmt
  
  if (proc_name == ens_master) then
   do k=0,ens_images-1 
    call mpi_wait(out%req_master(k),mpi_status_ignore,stat)
    l=0
    do j=1,ubound(out%out,2)
     do i=1,ubound(out%out,1)
      if (out%grid(i,j) == ens_master+k) then
       l=l+1
       out%out(i,j)=out%recv(k)%dat(l)
      end if
     end do
    end do
   end do
   call datawrite(fmt,filename,out%out)
  end if
  
 end subroutine
 
 subroutine datawrite(type,filename, variable)

  IMPLICIT NONE
 
  integer i, j
  integer :: type
  character(*), intent(in) :: filename
  real (kind=8), intent(in), dimension(:,:) :: variable
  character(32) :: format
  
  write (format,"(a2,i3,a8)") '( ', ubound(variable,1), 'd24.16 )'

  if (type == 1) then
   open(unit=10,file=trim(filename),status='unknown')
     do j=1,ubound(variable,2)
      write(10,trim(format)) (variable(i,j),i=1,ubound(variable,1))
     end do
   close(10)
  else
   open(unit=10,file=trim(filename),form='unformatted',status='unknown')
     do j=1,ubound(variable,2)
      write(10) (variable(i,j),i=1,ubound(variable,1))
     end do
   close(10)
  end if

  return

 end subroutine
 
 
 
 
 subroutine write()
  use writeout_grid
  use writeout_variables
  use global
  use solver_variables
  use variables
  use parallel
  use operate
 
  implicit none
   
  integer :: i
  
  
  
  do i=1,nz
   q_h(i)%z=Ax(Ay(fz%z+zeta(i)%z))*merge(0.0d0,1.0d0/h(i)%z,(h(i)%z == 0.0d0))
   q_z(i)%z=(fz%z+zeta(i)%z)*merge(0.0d0,1.0d0/h_z(i)%z,(h_z(i)%z == 0.0d0))
  end do
  
  call init_var_write(utau_out)
  call init_var_write(vtau_out)
  call init_var_write(pres_out)
  call init_var_write(div_u_out)
  call init_var_write(depth_out(0))
  do i=1,nz
   call init_var_write(depth_out(i))
   call init_var_write(h_out(i))
   call init_var_write(u_out(i))
   call init_var_write(v_out(i))
   call init_var_write(tendh_out(i))
   call init_var_write(tendu_out(i))
   call init_var_write(tendv_out(i))
   call init_var_write(smagu_out(i))
   call init_var_write(smagv_out(i))
   call init_var_write(ke_out(i))
   call init_var_write(ape_out(i))
   call init_var_write(q_h_out(i))
   call init_var_write(q_z_out(i))
   call init_var_write(zeta_out(i))
  end do
  
  
 
 
 end subroutine
 
 subroutine do_write(num)
  use writeout_variables
  use parallel
 
  implicit none
  
  integer, intent(in) :: num
  character(32) :: filename
  integer :: i
  logical :: dir_e
  
  
  if (proc_name == ens_master) then
   write (filename, "(a7,i4.4,a2)") './data/', ens_name, '/.'
   inquire( file=filename, exist=dir_e )
   do while (.not.(dir_e))   
    write (filename, "(a14,i4.4)") 'mkdir -p data/', ens_name
    call system(filename,stat)
    write (filename, "(a7,i4.4,a2)") './data/', ens_name, '/.'
    inquire( file=filename, exist=dir_e )
   end do
   write (filename, "(a7,i4.4,a1,i4.4,a2)") './data/', ens_name, '/', num, '/.'
   inquire( file=filename, exist=dir_e )
   do while (.not.(dir_e))   
    write (filename, "(a14,i4.4,a1,i4.4)") 'mkdir -p data/', ens_name, '/', num
    call system(filename)
    write (filename, "(a7,i4.4,a1,i4.4,a2)") './data/', ens_name, '/', num, '/.'
    inquire( file=filename, exist=dir_e )
   end do
   
   write (filename, "(a5,i4.4,a1,i4.4,a6,i1.1,a4)") 'data/', ens_name, '/', num, '/utau_', 0, '.dat'
   call end_var_write(utau_out,filename,0)
   
   write (filename, "(a5,i4.4,a1,i4.4,a6,i1.1,a4)") 'data/', ens_name, '/', num, '/vtau_', 0, '.dat'
   call end_var_write(vtau_out,filename,0)
   
   write (filename, "(a5,i4.4,a1,i4.4,a6,i1.1,a4)") 'data/', ens_name, '/', num, '/pres_', 0, '.dat'
   call end_var_write(pres_out,filename,0)

   write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/div_u_', 0, '.dat'
   call end_var_write(div_u_out,filename,0)

   write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/depth_', 0, '.dat'
   call end_var_write(depth_out(0),filename,0)
   
   do i=1,nz
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/depth_', i, '.dat'
    call end_var_write(depth_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a3,i1.1,a4)") 'data/', ens_name, '/', num, '/h_', i, '.dat'
    call end_var_write(h_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a3,i1.1,a4)") 'data/', ens_name, '/', num, '/u_', i, '.dat'
    call end_var_write(u_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a3,i1.1,a4)") 'data/', ens_name, '/', num, '/v_', i, '.dat'
    call end_var_write(v_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/tendh_', i, '.dat'
    call end_var_write(tendh_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/tendu_', i, '.dat'
    call end_var_write(tendu_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/tendv_', i, '.dat'
    call end_var_write(tendv_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/smagu_', i, '.dat'
    call end_var_write(smagu_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/smagv_', i, '.dat'
    call end_var_write(smagv_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/ke_', i, '.dat'
    call end_var_write(ke_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/ape_', i, '.dat'
    call end_var_write(ape_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/q_h_', i, '.dat'
    call end_var_write(q_h_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/q_z_', i, '.dat'
    call end_var_write(q_z_out(i),filename,0)
   
    write (filename, "(a5,i4.4,a1,i4.4,a6,i1.1,a4)") 'data/', ens_name, '/', num, '/zeta_', i, '.dat'
    call end_var_write(zeta_out(i),filename,0)
   
   end do
  end if
  
 end subroutine

END MODULE
  
