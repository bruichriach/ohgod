module writeout_grid
 use params

 implicit none

 integer, save :: hgrid_out(mx,my),ugrid_out(mx+1,my),vgrid_out(mx,my+1),zgrid_out(mx+1,my+1)
 
end module


module grid
 use params

 implicit none
 
 integer, dimension(0:mx+1,0:my+1) :: init_grid
 
 contains
 
 
 subroutine set_write_grid()
  use writeout_grid
  use parallel
  
  implicit none
  
  integer i, j, k
  
  hgrid_out=-1
  ugrid_out=-1
  vgrid_out=-1
  zgrid_out=-1
  do k=0,proc_num-1
   do j=1,my
    do i=1,mx
     if (hgrid_out(i,j) == -1) then
      if (init_grid(i,j) == k) then
       hgrid_out(i,j)=k
      end if
     end if
    end do
   end do
   do j=1,my
    do i=1,mx+1
     if (ugrid_out(i,j) == -1) then
      if ((init_grid(i,j) == k).or.(init_grid(i-1,j) == k)) then
       ugrid_out(i,j)=k
      end if
     end if
    end do
   end do
   do j=1,my+1
    do i=1,mx
     if (vgrid_out(i,j) == -1) then
      if ((init_grid(i,j) == k).or.(init_grid(i,j-1) == k)) then
       vgrid_out(i,j)=k
      end if
     end if
    end do
   end do
   do j=1,my+1
    do i=1,mx+1
     if (zgrid_out(i,j) == -1) then
      if ((init_grid(i,j) == k).or.(init_grid(i-1,j) == k).or. &
          (init_grid(i,j-1) == k).or.(init_grid(i-1,j-1) == k)) then
       zgrid_out(i,j)=k
      end if
     end if
    end do
   end do
  end do
  
  
 end subroutine
 
 subroutine point_at(a,b)
  use global
  
  implicit none
  
  real(kind=8), target, intent(in) :: a
  real(kind=8), pointer, intent(out) :: b
  
  b => a
  
 end subroutine
 
 subroutine create_grid()
  use optimise
  use parallel
  use system
  
  implicit none
  
  integer, dimension(0:mx+1,0:my+1) :: node_grid
  integer, dimension(0:mx+1,0:my+1) :: core_grid
  integer :: grid_core, grid_node
  integer :: core, node
  integer :: i,j,k,l,m
  logical :: old_grid
  
  
  inquire( file='init_grid.dat', exist=old_grid )
  if (old_grid) then
   open(unit=10,file='init_grid.dat',form='unformatted',status='unknown')
     do j=0,my+1
      read(10) (init_grid(i,j),i=0,mx+1)
     end do
   close(10)
   if (maxval(init_grid) /= ens_images-1) old_grid = .false.
  end if
  
  if (.not.(old_grid)) then
  
  
  if (max_core > proc_num) then
   core = proc_num
  else
   core = max_core
  end if
 
  node = proc_num/core
  
  if (ens_num <= node) then
   grid_node = node/ens_num
   grid_core = core
  else
   grid_node = 1
   grid_core = ens_images
  end if
  
  ens_master=ens_images*(proc_name/ens_images)
  
  
  if (grid_core*grid_node*ens_num /= proc_num) then
   print *, grid_core, grid_node, ens_num, "ensemble split error"
   stop
  end if
  
  node_grid=-1
  node_grid(1:mx,1:my)=0
!  do j=1,my
!   do i=1,mx
!    if ((abs(dble(i)-0.5d0-dble(mx)/2.0d0) <= dble(mx)/2.0d0 - 8.0d0).or.   &
!         (abs(dble(j)-0.5d0-dble(my)/2.0d0) <= 24.0d0)) then
!     node_grid(i,j)=0
!    end if
!   end do
!  end do
!  do j=1,my
!   do i=1,mx
!    if (((dble(i)-0.5d0)-dble(mx)/2.0d0)**2+((dble(j)-0.5d0)-dble(my)/2.0d0)**2  &
!         <= (dble(mx)/2.0d0)**2) then
!     node_grid(i,j)=0
!    end if
!   end do
!  end do

   call init_guess(node_grid, grid_node, node_grid )

!  l = maxval(count((node_grid(1:mx,1:my) /= -1),2))
!  m=0
!  do j=1,my
!   if (m /= l) then
!    do i=1,mx
!     if (node_grid(i,j) /= -1) then
!     node_grid(i,j) = (grid_node*m)/l
!     end if
!    end do
!    m=m+1
!   end if
!  end do
  
  call run_optimise(node_grid)
  
  call mpi_barrier(mpi_comm_world,stat)
  if (proc_name == proc_master) then
   print *, 'node optimisation complete    time:', cputime()
  end if
  call mpi_barrier(mpi_comm_world,stat)
  
  
  init_grid=-1
  
  do k=0,grid_node-1
   core_grid=-1
   
   call init_guess(merge(0,-1,node_grid == k), grid_core , core_grid )
 
   
!   l = count((count((node_grid(1:mx,1:my) == k),1) /= 0))
!   m=0
!   do j=1,my
!    if (m /= l) then
!     if (any(node_grid(:,j) == k)) then
!     do i=1,mx
!      if (node_grid(i,j) == k) then
!      core_grid(i,j) = (grid_core*m)/l
!      end if
!     end do
!     m=m+1
!     end if
!    end if
!   end do
   
   call run_optimise(core_grid)
   
   
   do j=1,my
    do i=1,mx
     if (node_grid(i,j) == k) then
      init_grid(i,j) = k*grid_core + core_grid(i,j)
     end if
    end do
   end do
  end do
  
  end if
  
  init_grid=merge(-1,ens_master+init_grid,(init_grid == -1))
  
  call mpi_barrier(mpi_comm_world,stat)
  
  
  nx = count((count((init_grid(1:mx,1:my) == proc_name),2) /= 0))
  ny = count((count((init_grid(1:mx,1:my) == proc_name),1) /= 0))
  

  
  if (proc_name == proc_master) then
   open(unit=10,file='init_grid.dat',form='unformatted',status='unknown')
     do j=0,my+1
      write(10) (init_grid(i,j),i=0,mx+1)
     end do
   close(10)
  end if
  
  call mpi_barrier(mpi_comm_world,stat)
  if (proc_name == proc_master) then
   print *, 'core optimisation complete    time:', cputime()
  end if
  call mpi_barrier(mpi_comm_world,stat)
   
  
 end subroutine
 
 subroutine init_guess(mask, num, grid )
  implicit none
 
  integer, intent(in) :: mask(0:,0:), num
  integer :: nx, ny, mx, my, lx, ly, k, m, n, i, j
  integer, intent(out) :: grid(0:ubound(mask,1),0:ubound(mask,2))
  real(kind=8) :: tmp(0:ubound(mask,1),0:ubound(mask,2),0:num-1)
  real(kind=8) :: x(0:num-1),y(0:num-1)
  
  mx=ubound(mask,1)-1
  my=ubound(mask,2)-1
  
  nx=count((count((mask /= -1),2) /= 0))
  ny=count((count((mask /= -1),1) /= 0))
  
  lx=maxval(minloc(mask ,1,(mask >= 0)))-1
  ly=maxval(minloc(mask ,2,(mask >= 0)))-1
  
  
  
  k=floor(sqrt(dble(nx*ny)/dble(num)))
  
  n=ceiling(dble(ny)/dble(k))
  m=ceiling(dble(nx)/dble(k))
  if (nx == max(nx,ny)) then
   do while (n*m < num)
    m=m+1
   end do
  else
   do while (n*m < num)
    n=n+1
   end do
  end if
   
  
  do i=0,num-1
   x(i)=dble(lx)+0.5d0*dble(nx)/dble(m)+(dble(nx)/dble(m))*dble(mod(i,m))
   y(i)=dble(ly)+0.5d0*dble(ny)/dble(n)+(dble(ny)/dble(n))*dble(i/m)
  end do
  
  do j=0,my+1
   do i=0,mx+1
    if (mask(i,j) /= -1) then
     do k=0,num-1
      tmp(i,j,k)=(dble(i)-0.5d0-x(k))**2+(dble(j)-0.5d0-y(k))**2
     end do
    end if
   end do
  end do
  
  grid=merge(minloc(tmp,3)-1,mask,(mask /= -1))
  
  
 
 
  end subroutine
 
 
 
 subroutine create_var(dat,x,y,vec)
  use global
  use parallel
  
  implicit none
  
  type(var), intent(out) :: dat
  integer, intent(in) :: x, y
  integer :: i, j, k, l, m
  logical, intent(in) :: vec
  
  if (.not.((abs(x)+abs(y) == 0).or.(abs(x)+abs(y) == 1))) then
   print *, 'variable type error'
   call mpi_finalize(stat)
   stop
  end if
  
  dat%nx = count((count(((init_grid(1:mx+x,1:my+y) == proc_name).or.   &
                         (init_grid(1-x:mx,1:my+y) == proc_name).or.   &
                         (init_grid(1:mx+x,1-y:my) == proc_name)),2) /= 0))
  dat%lx = 0
  do i=1,mx+x
   if (count((count(((init_grid(1:i,1:my+y) == proc_name).or.   &
                     (init_grid(1-x:i-x,1:my+y) == proc_name).or.   &
                     (init_grid(1:i,1-y:my) == proc_name)),2) /= 0)) == 0) dat%lx = dat%lx+1
  end do
  dat%ny = count((count(((init_grid(1:mx+x,1:my+y) == proc_name).or.   &
                         (init_grid(1-x:mx,1:my+y) == proc_name).or.   &
                         (init_grid(1:mx+x,1-y:my) == proc_name)),1) /= 0))
  dat%ly = 0
  do i=1,my+y
   if (count((count(((init_grid(1:mx+x,1:i) == proc_name).or.   &
                         (init_grid(1-x:mx,1:i) == proc_name).or.   &
                         (init_grid(1:mx+x,1-y:i-y) == proc_name)),1) /= 0)) == 0) dat%ly = dat%ly+1
  end do
  
  dat%tag=max_tag
  max_tag = max_tag+1
  
  if (proc_name == proc_master) print *, 'used tag #', dat%tag
  
  allocate(dat%z(dat%nx,dat%ny))
  allocate(dat%exp(dat%nx,dat%ny))
  allocate(dat%p(dat%nx,dat%ny))
  allocate(dat%bx(0:dat%nx,dat%ny))
  allocate(dat%ux(dat%nx+1,dat%ny))
  allocate(dat%by(dat%nx,0:dat%ny))
  allocate(dat%uy(dat%nx,dat%ny+1))
  allocate(dat%mpi(0:proc_num-1))
  
  do j=1,dat%ny
   do i=0,dat%nx
    call point_at(null_field,dat%bx(i,j)%z)  
   end do
  end do
  do j=1,dat%ny
   do i=1,dat%nx+1
    call point_at(null_field,dat%ux(i,j)%z)  
   end do
  end do
  do j=0,dat%ny
   do i=1,dat%nx
    call point_at(null_field,dat%by(i,j)%z)  
   end do
  end do
  do j=1,dat%ny+1
   do i=1,dat%nx
    call point_at(null_field,dat%uy(i,j)%z)  
   end do
  end do
  
  
  do j=1,dat%ny
   do i=1,dat%nx
    if ((init_grid(dat%lx+i,dat%ly+j) == proc_name).or.   &
        (init_grid(dat%lx+i-x,dat%ly+j) == proc_name).or.   &
        (init_grid(dat%lx+i,dat%ly+j-y) == proc_name)) then
     dat%p(i,j)%core = .true.
     dat%p(i,j)%i = dat%lx+i
     dat%p(i,j)%j = dat%ly+j
     dat%p(i,j)%flux = .false.
    else
     dat%p(i,j)%core = .false.
    end if
   end do
  end do
  
  
    dat%synced=.true.
  
  do k=0,proc_num-1
   if (k /= proc_name) then
   
    i=count(((init_grid(1:mx+x,2:my) == k).or.   &
             (init_grid(1-x:mx,2:my) == k)).and.  &
            ((init_grid(1:mx+x,1:my-1) /= k).and.   &
             (init_grid(1-x:mx,1:my-1) /= k)).and.  &
            ((max(init_grid(1:mx+x,1:my-1),   &
                  init_grid(1-x:mx,1:my-1)) == proc_name))) + &
        count(((init_grid(1:mx+x,1:my-1) == k).or.   &
               (init_grid(1-x:mx,1:my-1) == k)).and.  &
              ((init_grid(1:mx+x,2:my) /= k).and.   &
               (init_grid(1-x:mx,2:my) /= k)).and.  &
              ((max(init_grid(1:mx+x,2:my),   &
                    init_grid(1-x:mx,2:my)) == proc_name))) + &
        count(((init_grid(2:mx,1:my+y) == k).or.   &
               (init_grid(2:mx,1-y:my) == k)).and.  &
              ((init_grid(1:mx-1,1:my+y) /= k).and.   &
               (init_grid(1:mx-1,1-y:my) /= k)).and.  &
              ((max(init_grid(1:mx-1,1:my+y),  &
                    init_grid(1:mx-1,1-y:my)) == proc_name))) + &
        count(((init_grid(1:mx-1,1:my+y) == k).or.   &
               (init_grid(1:mx-1,1-y:my) == k)).and.  &
              ((init_grid(2:mx,1:my+y) /= k).and.   &
               (init_grid(2:mx,1-y:my) /= k)).and.  &
              ((max(init_grid(2:mx,1:my+y),  &
                    init_grid(2:mx,1-y:my)) == proc_name)))
    j=count(((init_grid(1:mx+x,2:my) == proc_name).or.   &
             (init_grid(1-x:mx,2:my) == proc_name)).and.  &
            ((init_grid(1:mx+x,1:my-1) /= proc_name).and.   &
             (init_grid(1-x:mx,1:my-1) /= proc_name)).and.  &
            ((max(init_grid(1:mx+x,1:my-1),   &
                  init_grid(1-x:mx,1:my-1)) == k))) + &
        count(((init_grid(1:mx+x,1:my-1) == proc_name).or.   &
               (init_grid(1-x:mx,1:my-1) == proc_name)).and.  &
              ((init_grid(1:mx+x,2:my) /= proc_name).and.   &
               (init_grid(1-x:mx,2:my) /= proc_name)).and.  &
              ((max(init_grid(1:mx+x,2:my),   &
                    init_grid(1-x:mx,2:my)) == k))) + &
        count(((init_grid(2:mx,1:my+y) == proc_name).or.   &
               (init_grid(2:mx,1-y:my) == proc_name)).and.  &
              ((init_grid(1:mx-1,1:my+y) /= proc_name).and.   &
               (init_grid(1:mx-1,1-y:my) /= proc_name)).and.  &
              ((max(init_grid(1:mx-1,1:my+y),  &
                    init_grid(1:mx-1,1-y:my)) == k))) + &
        count(((init_grid(1:mx-1,1:my+y) == proc_name).or.   &
               (init_grid(1:mx-1,1-y:my) == proc_name)).and.  &
              ((init_grid(2:mx,1:my+y) /= proc_name).and.   &
               (init_grid(2:mx,1-y:my) /= proc_name)).and.  &
              ((max(init_grid(2:mx,1:my+y),  &
                    init_grid(2:mx,1-y:my)) == k)))
    if (horiz_loop) then
     i=i+count(((init_grid(1,1:my+y) == k).or.   &
               (init_grid(1,1-y:my) == k)).and.  &
              ((max(init_grid(mx,1:my+y),  &
                    init_grid(mx,1-y:my)) == proc_name))) + &
        count(((init_grid(mx,1:my+y) == k).or.   &
               (init_grid(mx,1-y:my) == k)).and.  &
              ((max(init_grid(1,1:my+y),  &
                    init_grid(1,1-y:my)) == proc_name)))
     j=j+count(((init_grid(1,1:my+y) == proc_name).or.   &
               (init_grid(1,1-y:my) == proc_name)).and.  &
              ((max(init_grid(mx,1:my+y),  &
                    init_grid(mx,1-y:my)) == k))) + &
        count(((init_grid(mx,1:my+y) == proc_name).or.   &
               (init_grid(mx,1-y:my) == proc_name)).and.  &
              ((max(init_grid(1,1:my+y),  &
                    init_grid(1,1-y:my)) == k)))
    end if
    if (verti_loop) then
     i=i+count(((init_grid(1:mx+x,1) == k).or.   &
             (init_grid(1-x:mx,1) == k)).and.  &
            ((max(init_grid(1:mx+x,my),   &
                  init_grid(1-x:mx,my)) == proc_name))) + &
        count(((init_grid(1:mx+x,my) == k).or.   &
               (init_grid(1-x:mx,my) == k)).and.  &
              ((max(init_grid(1:mx+x,1),   &
                    init_grid(1-x:mx,1)) == proc_name)))
     j=j+count(((init_grid(1:mx+x,1) == proc_name).or.   &
             (init_grid(1-x:mx,1) == proc_name)).and.  &
            ((max(init_grid(1:mx+x,my),   &
                  init_grid(1-x:mx,my)) == k))) + &
        count(((init_grid(1:mx+x,my) == proc_name).or.   &
               (init_grid(1-x:mx,my) == proc_name)).and.  &
              ((max(init_grid(1:mx+x,1),   &
                    init_grid(1-x:mx,1)) == k)))
    end if
   
    if (i /= 0) then
     allocate(dat%mpi(k)%s(i))
     allocate(dat%mpi(k)%s_mpi(i))
     dat%mpi(k)%s_req=mpi_request_null
    end if
    if (j /= 0) then
     allocate(dat%mpi(k)%r(j))
     allocate(dat%mpi(k)%r_mpi(j))
     dat%mpi(k)%r_req=mpi_request_null
     if ((x == 1).or.(y == 1)) then
      dat%mpi(k)%r%vec=-1
     else
      dat%mpi(k)%r%vec=1
     end if
    end if
    
   else
   
   
    i=count(((init_grid(1:mx+x,2:my) == proc_name).or.   &
             (init_grid(1-x:mx,2:my) == proc_name)).and.  &
            ((init_grid(1:mx+x,1:my-1) /= proc_name).and.   &
             (init_grid(1-x:mx,1:my-1) /= proc_name)).and.  &
            ((max(init_grid(1:mx+x,1:my-1),   &
                  init_grid(1-x:mx,1:my-1)) == -1))) + &
        count(((init_grid(1:mx+x,1:my-1) == proc_name).or.   &
               (init_grid(1-x:mx,1:my-1) == proc_name)).and.  &
              ((init_grid(1:mx+x,2:my) /= proc_name).and.   &
               (init_grid(1-x:mx,2:my) /= proc_name)).and.  &
              ((max(init_grid(1:mx+x,2:my),   &
                    init_grid(1-x:mx,2:my)) == -1))) + &
        count(((init_grid(2:mx,1:my+y) == proc_name).or.   &
               (init_grid(2:mx,1-y:my) == proc_name)).and.  &
              ((init_grid(1:mx-1,1:my+y) /= proc_name).and.   &
               (init_grid(1:mx-1,1-y:my) /= proc_name)).and.  &
              ((max(init_grid(1:mx-1,1:my+y),  &
                    init_grid(1:mx-1,1-y:my)) == -1))) + &
        count(((init_grid(1:mx-1,1:my+y) == proc_name).or.   &
               (init_grid(1:mx-1,1-y:my) == proc_name)).and.  &
              ((init_grid(2:mx,1:my+y) /= proc_name).and.   &
               (init_grid(2:mx,1-y:my) /= proc_name)).and.  &
              ((max(init_grid(2:mx,1:my+y),  &
                    init_grid(2:mx,1-y:my)) == -1)))
    if (horiz_loop) then
     i=i+count(((init_grid(1,1:my+y) == proc_name).or.   &
               (init_grid(1,1-y:my) == proc_name)).and.  &
              ((max(init_grid(mx,1:my+y),  &
                    init_grid(mx,1-y:my)) == k))) + &
        count(((init_grid(mx,1:my+y) == proc_name).or.   &
               (init_grid(mx,1-y:my) == proc_name)).and.  &
              ((max(init_grid(1,1:my+y),  &
                    init_grid(1,1-y:my)) == k)))
    else
     i=i+count(((init_grid(1,1:my+y) == proc_name).or.   &
               (init_grid(1,1-y:my) == proc_name))) + &
        count(((init_grid(mx,1:my+y) == proc_name).or.   &
               (init_grid(mx,1-y:my) == proc_name)))
    end if
    if (verti_loop) then
     i=i+count(((init_grid(1:mx+x,1) == proc_name).or.   &
             (init_grid(1-x:mx,1) == proc_name)).and.  &
            ((max(init_grid(1:mx+x,my),   &
                  init_grid(1-x:mx,my)) == k))) + &
        count(((init_grid(1:mx+x,my) == proc_name).or.   &
               (init_grid(1-x:mx,my) == proc_name)).and.  &
              ((max(init_grid(1:mx+x,1),   &
                    init_grid(1-x:mx,1)) == k)))
    else
     i=i+count(((init_grid(1:mx+x,1) == proc_name).or.   &
             (init_grid(1-x:mx,1) == proc_name))) + &
        count(((init_grid(1:mx+x,my) == proc_name).or.   &
               (init_grid(1-x:mx,my) == proc_name)))
    end if
   
     
    if (i /= 0) then
     allocate(dat%mpi(k)%r(i))
     allocate(dat%mpi(k)%r_mpi(i))
     allocate(dat%mpi(k)%s(i))
     allocate(dat%mpi(k)%s_mpi(i))
     dat%mpi(k)%r_req=mpi_request_null
     dat%mpi(k)%s_req=mpi_request_null
     if (vec) then
      dat%mpi(k)%r%vec=-1
     else
      dat%mpi(k)%r%vec=1
     end if
    end if
   
   end if
  end do
  
  
  do j=1,my+y
   do i=1,mx+x
    if ((init_grid(i,j) == proc_name).or.   &
        (init_grid(i-x,j) == proc_name).or.   &
        (init_grid(i,j-y) == proc_name)) then
     call point_at(dat%z(i-dat%lx,j-dat%ly),dat%bx(i-dat%lx,j-dat%ly)%z)  
     call point_at(dat%z(i-dat%lx,j-dat%ly),dat%ux(i-dat%lx,j-dat%ly)%z) 
     call point_at(dat%z(i-dat%lx,j-dat%ly),dat%by(i-dat%lx,j-dat%ly)%z) 
     call point_at(dat%z(i-dat%lx,j-dat%ly),dat%uy(i-dat%lx,j-dat%ly)%z)
    end if
   end do
  end do
  
  
  do k=0,proc_num-1
   m=0 
   l=0
   
   if (k /= proc_name) then
   
! receive ----- from right
   
    do j=1,my+y
     do i=1,mx-1
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i,j-y) == proc_name)).and.  &
              ((init_grid(i+1,j) /= proc_name).and.   &
               (init_grid(i+1,j-y) /= proc_name)).and.  &
              ((max(init_grid(i+1,j),  &
                    init_grid(i+1,j-y)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%ux(i+x-dat%lx+1,j-dat%ly)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end do
     if (horiz_loop) then
      if (((init_grid(mx,j) == proc_name).or.   &
               (init_grid(mx,j-y) == proc_name)).and.  &
              ((max(init_grid(1,j),  &
                    init_grid(1,j-y)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%ux(mx+x-dat%lx+1,j-dat%ly)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end if
    end do
   
! receive ----- from left
   
    do j=1,my+y
     if (horiz_loop) then
      if (((init_grid(1,j) == proc_name).or.   &
               (init_grid(1,j-y) == proc_name)).and.  &
              ((max(init_grid(mx,j),  &
                    init_grid(mx,j-y)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%bx(-dat%lx,j-dat%ly)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end if
     do i=2,mx
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i,j-y) == proc_name)).and.  &
              ((init_grid(i-1,j) /= proc_name).and.   &
               (init_grid(i-1,j-y) /= proc_name)).and.  &
              ((max(init_grid(i-1,j),  &
                    init_grid(i-1,j-y)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%bx(i-dat%lx-1,j-dat%ly)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end do
    end do
    
! receive ----- from top

    do i=1,mx+x
     do j=1,my-1
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i-x,j) == proc_name)).and.  &
              ((init_grid(i,j+1) /= proc_name).and.   &
               (init_grid(i-x,j+1) /= proc_name)).and.  &
              ((max(init_grid(i,j+1),  &
                    init_grid(i-x,j+1)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%uy(i-dat%lx,j+y-dat%ly+1)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end do
     if (verti_loop) then
      if (((init_grid(i,my) == proc_name).or.   &
               (init_grid(i-x,my) == proc_name)).and.  &
              ((max(init_grid(i,1),  &
                    init_grid(i-x,1)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%uy(i-dat%lx,my+y-dat%ly+1)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end if
    end do
    
! receive ----- from bottom
   
    do i=1,mx+x
     if (verti_loop) then
      if (((init_grid(i,1) == proc_name).or.   &
               (init_grid(i-x,1) == proc_name)).and.  &
              ((max(init_grid(i,my),  &
                    init_grid(i-x,my)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%by(i-dat%lx,-dat%ly)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end if
     do j=2,my
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i-x,j) == proc_name)).and.  &
              ((init_grid(i,j-1) /= proc_name).and.   &
               (init_grid(i-x,j-1) /= proc_name)).and.  &
              ((max(init_grid(i,j-1),  &
                    init_grid(i-x,j-1)) == k))) then
       m=m+1
       call point_at(dat%mpi(k)%r(m)%z,dat%by(i-dat%lx,j-dat%ly-1)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     end do
    end do
    
! send ---- to left

    do j=1,my+y
     if (horiz_loop) then
      if (((init_grid(mx,j) == k).or.   &
               (init_grid(mx,j-y) == k)).and.  &
              ((max(init_grid(1,j),  &
                    init_grid(1,j-y)) == proc_name))) then
       l=l+1
       call point_at(dat%z(1+x-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end if
     do i=2,mx
      if (((init_grid(i-1,j) == k).or.   &
               (init_grid(i-1,j-y) == k)).and.  &
              ((init_grid(i,j) /= k).and.   &
               (init_grid(i,j-y) /= k)).and.  &
              ((max(init_grid(i,j),  &
                    init_grid(i,j-y)) == proc_name))) then
       l=l+1
       call point_at(dat%z(i+x-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end do
    end do

! send ---- to right

    do j=1,my+y
     do i=1,mx-1
      if (((init_grid(i+1,j) == k).or.   &
               (init_grid(i+1,j-y) == k)).and.  &
              ((init_grid(i,j) /= k).and.   &
               (init_grid(i,j-y) /= k)).and.  &
              ((max(init_grid(i,j),  &
                    init_grid(i,j-y)) == proc_name)))  then
       l=l+1
       call point_at(dat%z(i-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end do
     if (horiz_loop) then
      if (((init_grid(1,j) == k).or.   &
               (init_grid(1,j-y) == k)).and.  &
              ((max(init_grid(mx,j),  &
                    init_grid(mx,j-y)) == proc_name))) then
       l=l+1
       call point_at(dat%z(mx-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end if
    end do

! send ---- to bottom

    do i=1,mx+x
     if (verti_loop) then
      if (((init_grid(i,my) == k).or.   &
               (init_grid(i-x,my) == k)).and.  &
              ((max(init_grid(i,1),  &
                    init_grid(i-x,1)) == proc_name))) then
       l=l+1
       call point_at(dat%z(i-dat%lx,1+y-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end if
     do j=2,my
      if (((init_grid(i,j-1) == k).or.   &
               (init_grid(i-x,j-1) == k)).and.  &
              ((init_grid(i,j) /= k).and.   &
               (init_grid(i-x,j) /= k)).and.  &
              ((max(init_grid(i,j),  &
                    init_grid(i-x,j)) == proc_name))) then
       l=l+1
       call point_at(dat%z(i-dat%lx,j+y-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end do
    end do

! send ---- to top

    do i=1,mx+x
     do j=1,my-1
      if  (((init_grid(i,j+1) == k).or.   &
               (init_grid(i-x,j+1) == k)).and.  &
              ((init_grid(i,j) /= k).and.   &
               (init_grid(i-x,j) /= k)).and.  &
              ((max(init_grid(i,j),  &
                    init_grid(i-x,j)) == proc_name))) then
       l=l+1
       call point_at(dat%z(i-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end do
     if (verti_loop) then
      if (((init_grid(i,1) == k).or.   &
               (init_grid(i-x,1) == k)).and.  &
              ((max(init_grid(i,my),  &
                    init_grid(i-x,my)) == proc_name))) then
       l=l+1
       call point_at(dat%z(i-dat%lx,my-dat%ly),dat%mpi(k)%s(l)%z)
      end if
     end if
    end do
    
   else
    
  ! receive ----- from right
   
    do j=1,my+y
     do i=1,mx-1
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i,j-y) == proc_name)).and.  &
              ((init_grid(i+1,j) /= proc_name).and.   &
               (init_grid(i+1,j-y) /= proc_name)).and.  &
              ((max(init_grid(i+1,j),  &
                    init_grid(i+1,j-y)) == -1))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%ux(i+x-dat%lx+1,j-dat%ly)%z)
       call point_at(dat%z(i-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (x == 1) then
        dat%p(i-dat%lx+x,j-dat%ly)%flux=.true.
       end if
      end if
     end do
     if (horiz_loop) then
      if (((init_grid(mx,j) == proc_name).or.   &
               (init_grid(mx,j-y) == proc_name)).and.  &
              ((max(init_grid(1,j),  &
                    init_grid(1,j-y)) == k))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%ux(mx+x-dat%lx+1,j-dat%ly)%z)
       call point_at(dat%z(1+x-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     else
      if ((init_grid(mx,j) == proc_name).or.   &
               (init_grid(mx,j-y) == proc_name)) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%ux(mx+x-dat%lx+1,j-dat%ly)%z)
       call point_at(dat%z(mx-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (x == 1) then
        dat%p(mx-dat%lx+x,j-dat%ly)%flux=.true.
       end if
      end if
     end if
    end do
   
! receive ----- from left
   
    do j=1,my+y
     if (horiz_loop) then
      if (((init_grid(1,j) == proc_name).or.   &
               (init_grid(1,j-y) == proc_name)).and.  &
              ((max(init_grid(mx,j),  &
                    init_grid(mx,j-y)) == k))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%bx(-dat%lx,j-dat%ly)%z)
       call point_at(dat%z(mx-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     else
      if ((init_grid(1,j) == proc_name).or.   &
               (init_grid(1,j-y) == proc_name)) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%bx(-dat%lx,j-dat%ly)%z)
       call point_at(dat%z(1-dat%lx+x,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (x == 1) then
        dat%p(1-dat%lx,j-dat%ly)%flux=.true.
       end if
      end if
     end if
     do i=2,mx
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i,j-y) == proc_name)).and.  &
              ((init_grid(i-1,j) /= proc_name).and.   &
               (init_grid(i-1,j-y) /= proc_name)).and.  &
              ((max(init_grid(i-1,j),  &
                    init_grid(i-1,j-y)) == -1))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%bx(i-dat%lx-1,j-dat%ly)%z)
       call point_at(dat%z(i-dat%lx+x,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (x == 1) then
        dat%p(i-dat%lx,j-dat%ly)%flux=.true.
       end if
      end if
     end do
    end do
    
! receive ----- from top

    do i=1,mx+x
     do j=1,my-1
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i-x,j) == proc_name)).and.  &
              ((init_grid(i,j+1) /= proc_name).and.   &
               (init_grid(i-x,j+1) /= proc_name)).and.  &
              ((max(init_grid(i,j+1),  &
                    init_grid(i-x,j+1)) == -1))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%uy(i-dat%lx,j+y-dat%ly+1)%z)
       call point_at(dat%z(i-dat%lx,j-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (y == 1) then
        dat%p(i-dat%lx,j-dat%ly+y)%flux=.true.
       end if
      end if
     end do
     if (verti_loop) then
      if (((init_grid(i,my) == proc_name).or.   &
               (init_grid(i-x,my) == proc_name)).and.  &
              ((max(init_grid(i,1),  &
                    init_grid(i-x,1)) == k))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%uy(i-dat%lx,my+y-dat%ly+1)%z)
       call point_at(dat%z(i-dat%lx,1+y-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     else
      if ((init_grid(i,my) == proc_name).or.   &
               (init_grid(i-x,my) == proc_name)) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%uy(i-dat%lx,my+y-dat%ly+1)%z)
       call point_at(dat%z(i-dat%lx,my-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (y == 1) then
        dat%p(i-dat%lx,my-dat%ly+y)%flux=.true.
       end if
      end if
     end if
    end do
    
! receive ----- from bottom
   
    do i=1,mx+x
     if (verti_loop) then
      if (((init_grid(i,1) == proc_name).or.   &
               (init_grid(i-x,1) == proc_name)).and.  &
              ((max(init_grid(i,my),  &
                    init_grid(i-x,my)) == k))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%by(i-dat%lx,-dat%ly)%z)
       call point_at(dat%z(i-dat%lx,my-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=1
      end if
     else
      if ((init_grid(i,1) == proc_name).or.   &
               (init_grid(i-x,1) == proc_name)) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%by(i-dat%lx,-dat%ly)%z)
       call point_at(dat%z(i-dat%lx,1+y-dat%ly),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (y == 1) then
        dat%p(i-dat%lx,1-dat%ly)%flux=.true.
       end if
      end if
     end if
     do j=2,my
      if (((init_grid(i,j) == proc_name).or.   &
               (init_grid(i-x,j) == proc_name)).and.  &
              ((init_grid(i,j-1) /= proc_name).and.   &
               (init_grid(i-x,j-1) /= proc_name)).and.  &
              ((max(init_grid(i,j-1),  &
                    init_grid(i-x,j-1)) == -1))) then
       m=m+1
       l=l+1
       call point_at(dat%mpi(k)%r(m)%z,dat%by(i-dat%lx,j-dat%ly-1)%z)
       call point_at(dat%z(i-dat%lx,j-dat%ly+y),dat%mpi(k)%s(l)%z)
       dat%mpi(k)%r(m)%dir=-1
       if (y == 1) then
        dat%p(i-dat%lx,j-dat%ly)%flux=.true.
       end if
      end if
     end do
    end do
    
     
   end if
   
  end do
  
  
  dat%p%x=(dble(dat%p%i)-0.5d0-0.5d0*dble(x))*dx
  dat%p%y=(dble(dat%p%j)-0.5d0-0.5d0*dble(y))*dy
  do k=0,proc_num-1
   if (allocated(dat%mpi(k)%r)) then
    dat%mpi(k)%r%slip=dble(max(dat%mpi(k)%r%vec,dat%mpi(k)%r%dir))
   end if
  end do
  
  
   
 end subroutine
 
 
 

 
 
 
 
 subroutine create_hparam(dat)
  use global
  use parallel
  
  implicit none
  
  type(var), intent(out) :: dat
  integer :: i, j
  
  dat%nx = count((count((init_grid(1:mx,1:my) == proc_name),2) /= 0))
  dat%lx = 0
  do i=1,mx
   if (count((count((init_grid(1:i,1:my) == proc_name),2) /= 0)) == 0) dat%lx = dat%lx+1
  end do
  dat%ny = count((count((init_grid(1:mx,1:my) == proc_name),1) /= 0))
  dat%ly = 0
  do i=1,my
   if (count((count((init_grid(1:mx,1:i) == proc_name),1) /= 0)) == 0) dat%ly = dat%ly+1
  end do
    
  allocate(dat%z(dat%nx,dat%ny))
  allocate(dat%p(dat%nx,dat%ny))
  allocate(dat%exp(dat%nx,dat%ny))
  
  
  do j=1,dat%ny
   do i=1,dat%nx
    if (init_grid(dat%lx+i,dat%ly+j) == proc_name) then
     dat%p(i,j)%core = .true.
     dat%p(i,j)%i = dat%lx+i
     dat%p(i,j)%j = dat%ly+j
     dat%p(i,j)%flux = .false.
    else
     dat%p(i,j)%core = .false.
    end if
   end do
  end do
  
  
  dat%p%x=(dble(dat%p%i)-0.5d0)*dx
  dat%p%y=(dble(dat%p%j)-0.5d0)*dy
   
   
 end subroutine
 
 
 
 subroutine create_uparam(dat)
  use global
  use parallel
  
  implicit none
  
  type(var), intent(out) :: dat
  integer :: i, j
  
  dat%nx = count((count((init_grid(1:mx,1:my) == proc_name),2) /= 0)) + 1
  dat%lx = 0
  do i=1,mx
   if (count((count((init_grid(1:i,1:my) == proc_name),2) /= 0)) == 0) dat%lx = dat%lx+1
  end do
  dat%ny = count((count((init_grid(1:mx,1:my) == proc_name),1) /= 0))
  dat%ly = 0
  do i=1,my
   if (count((count((init_grid(1:mx,1:i) == proc_name),1) /= 0)) == 0) dat%ly = dat%ly+1
  end do
  
  
  allocate(dat%z(dat%nx,dat%ny))
  allocate(dat%p(dat%nx,dat%ny))
  allocate(dat%exp(dat%nx,dat%ny))
  
  
  do j=1,dat%ny
   do i=1,dat%nx
    if ((init_grid(dat%lx+i,dat%ly+j) == proc_name) .or.   &
         (init_grid(dat%lx+i-1,dat%ly+j) == proc_name)) then
     dat%p(i,j)%core = .true.
     dat%p(i,j)%i = dat%lx+i
     dat%p(i,j)%j = dat%ly+j
     dat%p(i,j)%flux = .false.
    else
     dat%p(i,j)%core = .false.
    end if
   end do
  end do
  
  
  dat%p%x=(dble(dat%p%i)-1.0d0)*dx
  dat%p%y=(dble(dat%p%j)-0.5d0)*dy
   
   
 end subroutine
 
 
 
 
 subroutine create_vparam(dat)
  use global
  use parallel
  
  implicit none
  
  type(var), intent(out) :: dat
  integer :: i, j
  
  dat%nx = count((count((init_grid(1:mx,1:my) == proc_name),2) /= 0))
  dat%lx = 0
  do i=1,mx
   if (count((count((init_grid(1:i,1:my) == proc_name),2) /= 0)) == 0) dat%lx = dat%lx+1
  end do
  dat%ny = count((count((init_grid(1:mx,1:my) == proc_name),1) /= 0)) + 1
  dat%ly = 0
  do i=1,my
   if (count((count((init_grid(1:mx,1:i) == proc_name),1) /= 0)) == 0) dat%ly = dat%ly+1
  end do
  
  
  allocate(dat%z(dat%nx,dat%ny))
  allocate(dat%p(dat%nx,dat%ny))
  allocate(dat%exp(dat%nx,dat%ny))
  
  
  do j=1,dat%ny
   do i=1,dat%nx
    if ((init_grid(dat%lx+i,dat%ly+j) == proc_name) .or.   &
         (init_grid(dat%lx+i,dat%ly+j-1) == proc_name)) then
     dat%p(i,j)%core = .true.
     dat%p(i,j)%i = dat%lx+i
     dat%p(i,j)%j = dat%ly+j
     dat%p(i,j)%flux = .false.
    else
     dat%p(i,j)%core = .false.
    end if
   end do
  end do
  
  
  
  dat%p%x=(dble(dat%p%i)-0.5d0)*dx
  dat%p%y=(dble(dat%p%j)-1.0d0)*dy
   
   
 end subroutine
 
 
 
  subroutine create_zparam(dat)
  use global
  use parallel
  
  implicit none
  
  type(var), intent(out) :: dat
  integer :: i, j
  
  dat%nx = count((count((init_grid(1:mx,1:my) == proc_name),2) /= 0)) + 1
  dat%lx = 0
  do i=1,mx
   if (count((count((init_grid(1:i,1:my) == proc_name),2) /= 0)) == 0) dat%lx = dat%lx+1
  end do
  dat%ny = count((count((init_grid(1:mx,1:my) == proc_name),1) /= 0)) + 1
  dat%ly = 0
  do i=1,my
   if (count((count((init_grid(1:mx,1:i) == proc_name),1) /= 0)) == 0) dat%ly = dat%ly+1
  end do
  
  
  allocate(dat%z(dat%nx,dat%ny))
  allocate(dat%p(dat%nx,dat%ny))
  allocate(dat%exp(dat%nx,dat%ny))
  
  
  do j=1,dat%ny
   do i=1,dat%nx
    if ((init_grid(dat%lx+i,dat%ly+j) == proc_name) .or.   &
         (init_grid(dat%lx+i,dat%ly+j-1) == proc_name) .or. &
         (init_grid(dat%lx+i-1,dat%ly+j) == proc_name) .or.   &
         (init_grid(dat%lx+i-1,dat%ly+j-1) == proc_name)) then
     dat%p(i,j)%core = .true.
     dat%p(i,j)%i = dat%lx+i
     dat%p(i,j)%j = dat%ly+j
     dat%p(i,j)%flux = .false.
    else
     dat%p(i,j)%core = .false.
    end if
   end do
  end do
  
  
  
  dat%p%x=(dble(dat%p%i)-1.0d0)*dx
  dat%p%y=(dble(dat%p%j)-1.0d0)*dy
   
   
 end subroutine
 
 subroutine timestepped(dat)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  
  allocate(dat%tend0(dat%nx,dat%ny))
  allocate(dat%tend1(dat%nx,dat%ny))
  allocate(dat%tend2(dat%nx,dat%ny))
  allocate(dat%tmp(dat%nx,dat%ny))
  
 end subroutine
 
 
 subroutine stepforward(dat)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  
  
   dat%tmp => dat%tend2
   dat%tend2 => dat%tend1
   dat%tend1 => dat%tend0
   dat%tend0 => dat%tmp
  
 end subroutine
 
 
 
 subroutine remove_param(dat)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  
  deallocate(dat%z)
  deallocate(dat%p)
 
 end subroutine
 
end module
