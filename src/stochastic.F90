module linklists
 use global

 implicit none
 
 type llist
  real(kind=8) :: rand(6)
  real(kind=8) :: timenow = 0.0d0
  real(kind=8) :: timeend
  logical :: first, last
  type(llist), pointer :: next
  type(llist), pointer :: tmp
  type(var) :: u, v
 end type
 
 type(llist), pointer :: stochwind
 
end module
  
module llist_ops
 use linklists
 
 implicit none
 
 contains
 
 subroutine create_link(linkedlist)
  use linklists
 
  implicit none
  
  type(llist), pointer :: linkedlist
  
  if (associated(linkedlist)) then
   do while (.not.(linkedlist%last))
    linkedlist => linkedlist%next
   end do
   linkedlist%tmp => linkedlist%next
   nullify(linkedlist%next)
   allocate(linkedlist%next)
   linkedlist%last=.false.
   linkedlist%next%next => linkedlist%tmp
   nullify(linkedlist%tmp)
   linkedlist => linkedlist%next
  else
   allocate(linkedlist)
   linkedlist%next => linkedlist
   linkedlist%first=.true.
  end if
  call random_number(linkedlist%rand)
  call random_number(linkedlist%timeend)
  linkedlist%last=.true.
  
 end subroutine
  
 subroutine remove_link(linkedlist)
  use linklists
  use grid
  
  implicit none
  
  type(llist), pointer :: linkedlist
  
  if (associated(linkedlist)) then
   if (linkedlist%next%timenow >= linkedlist%next%timeend) then
    call remove_param(linkedlist%next%u)
    call remove_param(linkedlist%next%v)
    if (.not.(associated(linkedlist,linkedlist%next))) then
     if (linkedlist%first) linkedlist%next%next%first=.true.
     if (linkedlist%next%last) linkedlist%last=.true.
     linkedlist%tmp => linkedlist%next%next
     nullify(linkedlist%next%next)
     deallocate(linkedlist%next)
     linkedlist%next => linkedlist%tmp
     nullify(linkedlist%tmp)
    else
     deallocate(linkedlist)
    end if
   else
    linkedlist => linkedlist%next
   end if
  end if
  
 end subroutine
 
 subroutine step_link(linkedlist,timestep)
  use linklists
  
  implicit none
  
  type(llist), pointer :: linkedlist
  real(kind=8), intent(in) :: timestep
  
  
 ! call remove_link(linkedlist)
  linkedlist => linkedlist%next
  linkedlist%timenow = linkedlist%timenow + timestep
  
 end subroutine
  
 subroutine set_wind(n)
  use global
  use linklists
  use params
  use params
  use params
  use variables
  use operate
  use grid
 
  implicit none
  
  integer :: i, j, k, l
  integer, intent(in) :: n
  real(kind=8) :: r

  if (associated(stochwind)) then
   do while (.not.(stochwind%next%last))
    call remove_link(stochwind)
   end do
   call remove_link(stochwind)
  end if
  if (mod(n,ceiling(1.0d0*8.64d4/(dt/f0))) == 0.0d0) then
   call create_link(stochwind)
   call create_uparam(stochwind%u)
   call create_vparam(stochwind%v)
   stochwind%timeend=6.048d5*(4.0d0+12.0d0*stochwind%timeend)*f0
   stochwind%rand(1)=sign(tau*0.5d0,stochwind%rand(1))+  &
              tau*0.5d0*(stochwind%rand(1))
   stochwind%rand(3)=x0*(stochwind%rand(3))
   stochwind%rand(4)=y0*(stochwind%rand(4))
   stochwind%rand(5)=2.0d0*(stochwind%rand(5)-0.5d0)
   stochwind%rand(5)=0.25d0*x0*(stochwind%rand(5))
   stochwind%rand(6)=2.0d0*(stochwind%rand(6)-0.5d0)
   stochwind%rand(6)=0.25d0*y0*(stochwind%rand(6))
   stochwind%rand(2)=(1.5d0+(x0/4.0d0-1.5d0)*(stochwind%rand(2)))
   do j=1,utau%ny
    do i=1,utau%nx
     r=radius(utau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0),lim(stochwind%rand(4)+stochwind%rand(6),y0))/stochwind%rand(2)
      stochwind%u%z(i,j)=-stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))   
          
     r=radius(utau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0),lim(stochwind%rand(4)-stochwind%rand(6),y0))/stochwind%rand(2)
      stochwind%u%z(i,j)=stochwind%u%z(i,j)+stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*           &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))
                      
    end do
   end do
   do j=1,vtau%ny
    do i=1,vtau%nx
     r=radius(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0),lim(stochwind%rand(4)+stochwind%rand(6),y0))/stochwind%rand(2)
     stochwind%v%z(i,j)=stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*           &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))

     r=radius(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0),lim(stochwind%rand(4)-stochwind%rand(6),y0))/stochwind%rand(2)
     
     stochwind%v%z(i,j)= stochwind%v%z(i,j)-stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0))/(r*stochwind%rand(2)),r == 0.0d0)    &
             *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*               &
           (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*      &
           (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))


    end do
   end do
  end if  

  utau%z=0.0d0
  vtau%z=0.0d0
  call step_link(stochwind,dt)
  utau%z=utau%z+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%u%z
  vtau%z=vtau%z+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%v%z
!  print *, ceiling(stochwind%rand), ceiling(stochwind%timenow), &
!        ceiling(stochwind%timeend)
  do while (.not.(stochwind%last))
   call step_link(stochwind,dt)
   utau%z=utau%z+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%u%z
   vtau%z=vtau%z+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%v%z
!   print *, ceiling(stochwind%rand), ceiling(stochwind%timenow), &
!        ceiling(stochwind%timeend)
  end do
  
 end subroutine
  
  
 subroutine output_llist(linkedlist,varname)
  use linklists
  use parallel
  
  implicit none
  
  type(llist), pointer :: linkedlist
  character(*), intent(in) :: varname 
  character(32) :: filename, foldername, fullname
  
  if (proc_name == ens_master) then
   write (foldername, "(a5,i4.4)") 'out/', ens_name
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
   if (associated(linkedlist)) then
    open(unit=10,file=adjustl(fullname),status='unknown',form='unformatted')
    
    linkedlist => linkedlist%next
    write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    do while (.not.(stochwind%last))
     linkedlist => linkedlist%next
     write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    end do
    
    close(10)
   end if
  end if
  
 end subroutine
  
  
 subroutine writeout_llist(num,linkedlist,varname)
  use linklists
  use parallel
  
  implicit none
  
  type(llist), pointer :: linkedlist
  character(*), intent(in) :: varname 
  character(32) :: filename, foldername, fullname
  integer, intent(in) :: num
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
 
   write (foldername, "(a5,i4.4,a1,i4.4)") 'data/', ens_name,'/', num
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
   if (associated(linkedlist)) then
    open(unit=10,file=adjustl(fullname),status='unknown',form='unformatted')
    
    linkedlist => linkedlist%next
    write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    do while (.not.(stochwind%last))
     linkedlist => linkedlist%next
     write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    end do
    
    close(10)
   end if
  end if
  
 end subroutine
   
  
 subroutine input_llist(linkedlist,varname)
  use linklists
  use variables
  use params
  use operate
  use grid
  use parallel

  
  implicit none
  
  type(llist), pointer :: linkedlist
  type(llist), pointer :: tmplist
  character(*), intent(in) :: varname 
  logical :: file_exist
  integer :: endoffile, i, j
  real(kind=8) :: r, maxum
  integer :: k 
  character(32) :: filename, foldername, fullname
  
   write (foldername, "(a5,i4.4)") 'old/', ens_name
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
  inquire( file=adjustl(fullname), exist=file_exist)
  
  if (file_exist) then
 
   open(unit=10,file=adjustl(filename),status='unknown',form='unformatted')
   allocate(tmplist)
   read(10,iostat=endoffile) tmplist%rand, tmplist%timenow, tmplist%timeend
   do while (endoffile == 0)
    if (.not.(associated(linkedlist))) then
     linkedlist => tmplist
     linkedlist%next => linkedlist
     linkedlist%first = .true.
     linkedlist%last = .true.
     nullify(tmplist)
    else
     do while (.not.(linkedlist%last))
      linkedlist => linkedlist%next
     end do      
     linkedlist%tmp => linkedlist%next
     linkedlist%next => tmplist
     linkedlist%next%last = .true.
     linkedlist%last = .false.
     linkedlist%next%next => linkedlist%tmp
     nullify(linkedlist%tmp)
     nullify(tmplist)
     linkedlist => linkedlist%next
   
    end if
    call create_uparam(stochwind%u)
    call create_vparam(stochwind%v)
    

   do j=1,utau%ny
    do i=1,utau%nx
     r=radius(utau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0),lim(stochwind%rand(4)+stochwind%rand(6),y0))/stochwind%rand(2)
      stochwind%u%z(i,j)=-stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))   
          
     r=radius(utau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0),lim(stochwind%rand(4)-stochwind%rand(6),y0))/stochwind%rand(2)
      stochwind%u%z(i,j)= stochwind%u%z(i,j)+stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*           &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))
                      
    end do
   end do
   do j=1,vtau%ny
    do i=1,vtau%nx
     r=radius(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0),lim(stochwind%rand(4)+stochwind%rand(6),y0))/stochwind%rand(2)
     stochwind%v%z(i,j)=stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0))/(r*stochwind%rand(2)),r == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau%p(i,j),lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*           &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau%p(i,j),lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))

     r=radius(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0),lim(stochwind%rand(4)-stochwind%rand(6),y0))/stochwind%rand(2)
     
     stochwind%v%z(i,j)= stochwind%v%z(i,j)-stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0))/(r*stochwind%rand(2)),r == 0.0d0)    &
             *(1.0d0/(2.0d0*pi*r))*(1.0d0-exp(-r**2))*exp(-r**2/2.0d0)*               &
           (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau%p(i,j),lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*      &
           (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau%p(i,j),lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))


    end do
   end do




    allocate(tmplist)
    read(10,iostat=endoffile) tmplist%rand, tmplist%timenow, tmplist%timeend
   end do 
   
   deallocate(tmplist)
    
  
   close(10)
   
  end if
  
 end subroutine
 
end module
  
