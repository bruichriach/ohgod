

module sync

#include "include.h"
 
 implicit none
 
 contains
 
 
 subroutine start_sync(dat)
  use global
  use parallel
  
  implicit none
 
  type(var), intent(inout) :: dat
  integer :: i,k
  
  call end_sync(dat)
  if (dat%synced) then
   do k=0,proc_num-1
    if (allocated(dat%mpi(k)%r_mpi)) then
     call mpi_wait(dat%mpi(k)%r_req,mpi_status_ignore,stat)
     call mpi_irecv(dat%mpi(k)%r_mpi(1),ubound(dat%mpi(k)%r_mpi,1),mpi_double_precision, &
        k,dat%tag,mpi_comm_world,dat%mpi(k)%r_req,stat)
    end if
    if (allocated(dat%mpi(k)%s_mpi)) then
     call mpi_wait(dat%mpi(k)%s_req,mpi_status_ignore,stat)
     do i=1,ubound(dat%mpi(k)%s_mpi,1)
      dat%mpi(k)%s_mpi(i)=dat%mpi(k)%s(i)%z
     end do
     call mpi_isend(dat%mpi(k)%s_mpi(1),ubound(dat%mpi(k)%s_mpi,1),mpi_double_precision, &
        k,dat%tag,mpi_comm_world,dat%mpi(k)%s_req,stat)
    end if
   end do
   dat%synced = .false.
  end if
  
   
 end subroutine
 
 subroutine end_sync(dat)
  use global
  use parallel
  
  implicit none
 
  type(var), intent(inout) :: dat
  integer :: i,k
 
  
  if (.not.(dat%synced)) then
   do k=0,proc_num-1
    if (allocated(dat%mpi(k)%r_mpi)) then
     call mpi_wait(dat%mpi(k)%r_req,mpi_status_ignore,stat)
     do i=1,ubound(dat%mpi(k)%r_mpi,1)
      dat%mpi(k)%r(i)%z=sign(dat%mpi(k)%r_mpi(i),dat%mpi(k)%r_mpi(i)*dat%mpi(k)%r(i)%slip)
     end do
    end if
   end do
   dat%synced = .true.
  end if
  
   
  
 end subroutine
 
 
  
  subroutine integer_allmax(in, out)
   use parallel
   
   implicit none
  
   integer, intent(in) :: in
   integer, intent(out) :: out
   integer :: out_tmp(0:ens_images-1), send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0
   
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             MPI_INTEGER,ens_master+i,40000,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             MPI_INTEGER,ens_master+i,40000,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=max(out,out_tmp(i))
    call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
   end do
   
  
  end subroutine


  subroutine real_allsum(in, out)
   use parallel
   use test

   implicit none
  
   real (kind=8), intent(in) :: in
   real (kind=8), intent(out) :: out
   real (kind=8) :: out_tmp(0:ens_images-1)
   integer :: send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0.0d0
  
   test_tag=mod(test_tag+1,1000)
 
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             MPI_DOUBLE_PRECISION,ens_master+mod(proc_name+i,ens_images),40000+test_tag,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             MPI_DOUBLE_PRECISION,ens_master+mod(proc_name+ens_images-i,ens_images),40000+test_tag,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=out+out_tmp(i)
   end do
   
   call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
 
 
  
  end subroutine
  
  
  
  subroutine real_allmax(in, out)
   use parallel
   
   implicit none
  
   real (kind=8), intent(in) :: in
   real (kind=8), intent(out) :: out
   real (kind=8) :: out_tmp(0:ens_images-1)
   integer :: send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0.0d0
   
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             MPI_DOUBLE_PRECISION,ens_master+i,40000,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             MPI_DOUBLE_PRECISION,ens_master+i,40000,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=max(out,out_tmp(i))
    call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
   end do
   
   
  
  end subroutine
 
end module



module system

 contains

 function print_time()
  implicit none
  
  real(kind=8) :: time
  integer :: seconds, minutes, hours, days
  character(33) :: print_time
  
  time=cputime()
  days=floor(time/8.64d4)
  time=time-dble(days)*8.64d4
  hours=floor(time/3.6d3)
  time=time-dble(hours)*3.6d3
  minutes=floor(time/6.0d1)
  time=time-dble(minutes)*6.0d1
  seconds=floor(time)
  
  write(print_time, "(i2.1,a5,i3.2,a5,i3.2,a6,i3.2,a6)") &
     days, ' dys,', &
     hours, ' hrs,', &
     minutes, ' mins,', &
     seconds, ' secs.'

 end function
 
 
 function modeltime(n)
  use params
   
  implicit none
   
  integer, intent(in) :: n
  real(kind=8) :: t
  integer :: yrs, days, hours
  character(26) :: modeltime
   
   t=(dt/f0)*dble(n)
   yrs=floor(t/3.15576e7)
   t=t-dble(yrs)*3.15576e7
   days=floor(t/8.64e4)
   t=t-dble(days)*8.64e4
   hours=floor(t/3.6e3)
   
  
  write(modeltime, "(i4.3,a5,i4.3,a5,i3.2,a5)") &
     yrs, ' yrs,', &
     days, ' dys,', &
     hours, ' hrs.'

  end function
 
 function cputime()
  implicit none
  
  real(kind=8) :: cputime
    
  call cpu_time(cputime)
  
 end function
 
 function update(seconds)
  implicit none
  
  logical :: update
  real(kind=8), save :: oldtime=-1.0d0
  real(kind=8), intent(in) :: seconds
  
  if (cputime() - oldtime >= seconds) then
   oldtime=cputime()
   update=.true.
  else
   update=.false.
  end if
   
 end function
 
 subroutine makeseed()
  use sync
 
  implicit none
  
  integer :: dt(8)
  integer :: t_tmp, t
  integer :: n, i
  integer, allocatable :: seed(:)
  
  call random_seed(size = n)
  allocate(seed(n))
  call date_and_time(values=dt)
  t_tmp = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
         + dt(2) * 31 * 24 * 60 * 60 * 1000 &
         + dt(3) * 24 * 60 * 60 * 1000 &
         + dt(5) * 60 * 60 * 1000 &
         + dt(6) * 60 * 1000 + dt(7) * 1000 &
         + dt(8)
  call integer_allmax(t_tmp, t)
  do i = 1, n
   seed(i) = lcg(t)
  end do
  
 end subroutine
  
 function lcg(s)
  implicit none
 
  integer :: lcg
  integer :: s
  
  if (s == 0) then
   s = 104729
  else
   s = mod(s, 65536)
  end if
   s = mod(s * 32768 , 65535)
   lcg = int(mod(s, int(huge(0), 8)), kind(0))
   
  end function lcg
  

end module


