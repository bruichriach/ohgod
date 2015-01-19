
module params

#include "include.h"

 implicit none
 
 real (kind=8), parameter :: pi=4.0d0*atan(1.0d0)
 real (kind=8), parameter :: omega=2.0d0*pi/8.64d4
 real (kind=8), parameter :: g=9.81d0
 
  ! Degree/radian conversion parameters
 real (kind=8), parameter :: deg2rad=pi/180.0d0, rad2deg=180.0d0/pi
 
 integer, parameter :: nz=3
 integer :: ii

#ifndef ALLOW_STATIC_LAYER
#ifdef ALLOW_RIGID_LID
 real (kind=8), parameter :: gp(nz-1) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=8), parameter :: ngp(nz) = (/ gp(1:nz-1) , g /)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#else
 real (kind=8), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3, 9.81d0 /)
 real (kind=8), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif
#else
 real (kind=8), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=8), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif

 real (kind=8), parameter :: f0=2.0d0*omega*sin(deg2rad*70.0d0)
 real (kind=8), parameter :: h0=1.0d3
#ifndef ALLOW_STATIC_LAYER
 real (kind=8), parameter :: ld=sqrt(sum(gp(1:nz-1))*h0)/f0
#else
 real (kind=8), parameter :: ld=sqrt(sum(gp(1:nz))*h0)/f0
#endif



! grid spacing and domain dimensions in terms of Ld
 
 real (kind=8), parameter :: dx = 0.1d0, dy = 0.1d0
 integer, parameter :: mx=10, my=15
 real (kind=8), parameter :: x0=dble(mx)*dx, y0=dble(my)*dy
 integer :: nx, ny
 integer, parameter :: slip=-1 ! -1 for no-slip, 1 for no grad
 logical, parameter :: horiz_loop=.true., verti_loop=.true.

! timestep in terms of 1/f

#ifdef ALLOW_RIGID_LID
 real (kind=8), parameter :: maxh=1.25d3!h0
 real (kind=8), parameter :: dt = (0.25d0*dx*ld*f0)/dsqrt(sum(gp)*maxh)
#else
 real (kind=8), parameter :: maxh=1.25d3*1.01d0
 real (kind=8), parameter :: dt = (0.25d0*dx*ld*f0)/dsqrt(sum(gp)*maxh)
#endif


 
 REAL (KIND=8), dimension(4), PARAMETER :: write_time=(/ 0.0d0 , 0.0d0 , 3.0d0 , 0.0d0 /) ! years, months, weeks, days    
 REAL (KIND=8), dimension(4), PARAMETER :: total_time=(/ 500.0d0 , 1.0d0 , 0.0d0 , 0.0d0 /)  ! years, months, weeks, days    

#ifdef DO_TIME_AVERAGE
 REAL (KIND=8), dimension(4), PARAMETER :: average_time=(/ 500.0d0 , 0.0d0 , 0.0d0 , 0.0d0 /)  ! years, months, weeks, days 
#endif
 
 INTEGER, PARAMETER :: nsteps=50000
 
 INTEGER, PARAMETER :: nstop=ceiling((3.15576d7*total_time(1)+   &
       2.592d6*total_time(2)+6.048d5*total_time(3)+8.64d4*total_time(4))/(dt/f0))
 INTEGER, PARAMETER :: wstep=floor((3.15576e7*write_time(1)+   &
       2.592d6*write_time(2)+6.048d5*write_time(3)+8.64d4*write_time(4))/(dt/f0))
#ifdef DO_TIME_AVERAGE
 INTEGER, PARAMETER :: nstep=floor((3.15576d7*average_time(1)+ &
       2.592d6*average_time(2)+6.048d5*average_time(3)+ &
        8.64d4*average_time(4))/((2*nsteps)*(dt/f0)))
 INTEGER, PARAMETER :: mbarno=ceiling(dble(nstep)*dble(2*nsteps)  &
        /(2*dble(wstep)))
 INTEGER, PARAMETER :: barno=ceiling(min(  &
       dble(nstep)*dble(2*nsteps) , &
       dble(wstep)*dble(nstop/wstep)-dble(nstep)*dble(2*nsteps))  &
        /(2*dble(wstep)))
 INTEGER, PARAMETER :: nbarno=(1+mbarno)
#endif

     
 REAL (KIND=8), PARAMETER :: bf=2.0d-2
 real (kind=8), parameter :: res=0.0d0/(f0**2*ld*h0) 
 REAL (KIND=8), PARAMETER :: tau=2.0e-1/(f0**2*ld*h0*1026.0d0)      
 real (kind=8), parameter :: cvar=(1.0d0/pi**2)    
 


end module


module parallel
 
 implicit none

 include 'mpif.h'

 integer, parameter :: max_core=24
 integer, parameter :: ens_num=1
 integer :: proc_name, proc_num, proc_master
 integer :: ens_name, ens_images, ens_master
 integer :: max_tag=0
 integer :: max_write_tag=0
 integer :: stat
 
end module
