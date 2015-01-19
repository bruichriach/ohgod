

module operate


 
 interface radius
  module procedure :: radiuss, radius2
 end interface
 
 interface x_dist
  module procedure :: x_dists, x_dist2
 end interface
 
 interface y_dist
  module procedure :: y_dists, y_dist2
 end interface
 
 
 contains
 
 pure function lim(x,x0)
  
         implicit none

         real(kind=8), intent(in) :: x, x0
         real(kind=8) :: lim

         lim=x
         do while ((lim < 0.0d0).or.(lim > x0))
          if (lim < 0.0d0) lim=lim+x0
          if (lim > x0) lim=lim-x0
         end do

 end function         
 
 pure function radiuss(p,x,y)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p
  real(kind=8), intent(in) :: x, y
  real(kind=8) :: x1, y1
  real(kind=8) :: radiuss
  
  x1=p%x-x
  y1=p%y-y
  if (horiz_loop) then
   if (x1 > x0/2.0d0) x1=x1-x0
   if (x1 < -x0/2.0d0) x1=x1+x0
  end if
  if (verti_loop) then
   if (y1 > y0/2.0d0) y1=y1-y0
   if (y1 < -y0/2.0d0) y1=y1+y0
  end if
  radiuss=sqrt(x1**2+y1**2)
  
 end function
 
 
 pure function radius2(p,x,y)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p(:,:)
  real(kind=8), intent(in) :: x, y
  integer :: i, j
  real(kind=8) :: radius2(ubound(p,1),ubound(p,2))
  
  do j=1,ubound(p,2)
   do i=1,ubound(p,1)
    radius2(i,j)=radiuss(p(i,j),x,y)
   end do
  end do
  
 end function
 
 pure function x_dists(p,x)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p
  real(kind=8), intent(in) :: x
  real(kind=8) :: x1
  real(kind=8) :: x_dists
  
  
  x1=p%x-x
  if (horiz_loop) then
   if (x1 > x0/2.0d0) x1=x1-x0
   if (x1 < -x0/2.0d0) x1=x1+x0
  end if
  x_dists=x1
  
  
 end function
 pure function x_dist2(p,x)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p(:,:)
  real(kind=8), intent(in) :: x
  integer :: i, j
  real(kind=8) :: x_dist2(ubound(p,1),ubound(p,2))
  
  do j=1,ubound(p,2)
   do i=1,ubound(p,1)
    x_dist2(i,j)=x_dists(p(i,j),x)
   end do
  end do
  
 end function
 
 pure function y_dists(p,y)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p
  real(kind=8), intent(in) :: y
  real(kind=8) :: y1
  real(kind=8) :: y_dists
  
  
  y1=p%y-y
  if (verti_loop) then
   if (y1 > y0/2.0d0) y1=y1-y0
   if (y1 < -y0/2.0d0) y1=y1+y0
  end if
  y_dists=y1
  
  
 end function
 
 pure function y_dist2(p,y)
  use global
  use params
  
  implicit none
  
  type(point), intent(in) :: p(:,:)
  real(kind=8), intent(in) :: y
  integer :: i, j
  real(kind=8) :: y_dist2(ubound(p,1),ubound(p,2))
  
  do j=1,ubound(p,2)
   do i=1,ubound(p,1)
    y_dist2(i,j)=y_dists(p(i,j),y)
   end do
  end do
  
 end function

 
 pure function gaussian(p,x,y,a)
  use global
  
  implicit none
  
  type(point), intent(in) :: p
  real(kind=8), intent(in) :: x, y, a
  real(kind=8) :: gaussian
  
  gaussian=exp(-radius(p,x,y)**2/a**2)
  
 end function
 
 
 !-------------------------------------------------------------------------!
 
 ! OPERATIONS
 
 !-------------------------------------------------------------------------!
 
 
 
 
 function sGx(dat)
  use global
  use params
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sGx(dat%nx+1,dat%ny)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny
   do i=1,dat%nx+1
    sGx(i,j)=(dat%ux(i,j)%z-dat%bx(i-1,j)%z)/dx
   end do
  end do
  
 end function
 
 function sGy(dat)
  use global
  use params
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sGy(dat%nx,dat%ny+1)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny+1
   do i=1,dat%nx
    sGy(i,j)=(dat%uy(i,j)%z-dat%by(i,j-1)%z)/dy
   end do
  end do
  
 end function
 
 
 function sAx(dat)
  use global
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sAx(dat%nx+1,dat%ny)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny
   do i=1,dat%nx+1
    sAx(i,j)=0.5d0*(dat%ux(i,j)%z+dat%bx(i-1,j)%z)
   end do
  end do
  
 end function
 
 function sAxy(dat)
  use global
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sAxy(dat%nx+1,dat%ny-1)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny-1
   do i=1,dat%nx+1
    sAxy(i,j)=0.25d0*(dat%ux(i,j+1)%z+dat%bx(i-1,j+1)%z+ &
            dat%ux(i,j)%z+dat%bx(i-1,j)%z)
   end do
  end do
  
 end function
 
 function sAy(dat)
  use global
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sAy(dat%nx,dat%ny+1)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny+1
   do i=1,dat%nx
    sAy(i,j)=0.5d0*(dat%uy(i,j)%z+dat%by(i,j-1)%z)
   end do
  end do
  
 end function
 
 
 function sAyx(dat)
  use global
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sAyx(dat%nx-1,dat%ny+1)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny+1
   do i=1,dat%nx-1
    sAyx(i,j)=0.25d0*(dat%uy(i+1,j)%z+dat%by(i+1,j-1)%z+ &
            dat%uy(i,j)%z+dat%by(i,j-1)%z)
   end do
  end do
  
 end function
 
 function sGxx(dat)
  use global
  use params
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sGxx(dat%nx,dat%ny)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny
   do i=1,dat%nx
    sGxx(i,j)=(dat%ux(i+1,j)%z-2.0d0*dat%z(i,j)+dat%bx(i-1,j)%z)/dx**2
   end do
  end do
  
 end function
 
 
 function sGyy(dat)
  use global
  use params
  use sync
  
  implicit none
 
  type(var), intent(inout) :: dat
  real(kind=8) :: sGyy(dat%nx,dat%ny)
  integer :: i,j
  
  call end_sync(dat)
  do j=1,dat%ny
   do i=1,dat%nx
    sGyy(i,j)=(dat%uy(i,j+1)%z-2.0d0*dat%z(i,j)+dat%by(i,j-1)%z)/dy**2
   end do
  end do
  
 end function
 
 function Gx(dat)
  use params
  
  implicit none
  
  real(kind=8), intent(in) :: dat(:,:)
  real(kind=8) :: Gx(ubound(dat,1)-1,ubound(dat,2))
  integer :: i,j
  
  do j=1,ubound(dat,2)
   do i=1,ubound(dat,1)-1
    Gx(i,j)=(dat(i+1,j)-dat(i,j))/dx
   end do
  end do
  
 end function
 
 function Gy(dat)
  use params
  
  implicit none
  
  real(kind=8), intent(in) :: dat(:,:)
  real(kind=8) :: Gy(ubound(dat,1),ubound(dat,2)-1)
  integer :: i,j
  
  do j=1,ubound(dat,2)-1
   do i=1,ubound(dat,1)
    Gy(i,j)=(dat(i,j+1)-dat(i,j))/dy
   end do
  end do
  
 end function
 
 function Ax(dat)
  
  implicit none
  
  real(kind=8), intent(in) :: dat(:,:)
  real(kind=8) :: Ax(ubound(dat,1)-1,ubound(dat,2))
  integer :: i,j
  
  do j=1,ubound(dat,2)
   do i=1,ubound(dat,1)-1
    Ax(i,j)=0.5d0*(dat(i+1,j)+dat(i,j))
   end do
  end do
  
 end function
 
 function Ay(dat)
  
  implicit none
  
  real(kind=8), intent(in) :: dat(:,:)
  real(kind=8) :: Ay(ubound(dat,1),ubound(dat,2)-1)
  integer :: i,j
  
  do j=1,ubound(dat,2)-1
   do i=1,ubound(dat,1)
    Ay(i,j)=0.5d0*(dat(i,j+1)+dat(i,j))
   end do
  end do
  
 end function
 
!-------------------------------------------------------------------------!

! FUNCTIONS

!---------------------------------------------------------------------------!


 subroutine thickness_correct(h_tmp,h,n,s)
  use global
  use params
  use grid
 
  type(var), intent(inout) :: h_tmp(nz), h(nz)
  type(var), intent(in) :: s   
  integer, intent(in) :: n
  real(kind=8) :: depth(ubound(s%z,1),ubound(s%z,2))
    
  depth=s%z
  do k=1,nz-1
   depth=depth+h_tmp(k)%z
  end do
  
!   if (n >= 3) then
!    h(nz)%tend0=h(nz)%tend0+(12.0d0/23.0d0)*(-depth-h_tmp(nz)%z)/dt
!   else
!    if (n == 2) then
!     h(nz)%tend0=h(nz)%tend0+(2.0d0/3.0d0)*(-depth-h_tmp(nz)%z)/dt
!    else
!     if (n == 1) then
!      h(nz)%tend1=h(nz)%tend1+(-depth-h_tmp(nz)%z)/dt
!     else
!      if (n ==0) then
!       h(nz)%tend0=h(nz)%tend0+2.0d0*(-depth-h_tmp(nz)%z)/dt
!      end if
!     end if
!    end if
!   end if
  
  h_tmp(nz)%z=-depth
  
 end subroutine
 
 
 
 
 subroutine calc_depth(thickness)
  use sync
  use global
  use variables
  use parallel
  use params


  implicit none
 
  type(var), intent(in) :: thickness(nz)    
  integer :: k, j, i
 
  
  depth(0)%z = s%z
  do k=1,nz
   depth(k)%z=(depth(k-1)%z+thickness(k)%z)
  end do
  
  
  return
  
 end subroutine
 
 
 subroutine montpot()
  use sync
  use global
  use variables
  use params
  
  implicit none
 
                             
  integer :: k
 
  
#ifdef ALLOW_RIGID_LID
  mont(nz)%z=0.0d0
#else
  mont(nz)%z=ngp(nz)*depth(nz)%z
#endif
  do k=nz-1,1,-1
   mont(k)%z=mont(k+1)%z+ngp(k)*depth(k)%z
  end do
  
  do k=1,nz
   call start_sync(mont(k))
  end do
  
  
  return
  
 end subroutine montpot
  
 
 
 
 
!--------------------------------------------------------------------------!

! TENDANCIES

!--------------------------------------------------------------------------!


 
 subroutine tend_h(n)
  use global
  use params
  use variables
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => h(i)%tend1
    else
     tend => h(i)%tend0
    end if
    tend=-Gx(u(i)%z*sAx(h(i)))-Gy(v(i)%z*sAy(h(i)))
   end do
   
   nullify(tend)
  
  
 end subroutine
 
 subroutine tend_u(n)
  use global
  use params
  use variables
  use parallel
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => u(i)%tend1
    else
     tend => u(i)%tend0
    end if
    tend=0.0d0    &
     +smagu(i)%z  &
     +Ay(fz%z+zeta(i)%z)   &
     *sAxy(v(i))    &
     -sGx(mont(i))    &
     -sGx(ke(i))     &
     +0.0d0
    if (i == nz) tend = tend + utau%z*merge(0.0d0,1.0d0/sAx(h(nz)),  &
         (sAx(h(nz)) == 0.0d0)) 
    if (i == 1) tend = tend - bfricu%z*merge(0.0d0,   &
          1.0d0/sAx(h(1)),(sAx(h(1)) == 0.0d0))
   end do
   
   nullify(tend)
   
  
 end subroutine
 
 subroutine tend_v(n)
  use global
  use params
  use variables
  use parallel
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => v(i)%tend1
    else
     tend => v(i)%tend0
    end if
    tend=0.0d0    &
     +smagv(i)%z   &
     -Ax(fz%z+zeta(i)%z)   &
     *sAyx(u(i))    &
     -sGy(mont(i))   &
     -sGy(ke(i))    &
     +0.0d0
    if (i == nz) tend = tend + vtau%z*merge(0.0d0,1.0d0/sAy(h(nz)),  &
         (sAy(h(nz)) == 0.0d0))
    if (i == 1) tend = tend - bfricv%z*merge(0.0d0,   &
          1.0d0/sAy(h(1)),(sAy(h(1)) == 0.0d0))
  end do
   
   nullify(tend)
  
  
 end subroutine
 
 
!--------------------------------------------------------------------------!

! TIMESTEPPING

!--------------------------------------------------------------------------!
 
 
 pure function fe(x)
  use params
  use global
  
  implicit none
 
  type(var), intent(in) :: x
  real (kind=8) :: fe(x%nx,x%ny)
  
  fe=x%z + 0.5d0*dt*x%tend0
  
  return
  
 end function fe
 
 
 
 
 
  

 
 pure function rk2(x)
  use params
  use global
  
  implicit none
  
  type(var), intent(in) :: x
  real (kind=8) :: rk2(x%nx,x%ny)
 
  rk2 = x%z + dt*x%tend1
  
  return
  
 end function rk2
 

 
 
 
 
 
 pure function ab2(x)
  use params
  use global
  
  implicit none
 
  type(var), intent(in) :: x
  real (kind=8) :: ab2(x%nx,x%ny)
  
 
  ab2=x%z+dt*(1.5d0*x%tend0-0.5d0*x%tend1)
 
  return
  
 end function ab2
 

 
 
 
 
 pure function ab3(x)
  use params
  use global
  
  implicit none
 
  type(var), intent(in) :: x
  real (kind=8) :: ab3(x%nx,x%ny)
 
   ab3=x%z+dt*((23.0d0/12.0d0)*x%tend0-(4.0d0/3.0d0)*x%tend1+                  &
       (5.0d0/12.0d0)*x%tend2)
 
  return
  
 end function ab3
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
end module
