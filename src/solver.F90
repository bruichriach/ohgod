module solver_variables
 use global


 type(var) :: pres, p, r
 type(var) :: thavx,thavy
 type(var) :: inthavx,inthavy
 type(var) :: inty
 type(var), allocatable :: y(:)
 real (kind=8), allocatable, dimension(:,:) :: ap, z

end module

module solver
 
 
 implicit none
 
 
 
 contains
 
 subroutine set_solver()
  use solver_variables
  use variables
  use params
  use operate
  use grid
  use sync
  use params
  
  implicit none
  
  integer :: i
 
  allocate(y(nz))
  call create_var(pres,0,0,.false.)
  call create_var(p,0,0,.false.)
  call create_var(r,0,0,.false.)
  do i=1,nz
   call create_hparam(y(i))
  end do
  call create_hparam(inty)
  call create_uparam(thavx)
  call create_vparam(thavy)
  call create_uparam(inthavx)
  call create_vparam(inthavy)
  
  
  pres%z=0.0d0
  call start_sync(pres)
  r%z=0.0d0
  call start_sync(r)
  p%z=0.0d0
  call start_sync(p)
  
  
  
  
  thavx%z=-sAx(s)
  thavy%z=-sAy(s)
       
   
  inthavx%z=merge(0.0d0,1.0d0/thavx%z,(thavx%z == 0.0d0))
  inthavy%z=merge(0.0d0,1.0d0/thavy%z,(thavy%z == 0.0d0))
  
  
 end subroutine
 
  
 
 
 
 subroutine surfpressure(ierr)
  use solver_variables
  use parallel
  use operate
  use sync
 
 
  implicit none

  integer :: n
  integer, intent(out) :: ierr
  real (kind=8) :: alpha, rdotold, rdotnew, rdotfirst, rdotfrac
  real (kind=8) :: pdotap, p_0
  real (kind=8) :: pdotap_tmp, rdotnew_tmp, p_0_tmp
   
  
   
  ierr = 0
  p_0=0
 ! p_0_tmp=sum(pres%z(1:pres%nx,1:pres%ny))
 ! call real_allsum(p_0_tmp, p_0)
 ! p_0=p_0/real(mx*my)
  
  
  
 ! pres%z%z=0.0d0!pres%z%z-p_0
 ! r%z%z=0.0d0
 ! p%z%z=0.0d0
  
  ap=Gx(thavx%z*sGx(pres))+   &
            Gy(thavy%z*sGy(pres))
            
  
 
 
  
  r%z =  inty%z
  
  call start_sync(r)
  
   
  rdotfrac = 1.0d-16
   
  
  
  z=Ax(inthavx%z*sAx(r))+Ay(inthavy%z*sAy(r))
  
   
 
    
  rdotnew_tmp = sum(z*r%z,r%p%core)
  call real_allsum(rdotnew_tmp, rdotold)
  rdotfirst = rdotfrac*rdotold

  !if (proc_name == proc_master) print *, rdotold


  r%z =  inty%z - ap

  call start_sync(r)




  z=Ax(inthavx%z*sAx(r))+Ay(inthavy%z*sAy(r))


  p%z = z

  call start_sync(p)

  rdotnew_tmp = sum(z*r%z,r%p%core)
  call real_allsum(rdotnew_tmp, rdotold)



  !if (proc_name == proc_master) print *, rdotold
 
 
  do n = 1, 1000
  
     
   if (abs(rdotold) /= 0.0) then
   
      
   ap=Gx(thavx%z*sGx(p))+Gy(thavy%z*sGy(p))
  
   
   
   pdotap_tmp=sum(p%z*ap,p%p%core)
   call real_allsum(pdotap_tmp, pdotap)
    
   
   alpha = rdotold/pdotap
   
     
   
     
   r%z=r%z-alpha*ap
   
   call start_sync(r)
   
   
  
   pres%z=pres%z+alpha*p%z
   
   
   
   
   
   z=Ax(inthavx%z*sAx(r))+Ay(inthavy%z*sAy(r))
   
   
      
   rdotnew_tmp = sum(z*r%z,r%p%core)
   call real_allsum(rdotnew_tmp, rdotnew)
        
        
  
   p%z=z+(rdotnew/rdotold)*p%z
   
   call start_sync(p)
    
   
   
   rdotold=rdotnew
 !  if (proc_name == proc_master) print *, n, rdotold, pdotap


  

   if (abs(rdotold) <=  rdotfirst) then
  
 
   !if (proc_name == ens_master) print *, ens_name, n, rdotold
    call start_sync(pres)
    
    
    return
   end if
   
   else
   
    if (proc_name == ens_master) then
     print *, 'exact solve'
     print *, ens_name
    end if
   
    call start_sync(pres)

    
    return
   
   end if
   
  end do
   
  if (proc_name == ens_master) then
   print *, 'no convergence'
   print *, rdotold, ens_name
  end if
  ierr=1
   
    call start_sync(pres)
   
  return
 
 end subroutine
 
 subroutine prescorrection(n)
  use solver_variables
  use global
  use operate
  use sync
  use params
  use variables
 
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tendu(:,:), tendv(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tendu => u(i)%tend1
     tendv => v(i)%tend1
    else
     tendu => u(i)%tend0
     tendv => v(i)%tend0
    end if
    tendu=tendu+sGx(pres)
    tendv=tendv+sGy(pres)
   end do
  
  
  
  return
 
 end subroutine
    
end module
 

