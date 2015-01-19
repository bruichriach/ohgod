program main
#include "include.h"
 use system
 use parallel
 use grid
 use variables
 use operate
 use sync
 use solver_variables
 use solver
 use test
 use writeout_variables
 use writeout
 use llist_ops
#ifdef DO_TIME_AVERAGE
 use timeav
#endif

 implicit none
 
 
 integer :: i, j, k, n, l, m

 call mpi_init(stat)
 
 call makeseed()
 
 
! parallelisation

 proc_master = 0
 call mpi_comm_rank(mpi_comm_world, proc_name, stat);
 call mpi_comm_size(mpi_comm_world, proc_num, stat);
 
 ens_images = proc_num/ens_num
 ens_name = (proc_name)/ens_images
 ens_master = (ens_name)*ens_images
 
 call create_grid()
  
 call set_write_grid()
 
 call random_number(null_field)
 print *, null_field
 
     
 call mpi_barrier(mpi_comm_world,stat)
 if (proc_name == proc_master) then
  print *, 'allocating arrays...          time:', cputime()
 end if
 call mpi_barrier(mpi_comm_world,stat)
 
 call create_var(s,0,0,.false.)
 do i=0,nz
  call create_hparam(depth(i))
  call create_var(true_p(i),0,0,.false.)
  call create_var(true_d(i),0,0,.false.)
 end do
 do i=1,nz
  call create_var(h(i),0,0,.false.)
  call timestepped(h(i))
  call create_var(h_tmp(i),0,0,.false.)
  call create_var(u(i),1,0,.true.)
  call timestepped(u(i))
  call create_uparam(u_tmp(i))
  call create_var(v(i),0,1,.true.)
  call timestepped(v(i))
  call create_vparam(v_tmp(i))
  
  call create_var(h_u(i),1,0,.false.)
  call create_var(h_v(i),0,1,.false.)
  call create_zparam(h_z(i))
  
  call create_var(ke(i),0,0,.false.)
  call create_hparam(ape(i))
  call create_zparam(zeta(i))
  
  call create_var(lapu(i),1,0,.true.)
  call create_var(lapv(i),0,1,.true.)
  
  call create_var(smag(i),0,0,.false.)
  call create_hparam(tension(i))
  call create_zparam(strain(i))
  
  call create_uparam(smagu(i))
  call create_vparam(smagv(i))
  
  
  call create_hparam(q_h(i))
  call create_zparam(q_z(i))
  
  
  call create_var(mont(i),0,0,.false.)
 end do
 
 
 call create_uparam(bfricu)
 call create_vparam(bfricv)
  
 
 call create_uparam(utau)
 call create_vparam(vtau)
  
  call create_zparam(fz)
  
!--------------------------------------------------------------------------------------------------------------!

! INITIALISE PARAMETERS

!--------------------------------------------------------------------------------------------------------------!


 
 ! -------------------------------------------------------------------------------------------------------------------------------- !

 ! topography
 
  
 
  do j=1,s%ny
   do i=1,s%nx
    if (s%p(i,j)%core) then
     s%z(i,j)=-1.0d0+(1.0d0/4.0d0)*   &
        (max(1.0d0-(y_dist(s%p(i,j),y0/2.0d0)/(y0/4.0d0))**2,0.0d0) - &
        max(1.0d0-(y_dist(s%p(i,j),0.0d0)/(y0/4.0d0))**2,0.0d0)) 
    else
     s%z(i,j)=0.0d0
    end if
   end do
  end do 
  
  call start_sync(s)


  fz%z=1.0d0
  
  
  depth(0)%exp=s%z
  
  depth(0)%z=depth(0)%exp
  do i=1,nz
   depth(i)%exp=((3.0d0/4.0d0)*dble(i-nz))/dble(nz)
   depth(i)%z=depth(i)%exp
  end do
  
 
  do i=1,nz
   h(i)%z=depth(i)%z-depth(i-1)%z
   h(i)%exp=depth(i)%exp-depth(i-1)%exp
  end do
  
  
  
#ifdef ALLOW_RIGID_LID
  mont(nz)%exp=0.0d0
#else
  mont(nz)%exp=ngp(nz)*depth(nz)%exp
#endif
  do k=nz-1,1,-1
   mont(k)%exp=mont(k+1)%exp+ngp(k)*depth(k)%exp
  end do
  
  
  do i=1,nz
   u(i)%z=0.0d0
   v(i)%z=0.0d0
  end do
  
  do i=1,nz
   call read_var(h(i),i,'h','old',0,0)
   call read_var(u(i),i,'u','old',1,0)
   call read_var(v(i),i,'v','old',0,1)
  end do

 
  do i=1,nz
   call start_sync(h(i))
   call start_sync(u(i))
   call start_sync(v(i))
  end do
    
    
  call calc_depth(h)
  
  call montpot()
  
  do i=1,nz
   ke(i)%z=0.5d0*(Ax(u(i)%z**2)+Ay(v(i)%z**2))
   h_u(i)%z=sAx(h(i))
   h_v(i)%z=sAy(h(i))
   call start_sync(h_u(i))
   call start_sync(h_v(i))
   call start_sync(ke(i))
  end do
  
   
  call set_solver()
  
  
!--------------------------------------------------------------------------!

! INIT WRITEOUTS

!--------------------------------------------------------------------------!
  
  
  call mpi_barrier(mpi_comm_world,stat)
  if (proc_name == proc_master) then
   print *, 'setting up output arrays...   time:', cputime()
  end if
  call mpi_barrier(mpi_comm_world,stat)
  
#ifdef DO_TIME_AVERAGE
  call init_diagnostics()
#endif
  
  
  call init_writeouts(utau%z,utau%p,utau_out,ugrid_out)
  call init_writeouts(vtau%z,vtau%p,vtau_out,vgrid_out)
  do i=0,nz
   call init_writeouts(depth(i)%z,depth(i)%p,depth_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(h(i)%z,h(i)%p,h_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(u(i)%z,u(i)%p,u_out(i),ugrid_out)
  end do
  do i=1,nz
   call init_writeouts(v(i)%z,v(i)%p,v_out(i),vgrid_out)
  end do
  do i=1,nz
   call init_writeouts(h(i)%tend0,h(i)%p,tendh_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(u(i)%tend0,u(i)%p,tendu_out(i),ugrid_out)
  end do
  do i=1,nz
   call init_writeouts(v(i)%tend0,v(i)%p,tendv_out(i),vgrid_out)
  end do
  do i=1,nz
   call init_writeouts(smagu(i)%z,smagu(i)%p,smagu_out(i),ugrid_out)
  end do
  do i=1,nz
   call init_writeouts(smagv(i)%z,smagv(i)%p,smagv_out(i),vgrid_out)
  end do
  do i=1,nz
   call init_writeouts(q_h(i)%z,q_h(i)%p,q_h_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(q_z(i)%z,q_z(i)%p,q_z_out(i),zgrid_out)
  end do
  do i=1,nz
   call init_writeouts(ke(i)%z,ke(i)%p,ke_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(ape(i)%z,ape(i)%p,ape_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(mont(i)%z,mont(i)%p,mont_out(i),hgrid_out)
  end do
  do i=1,nz
   call init_writeouts(zeta(i)%z,zeta(i)%p,zeta_out(i),zgrid_out)
  end do
  call init_writeouts(pres%z,pres%p,pres_out,hgrid_out)
  call init_writeouts(inty%z,inty%p,div_u_out,hgrid_out)
  
  
  call mpi_barrier(mpi_comm_world,stat)
  if (proc_name == proc_master) then
   print *, 'starting initial file write   time:', cputime()
  end if
  call mpi_barrier(mpi_comm_world,stat)
  
  call write_config(s,hgrid_out,'s')
  call write_config(fz,zgrid_out,'f')
  
  call write_main('in')
  
#ifdef ALLOW_STOCHASTIC_WIND
  call input_llist(stochwind,'stochwind')
#endif
 
  
!--------------------------------------------------------------------------------------------------------------!

! Start RUN

!--------------------------------------------------------------------------------------------------------------!

  
 do n=0,nstop
 
 
  if (proc_name == proc_master) then
   if (update(2.0d0)) then
    write(*,"(a33,a8,a26,a1)", ADVANCE = "NO") print_time(), ' Model: ', &
      modeltime(n), CHAR(13)
   end if
  end if 
  
  do i=1,nz
   h_z(i)%z=0.5d0*(sAx(h_v(i))+sAy(h_u(i)))
  end do
  
  if (n == 0) then
   do i=1,nz
    h(i)%tmp=h(i)%z
    u(i)%tmp=u(i)%z
    v(i)%tmp=v(i)%z
   end do
  end if
  
  call tend_h(n)
  
  do i=1,nz
   lapu(i)%z=sGxx(u(i))+Gy(sGy(u(i)))
   lapu(i)%z=merge(0.0d0,   &
      lapu(i)%z,lapu(i)%p%flux)
   call start_sync(lapu(i))
   lapv(i)%z=Gx(sGx(v(i)))+sGyy(v(i))
   lapv(i)%z=merge(0.0d0,   &
      lapv(i)%z,lapv(i)%p%flux)
   call start_sync(lapv(i))
  end do
  
 
#ifdef ALLOW_STOCHASTIC_WIND
  call set_wind(n)
#endif
  
  do i=1,nz
   strain(i)%z=sGx(v(i))+sGy(u(i))
   tension(i)%z=Gx(u(i)%z)-Gy(v(i)%z)
   smag(i)%z=max(dx,dy)**4*cvar*sqrt(tension(i)%z**2+Ax(Ay(strain(i)%z**2)))
   call start_sync(smag(i))
  end do
  
  
  
  do i=1,nz
   zeta(i)%z=sGx(v(i))-sGy(u(i))
  end do
  do i=1,nz
   smagu(i)%z=-sAx(smag(i))*(sGxx(lapu(i))+Gy(sGy(lapu(i))))
   smagu(i)%z=merge(0.0d0,   &
      smagu(i)%z,smagu(i)%p%flux)
  end do
  bfricu%z=bf*(sAx(ke(1))/(4.0d-4))*u(1)%z
  call tend_u(n)
  do i=1,nz
   smagv(i)%z=-sAy(smag(i))*(Gx(sGx(lapv(i)))+sGyy(lapv(i)))
   smagv(i)%z=merge(0.0d0,   &
      smagv(i)%z,smagv(i)%p%flux)
  end do
  bfricv%z=bf*(sAy(ke(1))/(4.0d-4))*v(1)%z
  call tend_v(n)


#ifdef ALLOW_RIGID_LID

  if (n >= 3) then
   do i=1,nz
    h_tmp(i)%z=ab3(h(i))
   end do
  else
   if (n == 2) then
    do i=1,nz
     h_tmp(i)%z=ab2(h(i))
    end do
   else
    if (n == 1) then
     do i=1,nz
      h(i)%z=h(i)%tmp
      h_tmp(i)%z=rk2(h(i))
     end do
    else
     if (n ==0) then
      do i=1,nz
       h_tmp(i)%z=fe(h(i))
      end do
     end if
    end if
   end if
  end if
  call thickness_correct(h_tmp,h,n,s)
   
  do i=1,nz
   call start_sync(h_tmp(i))
  end do

  if (n >= 3) then
   do i=1,nz
    u_tmp(i)%z=ab3(u(i))
    v_tmp(i)%z=ab3(v(i))
   end do
  else
   if (n == 2) then
    do i=1,nz
     u_tmp(i)%z=ab2(u(i))
     v_tmp(i)%z=ab2(v(i))
    end do
   else
    if (n == 1) then
     do i=1,nz
      u(i)%z=u(i)%tmp
      v(i)%z=v(i)%tmp
      u_tmp(i)%z=rk2(u(i))
      v_tmp(i)%z=rk2(v(i))
     end do
    else
     if (n ==0) then
      do i=1,nz
       u_tmp(i)%z=fe(u(i))
       v_tmp(i)%z=fe(v(i))
      end do
     end if
    end if
   end if
  end if
   
  
  
  thavx%z=-sAx(s)
  thavy%z=-sAy(s)
   
  inthavx%z=merge(0.0d0,1.0d0/thavx%z,(thavx%z == 0.0d0))
  inthavy%z=merge(0.0d0,1.0d0/thavy%z,(thavy%z == 0.0d0))
  
  
   inty%z=0.0d0
   do i=1,nz
    y(i)%z=-Gx(u_tmp(i)%z*sAx(h_tmp(i)))   &
           -Gy(v_tmp(i)%z*sAy(h_tmp(i)))
 !      inty%z=inty%z+y(i)%z
    if (n >= 3) then
     inty%z=inty%z+(12.0d0/23.0d0)*y(i)%z/dt
    else
     if (n == 2) then
      inty%z=inty%z+(2.0d0/3.0d0)*y(i)%z/dt
     else
      if (n == 1) then
       inty%z=inty%z+y(i)%z/dt
      else
       if (n ==0) then
        inty%z=inty%z+2.0d0*y(i)%z/dt
       end if
      end if
     end if
    end if
   end do
  
  
  
   
   call surfpressure(ierr)
   call prescorrection(n)

#endif


#ifdef ALLOW_RIGID_LID
  do i=1,nz
   ape(i)%z=(h(i)%z-h(i)%exp)*(-pres%z+mont(i)%z-mont(i)%exp)
  end do
#else
  do i=1,nz
   ape(i)%z=(h(i)%z-h(i)%exp)*(mont(i)%z-mont(i)%exp)
  end do
#endif

    
#ifdef DO_TIME_AVERAGE
  do m=1,2*barno
   do l=-nsteps,nsteps
    if (mod(n-(m-nbarno)*wstep-l*nstep,2*mbarno*wstep) == 0) then
     call integerate_flux(m)
    end if
   end do
  end do
  
  
  do m=1,2*barno
   if (mod(n-nstep*nsteps-(m-nbarno)*wstep,2*mbarno*wstep) == 0) then
    call write_diagnostics(m,n)
   end if
  end do
    
  
  do m=1,2*barno
   if (mod(n+nstep*nsteps-(m-nbarno)*wstep,2*mbarno*wstep) == 0) then
    call reset_flux(m)
    call integerate_flux(m)
   end if
  end do
#endif

  
  if (mod(n,wstep) == 0) then
   if (n /= 0) then
    call do_write(n/wstep-1)
    if (proc_name == ens_master) then
     print *
     print *, 'writeout done.     time:', cputime(), 'ens_name:', ens_name
    end if
   end if
   call write()
#ifdef ALLOW_STOCHASTIC_WIND
   call writeout_llist(n/wstep,stochwind,'stochwind')
#endif
  end if



  if (n >= 3) then
   do i=1,nz
#ifdef ALLOW_RIGID_LID
    h(i)%z=h_tmp(i)%z
#else
    h(i)%z=ab3(h(i))
#endif
   end do
  else
   if (n == 2) then
    do i=1,nz
#ifdef ALLOW_RIGID_LID
    h(i)%z=h_tmp(i)%z
#else
     h(i)%z=ab2(h(i))
#endif
    end do
   else
    if (n == 1) then
     do i=1,nz
#ifdef ALLOW_RIGID_LID
    h(i)%z=h_tmp(i)%z
#else
      h(i)%z=h(i)%tmp
      h(i)%z=rk2(h(i))
#endif
     end do
    else
     if (n ==0) then
      do i=1,nz
#ifdef ALLOW_RIGID_LID
    h(i)%z=h_tmp(i)%z
#else
       h(i)%z=fe(h(i))
#endif
      end do
     end if
    end if
   end if
  end if
   
  do i=1,nz
   call start_sync(h(i))
  end do
  
  
    
  call calc_depth(h)
  
  call montpot()
  
  
  
  
  
  
  
  
 
  do i=1,nz
   if (n >= 3) then
     u(i)%z=ab3(u(i))
     v(i)%z=ab3(v(i))
   else
    if (n == 2) then
     u(i)%z=ab2(u(i))
     v(i)%z=ab2(v(i))
    else
     if (n == 1) then
#ifndef ALLOW_RIGID_LID
      u(i)%z=u(i)%tmp
      v(i)%z=v(i)%tmp
#endif
      u(i)%z=rk2(u(i))
      v(i)%z=rk2(v(i))
     else
      if (n ==0) then
       u(i)%z=fe(u(i))
       v(i)%z=fe(v(i))
      end if
     end if
    end if
   end if
   call start_sync(u(i))
   call start_sync(v(i))
   ke(i)%z=0.5d0*(Ax(u(i)%z**2)+Ay(v(i)%z**2))
   call start_sync(ke(i))
   h_u(i)%z=sAx(h(i))
   h_v(i)%z=sAy(h(i))
   call start_sync(h_u(i))
   call start_sync(h_v(i))
  end do
  
   
     
  
  if (n /= 0) then
   do i=1,nz
    call stepforward(h(i))
    call stepforward(u(i))
    call stepforward(v(i))
   end do
  end if
  
 end do
 
  call do_write(nstop/wstep)
  
  
#ifdef DO_TIME_AVERAGE
  call write_diagnostics(1,(nstop/wstep+1)*wstep)
#endif
  
  call write_main('out')
  
  
 
#ifdef ALLOW_STOCHASTIC_WIND
 call output_llist(stochwind,'stochwind')
#endif
 
 call mpi_finalize(stat)
 
 
 
 
 
 
 
 
end program
