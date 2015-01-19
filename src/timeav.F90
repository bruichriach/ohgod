module timeav_variables

#include "include.h"
#ifdef DO_TIME_AVERAGE

 use writeout_variables
 use global
 use params
 
 implicit none
  
 type average
  logical :: h, u, v, z
  type(var) :: zf_h
  type(var) :: zf_u
  type(var) :: zf_v
  type(var) :: zf_z
  type(var) :: z_h(2*barno)
  type(var) :: z_u(2*barno)
  type(var) :: z_v(2*barno)
  type(var) :: z_z(2*barno)
 end type
  
 type fin_var
  type(var) :: z
  type(out_var) :: out
 end type
  
 type(average) :: sh(nz)
 type(average) :: sape(nz)
 type(average) :: sm(nz)
 type(average) :: shmx(nz), shmy(nz)
 type(average) :: smx(nz), smy(nz)
 type(average) :: shq(nz)
 type(average) :: shqq(nz)
 type(average) :: shuq(nz), shvq(nz)
 type(average) :: shu(nz), shv(nz)
 type(average) :: shuu(nz), shvv(nz), shuv(nz)
 type(average) :: sdd(0:nz)
 type(average) :: sd(0:nz)
 type(average) :: shlapu(nz), shlapv(nz)
 type(average) :: shlap2u(nz), shlap2v(nz)
 type(average) :: shsmagu(nz), shsmagv(nz)
 type(average) :: sutau, svtau
 type(average) :: sbfricu, sbfricv
 type(average) :: shtendu(nz), shtendv(nz)
 type(average) :: stendh(nz)
 
 type(fin_var) :: b_u(nz), b_v(nz)
 type(fin_var) :: a_h(nz), a_m(nz)
 type(fin_var) :: a_d(0:nz)
 type(fin_var) :: c_q(nz), b_q(nz)
 type(fin_var) :: b_tendu(nz), b_tendv(nz)
 type(fin_var) :: a_tendh(nz)
 type(fin_var) :: b_smagu(nz), b_smagv(nz)
 
 type(fin_var) :: a_utau, a_vtau
 type(fin_var) :: a_bfricu, a_bfricv
 
 type(fin_var) :: bq_uq(nz), bq_vq(nz)
 type(fin_var) :: bq_qq(nz)
 
 type(var) :: d_huv(nz)
 type(var) :: d_huu(nz)
 type(var) :: d_hvv(nz)
 type(var) :: d_dd(0:nz)
 type(var) :: d_hm(nz)
 
 type(fin_var) :: em(nz), en(nz), en_h(nz), ek(nz), ep(nz), ep2(nz), er(nz), es(nz)
 type(fin_var) :: eu(nz), ev(nz), ef(nz)
 type(fin_var) :: ereynu(nz), eformu(nz)
 type(fin_var) :: ereynv(nz), eformv(nz)
 

end module

module timeav
  
 implicit none
 
 contains
 
 subroutine init_avarage(dat,h,u,v,z,uvec,vvec)
  use timeav_variables
  use grid
  
  type(average), intent(inout) :: dat
  integer :: j
  logical, intent(in) :: uvec, vvec
  logical, intent(in) :: h, u, v, z
  
  dat%h=h
  dat%u=u
  dat%v=v
  dat%z=z
  
  if (dat%h) then
   do j=1,2*barno
    call create_hparam(dat%z_h(j))
   end do
   call create_var(dat%zf_h,0,0,.false.)
  end if
  if (dat%u) then
   do j=1,2*barno
    call create_uparam(dat%z_u(j))
   end do
   call create_var(dat%zf_u,1,0,uvec)
  end if
  if (dat%v) then
   do j=1,2*barno
    call create_vparam(dat%z_v(j))
   end do
   call create_var(dat%zf_v,0,1,vvec)
  end if
  if (dat%z) then
   do j=1,2*barno
    call create_zparam(dat%z_z(j))
   end do
   call create_zparam(dat%zf_z)
  end if
   
 end subroutine
  
 
 subroutine init_diagnostics()
  use parallel
  use timeav_variables
  use writeout
  use grid
  
  implicit none
  
  integer :: i, j
  
  do i=1,nz
   call init_avarage(sh(i),.true.,.true.,.true.,.true.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(sape(i),.true.,.false.,.false.,.false.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(sm(i),.true.,.false.,.false.,.false.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(shmx(i),.false.,.true.,.false.,.false.,   &
          .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shmy(i),.false.,.false.,.true.,.false.,   &
         .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(smx(i),.false.,.true.,.false.,.false.,   &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(smy(i),.false.,.false.,.true.,.false.,   &
        .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shq(i),.true.,.true.,.true.,.true.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(shqq(i),.true.,.false.,.false.,.true.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(shuq(i),.false.,.true.,.false.,.false.,  &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shvq(i),.false.,.false.,.true.,.false.,  &
         .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shu(i),.true.,.true.,.true.,.true.,   &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shv(i),.true.,.true.,.true.,.true.,   &
         .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shuu(i),.true.,.false.,.false.,.false.,   &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shvv(i),.true.,.false.,.false.,.false.,    &
         .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shuv(i),.true.,.false.,.false.,.true.,    &
         .true.,.true.)
  end do
  do i=1,nz
   call init_avarage(stendh(i),.true.,.false.,.false.,.false.,  &
         .false.,.false.)
  end do
  do i=1,nz
   call init_avarage(shtendu(i),.false.,.true.,.false.,.false.,   &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shtendv(i),.false.,.false.,.true.,.false.,   &
         .false.,.true.)
  end do
  
  
  do i=0,nz
   call init_avarage(sdd(i),.true.,.false.,.false.,.false.,   &
         .false.,.false.)
  end do
  do i=0,nz
   call init_avarage(sd(i),.true.,.true.,.true.,.false.,    &
          .false.,.false.)
  end do
  
  
  do i=1,nz
   call init_avarage(shlapu(i),.false.,.true.,.false.,.false.,   &
            .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shlapv(i),.false.,.false.,.true.,.false.,    &
          .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shlap2u(i),.false.,.true.,.false.,.false.,   &
         .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shlap2v(i),.false.,.false.,.true.,.false.,    &
           .false.,.true.)
  end do
  do i=1,nz
   call init_avarage(shsmagu(i),.false.,.true.,.false.,.false.,   &
           .true.,.false.)
  end do
  do i=1,nz
   call init_avarage(shsmagv(i),.false.,.false.,.true.,.false.,    &
           .false.,.true.)
  end do
  

  call init_avarage(sutau,.false.,.true.,.false.,.false.,   &
         .true.,.false.)
  
  call init_avarage(svtau,.false.,.false.,.true.,.false.,   &
         .false.,.true.)

  call init_avarage(sbfricu,.false.,.true.,.false.,.false.,   &
         .true.,.false.)
  
  call init_avarage(sbfricv,.false.,.false.,.true.,.false.,   &
         .false.,.true.)
         
         
  do i=0,nz
   call create_hparam(a_d(i)%z)
   call init_writeouts(a_d(i)%z%z,a_d(i)%z%p,a_d(i)%out,hgrid_out)
      
   call create_var(d_dd(i),0,0,.false.)
  end do
         
  do i=1,nz
   call create_hparam(a_h(i)%z)
   call init_writeouts(a_h(i)%z%z,a_h(i)%z%p,a_h(i)%out,hgrid_out)
   
   call create_hparam(a_m(i)%z)
   call init_writeouts(a_m(i)%z%z,a_m(i)%z%p,a_m(i)%out,hgrid_out)
   
   call create_var(b_u(i)%z,1,0,.true.)
   call init_writeouts(b_u(i)%z%z,b_u(i)%z%p,b_u(i)%out,ugrid_out)
   
   call create_var(b_v(i)%z,0,1,.true.)
   call init_writeouts(b_v(i)%z%z,b_v(i)%z%p,b_v(i)%out,vgrid_out)
   
   call create_zparam(c_q(i)%z)
   call init_writeouts(c_q(i)%z%z,c_q(i)%z%p,c_q(i)%out,zgrid_out)
   
   call create_zparam(b_q(i)%z)
   call init_writeouts(b_q(i)%z%z,b_q(i)%z%p,b_q(i)%out,zgrid_out)
   
   call create_zparam(d_huv(i))
   
   call create_var(d_huu(i),0,0,.false.)
   
   call create_var(d_hvv(i),0,0,.false.)
   
   call create_hparam(d_hm(i))
   
   call create_var(eu(i)%z,1,0,.false.)
   call init_writeouts(eu(i)%z%z,eu(i)%z%p,eu(i)%out,ugrid_out)
   
   call create_var(ev(i)%z,0,1,.false.)
   call init_writeouts(ev(i)%z%z,ev(i)%z%p,ev(i)%out,vgrid_out)
   
   call create_zparam(ef(i)%z)
   call init_writeouts(ef(i)%z%z,ef(i)%z%p,ef(i)%out,zgrid_out)
   
   call create_hparam(em(i)%z)
   call init_writeouts(em(i)%z%z,em(i)%z%p,em(i)%out,hgrid_out)
   
   call create_hparam(en_h(i)%z)
   call init_writeouts(en_h(i)%z%z,en_h(i)%z%p,en_h(i)%out,hgrid_out)
   
   call create_zparam(en(i)%z)
   call init_writeouts(en(i)%z%z,en(i)%z%p,en(i)%out,zgrid_out)
   
   call create_uparam(er(i)%z)
   call init_writeouts(er(i)%z%z,er(i)%z%p,er(i)%out,ugrid_out)
   
   call create_vparam(es(i)%z)
   call init_writeouts(es(i)%z%z,es(i)%z%p,es(i)%out,vgrid_out)
   
   call create_hparam(ep(i)%z)
   call init_writeouts(ep(i)%z%z,ep(i)%z%p,ep(i)%out,hgrid_out)
   
   call create_hparam(ep2(i)%z)
   call init_writeouts(ep2(i)%z%z,ep2(i)%z%p,ep2(i)%out,hgrid_out)
   
   call create_hparam(ek(i)%z)
   call init_writeouts(ek(i)%z%z,ek(i)%z%p,ek(i)%out,hgrid_out)
   
   call create_uparam(ereynu(i)%z)
   call init_writeouts(ereynu(i)%z%z,ereynu(i)%z%p,ereynu(i)%out,ugrid_out)
   
   call create_uparam(eformu(i)%z)
   call init_writeouts(eformu(i)%z%z,eformu(i)%z%p,eformu(i)%out,ugrid_out)
      
   call create_vparam(ereynv(i)%z)
   call init_writeouts(ereynv(i)%z%z,ereynv(i)%z%p,ereynv(i)%out,vgrid_out)
   
   call create_vparam(eformv(i)%z)
   call init_writeouts(eformv(i)%z%z,eformv(i)%z%p,eformv(i)%out,vgrid_out)
   
   
   
   call create_var(bq_uq(i)%z,1,0,.true.)
   call init_writeouts(bq_uq(i)%z%z,bq_uq(i)%z%p,bq_uq(i)%out,ugrid_out)
   
   call create_var(bq_vq(i)%z,0,1,.true.)
   call init_writeouts(bq_vq(i)%z%z,bq_vq(i)%z%p,bq_vq(i)%out,vgrid_out)
   
   call create_zparam(bq_qq(i)%z)
   call init_writeouts(bq_qq(i)%z%z,bq_qq(i)%z%p,bq_qq(i)%out,zgrid_out)
   
   
   call create_hparam(a_tendh(i)%z)
   call init_writeouts(a_tendh(i)%z%z,a_tendh(i)%z%p,a_tendh(i)%out,hgrid_out)
   
   call create_uparam(b_tendu(i)%z)
   call init_writeouts(b_tendu(i)%z%z,b_tendu(i)%z%p,b_tendu(i)%out,ugrid_out)
   
   call create_vparam(b_tendv(i)%z)
   call init_writeouts(b_tendv(i)%z%z,b_tendv(i)%z%p,b_tendv(i)%out,vgrid_out)
   
   call create_uparam(b_smagu(i)%z)
   call init_writeouts(b_smagu(i)%z%z,b_smagu(i)%z%p,b_smagu(i)%out,ugrid_out)
   
   call create_vparam(b_smagv(i)%z)
   call init_writeouts(b_smagv(i)%z%z,b_smagv(i)%z%p,b_smagv(i)%out,vgrid_out)
   
  end do
   
   call create_var(a_utau%z,1,0,.true.)
   call init_writeouts(a_utau%z%z,a_utau%z%p,a_utau%out,ugrid_out)
   
   call create_var(a_vtau%z,0,1,.true.)
   call init_writeouts(a_vtau%z%z,a_vtau%z%p,a_vtau%out,vgrid_out)
   
   call create_var(a_bfricu%z,1,0,.true.)
   call init_writeouts(a_bfricu%z%z,a_bfricu%z%p,a_bfricu%out,ugrid_out)
   
   call create_var(a_bfricv%z,0,1,.true.)
   call init_writeouts(a_bfricv%z%z,a_bfricv%z%p,a_bfricv%out,vgrid_out)
  
  
  
 end subroutine
 
 subroutine reset_sum(dat,m)
  use timeav_variables
  
  type(average), intent(inout) :: dat
  integer, intent(in) :: m
  
  if (dat%h) then
   dat%z_h(m)%z=0.0d0
  end if
  if (dat%u) then
   dat%z_u(m)%z=0.0d0
  end if
  if (dat%v) then
   dat%z_v(m)%z=0.0d0
  end if
  if (dat%z) then
   dat%z_z(m)%z=0.0d0
  end if
  
 end subroutine
 
 subroutine reset_flux(m)
  use timeav_variables
  
  implicit none
 
  integer, intent(in) :: m
  integer :: i
  
  do i=1,nz
   call reset_sum(sh(i),m)
   call reset_sum(sape(i),m)
   call reset_sum(sm(i),m)
   call reset_sum(shmx(i),m)
   call reset_sum(shmy(i),m)
   call reset_sum(smx(i),m)
   call reset_sum(smy(i),m)
   call reset_sum(shq(i),m)
   call reset_sum(shqq(i),m)
   call reset_sum(shuq(i),m)
   call reset_sum(shvq(i),m)
   call reset_sum(shu(i),m)
   call reset_sum(shv(i),m)
   call reset_sum(shuu(i),m)
   call reset_sum(shvv(i),m)
   call reset_sum(shuv(i),m)
   call reset_sum(shlapu(i),m)
   call reset_sum(shlapv(i),m)
   call reset_sum(shlap2u(i),m)
   call reset_sum(shlap2v(i),m)
   call reset_sum(shsmagu(i),m)
   call reset_sum(shsmagv(i),m)
   call reset_sum(stendh(i),m)
   call reset_sum(shtendu(i),m)
   call reset_sum(shtendv(i),m)
  end do
  do i=0,nz
   call reset_sum(sd(i),m)
   call reset_sum(sdd(i),m)
  end do
  call reset_sum(sutau,m)
  call reset_sum(svtau,m)
  call reset_sum(sbfricu,m)
  call reset_sum(sbfricv,m)
   
 end subroutine
 
 subroutine integerate_flux(m)
  use timeav_variables
  use variables
  use operate
  use sync
  use solver_variables
  use params
  
  implicit none
  
  integer, intent(in) :: m
  integer :: i, j

#ifdef ALLOW_RIGID_LID
   true_d(nz)%z=-pres%z/ngp(nz)
   true_p(nz)%z=0.0d0
#else
   true_d(nz)%z=depth(nz)%z
   true_p(nz)%z=0.0d0
#endif
  call start_sync(true_d(nz))
  call start_sync(true_p(nz))
  do i=nz-1,0,-1
   true_d(i)%z=depth(i)%z
   true_p(i)%z=0.0d0
   do j=nz,i+1,-1
    true_p(i)%z=true_p(i)%z+(true_d(j)%z-true_d(i)%z)*ngp(j)
   end do
   call start_sync(true_d(i))
   call start_sync(true_p(i))
  end do
  
  do i=1,nz
   sh(i)%z_h(m)%z=sh(i)%z_h(m)%z+h(i)%z
   sh(i)%z_u(m)%z=sh(i)%z_u(m)%z+h_u(i)%z
   sh(i)%z_v(m)%z=sh(i)%z_v(m)%z+h_v(i)%z
   sh(i)%z_z(m)%z=sh(i)%z_z(m)%z+h_z(i)%z
  end do
  
  do i=1,nz
   sape(i)%z_h(m)%z=sape(i)%z_h(m)%z+ape(i)%z
  end do
  
  do i=1,nz
   sm(i)%z_h(m)%z=sm(i)%z_h(m)%z+mont(i)%z-pres%z
  end do
  
  do i=1,nz
   shmx(i)%z_u(m)%z=shmx(i)%z_u(m)%z+sAx(h(i))*(sGx(mont(i))-sGx(pres))
  end do
  
  do i=1,nz
   shmy(i)%z_v(m)%z=shmy(i)%z_v(m)%z+sAy(h(i))*(sGy(mont(i))-sGy(pres))
  end do
  
  do i=1,nz
   smx(i)%z_u(m)%z=smx(i)%z_u(m)%z+sGx(mont(i))-sGx(pres)
  end do
  
  do i=1,nz
   smy(i)%z_v(m)%z=smy(i)%z_v(m)%z+sGy(mont(i))-sGy(pres)
  end do
  
  do i=1,nz
   shq(i)%z_h(m)%z=shq(i)%z_h(m)%z+Ax(Ay(fz%z+zeta(i)%z))
   shq(i)%z_u(m)%z=shq(i)%z_u(m)%z+Ay(fz%z+zeta(i)%z)
   shq(i)%z_v(m)%z=shq(i)%z_v(m)%z+Ax(fz%z+zeta(i)%z)
   shq(i)%z_z(m)%z=shq(i)%z_z(m)%z+fz%z+zeta(i)%z
  end do
  
  do i=1,nz
   shqq(i)%z_h(m)%z=shqq(i)%z_h(m)%z+Ax(Ay((fz%z+zeta(i)%z)**2))*  &
          merge(0.0d0,1.0d0/h(i)%z,(h(i)%z == 0.0d0))
   shqq(i)%z_z(m)%z=shqq(i)%z_z(m)%z+(fz%z+zeta(i)%z)**2*  &
          merge(0.0d0,1.0d0/h_z(i)%z,(h_z(i)%z == 0.0d0))
  end do
  
  do i=1,nz
   shuq(i)%z_u(m)%z=shuq(i)%z_u(m)%z+u(i)%z*Ay(fz%z+zeta(i)%z)
  end do
  
  do i=1,nz
   shvq(i)%z_v(m)%z=shvq(i)%z_v(m)%z+v(i)%z*Ax(fz%z+zeta(i)%z)
  end do
  
  do i=1,nz
   shu(i)%z_h(m)%z=shu(i)%z_h(m)%z+h(i)%z*Ax(u(i)%z)
   shu(i)%z_u(m)%z=shu(i)%z_u(m)%z+h_u(i)%z*u(i)%z
   shu(i)%z_v(m)%z=shu(i)%z_v(m)%z+h_v(i)%z*Ax(sAy(u(i)))
   shu(i)%z_z(m)%z=shu(i)%z_z(m)%z+h_z(i)%z*sAy(u(i))
  end do
  
  do i=1,nz
   shv(i)%z_h(m)%z=shv(i)%z_h(m)%z+h(i)%z*Ay(v(i)%z)
   shv(i)%z_u(m)%z=shv(i)%z_u(m)%z+h_u(i)%z*Ay(sAx(v(i)))
   shv(i)%z_v(m)%z=shv(i)%z_v(m)%z+h_v(i)%z*v(i)%z
   shv(i)%z_z(m)%z=shv(i)%z_z(m)%z+h_z(i)%z*sAx(v(i))
  end do
  
  
  do i=1,nz
   stendh(i)%z_h(m)%z=stendh(i)%z_h(m)%z+h(i)%tend0
  end do
  
  do i=1,nz
   shtendu(i)%z_u(m)%z=shtendu(i)%z_u(m)%z+h_u(i)%z*u(i)%tend0
  end do
  
  do i=1,nz
   shtendv(i)%z_v(m)%z=shtendv(i)%z_v(m)%z+h_v(i)%z*v(i)%tend0
  end do
  
  do i=1,nz
   shuu(i)%z_h(m)%z=shuu(i)%z_h(m)%z+h(i)%z*Ax(u(i)%z**2)
  end do
  
  do i=1,nz
   shvv(i)%z_h(m)%z=shvv(i)%z_h(m)%z+h(i)%z*Ay(v(i)%z**2)
  end do
  
  do i=1,nz
   shuv(i)%z_h(m)%z=shuv(i)%z_h(m)%z+h(i)%z*Ax(u(i)%z)*Ay(v(i)%z)
   shuv(i)%z_z(m)%z=shuv(i)%z_z(m)%z+h_z(i)%z*sAy(u(i))*sAx(v(i))
  end do
  
  do i=0,nz
   sd(i)%z_h(m)%z=sd(i)%z_h(m)%z+true_d(i)%z
   sd(i)%z_u(m)%z=sd(i)%z_u(m)%z+sAx(true_d(i))
   sd(i)%z_v(m)%z=sd(i)%z_v(m)%z+sAy(true_d(i))
  end do
  
  do i=0,nz
   sdd(i)%z_h(m)%z=sdd(i)%z_h(m)%z+true_d(i)%z**2
  end do
  
  do i=1,nz
   shlapu(i)%z_u(m)%z=shlapu(i)%z_u(m)%z+h_u(i)%z*lapu(i)%z
  end do
  
  do i=1,nz
   shlapv(i)%z_v(m)%z=shlapv(i)%z_v(m)%z+h_v(i)%z*lapv(i)%z
  end do
  
  do i=1,nz
   shlap2u(i)%z_u(m)%z=shlap2u(i)%z_u(m)%z+   &
       h_u(i)%z*(sGxx(lapu(i))+Gy(sGy(lapu(i))))
  end do
  
  do i=1,nz
   shlap2v(i)%z_v(m)%z=shlap2v(i)%z_v(m)%z+    &
       h_v(i)%z*(Gx(sGx(lapv(i)))+sGyy(lapv(i)))
  end do
  
  do i=1,nz
   shsmagu(i)%z_u(m)%z=shsmagu(i)%z_u(m)%z+   &
       h_u(i)%z*smagu(i)%z
  end do
  
  do i=1,nz
   shsmagv(i)%z_v(m)%z=shsmagv(i)%z_v(m)%z+    &
       h_v(i)%z*smagv(i)%z
  end do
  
  
   sutau%z_u(m)%z=sutau%z_u(m)%z+utau%z
  
   svtau%z_v(m)%z=svtau%z_v(m)%z+vtau%z
  
   sbfricu%z_u(m)%z=sbfricu%z_u(m)%z+bfricu%z
  
   sbfricv%z_v(m)%z=sbfricv%z_v(m)%z+bfricv%z
  
 end subroutine
  
 subroutine write_diagnostics(m,n)
  use timeav_variables
  use variables
  use writeout
  use sync
  use operate
  use parallel
  use params
  
  implicit none
  
  integer, intent(in) :: m, n
  character(32) :: filename
  integer :: i, num
  logical :: dir_e
  
  
  
  if (n >= wstep + 2*nsteps*nstep ) then
   num = (n-nsteps*nstep)/wstep-1
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
    
   
    do i=0,nz
    
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/a_d_', i, '.dat'
     call end_var_write(a_d(i)%out,filename,0)
     
    end do
   
    do i=1,nz
    
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/a_h_', i, '.dat'
     call end_var_write(a_h(i)%out,filename,0)
    
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/a_m_', i, '.dat'
     call end_var_write(a_m(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/b_u_', i, '.dat'
     call end_var_write(b_u(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/b_v_', i, '.dat'
     call end_var_write(b_v(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a9,i1.1,a4)") 'data/', ens_name, '/', num, '/a_tendh_', i, '.dat'
     call end_var_write(a_tendh(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a9,i1.1,a4)") 'data/', ens_name, '/', num, '/b_tendu_', i, '.dat'
     call end_var_write(b_tendu(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a9,i1.1,a4)") 'data/', ens_name, '/', num, '/b_tendv_', i, '.dat'
     call end_var_write(b_tendv(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a9,i1.1,a4)") 'data/', ens_name, '/', num, '/b_smagu_', i, '.dat'
     call end_var_write(b_smagu(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a9,i1.1,a4)") 'data/', ens_name, '/', num, '/b_smagv_', i, '.dat'
     call end_var_write(b_smagv(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/c_q_', i, '.dat'
     call end_var_write(c_q(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/b_q_', i, '.dat'
     call end_var_write(b_q(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/eu_', i, '.dat'
     call end_var_write(eu(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/ev_', i, '.dat'
     call end_var_write(ev(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/ef_', i, '.dat'
     call end_var_write(ef(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/ek_', i, '.dat'
     call end_var_write(ek(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/ep_', i, '.dat'
     call end_var_write(ep(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a5,i1.1,a4)") 'data/', ens_name, '/', num, '/ep2_', i, '.dat'
     call end_var_write(ep2(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/em_', i, '.dat'
     call end_var_write(em(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/en_', i, '.dat'
     call end_var_write(en(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a6,i1.1,a4)") 'data/', ens_name, '/', num, '/en_h_', i, '.dat'
     call end_var_write(en_h(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/er_', i, '.dat'
     call end_var_write(er(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a4,i1.1,a4)") 'data/', ens_name, '/', num, '/es_', i, '.dat'
     call end_var_write(es(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a8,i1.1,a4)") 'data/', ens_name, '/', num, '/ereynu_', i, '.dat'
     call end_var_write(ereynu(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a8,i1.1,a4)") 'data/', ens_name, '/', num, '/eformu_', i, '.dat'
     call end_var_write(eformu(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a8,i1.1,a4)") 'data/', ens_name, '/', num, '/ereynv_', i, '.dat'
     call end_var_write(ereynv(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a8,i1.1,a4)") 'data/', ens_name, '/', num, '/eformv_', i, '.dat'
     call end_var_write(eformv(i)%out,filename,0)
     
     
     write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/bq_uq_', i, '.dat'
     call end_var_write(bq_uq(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/bq_vq_', i, '.dat'
     call end_var_write(bq_vq(i)%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a7,i1.1,a4)") 'data/', ens_name, '/', num, '/bq_qq_', i, '.dat'
     call end_var_write(bq_qq(i)%out,filename,0)
     
     
     
    end do
    
    
     
     write (filename, "(a5,i4.4,a1,i4.4,a11)") 'data/', ens_name, '/', num, '/a_utau.dat'
     call end_var_write(a_utau%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a11)") 'data/', ens_name, '/', num, '/a_vtau.dat'
     call end_var_write(a_vtau%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a13)") 'data/', ens_name, '/', num, '/a_bfricu.dat'
     call end_var_write(a_bfricu%out,filename,0)
     
     write (filename, "(a5,i4.4,a1,i4.4,a13)") 'data/', ens_name, '/', num, '/a_bfricv.dat'
     call end_var_write(a_bfricv%out,filename,0)
     
        
   end if
   
  end if
   


  if (n >= 2*nsteps*nstep) then
  
  
   do i=1,nz
    sh(i)%zf_h%z=sh(i)%z_h(m)%z/dble(2*nsteps+1)
    sh(i)%zf_u%z=sh(i)%z_u(m)%z/dble(2*nsteps+1)    
    sh(i)%zf_v%z=sh(i)%z_v(m)%z/dble(2*nsteps+1)
    sh(i)%zf_z%z=sh(i)%z_z(m)%z/dble(2*nsteps+1)
    
    sape(i)%zf_h%z=sape(i)%z_h(m)%z/dble(2*nsteps+1)
    
    sm(i)%zf_h%z=sm(i)%z_h(m)%z/dble(2*nsteps+1)
    
    shmx(i)%zf_u%z=shmx(i)%z_u(m)%z/dble(2*nsteps+1)
   
    shmy(i)%zf_v%z=shmy(i)%z_v(m)%z/dble(2*nsteps+1)
   
    smx(i)%zf_u%z=smx(i)%z_u(m)%z/dble(2*nsteps+1)

    smy(i)%zf_v%z=smy(i)%z_v(m)%z/dble(2*nsteps+1)
          
    shq(i)%zf_h%z=shq(i)%z_h(m)%z/dble(2*nsteps+1)
    shq(i)%zf_u%z=shq(i)%z_u(m)%z/dble(2*nsteps+1)    
    shq(i)%zf_v%z=shq(i)%z_v(m)%z/dble(2*nsteps+1)
    shq(i)%zf_z%z=shq(i)%z_z(m)%z/dble(2*nsteps+1)
    
    shqq(i)%zf_h%z=shqq(i)%z_h(m)%z/dble(2*nsteps+1)
    shqq(i)%zf_z%z=shqq(i)%z_z(m)%z/dble(2*nsteps+1)
    
    shuq(i)%zf_u%z=shuq(i)%z_u(m)%z/dble(2*nsteps+1)  
       
    shvq(i)%zf_v%z=shvq(i)%z_v(m)%z/dble(2*nsteps+1)
    
    shu(i)%zf_h%z=shu(i)%z_h(m)%z/dble(2*nsteps+1)
    shu(i)%zf_u%z=shu(i)%z_u(m)%z/dble(2*nsteps+1)    
    shu(i)%zf_v%z=shu(i)%z_v(m)%z/dble(2*nsteps+1)
    shu(i)%zf_z%z=shu(i)%z_z(m)%z/dble(2*nsteps+1)
    
    shv(i)%zf_h%z=shv(i)%z_h(m)%z/dble(2*nsteps+1)
    shv(i)%zf_u%z=shv(i)%z_u(m)%z/dble(2*nsteps+1)    
    shv(i)%zf_v%z=shv(i)%z_v(m)%z/dble(2*nsteps+1)
    shv(i)%zf_z%z=shv(i)%z_z(m)%z/dble(2*nsteps+1)
    
    stendh(i)%zf_h%z=stendh(i)%z_h(m)%z/dble(2*nsteps+1)
    
    shtendu(i)%zf_u%z=shtendu(i)%z_u(m)%z/dble(2*nsteps+1)  
       
    shtendv(i)%zf_v%z=shtendv(i)%z_v(m)%z/dble(2*nsteps+1)
    
    shsmagu(i)%zf_u%z=shsmagu(i)%z_u(m)%z/dble(2*nsteps+1)  
       
    shsmagv(i)%zf_v%z=shsmagv(i)%z_v(m)%z/dble(2*nsteps+1)
    
    shuu(i)%zf_h%z=shuu(i)%z_h(m)%z/dble(2*nsteps+1)
    call start_sync(shuu(i)%zf_h)
    
    shvv(i)%zf_h%z=shvv(i)%z_h(m)%z/dble(2*nsteps+1)
    call start_sync(shvv(i)%zf_h)
    
    shuv(i)%zf_h%z=shuv(i)%z_h(m)%z/dble(2*nsteps+1)
    shuv(i)%zf_z%z=shuv(i)%z_z(m)%z/dble(2*nsteps+1)
    
   end do
   
   do i=0,nz
    sd(i)%zf_h%z=sd(i)%z_h(m)%z/dble(2*nsteps+1)
    sd(i)%zf_u%z=sd(i)%z_u(m)%z/dble(2*nsteps+1)
    sd(i)%zf_v%z=sd(i)%z_v(m)%z/dble(2*nsteps+1)
    
    sdd(i)%zf_h%z=sdd(i)%z_h(m)%z/dble(2*nsteps+1)
    call start_sync(sdd(i)%zf_h)
  end do
    
    sutau%zf_u%z=sutau%z_u(m)%z/dble(2*nsteps+1)    
    svtau%zf_v%z=svtau%z_v(m)%z/dble(2*nsteps+1)
    
    sbfricu%zf_u%z=sbfricu%z_u(m)%z/dble(2*nsteps+1)    
    sbfricv%zf_v%z=sbfricv%z_v(m)%z/dble(2*nsteps+1)
    
    do i=0,nz
     a_d(i)%z%z=sd(i)%zf_h%z
    end do
    
    do i=1,nz
     a_h(i)%z%z=sh(i)%zf_h%z
     a_m(i)%z%z=sm(i)%zf_h%z
     b_u(i)%z%z=shu(i)%zf_u%z/sh(i)%zf_u%z
     call start_sync(b_u(i)%z)
     b_v(i)%z%z=shv(i)%zf_v%z/sh(i)%zf_v%z
     call start_sync(b_v(i)%z)
     b_q(i)%z%z=shq(i)%zf_z%z/sh(i)%zf_z%z
     
     a_tendh(i)%z%z=stendh(i)%zf_h%z
     b_tendu(i)%z%z=shtendu(i)%zf_u%z/sh(i)%zf_u%z
     b_tendv(i)%z%z=shtendv(i)%zf_v%z/sh(i)%zf_v%z
     
     b_smagu(i)%z%z=shsmagu(i)%zf_u%z/sh(i)%zf_u%z
     b_smagv(i)%z%z=shsmagv(i)%zf_v%z/sh(i)%zf_v%z
     
     bq_uq(i)%z%z=shuq(i)%zf_u%z/sh(i)%zf_u%z-Ay(b_q(i)%z%z)*b_u(i)%z%z
     bq_vq(i)%z%z=shvq(i)%zf_v%z/sh(i)%zf_v%z-Ax(b_q(i)%z%z)*b_v(i)%z%z
     bq_qq(i)%z%z=shqq(i)%zf_z%z/sh(i)%zf_z%z-b_q(i)%z%z**2
     
     d_huv(i)%z=(shu(i)%zf_z%z*shv(i)%zf_z%z)/sh(i)%zf_z%z
     d_huu(i)%z=Ax(shu(i)%zf_u%z**2)/sh(i)%zf_h%z
     call start_sync(d_huu(i))
     d_hvv(i)%z=Ay(shv(i)%zf_v%z**2)/sh(i)%zf_h%z
     call start_sync(d_hvv(i))
    
     d_hm(i)%z=(sh(i)%zf_h%z-h(i)%exp)*(sm(i)%zf_h%z-mont(i)%exp)
    end do
     
    do i=0,nz
     d_dd(i)%z=sd(i)%zf_h%z**2
     call start_sync(d_dd(i))
     
    end do
    
    
     a_utau%z%z=sutau%zf_u%z
     a_vtau%z%z=svtau%zf_v%z
    
     a_bfricu%z%z=sbfricu%zf_u%z
     a_bfricv%z%z=sbfricv%zf_v%z
     
    
    
    
    
    do i=1,nz
     c_q(i)%z%z=(fz%z+sGx(b_v(i)%z)-sGy(b_u(i)%z))*  &
          merge(0.0d0,1.0d0/sh(i)%zf_z%z,(sh(i)%zf_z%z == 0.0d0))
          
     eu(i)%z%z=(0.0d0+  &
               sGx(shuu(i)%zf_h)-sGx(d_huu(i)) +  &
               Gy(shuv(i)%zf_z%z)-Gy(d_huv(i)%z) +  &
               shmx(i)%zf_u%z - sh(i)%zf_u%z*smx(i)%zf_u%z  + &
               0.0d0)/sh(i)%zf_u%z
     call start_sync(eu(i)%z)
          
     ev(i)%z%z=(0.0d0+  &
               sGy(shvv(i)%zf_h)-sGy(d_hvv(i))+  &
               Gx(shuv(i)%zf_z%z)-Gx(d_huv(i)%z) +  &
               shmy(i)%zf_v%z-sh(i)%zf_v%z*smy(i)%zf_v%z+ &
               0.0d0)/sh(i)%zf_v%z
     call start_sync(ev(i)%z)
     
     em(i)%z%z=0.5d0*(shvv(i)%zf_h%z/sh(i)%zf_h%z-Ay(b_v(i)%z%z**2)  -  &
            shuu(i)%zf_h%z/sh(i)%zf_h%z+Ax(b_u(i)%z%z**2))
            
     ek(i)%z%z=0.5d0*(shvv(i)%zf_h%z/sh(i)%zf_h%z-Ay(b_v(i)%z%z**2)  +  &
            shuu(i)%zf_h%z/sh(i)%zf_h%z-Ax(b_u(i)%z%z**2))
            
     en(i)%z%z=shuv(i)%zf_z%z/sh(i)%zf_z%z-  &
            sAx(b_v(i)%z)*sAy(b_u(i)%z)
            
     en_h(i)%z%z=shuv(i)%zf_h%z/sh(i)%zf_h%z-  &
            Ay(b_v(i)%z%z)*Ax(b_u(i)%z%z)
            
     ep(i)%z%z=0.5d0*ngp(i)*(sdd(i)%zf_h%z &
         -d_dd(i)%z)
            
     ep2(i)%z%z=0.5d0*(sape(i)%zf_h%z &
         -d_hm(i)%z)
               
               
     ereynu(i)%z%z=(sGx(shuu(i)%zf_h)-sGx(d_huu(i)) +  &
               Gy(shuv(i)%zf_z%z)-Gy(d_huv(i)%z)   &
               )/sh(i)%zf_u%z
          
     ereynv(i)%z%z=(0.0d0+  &
               sGy(shvv(i)%zf_h)-sGy(d_hvv(i))  +  &
               Gx(shuv(i)%zf_z%z)-Gx(d_huv(i)%z) +  &
               0.0d0)/sh(i)%zf_v%z
            
     eformu(i)%z%z=(shmx(i)%zf_u%z-sh(i)%zf_u%z*smx(i)%zf_u%z  +  &
               0.0d0)/sh(i)%zf_u%z
          
     eformv(i)%z%z=(0.0d0+  &
               shmy(i)%zf_v%z-sh(i)%zf_v%z*smy(i)%zf_v%z +   &
               0.0d0)/sh(i)%zf_v%z
     
    end do
    
    do i=1,nz
    
     er(i)%z%z=shmx(i)%zf_u%z-sh(i)%zf_u%z*smx(i)%zf_u%z
            
            
     es(i)%z%z=shmy(i)%zf_v%z-sh(i)%zf_v%z*smy(i)%zf_v%z
               
                   
    end do
    
    do i=1,nz
     ef(i)%z%z=(sGx(ev(i)%z)-sGy(eu(i)%z))/sh(i)%zf_z%z
    end do
    
    
    
    do i=0,nz
     call init_var_write(a_d(i)%out)
    end do
    do i=1,nz
     call init_var_write(a_h(i)%out)
     call init_var_write(a_m(i)%out)
     call init_var_write(b_u(i)%out)
     call init_var_write(b_v(i)%out)
     call init_var_write(a_tendh(i)%out)
     call init_var_write(b_tendu(i)%out)
     call init_var_write(b_tendv(i)%out)
     call init_var_write(b_smagu(i)%out)
     call init_var_write(b_smagv(i)%out)
     call init_var_write(c_q(i)%out)
     call init_var_write(b_q(i)%out)
     call init_var_write(eu(i)%out)
     call init_var_write(ev(i)%out)
     call init_var_write(ef(i)%out)
     call init_var_write(ek(i)%out)
     call init_var_write(em(i)%out)
     call init_var_write(en(i)%out)
     call init_var_write(en_h(i)%out)
     call init_var_write(ep(i)%out)
     call init_var_write(ep2(i)%out)
     call init_var_write(er(i)%out)
     call init_var_write(es(i)%out)
     call init_var_write(ereynu(i)%out)
     call init_var_write(ereynv(i)%out)
     call init_var_write(eformu(i)%out)
     call init_var_write(eformv(i)%out)
     
     call init_var_write(bq_uq(i)%out)
     call init_var_write(bq_vq(i)%out)
     call init_var_write(bq_qq(i)%out)
     
    end do
    
     call init_var_write(a_utau%out)
     call init_var_write(a_vtau%out)
     
     call init_var_write(a_bfricu%out)
     call init_var_write(a_bfricv%out)
     
     
  end if
   
  
 
 
 
 end subroutine
 
#endif
 
end module
