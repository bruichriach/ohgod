module global

 implicit none
 
 type point
  logical :: core
  real (kind = 8) :: x, y
  integer :: i, j
  logical :: flux
 end type
 
 type use
  real(kind=8) , pointer :: z
 end type
 
 type bound
  real (kind = 8) :: z
  integer :: vec, dir
  real(kind = 8) :: slip
 end type
 
 type mpi_sendrecv
  type(use), allocatable :: s(:)
  type(bound), allocatable :: r(:)
  real(kind = 8), allocatable :: s_mpi(:)
  real(kind = 8), allocatable :: r_mpi(:)
  integer :: s_req
  integer :: r_req
 end type
 
 type var
  type(point), pointer :: p(:,:)
  real(kind=8), allocatable :: z(:,:)
  real(kind=8), allocatable :: exp(:,:)
  real(kind=8), allocatable :: spng(:,:)
  real(kind=8), pointer :: tend0(:,:)
  real(kind=8), pointer :: tend1(:,:)
  real(kind=8), pointer :: tend2(:,:)
  real(kind=8), pointer :: tmp(:,:)
  logical, allocatable :: mask(:,:)
  type(use), allocatable :: bx(:,:), ux(:,:), by(:,:), uy(:,:)
  type(mpi_sendrecv), allocatable :: mpi(:)
  integer :: nx, ny
  integer :: lx, ly
  integer :: tag
  logical :: synced
 end type
 
 real(kind=8), target :: null_field
 
 
end module


module variables
 use params
 use global
 
 implicit none
 
 type(var) :: h(nz), u(nz), v(nz)
 type(var) :: h_z(nz), h_u(nz), h_v(nz)
 type(var) :: h1(nz), u1(nz), v1(nz)
 type(var) :: lapu(nz), lapv(nz)
 type(var) :: bfricu, bfricv
 type(var) :: mont(nz)
 type(var) :: true_d(0:nz), true_p(0:nz)
 type(var) :: smag(nz), smagu(nz), smagv(nz)
 type(var) :: ke(nz)
 type(var) :: ape(nz)
 
 type(var) :: s
 

 type(var) :: depth(0:nz)
 type(var) :: zeta(nz), q_h(nz), q_z(nz)
 type(var) :: h_tmp(nz), u_tmp(nz), v_tmp(nz)
 type(var) :: fz
 type(var) :: utau, vtau
 type(var) :: strain(nz), tension(nz)
 
end module


module test

 implicit none

 integer, save :: ierr, ierr_tmp, test_tag=0
 real(kind=8), save :: check(4)


end module
  
