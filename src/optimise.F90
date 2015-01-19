module optimise
 use params

 implicit none

 contains

  subroutine run_optimise(grid)
   integer, intent(inout) :: grid(0:mx+1,0:my+1)
   integer :: grid_old(0:mx+1,0:my+1)
   integer :: proc_num
   integer :: proc1,proc2
   real(kind=8):: energy_old, energy_old2
   real(kind=8) :: energy

   proc_num = maxval(grid)+1
   
   energy_old2=0.0
   energy=energy_metric(grid,proc_num)
   do while (energy /= energy_old2)
    grid_old=grid
    energy_old2=energy
    energy_old=0.0d0
    do while (energy /= energy_old)
     energy_old=energy
     do proc1=0,proc_num-1
      do proc2=0,proc_num-1
       if (proc1 /= proc2) then
        if (any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(0:mx-1,1:my) == proc2).or.    &
     (grid(2:mx+1,1:my) == proc2)))).or.    &
     any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(1:mx,0:my-1) == proc2).or.   &
     (grid(1:mx,2:my+1) == proc2))))) then
         call box_swap(grid,proc_num,proc1,proc2)
        end if
       end if
      end do
     end do
     energy=energy_metric(grid,proc_num)
     if (all(grid == grid_old)) exit
    end do
    energy_old=0.0d0
    do while (energy /= energy_old)
     grid_old=grid
     energy_old=energy
     do proc1=0,proc_num-1
      do proc2=0,proc_num-1
       if (proc1 /= proc2) then
        if (any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(0:mx-1,1:my) == proc2).or.    &
     (grid(2:mx+1,1:my) == proc2)))).or.    &
     any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(1:mx,0:my-1) == proc2).or.   &
     (grid(1:mx,2:my+1) == proc2))))) then
         call line_swap(grid,proc_num,proc1,proc2)
        end if
       end if
      end do
     end do
     energy=energy_metric(grid,proc_num)
     if (all(grid == grid_old)) exit
    end do
    energy_old=0.0d0
   ! do while (energy /= energy_old)
   !  do proc1=0,proc_num-1
   !   do proc2=0,proc_num-1
   !    if (proc1 /= proc2) then
   !     energy_old=energy
   !     call cell_swap(grid,proc_num,proc1,proc2)
   !     energy=energy_metric(grid,proc_num)
   !    end if
   !   end do
   !  end do
   ! end do
   end do

 end subroutine

 subroutine line_swap(grid,proc_num,proc1,proc2)
  integer, intent(inout) :: grid(0:mx+1,0:my+1)
  integer, intent(in) :: proc_num,proc1,proc2
  integer:: grid_tmp(0:mx+1,0:my+1)
  integer :: i,j,k
  real(kind=8):: energy_old
  real(kind=8) :: energy

  energy_old=0
  energy=energy_metric(grid,proc_num)
  do while (energy /= energy_old)
   energy_old=energy
   do i=1,ubound(grid,1)-1
    if (any((grid(i,:) == proc1) .and. (grid(i+1,:) == proc2))) then
     grid_tmp=grid
     do k=1,ubound(grid,2)
      if (((grid(i,k) == proc1) .and. (grid(i+1,k) == proc2)) .or. &
        ((grid(i,k) == -1) .and. (grid(i+1,k) == proc2))) grid_tmp(i+1,k)=proc1
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
     grid_tmp=grid
     do k=1,ubound(grid,2)
      if (((grid(i,k) == proc1) .and. (grid(i+1,k) == proc2)) .or. &
        ((grid(i,k) == proc1) .and. (grid(i+1,k) == -1))) grid_tmp(i,k)=proc2
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if

     grid_tmp=grid
     do k=1,ubound(grid,2)
      if (((grid(i,k) == proc1) .and. (grid(i+1,k) /= -1)) .or. &
        ((grid(i,k) == -1) .and. (grid(i+1,k) == proc2))) grid_tmp(i+1,k)=proc1
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
     grid_tmp=grid
     do k=1,ubound(grid,2)
      if (((grid(i,k) /= -1) .and. (grid(i+1,k) == proc2)) .or. &
        ((grid(i,k) == proc1) .and. (grid(i+1,k) == -1))) grid_tmp(i,k)=proc2
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
    end if
   end do



   do j=1,ubound(grid,2)-1
    if (any((grid(:,j) == proc1) .and. (grid(:,j+1) == proc2))) then
     grid_tmp=grid
     do k=1,ubound(grid,1)
      if (((grid(k,j) == proc1) .and. (grid(k,j+1) == proc2)) .or. &
        ((grid(k,j) == -1) .and. (grid(k,j+1) == proc2))) grid_tmp(k,j+1)=proc1
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
     grid_tmp=grid
     do k=1,ubound(grid,1)
      if (((grid(k,j) == proc1) .and. (grid(k,j+1) == proc2)) .or. &
        ((grid(k,j) == proc1) .and. (grid(k,j+1) == -1))) grid_tmp(k,j)=proc2
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if

     grid_tmp=grid
     do k=1,ubound(grid,1)
      if (((grid(k,j) == proc1) .and. (grid(k,j+1) /= -1)) .or. &
        ((grid(k,j) == -1) .and. (grid(k,j+1) == proc2))) grid_tmp(k,j+1)=proc1
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
     grid_tmp=grid
     do k=1,ubound(grid,1)
      if (((grid(k,j) /= -1) .and. (grid(k,j+1) == proc2)) .or. &
        ((grid(k,j) == proc1) .and. (grid(k,j+1) == -1))) grid_tmp(k,j)=proc2
     end do
     if (energy_metric(grid_tmp,proc_num) < energy) then
      grid=grid_tmp
      energy_old=energy
      energy=energy_metric(grid_tmp,proc_num)
     end if
    end if
   end do
  end do
  
  
 end subroutine

 subroutine cell_swap(grid,proc_num,proc1,proc2)
  integer, intent(inout) :: grid(0:mx+1,0:my+1)
  integer, intent(in) :: proc_num,proc1,proc2
  integer:: grid_tmp(0:mx+1,0:my+1)
  integer :: i,j
  real(kind=8):: energy_old
  real(kind=8) :: energy

  energy_old=0
  energy=energy_metric(grid,proc_num)
  do while (energy /= energy_old)
   energy_old=energy
      do j=1,ubound(grid,2)
       do i=1,ubound(grid,1)-1
        if ((grid(i,j) == proc1) .and. (grid(i+1,j) == proc2)) then
         grid_tmp=grid
         grid_tmp(i+1,j)=proc1
         if (energy_metric(grid_tmp,proc_num) <= energy) then
          grid=grid_tmp
          energy_old=energy
          energy=energy_metric(grid_tmp,proc_num)

         else
          grid_tmp(i,j)=proc2
          if (energy_metric(grid_tmp,proc_num) <= energy) then
           grid=grid_tmp
           energy_old=energy
           energy=energy_metric(grid_tmp,proc_num)

          end if
         end if
        end if
       end do
      end do
      do j=1,ubound(grid,2)-1
       do i=1,ubound(grid,1)
        if ((grid(i,j) == proc1) .and. (grid(i,j+1) == proc2)) then
         grid_tmp=grid
         grid_tmp(i,j+1)=proc1
         if (energy_metric(grid_tmp,proc_num) <= energy) then
          grid=grid_tmp
          energy_old=energy
          energy=energy_metric(grid_tmp,proc_num)

         else
          grid_tmp(i,j)=proc2
          if (energy_metric(grid_tmp,proc_num) <= energy) then
           grid=grid_tmp
           energy_old=energy
           energy=energy_metric(grid_tmp,proc_num)

          end if
         end if
        end if
       end do
      end do
  end do
  
 end subroutine

 subroutine box_swap(grid,proc_num,proc1,proc2)
  integer, intent(inout) :: grid(0:mx+1,0:my+1)
  integer, intent(in) :: proc_num,proc1,proc2
  integer:: grid_tmp(0:mx+1,0:my+1)
  integer :: i,j,k,l
  real(kind=8):: energy_old
  real(kind=8) :: energy

  energy_old=0
  energy=energy_metric(grid,proc_num)
  do while (energy /= energy_old)
   energy_old=energy
    grid_tmp=grid
     if (any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(0:mx-1,1:my) == proc2).or.    &
     (grid(2:mx+1,1:my) == proc2)))) ) then
     do l=1,count((count((grid == proc1) .or. (grid == proc2),1) /= 0))-1
      k=l
      do j=1,my
       if (any(grid(:,j) == proc1) .or. any(grid(:,j) == proc2)) then
        if (k /= 0) then
         do i=1,mx
          if ((grid(i,j) == proc1) .or. (grid(i,j) == proc2)) then
           grid_tmp(i,j)=proc1
          end if
         end do
         k=k-1
        else
         do i=1,mx
          if ((grid(i,j) == proc1) .or. (grid(i,j) == proc2)) then
           grid_tmp(i,j)=proc2
          end if
         end do
        end if
       end if
      end do 
      if (energy_metric(grid_tmp,proc_num) < energy) then
       grid=grid_tmp
       energy_old=energy
       energy=energy_metric(grid_tmp,proc_num)
      end if
     end do
    end if
    grid_tmp=grid
    if (any(((grid(1:mx,1:my) == proc1).and.  &
     ((grid(1:mx,0:my-1) == proc2).or.   &
     (grid(1:mx,2:my+1) == proc2))))) then
    do k=1,count((count((grid == proc1) .or. (grid == proc2),2) /= 0))-1
    l=k
    do i=1,mx
     if (any(grid(i,:) == proc1) .or. any(grid(i,:) == proc2)) then
      if (l /= 0) then
       do j=1,my
        if ((grid(i,j) == proc1) .or. (grid(i,j) == proc2)) then
         grid_tmp(i,j)=proc1
        end if
       end do
       l=l-1
      else
       do j=1,my
        if ((grid(i,j) == proc1) .or. (grid(i,j) == proc2)) then
         grid_tmp(i,j)=proc2
        end if
       end do
      end if
     end if
    end do
    if (energy_metric(grid_tmp,proc_num) < energy) then
     grid=grid_tmp
     energy_old=energy
     energy=energy_metric(grid_tmp,proc_num)
    end if
    end do
    end if
   end do
 end subroutine

 function energy_metric(grid,proc_num)
  integer, intent(in) :: grid(0:mx+1,0:my+1)
  integer, intent(in) :: proc_num
  real(kind=8) :: energy_metric
  real(kind=8) :: num_edges, edges(0:proc_num-1), mean_edges, &
                     var_edges
  real(kind=8) :: mean_area, area_varience, area
  integer :: k

  energy_metric = 0

  mean_area=count((grid(1:mx,1:my) /= -1))
  mean_area=mean_area/real(proc_num)
  area_varience=0 
  num_edges=0 
  do k=0,proc_num-1
   edges(k)=0
   edges(k)=edges(k)+count(((grid(1:mx,1:my) == k).and.  &
     (((grid(0:mx-1,1:my) /= k).and.(grid(0:mx-1,1:my) /= -1)).or.&
     ((grid(2:mx+1,1:my) /= k).and.(grid(2:mx+1,1:my) /= -1)))))
   edges(k)=edges(k)+count(((grid(1:mx,1:my) == k).and.  &
     (((grid(1:mx,0:my-1) /= k).and.(grid(1:mx,0:my-1) /= -1)).or.&
     ((grid(1:mx,2:my+1) /= k).and.(grid(1:mx,2:my+1) /= -1)))))
   area=count((grid(1:mx,1:my) == k))
   area_varience=area_varience+(area-mean_area)**2
  end do
  mean_edges=sum(edges**2)/real(proc_num)
  var_edges=sum((edges-mean_edges)**2)/real(proc_num)
  energy_metric=abs((mean_edges+0.0d0*sqrt(var_edges))/(mean_area-sqrt(area_varience/real(proc_num))))

  end function

end module
