module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
contains
 subroutine init_random_seed(seed)
  implicit none
  integer, intent(out) :: seed
  integer   day,hour,i4_huge,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
  doext: do
    if(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       cycle doext
    else
       exit doext
    end if
  enddo doext
  doext2: do
    if (1.0D+00<temp) then
       temp=temp-1.0D+00
       cycle doext2
    else
       exit doext2
    end if
  end do doext2
  seed=int(dble(i4_huge)*temp)
  if(seed == 0)       seed = 1
  if(seed == i4_huge) seed = seed-1
  return
 end subroutine init_random_seed
 integer function randint(i,j,seed)
  real               ::  a,b
  integer,intent(in) ::  i,j,seed
  a = real(i)
  b = real(j)
  randint=int(r4_uniform(a,b+1.0,seed))
 end function randint

 real function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed == 0) then
   write(*,'(b1)')'R4_UNIFORM - Fatal error!'
   write(*,'(b1)')'Input value of SEED = 0.'
   stop '[ERROR] Chiquitan chiquitintatan'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  return
 end function r4_uniform
end module

program init_renormalisation
 use iso_fortran_env
 use mod_random 
 implicit none
 real :: a(8)
 integer :: i,j,k,l,n,seed
 integer :: n_T=60
 integer :: n_average=0 
 real    :: x_cell(1:27),x_average,eta
 real    :: temperature=448.0     ! K
 real    :: k_B  = 8.61734E-5     ! eV/K
 real    :: beta, x, partition = 0.0
 real    :: probability_cell(1:27),free(1:27),free_min
 real    :: e_1,e_2
 call init_random_seed(seed)
 partition = 0.0
 beta = 1.0/(k_B*temperature)
 write(6,*)"Parameters:"
 read(5,*)  n_average 
 x_average=(n_average/real(n_T))
 do i=1,8
  read(5,*)  a(i)
  write(6,*) a(i)
 end do
 write(6,*)"======"
 do i=1,n_T
  x=(i/real(n_T))
  write(6,*)x,Free_Energy(a,x)
 end do
 write(6,*)"======"
 x_cell(1:27) = x_average
  do l=1,27
   e_1=Free_Energy(a,x_cell(l))
  end do
 write(6,'(27(i2,1x))') (int(n_T*x_cell( l )) , l=1,27)
 write(6,*)"======"
 call update_probability(x_cell,free,partition,probability_cell)
 do k=1,1000
  i=randint(1,27,seed) 
  j=randint(1,27,seed)
  do while ( i==j .or. x_cell(i) <= 0.0 .or. x_cell(j) <= 0.0 .or. x_cell(i) >= 1.0 .or. x_cell(j)>=1.0)
   i=randint(1,27,seed)
   j=randint(1,27,seed)
  end do
  e_1=Free_Energy(a,x_cell(i)) + Free_Energy(a,x_cell(j))
  n=randint(1,3,seed)
  x_cell(i)=x_cell(i)+n/60.0
  x_cell(j)=x_cell(j)-n/60.0
  e_2=Free_Energy(a,x_cell(i)) + Free_Energy(a,x_cell(j))
  if(e_1 >= e_2 ) then
   call update_probability(x_cell,free,partition,probability_cell)
   write(6,'(27(i2,1x),a1,1x,i2,1x,i2)') (int(n_T*x_cell( l )) ,l=1,27),':',int(n_t*averages()),int(sum(x_cell)*(n_T/27.0))
  else
   eta=R4_UNIFORM(0.0,1.0,SEED)
   if ( eta < exp(-( e_2 - e_1 )/(k_B*temperature)) ) then
    call update_probability(x_cell,free,partition,probability_cell)
    write(6,'(27(i2,1x),a1,1x,i2,1x,i2)') (int(n_T*x_cell( l )) ,l=1,27),':',int(n_t*averages()),int(sum(x_cell)*(n_T/27.0))
   else
   x_cell(i)=x_cell(i)-n/60.0
   x_cell(j)=x_cell(j)+n/60.0
   call update_probability(x_cell,free,partition,probability_cell)
   end if
  end if
 end do
 contains
 real function averages()
  implicit none
  integer :: i
  averages=0.0
  do i=1,27
   averages=averages+x_cell(i)*probability_cell(i)
  end do
  return
 end function averages
 subroutine update_probability(x_cell,free,partition,probability_cell)
 implicit none 
  integer          :: j
  real,intent(in)  :: x_cell(1:27)
  real,intent(out) :: probability_cell(1:27)
  real,intent(out) :: partition,free(1:27)
  real :: free_min = 0.0
  partition=0.0
  do j=1,27
   free(j)=Free_Energy(a,x_cell(j) )
  end do
  free_min=minval(free)
  do j =1,27
   partition=partition + exp(-beta*(Free_Energy(a,x_cell(j) )-free_min ))
  end do
  do j=1,27
   probability_cell(j)=exp(-beta*(Free_Energy(a,x_cell(j) ) - free_min))/partition
  end do
  return
 end subroutine 
 real function Free_Energy(a,x)
  implicit none
  real, intent(in) :: a(8)
  real,intent(in) :: x
  Free_Energy = a(1)*(x-a(2))**2 + a(3)*(x-a(4))**3 + &
                a(5)*(x-a(6))**4 + a(7)*(x-a(8))**5
  return
 end function Free_Energy
end program init_renormalisation
