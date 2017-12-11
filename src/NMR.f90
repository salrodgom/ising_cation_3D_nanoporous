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
  parameter (i4_huge=2147483277)
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
!
 real function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483277)
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
! main program
program NMR
 use iso_fortran_env
 implicit none
 integer                         :: n_Ge = 0
 integer,parameter               :: n_resonator_max = 26
 integer                         :: i,j,k,ierr=0
 integer                         :: u = 100
 integer                         :: n_lines = 0
 integer,parameter               :: n_configurations_max = 1000
 integer                         :: virtual_n_configurations = 0
 integer,parameter               :: virtual_n_configurations_max = 1000
 integer,parameter               :: n_cells = 4*4*4
 integer,parameter               :: initialisation_steps = 10000000
 real                            :: virtual_n_Ge(n_cells,virtual_n_configurations_max)
 logical                         :: Quasi_Random_Virtual_Structures_flag = .false.
 logical                         :: forcefield = .false., inquite_Geometry_information_debug = .false.
 logical                         :: inquire_CIFFiles_flag = .false., inquite_Geometry_information=.false.
! arguments in line:
 integer                         :: num_args = 0
 character(len=100),dimension(:), allocatable :: args
! allocatables
 type                            :: virtual_extension
  real                           :: n_Ge
 end type
 type  :: microstate
  integer                        :: id = 0
  integer                        :: n_configurations = 0
  real                           :: peso(1:n_configurations_max) = 0.0
  real                           :: partition_function = 0.0
  character(len=1024)            :: CIFFile(1:n_configurations_max) = " "
  integer                        :: resonator(1:n_configurations_max,0:n_resonator_max-1) = 0
  real                           :: T_population(1:n_configurations_max,1:5) = 0
  real                           :: extended_virtual_n_Ge = 0.0
  real                           :: extended_virtual_n_Ge_sigma = 0.0
  integer                        :: extended_n_configurations = 0
  type(virtual_extension)        :: extended_ensemble(1:n_cells,virtual_n_configurations_max)
! cell properties:
  real                           :: cell_parameter(6,0:n_configurations_max)
  real                           :: volume(0:n_configurations_max)
  real                           :: rv_matrix(3,3,1:n_configurations_max)
  real                           :: vr_matrix(3,3,1:n_configurations_max)
  real                           :: distance_TO(0:2,0:n_configurations_max)
  real                           :: distance_TT(0:2,0:n_configurations_max)
  real                           :: angle_OTO(0:2,0:n_configurations_max)
  real                           :: angle_TOT(0:2,0:n_configurations_max)
 end type                        
 type(microstate)                :: ensemble(0:60)
 real                            :: count_(0:60,0:n_resonator_max-1+5) = 0.0
 real                            :: NMR_(0:60,0:n_resonator_max-1+5) = 0.0
 real                            :: NMR_unify(0:60,1:4) = 0.0
 real                            :: suma,x(1:16)
 real                            :: cell_0(6),vr(3,3),rv(3,3)
 integer                         :: cont,resonator
 real                            :: r4_average(0:60,1:8)=0.0
 character(len=1024)             :: filename = " "
 character(len=1024)             :: line = " "
 logical                         :: lex
 ! Read Arguments in line
 num_args=command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 if(num_args >= 1)then
  do i=1,num_args
   select case(args(i))
    case('-h','--help')
     call print_help()
     stop
    case('-f','--forcefield')
     forcefield=.true.
    case('-e','--QSVS')
     Quasi_Random_Virtual_Structures_flag = .true.
    case('-g','--inquire-CIF')
     inquire_CIFFiles_flag = .true.
    case('-t','--inquire-topology')
     inquite_Geometry_information = .true.
    case('-x')
     inquite_Geometry_information_debug = .true.
   end select
  end do
 end if
! initialisation:
 ensemble(0)%n_configurations=1
 ensemble(0)%peso(1)=1.0
 ensemble(0)%CIFFile(1)='0_0.cif'
 ensemble(0)%resonator(1,0)=1
 ensemble(0)%T_population(1,1:5)=0
 ensemble(0)%extended_virtual_n_Ge=0.0
 ensemble(0)%extended_n_configurations = 1
 ensemble(0)%extended_ensemble(1:n_cells,1:1)%n_Ge = 0.0
 ensemble(0)%distance_TO(1,1)=1.620421
 ensemble(0)%distance_TO(2,1)=0.0
 ensemble(0)%distance_TT(1,1)=3.000330
 ensemble(0)%distance_TT(2,1)=0.0
 ensemble(0)%angle_OTO(1,1)=109.45
 ensemble(0)%angle_OTO(2,1)=0.0
 ensemble(0)%angle_TOT(1,1)=142.151852
 ensemble(0)%angle_TOT(2,1)=0.0
 ensemble(0)%T_population(1,1:5)=0.0
 !
 ensemble(60)%n_configurations=1
 ensemble(60)%peso(1)=1.0
 ensemble(60)%CIFFile(1)='60_0.cif'
 ensemble(60)%resonator(1,n_resonator_max-1)=1
 ensemble(60)%T_population(1,1:5)=1
 ensemble(60)%extended_virtual_n_Ge=60.0
 ensemble(60)%extended_n_configurations = 1
 ensemble(60)%extended_ensemble(1:n_cells,1:1)%n_Ge = 60.0
 ensemble(60)%distance_TO(1,0:)=1.722453
 ensemble(60)%distance_TO(2,0:)=1.724428
 ensemble(60)%distance_TT(1,0:)=3.214147
 ensemble(60)%distance_TT(2,0:)=3.214147
 ensemble(60)%angle_OTO(1,0:)=109.157674
 ensemble(60)%angle_OTO(2,0:)=109.157674
 ensemble(60)%angle_TOT(1,0:)=138.354665
 ensemble(60)%angle_TOT(2,0:)=138.354665
 ensemble(60)%T_population(1,1:5)=1.0
 !
 count_(0,0) = 6.0
 count_(60,n_resonator_max-1) = 6.0
 !
 if( inquire_CIFFiles_flag ) then
  ! 0
  call inquire_CIFFile( ensemble(0)%CIFFile(1), cell_0, rv, vr)
  ensemble(0)%cell_parameter(1:6,1) = cell_0(1:6)
  ensemble(0)%rv_matrix(1:3,1:3,1) = rv(1:3,1:3)
  ensemble(0)%vr_matrix(1:3,1:3,1) = vr(1:3,1:3)
  call cell(ensemble(0)%rv_matrix(1:3,1:3,1),ensemble(0)%vr_matrix(1:3,1:3,1),ensemble(0)%cell_parameter(1:6,1) )
  ensemble(0)%volume(1) = volume( ensemble(0)%rv_matrix(1:3,1:3,1)  )
  ! 60
  call inquire_CIFFile( ensemble(60)%CIFFile(1), cell_0, rv, vr)
  ensemble(60)%cell_parameter(1:6,1) = cell_0(1:6)
  ensemble(60)%rv_matrix(1:3,1:3,1) = rv(1:3,1:3)
  ensemble(60)%vr_matrix(1:3,1:3,1) = vr(1:3,1:3)
  call cell(ensemble(60)%rv_matrix(1:3,1:3,1),ensemble(60)%vr_matrix(1:3,1:3,1),ensemble(60)%cell_parameter(1:6,1) )
  ensemble(60)%volume(1) = volume( ensemble(60)%rv_matrix(1:3,1:3,1)  )
 end if
 !
 if( inquite_Geometry_information_debug ) then
  open(825,file="overall_input.txt",status='old',iostat=ierr)
  if(ierr/=0) stop
  open(826,file="germanium_input.txt",status='old',iostat=ierr)
  if(ierr/=0) stop
 end if
 ! ---------------------------------------------------------------------------------------------------------------------
 read_properties: do n_Ge=1,59
  n_lines=0
  open(u, file='inputs/analysis_'//trim(str(n_Ge))//'.txt',status='old',iostat=ierr)
  if(ierr/=0) then
   ensemble(n_Ge)%n_configurations=1
   ensemble(n_Ge)%id=n_Ge
   ensemble(n_Ge)%peso(1)=0
   ensemble(n_Ge)%resonator(1,0:n_resonator_max-1)=0
   ensemble(n_Ge)%T_population(1,1:5)=0
   ensemble(n_Ge)%CIFFile(1)=" "
   ensemble(n_Ge)%partition_function=0.0
   cycle read_properties
  end if
  if( inquite_Geometry_information ) then
   open(345, file='inputs/configuration_topology.'//trim(str(n_Ge))//'.txt',status='old',iostat=ierr)
   if(ierr/=0) stop "File not found."
  end if
  if ( inquite_Geometry_information_debug ) then
   read(825,*) j,ensemble(n_Ge)%distance_TO(1,0),ensemble(n_Ge)%angle_OTO(1,0),&
                 ensemble(n_Ge)%distance_TT(1,0),ensemble(n_Ge)%angle_TOT(1,0)
   read(826,*) j,ensemble(n_Ge)%distance_TO(2,0),ensemble(n_Ge)%angle_OTO(2,0),&
                 ensemble(n_Ge)%distance_TT(2,0),ensemble(n_Ge)%angle_TOT(2,0)
   do k=1,2
    do i=1,ensemble(n_Ge)%n_configurations
     ensemble(n_Ge)%distance_TO(k,i)=ensemble(n_Ge)%distance_TO(k,0)
     ensemble(n_Ge)%distance_TT(k,i)=ensemble(n_Ge)%distance_TT(k,0)
       ensemble(n_Ge)%angle_OTO(k,i)=  ensemble(n_Ge)%angle_OTO(k,0)
       ensemble(n_Ge)%angle_TOT(k,i)=  ensemble(n_Ge)%angle_TOT(k,0)
    end do
   end do
  end if
  inquire(u,EXIST=lex,name=filename)
  inquire_analysis: do
   read(u,'(a)',iostat=ierr) line
   if (ierr/=0) exit inquire_analysis
   n_lines=n_lines+1
  end do inquire_analysis
  rewind(u)
  ensemble(n_Ge)%n_configurations=n_lines-1
  ensemble(n_Ge)%id=n_Ge
  do i=1,ensemble(n_Ge)%n_configurations
   read(u,'(a)',iostat=ierr) line
   if (ierr/=0) exit
   read(line,*)ensemble(n_Ge)%peso(i),(ensemble(n_Ge)%resonator(i,j),j=0,n_resonator_max-1),&
    (ensemble(n_Ge)%T_population(i,j),j=1,5),ensemble(n_Ge)%CIFFile(i)
   if( inquire_CIFFiles_flag ) then
    vjw: do j=1,1024
     if(line(j:j)=="S") then
      cont=j
      exit vjw 
     end if
    end do vjw
    pql: do j=cont,1024
     if(line(j:j)==" ")then
      resonator=j
      exit pql
     end if
    end do pql
    read(line(cont:resonator),'(a)') ensemble(n_Ge)%CIFFile(i)
    ensemble(n_Ge)%CIFFile(i)="inputs/"//ensemble(n_Ge)%CIFFile(i)
    ensemble(n_Ge)%CIFFile(i)=adjustl( trim( ensemble(n_Ge)%CIFFile(i) ) )
    call inquire_CIFFile( ensemble(n_Ge)%CIFFile(i), cell_0, rv, vr)
    ensemble(n_Ge)%cell_parameter(1:6,i) = cell_0(1:6)
    ensemble(n_Ge)%rv_matrix(1:3,1:3,i) = rv(1:3,1:3)
    ensemble(n_Ge)%vr_matrix(1:3,1:3,i) = vr(1:3,1:3)
    call cell(ensemble(n_Ge)%rv_matrix(1:3,1:3,i),ensemble(n_Ge)%vr_matrix(1:3,1:3,i),ensemble(n_Ge)%cell_parameter(1:6,i) )
    ensemble(n_Ge)%volume(i) = volume( ensemble(n_Ge)%rv_matrix(1:3,1:3,i)  )
    call Symmetrising(ensemble(n_Ge)%cell_parameter(1:6,i) ,ensemble(n_Ge)%rv_matrix(1:3,1:3,i),   &
         ensemble(n_Ge)%vr_matrix(1:3,1:3,i),   &
         ensemble(n_Ge)%volume(i) )
   end if
   if( inquite_Geometry_information ) then
     read(345,'(a)') line
     if(scan(line,"nan")==0)then
      read(line,*) ( ensemble(n_Ge)%distance_TO(j,i),j=1,2), ( ensemble(n_Ge)%distance_TT(j,i),j=1,2), &
                   (   ensemble(n_Ge)%angle_OTO(j,i),j=1,2), (   ensemble(n_Ge)%angle_TOT(j,i),j=1,2)
     else
      write(6,'(a,a)')'Found NaN: ',line
      ensemble(n_Ge)%peso(i)=0.0 
      ensemble(n_Ge)%distance_TO(1:2,i)=0.0
      ensemble(n_Ge)%distance_TT(1:2,i)=0.0
      ensemble(n_Ge)%angle_OTO(1:2,i)=0.0
      ensemble(n_Ge)%angle_TOT(1:2,i)=0.0
     end if
   end if
  end do
  if(inquite_Geometry_information) close(345)
  ensemble(n_Ge)%partition_function=sum( ensemble(n_Ge)%peso( 1:ensemble(n_Ge)%n_configurations ) )
  do i=1,ensemble(n_Ge)%n_configurations
   do j=0,n_resonator_max-1+5
    !FORBIDEN RESONATORS
    if (j/=2.and.j/=6.and.j/=10.and.j/=16.and.j/=20) then
     if(forcefield)then
      if(ensemble(n_Ge)%partition_function/=0)then
       if(j<=n_resonator_max-1)then
        count_(n_Ge,j) = count_(n_Ge,j) + &
        ensemble(n_Ge)%resonator(i,j)*ensemble(n_Ge)%peso(i)/ensemble(n_Ge)%partition_function
       else
        count_(n_Ge,j) = count_(n_Ge,j) + &
        ensemble(n_Ge)%T_population(i,j-n_resonator_max+1)*ensemble(n_Ge)%peso(i)/ensemble(n_Ge)%partition_function
       end if
      end if
     else
      if(j<=n_resonator_max-1)then
       count_(n_Ge,j) = count_(n_Ge,j) + ensemble(n_Ge)%resonator(i,j)
      else  
       count_(n_Ge,j) = count_(n_Ge,j) + ensemble(n_Ge)%T_population(i,j-n_resonator_max+1)
      end if
     end if
    end if
   end do 
  end do
  write(6,'(a,1x,i4)')'N configurations:',ensemble(n_Ge)%n_configurations
  do i=1,ensemble(n_Ge)%n_configurations
   write(6,'(a,1x,a)')        'CIFFile:',ensemble(n_Ge)%CIFFile(i)(1:clen_trim(ensemble(n_Ge)%CIFFile(i)))
   if( inquire_CIFFiles_flag ) then
    write(6,'(a,1x,6(f14.7,1x))')    'Cell Parameters:',(ensemble(n_Ge)%cell_parameter(j,i),j=1,6)
    write(6,'(a,1x,f14.7)')    'Volume:',ensemble(n_Ge)%volume(i)
   end if
   if ( inquite_Geometry_information ) then
    write(6,*) (ensemble(n_Ge)%distance_TO(j,i),j=1,2)
    write(6,*) (ensemble(n_Ge)%distance_TT(j,i),j=1,2)
    write(6,*) (  ensemble(n_Ge)%angle_OTO(j,i),j=1,2)
    write(6,*) (  ensemble(n_Ge)%angle_TOT(j,i),j=1,2)
   end if
   if ( inquite_Geometry_information_debug ) then
    write(6,*) ensemble(n_Ge)%distance_TO(1,i),ensemble(n_Ge)%angle_OTO(1,i),&
               ensemble(n_Ge)%distance_TT(1,i),ensemble(n_Ge)%angle_TOT(1,i)
    write(6,*) ensemble(n_Ge)%distance_TO(2,i),ensemble(n_Ge)%angle_OTO(2,i),&
               ensemble(n_Ge)%distance_TT(2,i),ensemble(n_Ge)%angle_TOT(2,i) 
   end if
   write(6,'(a,1x,f14.7)')    'Peso:',ensemble(n_Ge)%peso(i)
   write(6,'(a,1x,f14.7)')    'Partition Function:',ensemble(n_Ge)%partition_function
   write(6,'(a,1x,26(i2,1x))')'Resonadores:',( ensemble(n_Ge)%resonator(i,j),j=0,n_resonator_max-1 )
   write(6,'(a,1x,5(i2,1x),1x,a,1x,i2)') 'Populations:',(int(ensemble(n_Ge)%T_population(i,j)),j=1,5),&
    ':',int(sum(ensemble(n_Ge)%T_population(i,1:5)))
  end do
  write(6,*)'===='
  close(u)
 end do read_properties
 if( inquite_Geometry_information_debug) then
  close(825)
  close(826)
 end if
 write(6,'(a)')'Counts:'
 do n_Ge=0,60
  write(6,'(26(f5.2,1x))')(count_(n_Ge,j), j=0,n_resonator_max-1 + 5 )
 end do
 write(6,*)'===='
 write(6,'(26(i4,1x))')( sum(int(count_(0:60,j)) ), j=0,n_resonator_max-1 + 5 )
 write(6,*)'===='
 !
 if( Quasi_Random_Virtual_Structures_flag.and.n_cells>=2 ) write(6,'(a)')"Quasi Random Virtual Structures"
 do n_Ge=1,59
  if( Quasi_Random_Virtual_Structures_flag .eqv. .false. .or. n_cells == 1 ) then
   ensemble(n_Ge)%extended_n_configurations = 1
   ensemble(n_Ge)%extended_virtual_n_Ge = real(n_Ge)
   ensemble(n_Ge)%extended_ensemble(1:n_cells,1)%n_Ge = real(n_Ge)
  else
   call Quasi_Random_Virtual_Structures(n_cells,virtual_n_configurations_max,n_Ge,virtual_n_configurations,virtual_n_Ge)
   ensemble(n_Ge)%extended_n_configurations = virtual_n_configurations
   ensemble(n_Ge)%extended_virtual_n_Ge=0.0
   do i = 1, ensemble(n_Ge)%extended_n_configurations
    do k = 1,n_cells
     ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge = virtual_n_Ge(k,i)
     ensemble(n_Ge)%extended_virtual_n_Ge=ensemble(n_Ge)%extended_virtual_n_Ge+&
      ( ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge/real(n_cells*ensemble(n_Ge)%extended_n_configurations) )
    end do
   end do
  end if
  write(6,'(i3,1x,f14.7,1x,i5)')n_Ge,ensemble(n_Ge)%extended_virtual_n_Ge,ensemble(n_Ge)%extended_n_configurations
 end do
 ensemble(n_Ge)%T_population(1:n_configurations_max,1:5) = 0.0
 do n_Ge=0,60
  do resonator=0,n_resonator_max-1+5
   NMR_(n_Ge,resonator)=0.0
   if (resonator/=2.and.resonator/=6.and.resonator/=10.and.resonator/=16.and.resonator/=20) then
   do i=1,ensemble(n_Ge)%extended_n_configurations
    do k = 1,n_cells   
     NMR_(n_Ge,resonator)=NMR_(n_Ge,resonator)+count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator)
     select case(resonator)
      case(0)
       NMR_unify(n_Ge,1)=NMR_unify(n_Ge,1)+count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator)
      case(1,4,5)
       NMR_unify(n_Ge,2)=NMR_unify(n_Ge,2)+count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator)
      case(3,7,15,19,23)
       NMR_unify(n_Ge,3)=NMR_unify(n_Ge,3)+count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator)
      case(11,18,24,25)
       NMR_unify(n_Ge,4)=NMR_unify(n_Ge,4)+count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator)
      case(n_resonator_max,n_resonator_max+1,n_resonator_max+2,n_resonator_max+3,n_resonator_max+4)
       ensemble(n_Ge)%T_population(1,resonator-n_resonator_max+1) = ensemble(n_Ge)%T_population(1,resonator-n_resonator_max+1) + &
        count_( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) ,resonator )
     end select
    end do
   end do
   end if
  end do
  suma=sum( NMR_(n_Ge,0:n_resonator_max-1))
  NMR_(n_Ge,0:n_resonator_max-1)=NMR_(n_Ge,0:n_resonator_max-1)/suma
  suma=sum(NMR_unify(n_Ge,1:4))
  NMR_unify(n_Ge,1:4)=NMR_unify(n_Ge,1:4)/suma
 end do
 open(u,file="populations.txt")
 do n_Ge=0,60
  write(u,'(f14.7,1x,5(f14.7,1x))') ensemble(n_Ge)%extended_virtual_n_Ge,&
  ( ensemble(n_Ge)%T_population(1,j)/real(n_cells*ensemble(n_Ge)%extended_n_configurations) ,j=1,5 )
 end do
 close(u)
 open(u,file="NMR.txt")
 do n_Ge=0,60
  write(u,'(f14.7,1x,27(f14.7,1x))') ensemble(n_Ge)%extended_virtual_n_Ge,&
   ( NMR_(n_Ge,resonator), resonator=0,n_resonator_max-1 ),sum( NMR_(n_Ge,0:n_resonator_max-1) )
 end do
 close(u)
 open(u,file="NMR_unify.txt")
 do n_Ge=0,60
  write(u,'(f14.7,1x,4(f14.7,1x))') ensemble(n_Ge)%extended_virtual_n_Ge,&
   ( NMR_unify(n_Ge,resonator), resonator=1,4 )
 end do
 close(u)
 !! cell averages
 do n_Ge=1,59
  do j=1,6
   call average(ensemble(n_Ge)%n_configurations,&
                ensemble(n_Ge)%cell_parameter(j,0:ensemble(n_Ge)%n_configurations),&
                ensemble(n_Ge)%peso )
  end do
  call average(ensemble(n_Ge)%n_configurations,ensemble(n_Ge)%volume,ensemble(n_Ge)%peso)
  do j=1,2
   call average( ensemble(n_Ge)%n_configurations,&
                 ensemble(n_Ge)%distance_TO(j,0:ensemble(n_Ge)%n_configurations),&
                 ensemble(n_Ge)%peso )
   call average( ensemble(n_Ge)%n_configurations,&
                 ensemble(n_Ge)%distance_TT(j,0:ensemble(n_Ge)%n_configurations),&
                 ensemble(n_Ge)%peso )
   call average( ensemble(n_Ge)%n_configurations,&
                 ensemble(n_Ge)%angle_OTO(j,0:ensemble(n_Ge)%n_configurations),&
                 ensemble(n_Ge)%peso )
   call average( ensemble(n_Ge)%n_configurations,&
                 ensemble(n_Ge)%angle_TOT(j,0:ensemble(n_Ge)%n_configurations),&
                 ensemble(n_Ge)%peso )
  end do
 end do
! 0
!_cell_length_a                       11.9021
!_cell_length_b                       11.9021
!_cell_length_c                       29.6920
! 60
!_cell_length_a                       12.663980
!_cell_length_b                       12.662056
!_cell_length_c                       30.972105
!ensemble(0)%cell_parameter(1,0:1)=11.9021
!ensemble(0)%cell_parameter(2,0:1)=11.9021
!ensemble(0)%cell_parameter(3,0:1)=29.6920
!ensemble(0)%cell_parameter(4,0:1)=90.0
!ensemble(0)%cell_parameter(5,0:1)=90.0
!ensemble(0)%cell_parameter(6,0:1)=120.0
!
!
!ensemble(60)%cell_parameter(1,0:1)=12.663
!ensemble(60)%cell_parameter(2,0:1)=12.663
!ensemble(60)%cell_parameter(3,0:1)=30.9721
!ensemble(60)%cell_parameter(4,0:1)=90.0
!ensemble(60)%cell_parameter(5,0:1)=90.0
!ensemble(60)%cell_parameter(6,0:1)=120.0
!
 if( inquire_CIFFiles_flag ) then
  ! 0
  call inquire_CIFFile( ensemble(0)%CIFFile(1), cell_0, rv, vr)
  ensemble(0)%cell_parameter(1:6,0) = cell_0(1:6)
  ensemble(0)%rv_matrix(1:3,1:3,1) = rv(1:3,1:3)
  ensemble(0)%vr_matrix(1:3,1:3,1) = vr(1:3,1:3)
  call cell(ensemble(0)%rv_matrix(1:3,1:3,1),ensemble(0)%vr_matrix(1:3,1:3,1),ensemble(0)%cell_parameter(1:6,0) )
  ensemble(0)%volume(0) = volume( ensemble(0)%rv_matrix(1:3,1:3,1)  )
  ! 60
  call inquire_CIFFile( ensemble(60)%CIFFile(1), cell_0, rv, vr)
  ensemble(60)%cell_parameter(1:6,0) = cell_0(1:6)
  ensemble(60)%rv_matrix(1:3,1:3,1) = rv(1:3,1:3)
  ensemble(60)%vr_matrix(1:3,1:3,1) = vr(1:3,1:3)
  call cell(ensemble(60)%rv_matrix(1:3,1:3,1),ensemble(60)%vr_matrix(1:3,1:3,1),ensemble(60)%cell_parameter(1:6,0) )
  ensemble(60)%volume(0) = volume( ensemble(60)%rv_matrix(1:3,1:3,1)  )
 end if

 ! QRVS averages:
 if( Quasi_Random_Virtual_Structures_flag.and.n_cells>=2 )then
  !ensemble(n_Ge)%distance_TO(,0)=ensemble(n_Ge)%distance_TO(j,0)
  ! properties: TO
  r4_average(0:60,1:8)=0.0
  do n_Ge = 0,60
   do j=1,8
    resonator=0
    suma=0.0
    do i=1,ensemble(n_Ge)%extended_n_configurations
     do k=1,n_cells
      resonator=resonator+1
      cont = int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge )
      if(cont==0)cont=1
      if(cont==60)cont=59
      kopwq: select case(j)
       case(1,2)
       suma=suma+ensemble(cont)%distance_TO(j  ,0)
       case(3,4)
       suma=suma+ensemble(cont)%distance_TT(j-2,0)
       case(5,6)
       suma=suma+  ensemble(cont)%angle_OTO(j-4,0)
       case(7,8)
       suma=suma+  ensemble(cont)%angle_TOT(j-6,0)
      end select kopwq
     end do
    end do 
    r4_average(n_Ge,j)=suma/real(resonator)
   end do 
  end do
  do n_Ge=0,60
   do j=1,8 
    kosa: select case(j)
     case(1,2)
     ensemble(n_Ge)%distance_TO(j,0)=r4_average(n_Ge,j)
     case(3,4)
     ensemble(n_Ge)%distance_TT(j-2,0)=r4_average(n_Ge,j)
     case(5,6)
     ensemble(n_Ge)%angle_OTO(j-4,0)=r4_average(n_Ge,j)
     case(7,8)
     ensemble(n_Ge)%angle_TOT(j-6,0)=r4_average(n_Ge,j)
    end select kosa
   end do
  end do
  ! cell parameters 
  r4_average(0:60,1:8)=0.0
  do n_Ge=0,60
   cps: do j=1,6
    resonator=0
    suma=0.0
    ik_extended: do i=1,ensemble(n_Ge)%extended_n_configurations
     do k = 1,n_cells
      resonator=resonator+1
      cont = int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge )
      suma=suma+ensemble(cont)%cell_parameter(j,0) 
     end do
    end do ik_extended
    r4_average(n_Ge,j)=suma/real(resonator)
   end do cps
  end do
  do n_Ge=0,60
   do j=1,6
    ensemble(n_Ge)%cell_parameter(j,0)=r4_average(n_Ge,j)
   end do
  end do
  ! volume:
  do n_Ge=0,60
   suma = 0.0
   do i=1,ensemble(n_Ge)%extended_n_configurations
    do k = 1,n_cells
     suma=suma+ensemble( int(ensemble(n_Ge)%extended_ensemble(k,i)%n_Ge) )%volume(0)
    end do
   end do
   r4_average(n_Ge,1)=suma/real(n_cells*ensemble(n_Ge)%extended_n_configurations)
  end do
  do n_Ge=0,60
   ensemble(n_Ge)%volume(0)=r4_average(n_Ge,1)
  end do
 end if
 open(u,file="cell.txt")
 open(789,file="properties.txt")
 do n_Ge=0,60
  write(u,*)ensemble(n_Ge)%extended_virtual_n_Ge,( ensemble(n_Ge)%cell_parameter(j,0),j=1,6 ), ensemble(n_Ge)%volume(0)
  write(789,*)ensemble(n_Ge)%extended_virtual_n_Ge,(ensemble(n_Ge)%distance_TO(j,0),j=1,2),&
   (ensemble(n_Ge)%distance_TT(j,0),j=1,2),(ensemble(n_Ge)%angle_OTO(j,0),j=1,2),(ensemble(n_Ge)%angle_TOT(j,0),j=1,2) 
 end do
 close(u)
 close(789)
 stop
 !
 contains
!
 subroutine average(n,data_,peso)
  implicit none
  integer,intent(in) :: n
  real,intent(inout) :: data_(0:n)
  real,intent(in)    :: peso(1:n)
  real               :: suma = 0.0,w(1:n)
  integer            :: i
  !
  !check: data
  w(1:n)=peso(1:n)
  suma=0.0
  do i=1,n
   if(ISNAN(data_(i)))then
    w(i)=0.0
   end if
   suma=suma+w(i)
  end do
  w(1:n)=w(1:n)/suma
  suma=0.0
  do i=1,n
   if(ISNAN(data_(i)).eqv..false.) suma=suma+data_(i)*w(i)
  end do
  data_(0)=suma
  return
 end subroutine average
!
 character(len=20) function str(k)
 !Convert an integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
 end function str
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 subroutine symmetrising(cell_0,rv,vr,vol)
  implicit none
  real,intent(inout) :: cell_0(1:6)
  real,intent(inout) :: rv(3,3),vr(3,3)
  real,intent(in)    :: vol
  cell_0(4)=90.0
  cell_0(5)=90.0
  cell_0(6)=120.0
  cell_0(1) = extract_a(cell_0,vol)
  cell_0(2) = cell_0(1)
  call cell(rv,vr,cell_0)
  return
 end subroutine symmetrising
!
 real function csc_deg(x)
  implicit none
  real :: x
  csc_deg = -(2*sin(x))/(cos(2*x)-1.0)
  return
 end function csc_deg
 real function cot_deg(x)
  implicit none
  real :: x
  cot_deg = -(sin(2*x))/(cos(2*x)-1.0)
  return
 end function cot_deg
 !
 real function extract_a(cell_0,v)
  implicit none
  real :: cell_0(6),v
  real :: pi,DEGTORAD
  real :: alp,bet,gam
  real :: cosa,cosb,cosg
  real :: sing
  pi = ACOS(-1.0)
  DEGTORAD=pi/180.0
  IF(cell_0(4) == 90.0) THEN
    cosa = 0.0
  ELSE
    ALP=cell_0(4)*degtorad
    COSA=cos(ALP)
  ENDIF
  IF(cell_0(5) == 90.0) THEN
    cosb = 0.0
  ELSE
    bet = cell_0(5)*degtorad
    cosb = cos(bet)
  ENDIF
  IF(cell_0(6) == 90.0) then
    sing = 1.0
    cosg = 0.0
  ELSE
    gam = cell_0(6)*degtorad
    sing = sin(gam)
    cosg = cos(gam)
  ENDIF
  extract_a = (sqrt(v)*sqrt(csc_deg(gam)))/( sqrt(cell_0(3))*( &
   1.0-cosb**2-(cosb**2)*(cot_deg(gam)**2)+&
   2.0*cosa*cosb*cot_deg(gam)*csc_deg(gam)-(cosa**2)*(csc_deg(gam)**2))**(1.0/4.0))
 ! Mathematica Solution: 
 ! Sqrt[Csc[\[Gamma]]])/(Sqrt[
 !   c] (1 - Cos[\[Beta]]^2 - Cos[\[Beta]]^2 Cot[\[Gamma]]^2 + 
 !     2 Cos[\[Alpha]] Cos[\[Beta]] Cot[\[Gamma]] Csc[\[Gamma]] - 
 !     Cos[\[Alpha]]^2 Csc[\[Gamma]]^2)^(1/4))
 ! extract_a = (sqrt(v)*sqrt(csc_deg(gam))/ &
 !  ( cell_0(3))*(1.0-cosb**2-(cosb**2)*(cot_deg(gam)**2)+&
 !  2.0*cosa*cosb*cot_deg(gam)*csc_deg(gam)-(cosa**2)*(csc_deg(gam)**2))**(1.0/4.0)) )
 !extract_a= (V*csc_deg(gam))/(cell_0(2)*cell_0(3)*Sqrt( &
 !  1.0-cosb**2-(cosb**2)*(cot_deg(gam)**2)+&
 !  2.0*cosa*cosb*cot_deg(gam)*csc_deg(gam)-(cosa**2)*(csc_deg(gam)**2))**(1.0/4.0))
  return
 end function extract_a
!
 subroutine Quasi_Random_Virtual_Structures(n_cells,n_virtual_systems_max,n_average,n_virtual_systems,virtual_n_Ge)
  use mod_random 
  implicit none
  integer,intent(in)  :: n_cells
  integer,intent(in)  :: n_average
  integer,intent(in)  :: n_virtual_systems_max
  integer,intent(out) :: n_virtual_systems
  real,intent(out)    :: virtual_n_Ge(n_cells,n_virtual_systems_max)
  real    :: a(8)
  integer :: i,j,k,l,n,seed
  integer :: n_T=60
  real    :: x_cell(1:n_cells),x_average,eta
  real    :: temperature=448.0     ! K
  real    :: k_B  = 8.61734E-5     ! eV/K
  real    :: beta, x, partition = 0.0
  real    :: probability_cell(n_cells),free(n_cells)
  real    :: e_1,e_2
  integer :: production_steps
  !
  n_virtual_systems = 0
  production_steps = n_virtual_systems_max
  !
  a(1) = -20.8423
  a(2) = -3.23967
  a(3) =  2.37499
  a(4) = -4.32484
  a(5) =  1.99623
  a(6) =  0.751645
  a(7) = -9.70216
  a(8) =  0.976337
  !
  virtual_n_Ge(1:n_cells,1:n_virtual_systems_max) = real(n_average)
  !
  call init_random_seed(seed)
  partition = 0.0
  beta = 1.0/(k_B*temperature)
  x_average=(n_average/real(n_T))
  do i=1,n_T
   x=(i/real(n_T))
  end do
  x_cell(1:n_cells) = x_average
   do l=1,n_cells
    e_1=Free_Energy(a,x_cell(l))
   end do
  call update_probability(a,beta,n_cells,x_cell,free,partition,probability_cell)
  search_configurations: do k=1,initialisation_steps+production_steps
   i=randint(1,n_cells,seed) 
   j=randint(1,n_cells,seed)
   do while ( i==j )
    i=randint(1,n_cells,seed)
    j=randint(1,n_cells,seed)
   end do
   e_1=Free_Energy(a,x_cell(i)) + Free_Energy(a,x_cell(j))
   n=randint(1,3,seed)
   x_cell(i)=x_cell(i)+n/60.0
   x_cell(j)=x_cell(j)-n/60.0
   e_2=Free_Energy(a,x_cell(i)) + Free_Energy(a,x_cell(j))
   if ( e_1 >= e_2 ) then
    call update_probability(a,beta,n_cells,x_cell,free,partition,probability_cell)
     if(k>initialisation_steps.and.n_virtual_systems_max>n_virtual_systems)then
     n_virtual_systems=n_virtual_systems+1
     do i=1,n_cells
      virtual_n_Ge(i,n_virtual_systems)=n_T*x_cell( i )
     end do
     end if
   else
    eta=R4_UNIFORM(0.0,1.0,SEED)
     if ( eta < exp(-( e_2 - e_1 )/(k_B*temperature)) ) then
      call update_probability(a,beta,n_cells,x_cell,free,partition,probability_cell)
      if(k>initialisation_steps.and.n_virtual_systems_max>n_virtual_systems)then
      n_virtual_systems=n_virtual_systems+1
      do i=1,n_cells
       virtual_n_Ge(i,n_virtual_systems)=n_T*x_cell( i )
      end do
     end if
    else
    x_cell(i)=x_cell(i)-n/60.0
    x_cell(j)=x_cell(j)+n/60.0
    call update_probability(a,beta,n_cells,x_cell,free,partition,probability_cell)
    end if
   end if
   if (n_virtual_systems_max<=n_virtual_systems) exit search_configurations
  end do search_configurations
  return
 end subroutine Quasi_Random_Virtual_Structures
 subroutine update_probability(a,beta,n_cells,x_cell,free,partition,probability_cell)
 implicit none 
  integer           :: j
  real,intent(in)   :: a(8),beta
  integer,intent(in):: n_cells
  real,intent(in)   :: x_cell(n_cells)
  real,intent(out)  :: probability_cell(n_cells)
  real,intent(out)  :: partition,free(n_cells)
  real :: free_min = 0.0
  partition=0.0
  do j=1,n_cells
   free(j)=Free_Energy(a,x_cell(j) )
  end do
  free_min=minval(free)
  do j =1,n_cells
   partition=partition + exp(-beta*(Free_Energy(a,x_cell(j) )-free_min ))
  end do
  do j=1,n_cells
   probability_cell(j)=exp(-beta*(Free_Energy(a,x_cell(j) ) - free_min))/partition
  end do
  return
 end subroutine 
 real function Free_Energy(a,x)
  implicit none
  real, intent(in)   :: a(8)
  real, intent(in)   :: x
  real,parameter     :: muchisimo = 1e12
  if (x<0.or.x>1) then
   Free_Energy = muchisimo
  else
   if(forcefield)then
    Free_Energy = a(1)*(x-a(2))**2 + a(3)*(x-a(4))**3 + &
                a(5)*(x-a(6))**4 + a(7)*(x-a(8))**5
   else
    Free_Energy = 0.0
   end if
  end if
  return
 end function Free_Energy
 !
 subroutine inquire_CIFFile(CIFFilename,cell_0,rv,vr)
  implicit none
  character(len=1024),intent(in) :: CIFFilename
  character(len=80)              :: string_stop_head= "_atom_site_charge"
  integer                        :: uu = 120
  character(len=20)              :: spam
  real,intent(out)               :: cell_0(6)
  real,intent(out)               :: rv(3,3),vr(3,3)
  open(uu,file=CIFfilename,status='old',iostat=ierr)
  if(ierr/=0) then
   write(6,'(a,1x,a)')'[ERROR] opening CIF file', CIFfilename
   stop
  end if
  read_cif: do
   read(uu,'(a)',iostat=ierr) line
   if(ierr/=0) exit read_cif
   if(line(1:14)=="_cell_length_a")then
    read(line,*)spam,cell_0(1)
    cycle read_cif
   end if
   if(line(1:14)=="_cell_length_b")then
    read(line,*)spam,cell_0(2)
    cycle read_cif
   end if
   if(line(1:14)=="_cell_length_c")then
    read(line,*)spam,cell_0(3)
    cycle read_cif
   end if
   if(line(1:17)=="_cell_angle_alpha")then
    read(line,*)spam,cell_0(4)
    cycle read_cif
   end if
   if(line(1:16)=="_cell_angle_beta")then
    read(line,*)spam,cell_0(5)
    cycle read_cif
   end if
   if(line(1:17)=="_cell_angle_gamma")then
    read(line,*)spam,cell_0(6)
    cycle read_cif
   end if
   if(line(1:)==string_stop_head) exit read_cif
  end do read_cif
  call cell(rv,vr,cell_0)
  return
 end subroutine inquire_CIFFile
 !
 SUBROUTINE cell(rv,vr,cell_0)
  implicit none
  real, intent(in)  :: cell_0(6)
  real, intent(out) :: rv(3,3),vr(3,3)
  real, parameter   :: pi = ACOS(-1.0)
  real :: alp,bet
  real :: cosa,cosb,cosg
  real :: gam,sing
  real :: DEGTORAD
  DEGTORAD=pi/180.0
  IF(cell_0(4) == 90.0) THEN
    cosa = 0.0
  ELSE
    ALP=cell_0(4)*degtorad
    COSA=cos(ALP)
  ENDIF
  IF(cell_0(5) == 90.0) THEN
    cosb = 0.0
  ELSE
    bet = cell_0(5)*degtorad
    cosb = cos(bet)
  ENDIF
  IF(cell_0(6) == 90.0) then
    sing = 1.0
    cosg = 0.0
  ELSE
    gam = cell_0(6)*degtorad
    sing = sin(gam)
    cosg = cos(gam)
  ENDIF
  rv(1,1) = cell_0(1)
  rv(1,2) = cell_0(2)*cosg
  rv(1,3) = cell_0(3)*cosb
  rv(2,1) = 0.0
  rv(2,2) = cell_0(2)*sing
  rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
  rv(3,1) = 0.0
  rv(3,2) = 0.0
  rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
  call inverse(rv,vr,3)
  !WRITE(6,'(a)') 'Cell:'
  !WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
  !WRITE(6,'(a)')'Linear Transformation Operator:'
  !DO i=1,3
  ! WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
  !ENDDO
  !WRITE(6,'(a)')'----------------------------------------'
  !WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
  !DO i=1,3
  ! WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
  !ENDDO
  !WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0)
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
 real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  RETURN
 end function
 !
 subroutine print_help()
    print '(a)', '  -h, --help, Print usage information and exit'
    print '(a)', '  -f, --forcefield, Use free energy for Average Calculations'
    print '(a)', '  -e, --QSVS, Quasi-Random-Virtual Structure'
    print '(a)', '  -g, --inquire-CIF'
    print '(a)', 'If you need a bigger system change n_cells or initialisation_steps'
 end subroutine print_help
end program NMR
