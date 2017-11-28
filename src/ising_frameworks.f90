program main
 implicit none
 integer                    :: op1,ii,deg,configuration_file=456
 integer                    :: i,j,k,ijk,l,m,n,lmn,jj,kk,ll,mm,nn,o
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer                    :: n_Ge = 12
 integer                    :: MC_steps = 100
 integer,parameter          :: n_configurations = 0
 character(len=120)         :: file_name,line
 real,parameter             :: temperature = 2000.0 !448.0
 integer, parameter         :: NOPMAX=10000
 integer                    :: k_max_1,k_max_2,k_max_3,k_max_4
 real                       :: ener_0=-7745.86721305, epsilon_, energy_123,energy
 real                       :: ener_end_4 = 0 !-7716.685
 !
 character(len=4)   :: lab(n_atoms)
 real                       :: ener_1(n_atoms) = 0.0
 real                       :: ener_2(n_atoms,n_atoms) = 0.0
 real                       :: ener_3(n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: ener_4(n_atoms,n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: minener_4(1:4),choose_one
 real                       :: cell_0(1:6)
 logical                    :: MC_flag = .false.,no_presente=.false.,restart_file=.false.
 logical                    :: debug = .false.
 integer,allocatable        :: RestartPosition(:,:)
 real,dimension(NOPMAX,3,3) :: mgroup1
 real,dimension(NOPMAX,3)   :: vgroup1
 real,dimension(0:n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(0:n_configurations,n_atoms,2) :: label
 character(len=4),dimension(n_atoms) :: label_debug
! arguments in line:
 integer                                      :: num_args = 0
 character(len=100),dimension(:), allocatable :: args
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
    case('-n','--n_Ge')
     read(args(i+1),*) n_Ge
    case('-mc','--MC_Cycles')
     read(args(i+1),*) MC_steps
     MC_flag = .true.
    case('-r','--restart')
     restart_file=.true.
    case('-w4b','--without-four-body')
     no_presente=.true.
   end select
  end do
 end if
 if(restart_file) then
  allocate(RestartPosition(n_Ge,3))
 end if
 write(6,*) '# substitutions Si/Ge: ', n_Ge
 do ii=0,n_configurations
  !write(6,'(i4.4)') ii
  write (file_name, '( "src/c", I4.4, ".gin" )' )0
  !gulp-1-subs/c0000.gin
  write(6,*)'Configuration',ii, file_name
  write(6,*)'==============='
  open(unit=configuration_file,file=file_name,status='old',iostat=err_apertura)
  if(err_apertura/=0) stop "[ERROR] reading file // file_name // "
  reading_file: do
   READ (configuration_file,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit reading_file
   if(line(1:5)=='cell ') then
    write(6,*)'Cell parameters'
    read(configuration_file,'(A)',iostat=err_apertura) line
    if ( err_apertura /= 0 ) exit reading_file
    read(line,*) ( cell_0(j),j=1,6 )
    write(6,*) ( cell_0(j),j=1,6 )
   end if
   if(line(1:5)==' frac') then
    write(6,*)'Atom coordinates'
    do i=1,n_atoms
     !Si   core    0.0298000    0.8771000    0.4669700
     read (configuration_file,'(A)',iostat=err_apertura) line
     if ( err_apertura /= 0 ) exit reading_file
     if(line(1:7)=='species') exit reading_file
     read(line,*)(label(ii,i,j),j=1,2),(cryst_coor(ii,i,j),j=1,3)
     write(6,'(a4,1x,a4)') (label(ii,i,j),j=1,2)
     write(6,*) (cryst_coor(ii,i,j),j=1,3)
    end do
   end if
  end do reading_file
  close(configuration_file)
 end do 
!
 write(6,*)'Read relation between unique and dependent configurations'
! {{{ 
 write(6,*)'One substitution:'
 k_max_1 = 0
 open(unit=111,file='gulp-1-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=112,file='gulp-1-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=113,file='gulp-1-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read(111,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
!
 read(112,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_1 !<-N.configuraciones totales
 ener_1(1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop
!
 if(n_Ge>=1)then
 read_matrix_1: do
  read(111,'(a)',iostat=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_1
  read(line,*)i,deg,k
  read(113,*,iostat=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_1
  ener_1(k)=epsilon_-ener_0
  write(6,*)k,ener_1(k) !/real(deg)
  do j=2,deg
   READ (112,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_1
   read(line,*)l,m,i,n
   ener_1(i)=ener_1(k)
   write(6,*)k,i,ener_1(i) !/real(deg)
  end do
 end do read_matrix_1
 end if
 close(111)
 close(112)
 close(113)
! }}}
 write(6,*)'Two substitutions:'
 k_max_2 = 0
 open(unit=121,file='gulp-2-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=122,file='gulp-2-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=123,file='gulp-2-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"

 read(121,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(122,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_2 !<-N.configuraciones totales
 ener_2(1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop
 if(n_Ge>=2)then
 read_matrix_2: do
  read(121,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_2
  read(line,*)k,deg,i,j
  READ(123,*,IOSTAT=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_2
  if(n_Ge<=2)then
   choose_one = 1.0
  else
   choose_one = real(n_Ge-1)
  end if
  epsilon_ = (epsilon_-ener_0-ener_1(i)-ener_1(j))/choose_one
   !real(choose(n_Ge,2))
  ener_2(i,j)= epsilon_
  ener_2(j,i)= ener_2(i,j)
  write(6,*)i,j,ener_2(i,j)
  do k=2,deg
   read(122,'(A)',iostat=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_2
   ii=ii+1
   read(line,*)jj,kk,l,m,n
   ener_2(l,m) = epsilon_
   ener_2(m,l) = epsilon_
   write(6,*)i,j,ener_2(l,m)
  end do
 end do read_matrix_2
 end if
 close(121)
 close(122)
 close(123)
!!!
 write(6,*)'3 substitutions:'
 k_max_3 = 0
 open(unit=131,file='gulp-3-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=132,file='gulp-3-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=133,file='gulp-3-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read(131,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(132,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_3 !<-N.configuraciones totales
 ener_3(1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop 
 if( n_Ge>=3)then
 read_matrix_3: do
  read(131,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_3
  read(line,*)ijk,deg,i,j,k
  read(133,*,IOSTAT=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_3
  if(n_Ge<=3)then
   choose_one=1.0
  else
   choose_one=(n_Ge-2)*(n_Ge-3)/2.0
  end if
  epsilon_ = (epsilon_-ener_0-ener_1(i)-ener_1(k)-ener_1(j)-ener_2(i,j)-&
   ener_2(j,k)-ener_2(i,k))/choose_one
  !real(choose(n_Ge,3))
  ener_3(i,j,k) = epsilon_
  ener_3(j,i,k) = epsilon_
  ener_3(i,k,j) = epsilon_
  ener_3(k,j,i) = epsilon_
  ener_3(k,i,j) = epsilon_
  ener_3(j,k,i) = epsilon_
  !deg_3(i,j,k) = deg
  !deg_3(j,i,k) = deg
  !deg_3(i,k,j) = deg
  !deg_3(k,j,i) = deg
  !deg_3(k,i,j) = deg
  !deg_3(j,k,i) = deg
  write(6,*)i,j,k,ener_3(i,j,k),choose_one
  do k=2,deg
   READ (132,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_3
   read(line,*)jj,kk,l,m,n,ii
   ener_3(l,m,n) = epsilon_
   ener_3(m,l,n) = epsilon_
   ener_3(l,n,m) = epsilon_
   ener_3(n,m,l) = epsilon_
   ener_3(m,n,l) = epsilon_
   ener_3(n,l,m) = epsilon_
   !deg_3(l,m,n) = deg
   !deg_3(m,l,n) = deg
   !deg_3(l,n,m) = deg
   !deg_3(n,m,l) = deg
   !deg_3(m,n,l) = deg
   !deg_3(n,l,m) = deg
   write(6,*)i,j,k,ener_3(l,m,n),choose_one
  end do
 end do read_matrix_3
 end if
 close(131)
 close(132)
 close(133)
 write(6,*)'4 substitutions:'
 k_max_4 = 0
 open(unit=141,file='gulp-4-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=142,file='gulp-4-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=143,file='gulp-4-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read(141,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(142,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_4 !<-N.configuraciones totales
 ener_4(1:n_atoms,1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop 
 if (n_Ge>=4.or.no_presente.eqv..false.)then
 read_matrix_4: do
  read(141,'(A)',IOSTAT=err_apertura) line
  if( err_apertura /= 0 ) exit read_matrix_4
  read(line,*)ijk,deg,i,j,k,l
  no_presente=.true.
  check_energy_4: do
   read(143,*,iostat=err_apertura) epsilon_,jj
   if( err_apertura /= 0 )then
    !rewind(143)
    no_presente=.true.
    exit check_energy_4
   end if
   if(jj==ijk)then
    no_presente = .false.
    if(n_Ge<=4)then
     choose_one=1.0
    else
     choose_one = (n_Ge-3)*(n_Ge-4)/2.0
    end if
    epsilon_ = (epsilon_-ener_0-ener_1(i)-ener_1(j)-ener_1(k)-ener_1(l)-&
     ener_2(i,j)-ener_2(j,k)-ener_2(k,l)-ener_2(i,l)-ener_2(j,l)-ener_2(k,l)-&
     ener_3(i,j,k)-ener_3(i,k,j)-ener_3(j,i,k)-ener_3(j,k,i)-ener_3(k,i,j)-ener_3(k,j,i)-&
     ener_3(i,j,l)-ener_3(i,l,j)-ener_3(j,i,l)-ener_3(j,l,i)-ener_3(l,i,j)-ener_3(l,j,i)-&
     ener_3(i,l,k)-ener_3(i,k,l)-ener_3(l,i,k)-ener_3(l,k,i)-ener_3(l,k,i)-ener_3(k,i,l)-&
     ener_3(k,l,i)-ener_3(l,j,k)-ener_3(l,k,j)-ener_3(j,l,k)-ener_3(j,k,l)-ener_3(k,l,j))/choose_one
    write(6,*)i,j,k,l,epsilon_,choose_one
    exit check_energy_4
   end if
  end do check_energy_4
  rewind(143) 
  if(no_presente)then
   epsilon_=ener_end_4
   !do i=1,n_atoms
   ! lab(i) = label(ii,i,1)
   !end do
   !epsilon_ = energy(n_atoms,n_T_atoms,lab,ener_0,ener_1,ener_2,ener_3,ener_4)
   !
   epsilon_ = (epsilon_-ener_0-ener_1(i)-ener_1(k)-ener_1(j)-ener_1(l)-ener_2(i,j)-ener_2(j,k)-&
     ener_2(i,k)-ener_2(i,l)-ener_2(j,l)-ener_2(k,l)-&
     ener_3(i,j,k)-ener_3(i,k,j)-ener_3(j,i,k)-ener_3(j,k,i)-ener_3(k,i,j)-ener_3(k,j,i)-&
     ener_3(i,j,l)-ener_3(i,l,j)-ener_3(j,i,l)-ener_3(j,l,i)-ener_3(l,i,j)-ener_3(l,j,i)-&
     ener_3(i,l,k)-ener_3(i,k,l)-ener_3(l,i,k)-ener_3(l,k,i)-ener_3(l,k,i)-ener_3(k,i,l)-&
     ener_3(k,l,i)-ener_3(l,j,k)-ener_3(l,k,j)-ener_3(j,l,k)-ener_3(j,k,l)-ener_3(k,l,j))/choose_one
   !epsilon_ = 0.0
  end if
  if( epsilon_ == 0 .and. ener_4(i,j,k,l) /= 0 ) stop
  ener_4(i,j,k,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(i,j,l,k) /= 0 ) stop
  ener_4(i,j,l,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(i,k,l,j) /= 0 ) stop
  ener_4(i,k,l,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(i,k,j,l) /= 0 ) stop
  ener_4(i,k,j,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(i,l,j,k) /= 0 ) stop
  ener_4(i,l,j,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(i,l,k,j) /= 0 ) stop
  ener_4(i,l,k,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,i,k,l) /= 0 ) stop
  ener_4(j,i,k,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,i,l,k) /= 0 ) stop
  ener_4(j,i,l,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,k,i,l) /= 0 ) stop
  ener_4(j,k,i,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,k,l,i) /= 0 ) stop
  ener_4(j,k,l,i) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,l,k,i) /= 0 ) stop
  ener_4(j,l,k,i) = epsilon_
  if( epsilon_ == 0 .and. ener_4(j,l,i,k) /= 0 ) stop
  ener_4(j,l,i,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,i,j,l) /= 0 ) stop
  ener_4(k,i,j,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,i,l,j) /= 0 ) stop
  ener_4(k,i,l,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,j,i,l) /= 0 ) stop
  ener_4(k,j,i,l) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,j,l,i) /= 0 ) stop
  ener_4(k,j,l,i) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,l,i,j) /= 0 ) stop
  ener_4(k,l,i,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(k,l,j,i) /= 0 ) stop
  ener_4(k,l,j,i) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,i,j,k) /= 0 ) stop
  ener_4(l,i,j,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,i,k,j) /= 0 ) stop
  ener_4(l,i,k,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,j,i,k) /= 0 ) stop
  ener_4(l,j,i,k) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,j,k,i) /= 0 ) stop
  ener_4(l,j,k,i) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,k,i,j) /= 0 ) stop
  ener_4(l,k,i,j) = epsilon_
  if( epsilon_ == 0 .and. ener_4(l,k,j,i) /= 0 ) stop
  ener_4(l,k,j,i) = epsilon_
  !deg_4(i,j,k,l) = deg
  !deg_4(i,j,l,k) = deg
  !deg_4(i,k,l,j) = deg
  !deg_4(i,k,j,l) = deg
  !deg_4(i,l,j,k) = deg
  !deg_4(i,l,k,j) = deg
  !deg_4(j,i,k,l) = deg
  !deg_4(j,i,l,k) = deg
  !deg_4(j,k,i,l) = deg
  !deg_4(j,k,l,i) = deg
  !deg_4(j,l,k,i) = deg
  !deg_4(j,l,i,k) = deg
  !deg_4(k,i,j,l) = deg
  !deg_4(k,i,l,j) = deg
  !deg_4(k,j,i,l) = deg
  !deg_4(k,j,l,i) = deg
  !deg_4(k,l,i,j) = deg
  !deg_4(k,l,j,i) = deg
  !deg_4(l,i,j,k) = deg
  !deg_4(l,i,k,j) = deg
  !deg_4(l,j,i,k) = deg
  !deg_4(l,j,k,i) = deg
  !deg_4(l,k,i,j) = deg
  !deg_4(l,k,j,i) = deg
  write(6,*)i,j,k,l,ener_4(i,j,k,l)
  ll=l
  do kk=2,deg
   READ (142,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_4
   read(line,*)jj,lmn,l,m,n,o,ii
!   24 permutaciones
    ener_4(l,m,n,o) = epsilon_
    ener_4(l,m,o,n) = epsilon_
    ener_4(l,n,o,m) = epsilon_
    ener_4(l,n,m,o) = epsilon_
    ener_4(l,o,m,n) = epsilon_
    ener_4(l,o,n,m) = epsilon_
    ener_4(m,l,n,o) = epsilon_
    ener_4(m,l,o,n) = epsilon_
    ener_4(m,n,l,o) = epsilon_
    ener_4(m,n,o,l) = epsilon_
    ener_4(m,o,n,l) = epsilon_
    ener_4(m,o,l,n) = epsilon_
    ener_4(n,l,m,o) = epsilon_
    ener_4(n,l,o,m) = epsilon_
    ener_4(n,m,l,o) = epsilon_
    ener_4(n,m,o,l) = epsilon_
    ener_4(n,o,l,m) = epsilon_
    ener_4(n,o,m,l) = epsilon_
    ener_4(o,l,m,n) = epsilon_
    ener_4(o,l,n,m) = epsilon_
    ener_4(o,m,l,n) = epsilon_
    ener_4(o,m,n,l) = epsilon_
    ener_4(o,n,l,m) = epsilon_
    ener_4(o,n,m,l) = epsilon_
!
    !deg_4(l,m,n,o) = deg
    !deg_4(l,m,o,n) = deg
    !deg_4(l,n,o,m) = deg
    !deg_4(l,n,m,o) = deg
    !deg_4(l,o,m,n) = deg
    !deg_4(l,o,n,m) = deg
    !deg_4(m,l,n,o) = deg
    !deg_4(m,l,o,n) = deg
    !deg_4(m,n,l,o) = deg
    !deg_4(m,n,o,l) = deg
    !deg_4(m,o,n,l) = deg
    !deg_4(m,o,l,n) = deg
    !deg_4(n,l,m,o) = deg
    !deg_4(n,l,o,m) = deg
    !deg_4(n,m,l,o) = deg
    !deg_4(n,m,o,l) = deg
    !deg_4(n,o,l,m) = deg
    !deg_4(n,o,m,l) = deg
    !deg_4(o,l,m,n) = deg
    !deg_4(o,l,n,m) = deg
    !deg_4(o,m,l,n) = deg
    !deg_4(o,m,n,l) = deg
    !deg_4(o,n,l,m) = deg
    !deg_4(o,n,m,l) = deg
    write(6,*)i,j,k,ll,ener_4(l,m,n,o)
  end do
 end do read_matrix_4
 end if
 close(141)
 close(142)
 close(143)
! {{{
 if(debug)then
  write(6,*) 'debugging'
  open(unit=111,file='gulp-4-subs/ENERGIES',status='old',iostat=err_apertura)
  if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
  open(222,file="gulp-4-subs/OUTSOD_debug")
  label_debug(1:n_atoms)="Si  "
  ccc: do
   read(111,*,iostat=err_apertura) epsilon_,m
   if( err_apertura /= 0 ) exit ccc
   !write(6,*) epsilon_,m
   open(unit=112,file='gulp-4-subs/OUTSOD',status='old',iostat=err_apertura)
   read(112,*,iostat=err_apertura) line
   ddd: do
    read(112,*,iostat=err_apertura) ijk,deg,i,j,k,l
    if( err_apertura /= 0 ) exit ddd
    if(ijk==m)then
    label_debug(i)="Ge  "
    label_debug(j)="Ge  "
    label_debug(k)="Ge  "
    label_debug(l)="Ge  "
    !
     choose_one = ener_0+ener_1(i)+ener_1(j)+ener_1(k)+ener_1(l)+ener_2(i,j)+ener_2(j,k)+&
       ener_2(i,k)+  ener_2(i,l)+  ener_2(j,l)+  ener_2(k,l)+&
     ener_3(i,j,k)+ener_3(i,k,j)+ener_3(j,i,k)+ener_3(j,k,i)+ener_3(k,i,j)+ener_3(k,j,i)+&
     ener_3(i,j,l)+ener_3(i,l,j)+ener_3(j,i,l)+ener_3(j,l,i)+ener_3(l,i,j)+ener_3(l,j,i)+&
     ener_3(i,l,k)+ener_3(i,k,l)+ener_3(l,i,k)+ener_3(l,k,i)+ener_3(l,k,i)+ener_3(k,i,l)+&
     ener_3(k,l,i)+ener_3(l,j,k)+ener_3(l,k,j)+ener_3(j,l,k)+ener_3(j,k,l)+ener_3(k,l,j)+&
     ener_3(j,l,k)
     write(222,*)epsilon_,energy(n_atoms,n_T_atoms,label_debug,ener_0,ener_1,ener_2,ener_3,ener_4),choose_one
     label_debug(1:n_atoms)="Si  "
    end if
   end do ddd
   rewind(112)
   close(112)
  end do ccc
  close(111)
  close(222)
 end if
!
 if(restart_file)then
  ! RestartPosition(i,ii),i=1,n_Ge )
  do j=1,n_configurations
   write(6,*)"======================"
   write(6,*)"Selected configuration:"
   read(5,*) ( RestartPosition(i,j),i=1,n_Ge )
   write(6,*)"Configuration", j
   write(6,*)( RestartPosition(i,j),i=1,n_Ge )
   do i=1,n_Ge
    label(0,RestartPosition(i,j),1)="Ge  "
   end do
   call geometrical_properties(n_atoms,n_T_atoms,n_Ge,cell_0,cryst_coor,n_configurations,&
       label,ener_0,ener_1,ener_2,ener_3,ener_4,0) 
  end do
   do i=1,n_Ge
    label(0,RestartPosition(i,j),1)="Si  "
   end do
 else
 if(MC_flag)then
    call farthrest_nodes_subtitutions(n_atoms,n_T_atoms,n_Ge,ener_0,ener_1,&
     ener_2,ener_3,ener_4,cell_0,cryst_coor,&
     n_configurations,label,MC_steps,temperature,0)
 else
  do ii=1,n_configurations
   call geometrical_properties(n_atoms,n_T_atoms,n_Ge,cell_0,cryst_coor,n_configurations,&
       label,ener_0,ener_1,ener_2,ener_3,ener_4,ii)
  end do
 end if
 end if
 write(6,'(a)')'Finish simulation'
 stop
contains
  real (kind=4) function choose (n, k)
    implicit none
    integer, intent (in) :: n,k
    choose = factorial(n) / (factorial(k) * factorial(n - k))
    if (choose<1) choose=1
    return
  end function choose
  real (kind=4) function factorial(nn)
   implicit none
   integer, intent(in) :: nn
   integer             :: j
   factorial = 1
   if(nn<=0) then
    factorial = 1
   else
    do j=1,nn
     factorial = factorial*j
    end do
   end if
   return
  end function factorial
 subroutine print_help()
    print '(a)', '  -h,  --help, Print usage information and exit'
    print '(a)', '  -n,  --n_Ge, Use free energy for Average Calculations'
    print '(a)', '  -mc, --MC_Cycles, Quasi-Random-Virtual Structure'
    print '(a)', '  -r,  --restart, From Restart File'
    print '(a)', 'If you need a bigger system change n_cells or initialisation_steps'
 end subroutine print_help
end program main
