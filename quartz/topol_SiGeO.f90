program topol
 implicit none
 integer                      :: i,j,k,ierr = 0
 integer                      :: n_atoms = 0
 real,allocatable             :: xcryst(:,:),xcarte(:,:),dist(:,:)
 real                         :: r1(3),r2(3),cell_0(6),rv(3,3),vr(3,3),r
 integer,allocatable          :: adj(:,:)
 character(len=4),allocatable :: label(:,:)
 character(len=80)            :: line
!
! open gulp file:
 fileopen_gin: if(ierr==0)then
  do1: DO
    READ(5,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT do1
    IF(line(1:4)=='cell') EXIT do1
  ENDDO do1
  READ(5,*)(cell_0(k), k=1,6 )
  CALL cell(rv,vr,cell_0)
  READ(5,'(A)',IOSTAT=IERR) line ! frac
  DO
    READ(5,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    IF ((line(6:9)=='core'.or.line(6:9)=='shel')) n_atoms = n_atoms + 1
!species
    if(line(1:7)=="species") exit
  END DO
 end if fileopen_gin
 rewind(5)
!allocate
 allocate(xcryst(1:n_atoms,0:3),stat=ierr)
 allocate(label(1:n_atoms,1:2),stat=ierr)
 ALLOCATE(dist(n_atoms,n_atoms),STAT=IERR)
 allocate(adj(n_atoms,n_atoms),STAT=IERR)
!
 if(ierr/=0) stop '[error] allocate variable'
 do2: do
  READ(5,'(A)',IOSTAT=IERR) line
  IF( IERR /= 0 ) EXIT do2
  IF (line(1:4)=='frac'.or.line(1:5)==' frac') EXIT do2
 end do do2
 do k=1,n_atoms
  READ(5,'(A)',IOSTAT=IERR) line
  read(line,*) label(k,1),label(k,2),(xcryst(k,j),j=1,3)
  !READ(line(1:4),*)  label(k,1)
  !READ(line(6:9),*)  label(k,2)
!  Coordinates (gulp code: wcoord.f90 )
  !READ(line(13:22),*) xcryst(k,1)
  !READ(line(26:35),*) xcryst(k,2)
  !READ(line(39:48),*) xcryst(k,3)
  DO j=1,3
   xcryst(k,j)=MOD(xcryst(k,j)+100.0,1.0)
  END DO
  WRITE(6,*)label(k,1),label(k,2),(xcryst(k,j),j=1,3)
!
 end do
! distance matrix:
! ===============
  DO i=1,n_atoms
     DO j=i,n_atoms
        IF(i==j)dist(i,j)=0.0
        FORALL ( k=1:3 )
         r1(k)= xcryst(j,k)
         r2(k)= xcryst(i,k)
        END FORALL
        CALL make_distances(cell_0,r1,r2,rv,r)
        dist(i,j)=r
        dist(j,i)=dist(i,j)
        if( r >= 1.0 .and. r < 2.0 ) then
         adj(i,j)=1.0
         adj(j,i)=adj(i,j)
        end if
     END DO
  END DO
! labeleo:
  do i=1,n_atoms
   if(label(i,1)=="O   ")then 
    k=0
    do j=1,n_atoms
     if(adj(i,j)==1.0.and.label(j,1)=="Si  ") k=k+1
     if(adj(i,j)==1.0.and.label(j,1)=="Ge  ") k=k+2
    end do
    select case (k)
     case (2)
      write(6,*)' Si-O-Si', k
      ! Si-O-Si (+1,+1)
      label(i,1)="O2  "
     case (3)
      write(6,*)' Si-O-Ge', k
      ! Si-O-Ge (+1,+2) 
      label(i,1)="O3  "
     case (4)
   !   write(6,*)' Ge-O-Ge', k
      ! Ge-O-Ge (+2,+2)
      label(i,1)="O1  "
     case default
      label(i,1)=label(i,1)
    end select
   end if
  end do
!
  conectivity: DO i=1,n_atoms
     k=0
     DO j=1,n_atoms
      k = k + adj(i,j)
     ENDDO
     xcryst(i,0) = k
  END DO conectivity
  write(6,'(a)')'========================================'
  WRITE(6,'(1000(I1))')(int(xcryst(k,0)),k=1,n_atoms) ! }}
  write(6,'(a)')'========================================'
  call write_gin(cell_0,xcryst,n_atoms,label)
! 
  
  STOP ':P'
  CONTAINS
  SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
   IMPLICIT NONE
   REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
   REAL,    intent(out) :: dist
  !REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
   REAL                 :: d_image(1:125),image(3,125)
   INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
   REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
   k=0
   do l=-2,2
    do m=-2,2
      do n=-2,2
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
    enddo
   enddo
   dist=MINVAL(d_image)
   RETURN
  END SUBROUTINE
  REAL FUNCTION DISTANCE(atom,ouratom,rv)
   IMPLICIT NONE
   INTEGER :: j
   REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
   REAL :: rv(3,3)
   FORALL ( j=1:3 )
    o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
    o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
    dist(j) = o_ouratom(j) - o_atom(j)
   END FORALL
   DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
  END FUNCTION
 SUBROUTINE cell(rv,vr,cell_0)
! ======================
! GULP
! ======================
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
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
 RETURN
 END SUBROUTINE cell
 SUBROUTINE uncell(rv,cell_0)
! ======================
! GULP
! ======================
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
 SUBROUTINE inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!==========================================================
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
 SUBROUTINE write_gin(cell_0,xcryst,n_atoms,label)
  ! Catlow potential
  IMPLICIT NONE
  INTEGER,          intent(in) :: n_atoms
  CHARACTER (LEN=4),intent(in) :: label(1:n_atoms,2)
  REAL,             intent(in) :: xcryst(1:n_atoms,0:3),cell_0(1:6)
  OPEN(999,file='out.gin')
  WRITE(999,'(A)')'opti conv qok'
  write(999,'(a)')'stepmx     0.01'
  write(999,'(a)')'maxcyc opt 5000'
  WRITE(999,'(A)')'cell'
  WRITE(999,'(6f10.5)') (cell_0(j) , j=1,6)
  WRITE(999,'(A)')'fractional'
  do i=1,n_atoms
   WRITE(999,'(a,1x,a,1x,3f10.5)')label(i,1),label(i,2),xcryst(i,1),xcryst(i,2),xcryst(i,3)
  enddo
  write(999,'(a)')'species'
  write(999,'(a)')'Si    core  4.00000'
  write(999,'(a)')'O1    core  1.733957'
  write(999,'(a)')'O1    shel -3.733957'
  write(999,'(a)')'O2    core  0.870733'
  write(999,'(a)')'O2    shel -2.870733'
  write(999,'(a)')'O3    core  1.330431'
  write(999,'(a)')'O3    shel -3.330431'
  write(999,'(a)')'F     core  0.56'
  write(999,'(a)')'F     shel -1.56'
  write(999,'(a)')'Ge    core  4.0'
  write(999,'(a)')'buck'
  write(999,'(a)')'Si core    O1 shel     1315.2478 0.317759 10.141118 0.0 16.0'
  write(999,'(a)')'Si core    O2 shel     1315.2478 0.317759 10.141118 0.0 16.0'
  write(999,'(a)')'Si core    O3 shel     1315.2478 0.317759 10.141118 0.0 16.0'
  write(999,'(a)')'Ge core    O1 shel     1497.3996 0.325646 16.808599 0.0 16.0'
  write(999,'(a)')'Ge core    O2 shel     1497.3996 0.325646 16.808599 0.0 16.0'
  write(999,'(a)')'Ge core    O3 shel     1497.3996 0.325646 16.808599 0.0 16.0'
  write(999,'(a)')'O1 shel    O1 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'O1 shel    O2 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'O1 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'O2 shel    O2 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'O2 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'O3 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0'
  write(999,'(a)')'Si core     F shel     976.82887 0.282000 0.0       0.0 16.0'
  write(999,'(a)')'Ge core     F shel     681.47288 0.320000 0.0       0.0 16.0'
  write(999,'(a)')'F  shel     F shel     540.39761 0.262490 0.0       0.0 16.0'
  write(999,'(a)')'O1 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0'
  write(999,'(a)')'O2 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0'
  write(999,'(a)')'O3 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0'
  write(999,'(a)')'spring'
  write(999,'(a)')'O2     75.969800'
  write(999,'(a)')'O3    128.142790'
  write(999,'(a)')'O1    180.315770'
  write(999,'(a)')'F      33.452757'
  write(999,'(a)')'three'
  write(999,'(a)')'Si core O1 shel O1 shel 1.2614 109.47 1.9 1.9 3.5'
  write(999,'(a)')'Si core O1 shel O2 shel 1.2614 109.47 1.9 1.9 3.5'
  write(999,'(a)')'Si core O2 shel O2 shel 1.2614 109.47 1.9 1.9 3.5'
  write(999,'(a)')'cutp             16.0 1.0'
  write(999,'(a)')'cuts 0.8'
  write(999,'(a)')'output cif out.cif'
  write(999,'(a)')'dump every 1     out.res'
  RETURN
 END SUBROUTINE write_gin
end program topol
