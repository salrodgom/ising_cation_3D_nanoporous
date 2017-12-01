program main
 implicit none
 integer               :: i,j,k,deg,paso,ierr
 real                  :: energy, oto,ot,factor=0.0, sumfactor=0.0
 real                  :: energy0=99999999999.0
 real                  :: energy_average = 0.0
 real                  :: oto_average = 0.0
 real                  :: ot_average  = 0.0
 real                  :: tt,tot
 real                  :: tt_average  = 0.0
 real                  :: tot_average = 0.0
 real                  :: rv(3,3),vr(3,3),cell_0(6),volume_0=0.0
 real                  :: cell_1_tmp = 0.0
 real                  :: cell_1_average = 0.0
 real                  :: cell_1_average_qrt=0.0
 real                  :: cell_3_average = 0.0
 real,parameter        :: k_b = 8.6173303e-5
 real,parameter        :: temperature = 450.0
 real                  :: beta = 1.0/(k_b*temperature)
 character(len=40)     :: filename
 read_file: do
  read(5,*,iostat=ierr)i,j,energy,(cell_0(k),k=1,6),ot,oto,tt,tot,filename
  if (ierr/=0) exit read_file
  if(energy<=energy0)energy0=energy
 end do read_file
 rewind(5)
 read_file2: do
  read(5,*,iostat=ierr)i,j,energy,(cell_0(k),k=1,6),ot,oto,tt,tot,filename
  if (ierr/=0) exit read_file2
  factor=exp(-beta*(energy-energy0))
  sumfactor=sumfactor+factor
 end do read_file2
 !write(6,*) sumfactor,energy0
 rewind(5)
 read_again: do
  read(5,*,iostat=ierr)i,j,energy,(cell_0(k),k=1,6),ot,oto,tt,tot,filename
  if (ierr/=0) exit read_again
  call cell(rv,vr,cell_0)
  volume_0= volume(rv)
  cell_0(4)=90.0
  cell_0(5)=90.0
  cell_0(6)=120.0
  cell_1_tmp= extract_a(cell_0,volume_0)
  call cell(rv,vr,cell_0)
  factor=exp(-beta*(energy-energy0))/sumfactor
  tt_average = tt_average + tt*factor
  tot_average = tot_average + tot*factor
  ot_average = ot_average + ot*factor
  oto_average = oto_average + oto*factor
  energy_average= energy_average+energy*factor
  cell_1_average= cell_1_average+cell_1_tmp*factor
  cell_3_average= cell_3_average+cell_0(3)*factor
 end do read_again
 rewind(5)
 read_again2: do
  read(5,*,iostat=ierr)i,j,energy,(cell_0(k),k=1,6),filename
  if (ierr/=0) exit read_again2
  call cell(rv,vr,cell_0)
  cell_0(4)=90.0
  cell_0(5)=90.0
  cell_0(6)=120.0
  cell_1_tmp= extract_a(cell_0,volume_0)
  !cell_1_tmp=(0.09643*volume_0)/cell_0(3)
  call cell(rv,vr,cell_0)
  factor=exp(-beta*(energy-energy0))/sumfactor
  !cell_1_average_qrt=cell_1_average_qrt+((0.5*(cell_0(1)+cell_0(2))-cell_1_average)**2)*factor
  cell_1_average_qrt=cell_1_average_qrt+((cell_1_tmp-cell_1_average)**2)*factor
 end do read_again2
 cell_1_average_qrt=sqrt(cell_1_average_qrt)
 write(6,*)i,energy_average,cell_1_average,cell_1_average_qrt,cell_3_average,volume(rv),&
           ot_average,oto_average,tt_average,tot_average
contains
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
real function extract_a(cell_0,v)
 implicit none
 real :: cell_0(6),v
 real :: pi,DEGTORAD
 real :: alp,bet,gam
 real :: cosa,cosb,cosg
 real :: sina,sinb,sing
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
SUBROUTINE cell(rv,vr,cell_0)
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
! print*,'Cell:'
! WRITE(*,'(6F14.7)')( cell_0(j), j=1,6 )
! print*,'Box:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
! ENDDO
! WRITE(*,*)'----------------------------------------'
! WRITE(*,*)'bOX:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
! ENDDO
 RETURN
END SUBROUTINE cell
!
SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  REAL :: radtodeg,PI
  PI=ACOS(-1.0)
  radtodeg=180.0/PI
!
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
!  Avoid round off errors for 90.0 and 120.0 degrees
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
!
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
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
!===========================================================
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
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
!
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
end program main
