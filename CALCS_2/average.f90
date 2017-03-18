program ave
 implicit none
 real    :: ener(0:60), ener_0, free, entropy(1:59),entro
 real    :: temperature =2000 
 integer :: i,j,k
 open(111,file='energies')
 open(222,file='entropy')
 ener = 0
 do
  read(111,*,iostat=k) j, ener_0
  if (k/=0) exit 
  ener(j)=ener(j)+ener_0/4.0
 end do
 ener(0)=4*ener(0)
 ener(60)=4*ener(60)
 do i =0,60
  if(i/=0.and.i/=60) then
   read(222,*) j,entro
   entropy(j)=entro
  end if
  write(6,*)i,ener(i),entropy(i),ener(i)-entropy(i)*temperature
 end do
end program
