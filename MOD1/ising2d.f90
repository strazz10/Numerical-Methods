program ising2d
implicit none 

integer, parameter :: L = 4, N = 10**8   !!lunghezza reticolo, numero di campionamenti
integer :: i, j, x, y, xtemp, ytemp, lattice(L,L), iteration, temp, s
real :: beta	          !!temperatura inversa
real :: aux, acceptance, r
double precision :: totalE, totalMagn

print *, '###########################################################' !!bellurie
print *, 'Simulation of 2d Ising model, periodic boundary conditions;'
print *, 'Lattice lenght:', L
write(6,"(' Number of iterations:', es10.1)") real(N)
print *, '###########################################################'
1 print *, 'Input beta g.e. 0 to perform a measurement:'
read *, beta

if (beta<0) then
  print *, ''
  print *, '##################################'
  print *, 'Error: beta must be non negative'
  print *, '##################################'
  print *, ''
  goto 1
end if

open (unit = 20, file = 'ising.dat', status = 'unknown')

lattice = 0          !!inizializzo il reticolo
do i=1,L
  do j=1,L
  lattice(i,j) = 1
  enddo
enddo

acceptance = 0       !!seleziono un punto casuale
do iteration=1,N
  do x=1,L
    do y=1,L
    call random_number(aux)
    xtemp = int(aux*L)+1
    call random_number(aux)
    ytemp = int(aux*L)+1
               
    s = 0             !!calcolo la somma dei vicini
    temp = mod(xtemp+1,L)    
    if (temp == 0) then 
    temp = xtemp+1
    end if
    s = s + lattice(temp,ytemp)

    temp = xtemp-1
    if (temp < 1) then 
    temp = L
    end if
    s = s + lattice(temp,ytemp)

    temp = mod(ytemp+1,L)
    if (temp == 0) then 
    temp = ytemp+1
    end if
    s = s + lattice(xtemp,temp)

    temp = ytemp-1
    if (temp < 1) then 
    temp = L
    end if
    s = s + lattice(xtemp,temp)
    s = s*lattice(xtemp,ytemp)
                    
    r = 0                     !!step metropolis
    if (s < 0) then
    lattice(xtemp,ytemp) = -lattice(xtemp,ytemp)
    acceptance = acceptance + 1
    else 
      call random_number(r)
      if (r <= exp(-beta*2*s)) then
      lattice(xtemp,ytemp) = -lattice(xtemp,ytemp)
      acceptance = acceptance + 1
      end if
    end if
    enddo
  enddo
  write(20,*) totalE(lattice, L), totalMagn(lattice, L)         !!misura dopo aver iterato su tutto il reticolo
enddo

print *, ''
print *, 'beta=', beta
print *, 'acceptance=', acceptance/(N*L*L)
print *, ''

close(20)

end program ising2d	

function totalE(lattice, L) result(E)        !!energia per unità di L*L
implicit none
double precision :: E
integer :: x, y, lattice(L,L), L, mog, temp

temp = 0
mog = 0
do x=1,L
  do y=1,L
  mog = mod(x+1, L)
  if (mog==0) then
  mog = x+1
  end if
  temp = temp - lattice(mog,y)*lattice(x,y)
  
  mog = mod(y+1, L)
  if (mog==0) then
  mog = y+1
  end if
  temp = temp - lattice(x,mog)*lattice(x,y)
  
  enddo
enddo
E = dble(temp)/dble(L*L)

end function totalE

function totalMagn(lattice, L) result(m)     !!magnetizzazione
implicit none
double precision :: m
integer :: L, x, y, lattice(L,L), tmp

tmp = 0
do x=1,L
  do y=1,L
  tmp = tmp + lattice(x,y)
  enddo
enddo
m = dble(tmp)/dble(L*L)

end function totalMagn



