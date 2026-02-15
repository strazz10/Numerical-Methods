program ising_non_square
implicit none 

integer, parameter :: n_horiz = 2, n_vert = 2, N = 10**8                !!numero celle orizzontali e verticali, numero di campionamenti
integer, parameter :: rows = 4+3*(n_vert-1)+2, columns = 4+3*(n_horiz-1)+2  !!numero di elementi della matrice, ogni esagono sono 4 punti ma uno dopo il primo esagono non è indipendente, i +2 servono per le condizioni al bordo aperte
integer :: i, j, iteration, rowtemp, coltemp, s, norm
real :: beta	          !!temperatura inversa
real :: aux, acceptance, r
integer :: lattice(rows, columns) 
double precision :: totalE, totalMagn

print *, '###########################################################'  !!bellurie
print *, 'Simulation of 2d Ising model, non square lattice (hexagonal), open boundary conditions;'
print *, 'Number of horizontal, vertical cells:', n_horiz, n_vert
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

open (unit = 20, file = 'ising_hex.dat', status = 'unknown')

  print *, ''
  print *, '##################################'
  print *, 'Initial Lattice'
  print *, '##################################'
  print *, ''

lattice = 0          !!inizializzo il reticolo e lo stampo
do j=1,columns
   lattice(1,j) = 0
   lattice(rows,j) = 0
enddo

print*, lattice(1,:)
do i=2,rows-1
   if (i==2 .or. mod(i-1,3)==1) then
   do j=2,columns-1 
      if (j/=2 .and. mod(j-1,3)/=1) then
      lattice(i,j) = 1
      end if
   enddo
   else 
   do j=2,columns-1
      if (j==2 .or. mod(j-1,3)==1) then
      lattice(i,j) = 1
      end if
   enddo
   end if
   print *, lattice(i,:)
enddo
print*, lattice(rows,:)         !!fino a qui ok

norm = 0                        !!fattore di normalizzazione che conta quanti siti sono effettivamente occupati
do i=1,rows
   do j=1,columns
   norm = norm + lattice(i,j)
   enddo
enddo

acceptance = 0       !!seleziono un punto casuale, posso saltare le prime e ultime righe/colonne
do iteration=1,N
   do i=2,rows-1  
      do j=2,columns-1
   2  call random_number(aux)
      coltemp = int(aux*(columns))+1
      call random_number(aux)
      rowtemp = int(aux*(rows))+1
      
      if (lattice(rowtemp, coltemp)==0) then   !!tengo solo i punti del reticolo effettivo (sicuramente lentissimo)
      goto 2
      end if
     
      s = 0         !!somma sui vicini
      
      if (rowtemp-1==1 .or. mod(rowtemp-1,3)==1) then
      s = s + lattice(rowtemp, coltemp+1)+lattice(rowtemp, coltemp-1)
      s = s + lattice(rowtemp+1, coltemp-1)+lattice(rowtemp-1, coltemp-1)
      s = s + lattice(rowtemp-1, coltemp+1)+lattice(rowtemp+1, coltemp+1)
      else
      s = s + lattice(rowtemp+1, coltemp)+lattice(rowtemp-1, coltemp)
      s = s + lattice(rowtemp+1, coltemp-1)+lattice(rowtemp-1, coltemp-1)
      s = s + lattice(rowtemp-1, coltemp+1)+lattice(rowtemp+1, coltemp+1)
      end if

      s = s*lattice(rowtemp,coltemp)
      
      
      r = 0                     !!step metropolis
      if (s < 0) then
      lattice(rowtemp,coltemp) = -lattice(rowtemp,coltemp)
      acceptance = acceptance + 1
      else 
      call random_number(r)
      if (r <= exp(-beta*2*s)) then
      lattice(rowtemp,coltemp) = -lattice(rowtemp,coltemp)
      acceptance = acceptance + 1
      end if
      end if
      enddo
   enddo
write(20,*) totalE(lattice, rows, columns, norm), totalMagn(lattice, rows, columns, norm)         !!misura dopo aver iterato su tutto il reticolo
enddo 

print *, ''
print *, 'beta=', beta
print *, 'acceptance=', acceptance
print *, ''

close(20)

end program ising_non_square

function totalE(lattice, rows, columns, norm) result(E) !!definibile come E, E*(norm/rows*cols), possibilmente sbagliata
implicit none
double precision :: E
integer :: i, j, rows, columns, lattice(rows,columns), temp, norm
temp = 0
do i=2,rows-1
   do j=2,columns-1
   temp = temp - lattice(i,j)*lattice(i,j+1)
   temp = temp - lattice(i,j)*lattice(i+1,j)
   temp = temp - lattice(i,j)*lattice(i+1,j+1)
   temp = temp - lattice(i,j)*lattice(i+1,j-1)
   enddo
enddo
E = dble(temp)/dble(norm)
end function totalE

function totalMagn(lattice, rows, columns, norm) result(m)  !!definibile come m, m*(norm/rows*cols)
implicit none
double precision :: m
integer :: i, j, rows, columns, lattice(rows,columns), temp, norm
temp = 0
do i=2,rows-1
   do j=2,columns-1
   temp = temp + lattice(i,j)
   enddo
enddo
m = dble(temp)/dble(norm)
end function totalMagn









