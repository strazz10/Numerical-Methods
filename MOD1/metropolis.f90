function test_function(x) result(y)
implicit none
real,intent(in) :: x
real :: y, mu, sigma
real,parameter :: pi = 3.14159
mu = 0
sigma = 1
y = 1/(sqrt(2*pi*sigma**2)) * exp(-((x-mu)**2)/(2*sigma**2))
end function test_function

!!!!!!!!!!!!!!!!!!!!!!!!
!!funzioni che danno la probabilità di rifiuto: vanno confrontate con una distribuzione uniforme tra (0,1) che dà effettivamente il risultato SI/NO
function F1(x) result(y) 
implicit none
real,intent(in) :: x
real :: y
y = min(1.0,x)
end function F1

!!!!!!!!!!!!!!!!!!!!!!!!

function F2(x) result(y)
implicit none
real,intent(in) :: x
real :: y
y = x/1+x
end function F2

!!!!!!!!!!!!!!!!!!!!!!!!

subroutine metropolis(n)
implicit none
real, parameter :: delta = 2
real :: state, state_temp, r, interval, F1, F2, test_function
integer :: n, i, acceptance

open (unit = 20, file = 'metropolis_test.dat', status = 'unknown')
write (20,*) '#sampled pdf using metropolis algorithm, iterations =', n

state = 0
state_temp = 0
r = 0
interval = 0
acceptance = 0

call random_number(state) !stato iniziale

do i=1,n
call random_number(interval)
interval = (interval-0.5)*2*delta !intervallo simmetrico
state_temp = state + interval
  if (F1(test_function(state_temp)/test_function(state)) >= 1) then 
  state = state_temp
  acceptance = acceptance+1
  else 
      call random_number(r)
      if (r <= F1(test_function(state_temp)/test_function(state))) then
      state = state_temp
      acceptance = acceptance+1
      end if
  end if
write (20,*) i, state
end do
write (20,*) '#Accettanza = ', acceptance/n

close (20)
return(0)
end subroutine metropolis
!!!!!!!!!!!!!!!!!!!!!!!!

program main
implicit none
integer, parameter :: n = 1000000 !numero di iterazioni

call metropolis(n)
end program main






