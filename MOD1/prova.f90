FUNCTION TEST_FUNCTION(ii) RESULT(jj) !! funzione da integrare o altro
IMPLICIT NONE 
REAL, INTENT(IN) :: ii
REAL :: jj

jj = SQRT(1-(ii)**2) 

END FUNCTION TEST_FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MC_INTEGRAL() !! integrali con estimatori puntuali
IMPLICIT NONE 
INTEGER, PARAMETER :: nn = 1000 !! numero di campionamenti 
DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932563353534d0
INTEGER :: ii
REAL :: SET1(nn), SET2(nn), SET3(nn)
REAL :: ESTIMATE, STDDEV, TEST_FUNCTION

SET1 = 0
CALL RANDOM_NUMBER(SET1)

DO ii=1,nn
SET2(ii) = TEST_FUNCTION(SET1(ii))
SET3(ii) = TEST_FUNCTION(SET1(ii))**2
ENDDO

ESTIMATE = 1/REAL(nn)*(SUM(SET2))
STDDEV = SQRT(1/REAL(nn)*SUM(SET3) - ESTIMATE**2)*(1/SQRT(REAL(nn)))

PRINT *, 'MC INTEGRATION RESULT: ', ESTIMATE
PRINT *, 'MC STD DEV: ', STDDEV
PRINT *, 'EXACT INTEGRATION RESULT: ', pi/4.d0

END SUBROUTINE MC_INTEGRAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine trapezoid_integral() !!metodo "a la andrea sommo"
integer, parameter :: N = 1000
integer :: ii
real, parameter :: interval_low = 0, interval_high = 1 !!estremi di integrazione
real :: delta_x, integral

delta_x = (interval_high-interval_low)/N
integral = 0

do ii=0,N-1 !!forse più efficiente inizializzare prima un array su x, boh
integral = integral + 0.5*delta_x*(test_function(interval_low + ii*delta_x) + test_function(interval_low + (ii+1)*delta_x))
enddo

print *, 'TRAPEZOID INTEGRATION RESULT: ', integral
end subroutine trapezoid_integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MC_GAUSSIAN_GENERATOR() !! trasformata di Box-Muller
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932563353534d0
DOUBLE PRECISION, PARAMETER :: MU = 0.d0, SIGMA = 1 !! media e varianza della distribuzione
INTEGER, PARAMETER:: HOW_MANY = 10000 !! quanti ne vuoi
INTEGER :: ii, jj
DOUBLE PRECISION :: PHI(HOW_MANY), RADIUS(HOW_MANY), X(HOW_MANY), Y(HOW_MANY), Y1(HOW_MANY), Y2(HOW_MANY)

DOUBLE PRECISION, PARAMETER :: BIN_WIDTH = 0.2
INTEGER, PARAMETER :: HISTO_LENGHT = 32 !possibilmente pari
REAL :: HISTO_X(HISTO_LENGHT)
INTEGER :: HISTO_Y(HISTO_LENGHT), kk, ll, nn

X = 0
Y = 0
CALL RANDOM_NUMBER(X)
CALL RANDOM_NUMBER(Y)
RADIUS = 0
PHI = 0
DO ii=1,HOW_MANY
RADIUS(ii) = SQRT(-2*LOG(X(ii)))*SIGMA
PHI(ii) = 2*pi*Y(ii)
ENDDO

Y1 = 0
Y2 = 0
DO jj=1,HOW_MANY
Y1(jj) = RADIUS(jj)*COS(PHI(jj)) + MU
Y2(jj) = RADIUS(jj)*SIN(PHI(jj)) + MU
ENDDO

!PRINT *, 'TWO ARRAYS OF GAUSSIANLY DISTRIBUTED NUMBERS WITH MU,SIGMA:', REAL(MU), REAL(SIGMA) !!da commentare/no se vuoi effettivamente i numeri
!PRINT *, Y1
!PRINT *, Y2

!!plotting part

HISTO_X = 0
HISTO_Y = 0

DO kk=1,HOW_MANY

 DO ll=0,INT(HISTO_LENGHT/2)-1
 IF (ll < Y1(kk)/(BIN_WIDTH) .AND. Y1(kk)/(BIN_WIDTH) < ll+1) THEN
 HISTO_Y(ll+INT(HISTO_LENGHT/2)+1) = HISTO_Y(ll+INT(HISTO_LENGHT/2)+1)+1 !numeri positivi in 2° metà posti
 ELSEIF (-ll-1 < Y1(kk)/(BIN_WIDTH) .AND. Y1(kk)/(BIN_WIDTH) < -ll) THEN
 HISTO_Y(INT(HISTO_LENGHT/2)-ll) = HISTO_Y(INT(HISTO_LENGHT/2)-ll)+1 !numeri negativi 1° metà posti
 ENDIF
	
 IF (ll < Y2(kk)/(BIN_WIDTH) .AND. Y2(kk)/(BIN_WIDTH) < ll+1) THEN
 HISTO_Y(ll+INT(HISTO_LENGHT/2)+1) = HISTO_Y(ll+INT(HISTO_LENGHT/2)+1)+1 !numeri positivi in 2° metà posti
 ELSEIF (-ll-1 < Y2(kk)/(BIN_WIDTH) .AND. Y2(kk)/(BIN_WIDTH) < -ll) THEN
 HISTO_Y(INT(HISTO_LENGHT/2)-ll) = HISTO_Y(INT(HISTO_LENGHT/2)-ll)+1 !numeri negativi 1° metà posti
 ENDIF
 ENDDO

ENDDO

DO nn=1,INT(HISTO_LENGHT/2)
HISTO_X(nn) = (-BIN_WIDTH*((HISTO_LENGHT/2)+1-nn)) + MU
ENDDO
DO nn=INT(HISTO_LENGHT/2)+1,INT(HISTO_LENGHT) 
HISTO_X(nn) = (BIN_WIDTH*(nn-(HISTO_LENGHT/2))) + MU
ENDDO

OPEN (Unit = 20, file = 'gaussdata.dat', status = 'unknown')
WRITE (20,*) '#GENERATED POINTS N_TOT, SIGMA, MU', 2*HOW_MANY, REAL(SIGMA), REAL(MU)
WRITE (20,*) '#BIN          NUMBER PER BIN'
DO nn=1,HISTO_LENGHT 
PRINT *, HISTO_X(nn), HISTO_Y(nn)
WRITE (20,*) HISTO_X(nn), HISTO_Y(nn)
ENDDO
CLOSE (20)

END SUBROUTINE MC_GAUSSIAN_GENERATOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
IMPLICIT NONE 

CALL MC_INTEGRAL()
CALL TRAPEZOID_INTEGRAL()
!CALL MC_GAUSSIAN_GENERATOR()

STOP
END PROGRAM MAIN





