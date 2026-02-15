program clock_model
    implicit none                            !!mettendo tutto in double viene più preciso (male non ci fa, per esperienza)
    integer, parameter :: L = 40        
    integer, parameter :: N = 2*10**6      
    integer, parameter :: q = 4              !!parametro del modello
    integer :: i, j, step, new_state, ix, iy, k
    real(8) :: theta(L, L), theta_new(L, L)
    real(8) :: energy, energy_change, magn, rand_temp
    real(8) :: beta 
                          
                          
    call random_seed()               !!inizializzo il generatore e il reticolo
    read*, beta
    call random_angles(theta)
    print*, L*L, beta

    do step = 1, N
        do i = 1, L
            do j =1, L           
                !!do k = 1, L             linee per debugging
                !!print*, theta(k,:)
                !!end do
           
                call random_number(rand_temp)
                if (rand_temp < 0.5) then
                    new_state = int(mod(theta(i, j) + 1.0, real(q, kind=8)))  !! +-1
                else 
                     if (theta(i, j) > 0.9) then
                         new_state = int(theta(i, j) - 1.0)
                     else 
                         new_state = int(q - 1.0)
                     end if
                end if
            
                !!print*, ''
				!!print*, i, j, new_state
				!!print*, '' 
            
				energy_change = calculate_energy_change(theta, i, j, new_state)

				if (energy_change < 0.0) then   !!metro
					theta(i, j) = real(new_state, kind=8)
				else
					call random_number(rand_temp)
					if (rand_temp < exp(-energy_change * beta)) then
						theta(i, j) = real(new_state, kind=8) 
					end if
				end if
			end do
        end do
        energy = calculate_total_energy(theta)
        magn = calculate_magn(theta)
        print*, energy, magn
        !!print*, ''

    end do

contains

    subroutine random_angles(theta)      !!inizializza una matrice con valori casuali da 0 a q-1
        real(8), dimension(L, L) :: theta
        real(8) :: rand_temp
        integer :: i, j

        do i = 1, L
            do j = 1, L
                call random_number(rand_temp)
                theta(i, j) = int(rand_temp * q)
            end do
        end do
    end subroutine random_angles

    function calculate_total_energy(theta) result(total_energy)
        real(8), dimension(L, L) :: theta
        real(8) :: total_energy
        integer :: i, j
        integer :: neighbor_x, neighbor_y
        real(8) :: diff, cos_diff

        total_energy = 0.0
        do i = 1, L
            do j = 1, L
                neighbor_x = mod(i, L) + 1
                neighbor_y = mod(j, L) + 1

                diff = mod(real(theta(i, j) - theta(neighbor_x, j), kind=8), real(q, kind=8))
                cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
                total_energy = total_energy - cos_diff
                
                diff = mod(real(theta(i, j) - theta(i, neighbor_y), kind=8), real(q, kind=8))
                cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
                total_energy = total_energy - cos_diff
            end do
        end do
        total_energy = total_energy / (L * L)
    end function calculate_total_energy

    function calculate_energy_change(theta, ix, iy, new_state) result(energy_change)  !!funzione per il dE
        real(8), dimension(L, L) :: theta
        integer :: ix, iy, new_state
        real(8) :: energy_change
        real(8) :: diff, cos_diff
        integer :: neighbor_x, neighbor_y

        energy_change = 0.0

        neighbor_x = mod(ix, L) + 1
        neighbor_y = mod(iy, L) + 1

        diff = mod(real(new_state - theta(neighbor_x, iy), kind=8), real(q, kind=8))
        cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
        energy_change = energy_change - cos_diff + cos(2.0d0 * 3.141592653589793d0 * & 
        & mod(real(theta(ix, iy) - theta(neighbor_x, iy), kind=8), real(q, kind=8)) / real(q, kind=8))

        diff = mod(real(new_state - theta(ix, neighbor_y), kind=8), real(q, kind=8))
        cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
        energy_change = energy_change - cos_diff + cos(2.0d0 * 3.141592653589793d0 * & 
        & mod(real(theta(ix, iy) - theta(ix, neighbor_y), kind=8), real(q, kind=8)) / real(q, kind=8))
        
        if (ix > 1.9) then
            neighbor_x = mod(ix, L) - 1  
        else 
            neighbor_x = L 
        end if
        
        diff = mod(real(new_state - theta(neighbor_x, iy), kind=8), real(q, kind=8))
        cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
        energy_change = energy_change - cos_diff + cos(2.0d0 * 3.141592653589793d0 * & 
        & mod(real(theta(ix, iy) - theta(neighbor_x, iy), kind=8), real(q, kind=8)) / real(q, kind=8))
        
        if (iy > 1.9) then
            neighbor_y = mod(iy, L) - 1  
        else 
            neighbor_y = L 
        end if
        
        diff = mod(real(new_state - theta(ix, neighbor_y), kind=8), real(q, kind=8))
        cos_diff = cos(2.0d0 * 3.141592653589793d0 * diff / real(q, kind=8))
        energy_change = energy_change - cos_diff + cos(2.0d0 * 3.141592653589793d0 * & 
        & mod(real(theta(ix, iy) - theta(ix, neighbor_y), kind=8), real(q, kind=8)) / real(q, kind=8))
        
    end function calculate_energy_change

    function calculate_magn(theta) result(mag)
        real(8), dimension(L, L) :: theta
        real(8) :: mag, temp1, temp2
        integer :: i, j

        mag = 0.0
        temp1 = 0.0
        temp2 = 0.0
        do i = 1, L
            do j = 1, L
                temp1 = temp1 + & 
                cos(2.0d0 * 3.141592653589793d0 * real(theta(i, j), kind=8) / real(q, kind=8))
                temp2 = temp2 + &
                sin(2.0d0 * 3.141592653589793d0 * real(theta(i, j), kind=8) / real(q, kind=8))
            end do
        end do
        temp1 = temp1 / (L * L)
        temp2 = temp2 / (L * L)
        mag = sqrt(temp1*temp1 + temp2*temp2)
    end function calculate_magn

end program clock_model
