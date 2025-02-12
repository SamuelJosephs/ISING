module RGFlow
        use omp_lib 
        use ChainMesh, only : distance, chainMesh_t  
        implicit none 


        public :: compute_correlation_function 

        contains 

                ! This method calculates the two point correlation function for the spins in a Metropolic alsogirhtm Monte Carlo
                ! simulation. 
                ! Within this formalism, < sigma_x sigma_y > = Tr{p sigma_x sigma_y}. For a thermal state the density matrix is
                ! e^{-beta H}. The metropolis algorithm effectively importance samples this distribution, reducing the problem to
                ! finding the mean. Here we compute the correlation function as a function of the distance between lattice sites,
                ! this is allowed due to our periodic boundary conditions.
                subroutine compute_correlation_function(chainMesh, max_distance_in, correlation) 
                        type(ChainMesh_t), intent(in) :: chainMesh 
                        real, intent(in) :: max_distance_in 
                        real(kind=8), allocatable,  intent(out) :: correlation(:) 
                        integer, allocatable :: counts(:) 
                        integer :: i, j, numAtoms, numbins, bin 
                        real :: r, max_distance 
                        
                        if (.not. allocated(correlation)) then 
                                numbins = ceiling(max_distance_in / chainMesh%latticeParameter)
                                numAtoms = size(chainMesh%atoms)         
                                allocate(correlation(numbins))
                                max_distance = max_distance_in 

                        else 
                                numbins = size(correlation)
                                max_distance = numbins * chainMesh%latticeParameter
                                if (abs(max_distance - max_distance_in) > 1e-6) then 
                                        print *, "WARNING: subroutine compute_correlation_function recieved an already allocated & 
                                                correlation vector, overiding the provided max_distance."
                                end if 
                                
                        end if 
                        allocate(counts(numbins))
                        correlation = 0.0
                        counts = 0

                        !$omp parallel do 
                        do i = 1, size(chainMesh%atoms)
                                do j = i, size(chainMesh%atoms)
                                        ! Compute the distance between atoms i and j 
                                        r = distance(chainMesh, i, j)
                                        if (r < max_distance) then 
                                                bin = int(floor(r / chainMesh%latticeParameter)) + 1 ! Fortran has 1 based indexing 
                                                correlation(bin) = correlation(bin) + & 
                                                        chainMesh%atoms(i)%AtomParameters(1)*chainMesh%atoms(j)%AtomParameters(1)
                                                counts(bin) = counts(bin) + 1 
                                        end if 

                                end do 
                        end do 
                        !$omp end parallel do 
                        where (counts > 0) 
                                correlation = correlation / counts 
                        end where 
                print *, "correlation = ", correlation
                print *, "Size of correlation array = ", size(correlation)
                end subroutine compute_correlation_function


end module RGFlow
