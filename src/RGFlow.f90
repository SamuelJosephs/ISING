module RGFlow
        use omp_lib 
        use ChainMesh, only : distance, chainMesh_t, compute_unique_lattice_vectors
        use vecNd 
        implicit none 


        public :: compute_correlation_function 

        contains 

                ! This function turns a distance into an array index, or "bin". 
                ! We do this by recognising that each lattice distance is an integer multiple of a lattice vector G, the set of
                ! which are passed into the function. For n lattice vectors, there will be an integer multiple of n lattice
                ! distances below the maximum distance. As such the bin position may be calculated with modular arithmatic.
                ! We are using a block based indexing scheme, where we realise that all lattice distances that may occur are grouped
                ! into blocks. We first find the block offset with max_distance / |G_d|, then find the offset into the block using
                ! modulo(max_distance,abs(G_d)).
                function distance_to_bin(chainMesh, max_distance_in,distance,lattice_vector_array) result(bin)
                        type(ChainMesh_t), intent(in) :: chainMesh 
                        real(kind=8), intent(in) :: max_distance_in, distance 
                        type(vecNd_t), intent(in) :: lattice_vector_array(:)
                        real(kind=8) :: ratio, little_n 
                        integer :: numDistances, i, bin, G_d, N 
                        
                        if (size(lattice_vector_array) == 0) error stop "lattice_vector_array is empty"

                        if (abs(distance) < 1e-8) then 
                                bin = 1
                                return 
                        end if 
                        numDistances = 0 ! We only care about non zero distances for this calculation
                        G_d = -1 
                        do i = 1,size(lattice_vector_array)
                                if (abs(lattice_vector_array(i)) == 0.0) cycle

                                numDistances = numDistances + floor(max_distance_in / abs(lattice_vector_array(i)))
                                !print *, "distance = ", distance, "abs(G) = ", abs(lattice_vector_array(i))
                                !print *, "distance % abs(G) = ", modulo(distance,abs(lattice_vector_array(i)))
                                ratio = distance / abs(lattice_vector_array(i))
                                if (abs(ratio - nint(ratio)) < 1e-4) G_d = i 
                                !if (abs(modulo(distance,abs(lattice_vector_array(i)))) < 1e-3) then 
                                !        G_d = i 
                                !end if 
                        end do 

                        if (G_d == -1) then
                                print *, "########################################################"
                                print *, "Distance = ", distance
                                print *, "available lattice vectors = "
                                do i = 1,size(lattice_vector_array)
                                        print *, lattice_vector_array(i)%coords 
                                end do 
                                print *, "available lattice vector magnitudes = ", abs(lattice_vector_array)
                                print *, "available modulo = ", modulo(distance, abs(lattice_vector_array))
                                error stop "Could not find matching lattice vector for distance"
                        end if 
                        !bin = int(floor(max_distance_in / distance)) + floor(modulo(max_distance_in,abs(lattice_vector_array(G_d)))) + 2
                        N = size(lattice_vector_array) - 1 
                        little_n = distance / abs(lattice_vector_array(G_d))

                        bin = n*N + modulo(numDistances,nint(abs(lattice_vector_array(G_d)))) + 2 
                        print *, "bin = ", n*N, "Step into bin = ", modulo(numDistances,nint(abs(lattice_vector_array(G_d))))  
                        !print *, "bin = ",int(floor(max_distance_in / distance)), "step into bin = ", floor(modulo(max_distance_in,abs(lattice_vector_array(G_d))))  
                        ! 1 based indexing and +1 to avoid the bin meant for 0 distance gives a displacement of +2 

                        

                end function distance_to_bin 


                ! This method calculates the two point correlation function for the spins in a Metropolic alsogirhtm Monte Carlo
                ! simulation. 
                ! Within this formalism, < sigma_x sigma_y > = Tr{p sigma_x sigma_y}. For a thermal state the density matrix is
                ! e^{-beta H}. The metropolis algorithm effectively importance samples this distribution, reducing the problem to
                ! finding the mean. Here we compute the correlation function as a function of the distance between lattice sites,
                ! this is allowed due to our periodic boundary conditions.
                subroutine compute_correlation_function(chainMesh, max_distance_in, correlations) 
                        implicit none 
                        type(ChainMesh_t), intent(in) :: chainMesh 
                        real(kind=8), intent(in) :: max_distance_in 
                        real(kind=8), allocatable,  intent(out) :: correlations(:) 
                        integer, allocatable :: counts(:) 
                        integer :: i, j,  bin
                        integer(kind=8) :: numDistances 
                        real(kind=8) :: r
                        type(vecNd_t), allocatable :: lattice_vector_array(:)
                        real(kind=8) :: correlation 
                        integer(kind=8) :: temp

                        lattice_vector_array = compute_unique_lattice_vectors(chainMesh)
                        
                        do i = 1,size(lattice_vector_array)
                                print *, lattice_vector_array(i)%coords
                        end do 

                        if (size(lattice_vector_array) == 0) error stop "Error: lattice_vector_array is empty."
                        ! Compute the number of discrete dinstances within max_distance 
                        numDistances = 0_8 
                        do i = 1,size(lattice_vector_array)
                                if (abs(lattice_vector_array(i)) == 0) cycle 
                                numDistances = numDistances + floor(max_distance_in / abs(lattice_vector_array(i)))
                        end do
                        print *, "NumDistances = ", numDistances
                        allocate(correlations(numDistances + 2)) ! +1 for the zero bin 
                        allocate(counts(numDistances + 2))
                        correlations = 0.0 
                        counts = 0 
                        do i = 1,size(chainMesh%atoms)
                                do j = i, size(chainMesh%atoms)
                                        r = distance(chainMesh,i,j)
                                        if (r > max_distance_in) cycle 
                                        bin = distance_to_bin(chainMesh,max_distance_in,r,lattice_vector_array)
                                        print *, "numDistances = ", numDistances, "bin = ", bin
                                        
                                        correlations(bin) = correlations(bin) + & 
                                                chainMesh%atoms(i)%AtomParameters(1)*chainMesh%atoms(j)%AtomParameters(1) 
                                        counts(bin) = counts(bin) + 1 

                                end do 
                        end do 

                        do i = 1,size(correlations)
                                correlations(i) = correlations(i) / counts(i)
                        end do 
                        
                print *, "correlation = ", correlations 
                print *, "Size of correlation array = ", size(correlations)
                end subroutine compute_correlation_function


end module RGFlow
