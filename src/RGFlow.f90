module RGFlow
        use omp_lib 
        use ChainMesh, only : distance, chainMesh_t
        use vecNd 
        implicit none 


        public :: compute_correlation_function, write_correlations_to_file 

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
                
                function find_index(array, num) result(ind)
                        real(kind=8), intent(in) :: array(:)
                        real(kind=8), intent(in) ::  num 
                        integer :: ind, i 
                        ind = -1 
                        do i = 1,size(array)
                                if (abs(array(i) -  num)< 1e-4) then 
                                        ind = i 
                                        return 
                                end if 
                        end do 
                 end function find_index
                        

                ! This method calculates the two point correlation function for the spins in a Metropolic alsogirhtm Monte Carlo
                ! simulation. 
                ! Within this formalism, < sigma_x sigma_y > = Tr{p sigma_x sigma_y}. For a thermal state the density matrix is
                ! e^{-beta H}. The metropolis algorithm effectively importance samples this distribution, reducing the problem to
                ! finding the mean. Here we compute the correlation function as a function of the distance between lattice sites,
                ! this is allowed due to our periodic boundary conditions.
                subroutine compute_correlation_function(chainMesh, max_distance_in, correlations,positions) 
                        implicit none 
                        type(ChainMesh_t), intent(in) :: chainMesh 
                        real(kind=8), intent(in) :: max_distance_in 
                        real(kind=8), allocatable,  intent(out) :: correlations(:) 
                        integer, allocatable :: counts(:) 
                        real(kind=8), allocatable, intent(out) :: positions(:) 
                        integer :: i, j, ind  
                        real(kind=8) :: r 
                        real(kind=8), parameter :: zero = 0.0
                        allocate(correlations(0))
                        allocate(positions(0))
                        allocate(counts(0))
                        correlations = 0.0 
                        positions = 0.0 
                        counts = 0 
                        do i = 1,size(chainMesh%atoms)
                                ! Calculate correlation for atom i as a function of distance 
                                do j = 1, size(chainMesh%atoms)
                                r = distance(chainMesh,i,j) 
                                ! If r is in the positions, update the correlations appropriately, if not add it 
                                ind  = find_index(positions,r)
                                if (ind == -1) then 
                                        ! add position to the positions array 
                                        positions = [positions, r]
                                        correlations = [correlations, zero]
                                        counts = [counts, 0]
                                        ind = size(positions) !The index is not the last element 
                                end if 
                                correlations(ind)  = correlations(ind) + chainMesh%atoms(i)%AtomParameters(1)*chainMesh%atoms(j)%AtomParameters(1)
                                counts(ind) = counts(ind) + 1 
                                end do 
                        end do 
                        do i = 1, size(correlations)
                                correlations(i) = correlations(i) / counts(i)
                        end do 
                end subroutine compute_correlation_function

                subroutine write_correlations_to_file(correlations, positions, fileName)
                    implicit none
                    real(kind=8), intent(in) :: correlations(:)
                    real(kind=8), intent(in) :: positions(:)
                    character(len=*), intent(in) :: fileName
                    integer :: i, fileUnit

                    open(newunit=fileUnit,file=fileName,status="replace",action="write")

                    do i = 1, size(positions)
                        write(fileUnit, '(F12.6, 2X, F12.6)') positions(i), correlations(i)
                    end do

                    close(fileUnit)



                end subroutine write_correlations_to_file



end module RGFlow
