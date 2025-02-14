

module EnergyMin
        use Atom 
        use ChainMesh 
        use ChainMeshCell
        use Rand, only: makeRandom, algor_uniform_random, random
        use cube_partition, only: bounds_type, partition_unit_cube, dp 
        use omp_lib 
        implicit none 
        contains 

function AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex) result(Energy)
    type(chainMesh_t), intent(inout), target :: chainMesh 
    real(kind=8), intent(in) :: sigma
    type(Atom_t), pointer, intent(inout) :: atom
    integer, intent(in) :: threadID, atomIndex
    real(kind=8) :: Energy, spin
    type(Atom_t), pointer :: atoms(:)
    integer :: j, neighbor_idx
    if (.not. associated(atom)) then 
        error stop "atom is not a valid pointer, critical error: aborting"
    end if

    atoms => chainMesh%atoms 
    Energy = 0.0d0
    
    do j = 1, size(atom%NeighborList)
        ! Get the neighbor's index first
        neighbor_idx = atom%NeighborList(j)
        ! Then use it to access the neighbor's parameters
        !$omp atomic read 
        spin = atoms(neighbor_idx)%AtomParameters(1) 
!        Energy = Energy - dble(sigma) * dble(atom%AtomParameters(1)) * &
                !dble(atoms(neighbor_idx)%AtomParameters(1))
        Energy = Energy - dble(sigma) * dble(atom%AtomParameters(1)) * &
                dble(spin)
    end do 
!     print *, "Thread ", threadID, "Calculating energy of atomIndex ", atomIndex, " Energy = ", Energy

end function AtomEnergy 

                subroutine coordsFromThreadIndex(Ind,b,c, i,j,k) ! factored Ind into a*b*c. 
                        integer, intent(in) :: Ind,b,c 
                        integer, intent(out) :: i,j,k
                        integer :: remainder 
                        real(kind=dp), dimension(3) :: min_point, max_point
                       ! Ind = i (bc) + j(c) + k. This is just the formula for stride. This can be solved for i,j,k using the
                        ! following: 
                        i = Ind/(b*c)
                        remainder = modulo(Ind,b*c)
                        j = remainder / c 
                        k = modulo(remainder,c)
                       ! i,j,k range from (0,a-1), (0,b-1), and (0,c-1) respectively
                end subroutine coordsFromThreadIndex 

               ! Do the metropolis algorithm for the chain mesh 
                subroutine Metropolis(chainMesh, sigma, beta, nsteps)
                        implicit none
                        type(chainMesh_t), intent(inout), target :: chainMesh
                        real(kind=8), intent(in) :: sigma 
                        real(kind=8), intent(in) :: beta !1 / k_B T 
                        integer, intent(in) :: nsteps
                        integer :: threadID, numThreads, binsPerSide, a,b,c 
                        type(bounds_type), allocatable :: thread_bounds(:)
                        integer :: i, imin, imax, jmin, jmax, kmin, kmax 
                        type(bounds_type) :: thread_bound
                        logical :: critical_error, success
                        type(random) :: random
                        real(kind=dp) :: randomNumber, randi, randj, randk
                        real(kind=dp), dimension(3) :: min_point, max_point
                        integer :: N, randAtom, itemp, atomIndex, cell_index, num_atoms_in_cell, first_atom_index
                        integer :: tempint
                        type(Atom_t), pointer :: atom, tempAtomPointer 
                        real(kind=8) :: Energy, NewEnergy, Z, prob  
                        integer(kind=OMP_LOCK_KIND), allocatable :: lockArray(:), cell_lock_array(:)
                        real, allocatable :: atomParameterArray(:)
                        integer :: clock_count, clock_rate, attempt, j 
                        type(Atom_t), target :: tempAtom 
                        integer, dimension(27) :: cell_neighbor_array, sorted_cells 
                        call system_clock(clock_count,clock_rate)
                        allocate(lockArray(chainMesh%numAtoms))
                        allocate(cell_lock_array(chainMesh%numChainMeshCells))
                        allocate(atomParameterArray(size(chainMesh%atoms(1)%AtomParameters)))
                        do i=1,size(lockArray)
                            call omp_init_lock(lockArray(i))
                        end do
                        do i=1,size(cell_lock_array)
                           
                            call omp_init_lock(cell_lock_array(i))
                        end do 
                        ! Algorithm:
                        ! Each thread will get it's thread ID and figure out how many threads there are 
                        ! We will factor the number of threads into a,b,c such that a*b*c = num threads 
                        ! From this we will assign a range of chain mesh cells to each thread 
                        ! Each thread will do Metropolis on each of it's chain mesh cells 
                        critical_error = .false.
                        print *, "Entering parralell section"
                        !$omp parallel &
                        !$omp shared(chainMesh, critical_error, lockArray, beta, nsteps, sigma,cell_lock_array) &
                        !$omp private(threadID, numThreads, a, b, c, binsPerSide, thread_bounds, & 
                        !$omp         random, N, randAtom, itemp, i, imin, imax, jmin, jmax, &
                        !$omp         kmin, kmax, atomIndex, cell_index, Energy, NewEnergy, Z, prob, &
                        !$omp         randi, randj, randk, randomNumber, min_point, max_point, tempint, atom, tempAtom, &
                        !$omp         tempAtomPointer, atomParameterArray,cell_neighbor_array, sorted_cells, attempt,j,success) & 
                        !$omp default(private)
                       block 
                                !binsPerSide = chainMesh%numChainMeshCellsPerSide
                                N = chainMesh%numChainMeshCellsPerSide
                                threadID = omp_get_thread_num()
                                numThreads = omp_get_num_threads()
                                ! Compute which range of x,y,z for atoms within which the thread will randomly sample
                                thread_bounds = partition_unit_cube(numThreads)
                                ! if (threadID == 0) then !testing that it works 
                                !         print *, "num threads: ", numThreads 
                                !         do i=1,size(thread_bounds)
                                !                 print *, "i = ", i 
                                !                 print *, "lower bound: ", thread_bounds(i)%min_point
                                !                 print *, "upper bound: ", thread_bounds(i)%max_point 
                                                
                                !         end do
                                ! end if
                                thread_bound = thread_bounds(threadID + 1) ! +1 due to Fortran indexing
                                min_point = thread_bound%min_point
                                max_point = thread_bound%max_point 
                                imin = int(min_point(1)*N)
                                imax = int(max_point(1)*N)
                                !if (imax == N) print *, "imax == N"
                                !if (imin == 0) print *, "imin == 0"
                                jmin = int(min_point(2)*N)
                                jmax = int(max_point(2)*N)
                               !if (jmax == N) print *, "jmax == N"
                               !if (jmin == 0) print *, "jmin == 0"
                                kmin = int(min_point(3)*N)
                                kmax = int(max_point(3)*N)
                               !if (kmax == N) print *, "kmax == N" 
                               !if (kmin == 0) print *, "kmin == 0"
                                ! Check for sensible values
                                if (imin == imax) then 
                                        print *, "imin == imax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                if (jmin == jmax) then 
                                        print *, "jmin == jmax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                if (kmin == kmax) then 
                                        print *, "kmin == kmax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                ! Read only so no race conditions within this barrier 
                                if (critical_error) then 
                                        error stop "Critical Error detected, ending execution on all threads"
                                end if 

                                ! Avoid overlaps between threads 
                                if (imin /= 0 ) then 
                                        imin = imin + 1
                                end if 
                                if (jmin /= 0 ) then 
                                        jmin = jmin + 1
                                end if 
                                if (kmin /= 0) then 
                                        kmin = kmin + 1 
                                end if 
                                ! Now need to select a mesh cell i,j,k coordinate within the threads assigned bounds, then select an
                                ! atom within it 
                                ! To do this we will need the lobotomised CASTEP random number generator as the standard library one
                                ! is probably not thread safe (better safe than sorry).
                                ! This is a bitshift random number generator which returns numbers in the range [0,1], we will also need to warm up
                                ! the generator otherwise threads will make very similar choices for the first few steps.
                                ! We will use the threadID as the seed 
                                !print *, "Thread ", threadID, " is about to make a random number"
                                random=makeRandom(threadID + clock_count)
                                !print *, "Thread ", threadID, " has made a  random number"

                                !Warm up the generator 
                                do i=1,30
                                    randomNumber = algor_uniform_random(random)
                                end do
                                !print *, "Thread ", threadID, " has warmed up the random number generator"

                                do itemp=1,nsteps
                                        !$omp barrier
                                        !Randomly select i,j, and k 
                                        !randi = imin + floor(algor_uniform_random(random)*real((imax - imin)))
                                        !randj = jmin + floor(algor_uniform_random(random)*real((jmax - jmin)))
                                        !randk = kmin + floor(algor_uniform_random(random)*real((kmax - kmin)))
                                        !if (int(randi) == N - 1) print *, "randi == N - 1"
                                        !if (int(randi) == 0) print *, "randi == 0"
                                        !Now randomly select atom in this chain mesh cell
                                        !claude: 
                                       ! cell_index = nint(randi*N**2) + nint(randj*N) + nint(randk) + 1
                                        !cell_index = int(randi)*N**2 + int(randj)*N + int(randk) + 1 TODO: This is a temporary
                                        !exploration using an alternative sampling algorithm 
                                        randi = floor(algor_uniform_random(random) * N) 
                                        randj = floor(algor_uniform_random(random) * N)
                                        randk = floor(algor_uniform_random(random) * N)
                                        cell_index = int(randi)*N*N + int(randj)*N + int(randk) + 1
                                        !print *, "Thread ", threadID, " Is trying to acquire the lock of cell number ", cell_index

                                        call omp_set_lock(cell_lock_array(cell_index))
                                        !print *, "Thread ", threadID, "Has acquired the lock of cell: ", cell_index 
                                        !print *, "cell_index = ", cell_index, " out of ", N**3 
                                        !print *, "Thread ", threadID, "is setting it's locks: ", cell_neighbor_array
        
                                        num_atoms_in_cell = chainMesh%chainMeshCells(cell_index)%NumAtomsPerUnitCell
                                        randAtom = floor(algor_uniform_random(random) * real(num_atoms_in_cell))
                                        !randAtom = floor(algor_uniform_random(random)*real(chainMesh%chainMeshCells(nint(randi*N**2) + &
                                                !nint(randj*N) +  nint(randk) + 1 )%NumAtomsPerUnitCell))                                               ! randomly 
                                        !claude: 

                                        !cell_index = int(randi*N**2) + int(randj*N) + int(randk) + 1
!                                        cell_index = int(randi)*N**2 + int(randj)*N + int(randk) + 1 
                                        first_atom_index = chainMesh%chainMeshCells(cell_index)%firstAtomInMeshCell
                                        atom => chainMesh%atoms(first_atom_index)
                                        !print *, "Thread ", threadID, "has initialised atom"
                                        !atom = chainMesh%atoms(chainMesh%chainMeshCells(int(randi*N**2) + int(randj*N) + int(randk) + 1 &
                                        !)%firstAtomInMeshCell)
                                        !atomIndex = chainMesh%chainMeshCells(int(randi*N**2) + int(randj*N) + int(randk) + 1 &
                                        !)%firstAtomInMeshCell  
                                        atomIndex = first_atom_index
                                        if (atom%nextAtom == -1) then !DEBUG  
                                                print *, "Atom has next atom -1, critical error: aborting"
                                                error stop "Critical Error"
                                        end if
                                        !TODO: Do metropolis step then put in loop
                                        ! Check that all of the first elements of the chainMeshCell's first atoms are initialised:
                                        do i=1,size(chainMesh%chainMeshCells) !DEBUG 
                                            if (chainMesh%chainMeshCells(1)%firstAtomInMeshCell == -1) then 
                                                    print *, "chainMeshCell ", i, " is not properly initialised"
                                                    error stop "Critical error, aborting"
                                            end if 
                                        end do 
                                        do i = 1,randAtom ! Don't need + 1 as if randAtom=0 the od loop won't run, if randAtom=1 it will run
                                                          !once
                                            !print *, "i = ", i, "atom%nextatom = ", atom%nextAtom, " randAtom = ", randAtom
                                            atomIndex = atom%NextAtom !This should never be -1, if it is we would get an out of
                                                                      ! Out of bounds memory accsess.
                                            
                                            atom => chainMesh%atoms(atom%nextAtom)
                                            
                                        end do
                                        ! print *, "Thread ", threadID , " has selected atom with index", atomIndex
                                        !print*, "Thread ", threadID, "Has randomly selected atom ", randAtom, " In it's unit cell: ", &
                                        !int(randi), int(randj), int(randk)
                                        !TODO: Finish metropolis on these randomly selected atoms, maybe make atom a pointer such that we
                                        !can change it's value within the array
                                        if (.not. associated(atom)) then 
                                                print *, "Thread ", threadID, " has a not valid atom pointer with index ", atomIndex
                                                error stop "Threadm atom not a valid pointer"
                                        end if 
                                        !print *, "Thread ", threadID, " is passing atomIndex ", atomIndex
                                        Energy = AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex)
                                        ! print *, "Thread ", threadID, " is about to index and write to atomParameters of atom "
                                        ! print *, atomIndex

                                        ! atom%atomParameters(1) = - atom%atomParameters(1) ! TODO:Atom is a pointer so should make a copy of atomparameters(1)
                                        atomParameterArray = - atom%atomParameters !This is not being used but should be
                                        
                                        ! print *, "Thread ", threadID, " has indexed atomParameters"
                                        ! In order to avoid race conditions while limiting the use of locks, a copy of the atom
                                        ! being considered must be made. Annoyingly this 
                                        !tempAtom = makeAtom(atom%x,atom%y,atom%z,atomParameterArray,atom%numAtomParameters,atom%nextAtom)
                                        tempAtom = atomDeepCopy(atom)  
                                        tempAtomPointer => tempAtom
                                        tempAtom%atomParameters = -atom%atomParameters
                                        NewEnergy = AtomEnergy(chainMesh, sigma, tempAtomPointer,threadID,atomIndex)
                                        ! tempAtom should be automatically deallocated at the end of the scope
                                        !print *, "Thread ", threadID, " has obtained a probability"
                                        
                                        if (NewEnergy < Energy) then 
                                                Z = 1.0 
                                        else 
                                                Z = exp(- beta * ((NewEnergy - Energy)/sigma))
                                        end if 
                                        prob = algor_uniform_random(random)
                                        ! print *, "Thread ", threadID, " has obtained a probability"
                                        ! tempint = atomIndex
                                        ! print *, "Thread ", threadID, " has selected atomIndex to be ",tempint
                                        ! Barrier here to prevent race conditions from threads calculating energies (Reading) and
                                        ! updating the mesh (writing)
                                        !print *, "Debug: thread ", threadID, " Z = ", Z, " prob = ", prob, " atomIndex = ",atomIndex, "New Energy: ", NewEnergy, " Old Energy ", Energy
                                        if (Z >= prob) then 
                                                !call omp_set_lock(lockArray(atomIndex))
                                                do i = 1,size(atom%atomParameters)
                                                !     print *, "Thread ", threadID, " is about to write to atom"
                                                    chainMesh%atoms(atomIndex)%atomParameters(i) = - atom%atomParameters(i) 
                                                    !print *, "Thread ", threadID, " Writing to mesh with value: ", chainMesh%atoms(atomIndex)%atomParameters(i)
                                                end do 
                                                !call omp_unset_lock(lockArray(atomIndex))
                                                !print *, "Thread ", threadID, "Has unset the lock"
                                                ! print *, "Thread ", threadID, " has updated the mesh for atom index", atomIndex

                                        end if 
                                        call omp_unset_lock(cell_lock_array(cell_index))
                                        !print *, "Thread ", threadID, " has released the lock of cell number ", cell_index
                                        !$omp barrier
 
                                        !print *, "######################################################"
                                end do 

                      end block 

                     !$omp end parallel 
                     do i = 1,size(lockArray)
                        call omp_destroy_lock(lockArray(i))
                     end do 

                     do i=1,size(cell_lock_array)
                        call omp_destroy_lock(cell_lock_array(i))
                     end do 
            end subroutine Metropolis 


        subroutine Metropolis2(chainMesh, sigma, beta, nsteps)
                ! Each thread will randomly select an atom nsteps times and determine whether to flip the spin 
                type(ChainMesh_t), intent(inout), target :: chainMesh 
                real(kind=8), intent(in) :: sigma, beta 
                integer, intent(in) :: nsteps 
                integer(kind=OMP_LOCK_KIND), allocatable :: atom_lock_array(:)
                integer :: i, i2, time, threadID, atomIndex
                type(random) :: rand 
                real(kind=8) :: rand_num, oldEnergy, newEnergy, prob, Z
                type(Atom_t), pointer :: atom 
                allocate(atom_lock_array(size(chainMesh%atoms)))
                print *, "Atom lock array allocated with size:", size(atom_lock_array)
                do i = 1, size(atom_lock_array)
                        call omp_init_lock(atom_lock_array(i))
                end do
                print *, "Atom locks initialised"
                 
                !$omp parallel &
                !$omp shared(atom_lock_array,chainMesh,sigma, beta, nsteps) &
                !$omp private(rand,i,rand_num, time,threadID,atomIndex, oldEnergy, newEnergy,prob, Z, i2, atom) 
                block 
                        threadID = omp_get_thread_num()
                        call system_clock(time)
                        rand = makeRandom(time*threadID + modulo(threadID,time))
                        do i = 1,100
                                rand_num = algor_uniform_random(rand) ! Warm up the generator
                        end do 
                        do i = 1,nsteps
                                rand_num = algor_uniform_random(rand)
                                atomIndex = int(rand_num*dble(size(chainMesh%atoms)) - 1) + 1
                                if (atomIndex > size(atom_lock_array)) print *, "Atom Index = " ,atomIndex, "size(atom_lock_array)"&
                                        , size(atom_lock_array)
                                if (atomIndex < 1) print *, "atomIndex = ", atomIndex
                                atom => chainMesh%atoms(atomIndex) !This is not really needed but to satisfy the function here it is 
                                call omp_set_lock(atom_lock_array(atomIndex))

                                        oldEnergy = AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex)  
                                        atom%AtomParameters(1) = - atom%AtomParameters(1)
                                        newEnergy = AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex)
                                        atom%AtomParameters(1) = - atom%AtomParameters(1) 
                                        
                                        if (newEnergy < oldEnergy) then 
                                                Z = 1.0 
                                        else 
                                                Z = exp(- beta * ((newEnergy - oldEnergy)/sigma))
                                        end if 
                                        prob = algor_uniform_random(rand)
                                        if (Z >= prob) then 
                                                do i2 = 1,size(atom%atomParameters)
                                                    chainMesh%atoms(atomIndex)%atomParameters(i2) = - atom%atomParameters(i2) 
                                                end do 
                                        end if 
                                        call omp_unset_lock(atom_lock_array(atomIndex))


                        end do 

                end block 
                !$omp end parallel
                do i = 1,size(atom_lock_array)
                        call omp_destroy_lock(atom_lock_array(i))
                end do 
                print *, "Destroyed locks"
        end subroutine Metropolis2
subroutine WriteMagnetization(chainMesh, filename)
    implicit none
    type(ChainMesh_t), intent(in) :: chainMesh
    character(len=*), intent(in) :: filename
    real(kind = 8) :: total_mag, mag_per_atom
    integer :: i, unit_num
    
    ! Calculate average magnetization
    total_mag = 0.0
    do i = 1, size(chainMesh%atoms)
        total_mag = total_mag + chainMesh%atoms(i)%AtomParameters(1)
    end do
    mag_per_atom = total_mag / real(size(chainMesh%atoms))
    
    ! Open file in append mode
    open(newunit=unit_num, file=trim(filename), position='append', action='write')
    
    ! Write magnetization per atom to new line
    write(unit_num, *) mag_per_atom
    
    close(unit_num)
    
end subroutine WriteMagnetization


subroutine write_array_to_file(array, filename, iostat)
    implicit none
    
    ! Arguments
    real(kind=8), intent(in) :: array(:)          ! Input array (assumed-shape array)
    character(len=*), intent(in) :: filename  ! Name of output file
    integer, intent(out) :: iostat        ! Status indicator
    
    ! Local variables
    integer :: unit_num, i
    
    ! Open the file
    open(newunit=unit_num, file=trim(filename), status='replace', action='write', iostat=iostat)
    if (iostat /= 0) return
    
    ! Write array elements one per line
    do i = 1, size(array)
        write(unit_num, *, iostat=iostat) array(i)
        if (iostat /= 0) then
            close(unit_num)
            return
        end if
    end do
    
    ! Close the file
    close(unit_num)
end subroutine write_array_to_file
end module EnergyMin 
