
module EnergyMin
        use Atom 
        use ChainMesh 
        use ChainMeshCell
        use Rand, only: makeRandom, algor_uniform_random, random, Normal
        use cube_partition, only: bounds_type, partition_unit_cube, dp 
        use constants, only: gyromagnetic_ratio, Kb, Bohr_magneton
        use reciprocal_space_processes
        use omp_lib 
        implicit none 
        contains 


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

        subroutine UniformRandomSpin(vec3d, rand)
                type(vecNd_t), intent(inout) :: vec3d 
                type(random), intent(inout) :: rand
                real(kind=8), parameter :: pi = 3.14159265358979323846
                integer :: i 
                real(kind=8) :: theta,phi, u 

                u = algor_uniform_random(rand) 
                u = 2.0_8*u - 1 
                theta = acos(u) 
                phi = algor_uniform_random(rand) 
                phi = 2.0_8*pi*phi 

                vec3d%coords(1) = cos(theta)*sin(phi)
                vec3d%coords(2) = sin(theta)*sin(phi)
                vec3d%coords(3) = cos(theta)

        end subroutine UniformRandomSpin 

        subroutine UniformRandomInSphere(vec3d, rand, R)
                type(vecNd_t), intent(inout)   :: vec3d
                type(random), intent(inout)  :: rand
                real(kind=8), intent(in)     :: R
                real(kind=8) :: theta, phi, u, r_lower
                real(kind=8), parameter :: pi = 3.14159265358979323846
                ! First pick a uniformly random direction on the unit sphere
                u = 2.0_8*algor_uniform_random(rand) - 1.0_8
                theta = acos(u)
                phi = 2.0_8*pi*algor_uniform_random(rand)

                ! Now pick radius with the r^2 weight built in
                u = algor_uniform_random(rand)
                r_lower = R * u**(1.0_8/3.0_8)

                ! Combine into a point in the ball
                vec3d%coords(1) = r_lower * sin(theta)*cos(phi)
                vec3d%coords(2) = r_lower * sin(theta)*sin(phi)
                vec3d%coords(3) = r_lower * cos(theta)
        end subroutine UniformRandomInSphere



        subroutine GaussianStep(vec3d, rand, beta)
                type(vecNd_t), intent(inout) :: vec3d 
                type(random), intent(inout) :: rand 
                real(kind=8), intent(in) :: beta
                real(kind=8) :: T 
                integer :: i
                real(kind=8) :: stddev
                stddev = (2.0_8/25.0_8)*(1.0_8 / (beta*Bohr_magneton))**(0.2)
                do i = 1,size(vec3d)
                        vec3d%coords(i) = vec3d%coords(i) + Normal(rand,0.0_8, stddev)
                end do 
                vec3d = vec3d / abs(vec3d)
        end subroutine GaussianStep 

        subroutine calculateHeisenbergEnergy(chainMesh,atomIndex, J, J_prime, Dz, Dz_prime ,B , & 
                        lockArray, S_proposed, oldEnergy, newEnergy, demagnetisation_array, demag)
                use iso_fortran_env, only: stderr => error_unit
                implicit none
                type(ChainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: atomIndex 
                real(kind=8), intent(in) :: J, J_prime, Dz,Dz_prime, B 
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
                type(vecNd_t), intent(in) :: S_proposed 
                real(kind=8), intent(inout) ::  oldEnergy, newEnergy 
                real(kind=8), dimension(:,:), intent(in) :: demagnetisation_array
                logical, optional, intent(in) :: demag 

                real(kind = 8) :: Energy, x,y,z, tempEnergy, Hx, Hy, Hz, MdotH, MdotH_proposed
                integer :: i, atomIndexTemp, dim_i, nn , nn_index  
                type(vecNd_t) :: S, S_prime, atomPos1, atomPos2, tempVec,r, D, D_prime
                logical :: calculate_demag 

                calculate_demag = .True.
                if (present(demag)) calculate_demag = demag
                !call OMP_SET_LOCK(lockArray(atomIndex))
                oldEnergy = 0.0_8
                newEnergy = 0.0_8
                S = makeVecNdCheck(S, dble(chainMesh%atomSpins(atomIndex,:)))
                x = chainMesh%atoms(atomIndex)%x
                y = chainMesh%atoms(atomIndex)%y 
                z = chainMesh%atoms(atomIndex)%z 
                atomPos1 = makeVecNd([x,y,z]) 
                tempVec = Dz*(chainMesh%c_vec / abs(chainMesh%c_vec))
                Hx = demagnetisation_array(atomIndex,1)
                Hy = demagnetisation_array(atomIndex,2)
                Hz = demagnetisation_array(atomIndex,3)
                MdotH = S%coords(1)* Hx + S%coords(2)* Hy + S%coords(3)* Hz 
                MdotH_proposed = S_proposed%coords(1)* Hx + S_proposed%coords(2)* Hy + S_proposed%coords(3)* Hz


                do i = 1,size(chainMesh%atomShells(atomIndex,1)%NNList)
                     
                        atomIndexTemp = chainMesh%atomShells(atomIndex,1)%NNList(i)

                        if (atomIndexTemp == atomIndex) error stop "Encountered self interaction"
                        call OMP_SET_LOCK(lockArray(atomIndexTemp))
                                S_prime = makeVecNdCheck(S_prime,dble(chainMesh%atomSpins(atomIndexTemp,:)))
                        call OMP_UNSET_LOCK(lockArray(atomIndexTemp))
                        x = chainMesh%atoms(atomIndexTemp)%x
                        y = chainMesh%atoms(atomIndexTemp)%y 
                        z = chainMesh%atoms(atomIndexTemp)%z
                        atomPos2 = makeVecNdCheck(atomPos2, [x,y,z])
                        !r = atomPos1 - atomPos2
                        call distance_points_vec(chainMesh,atomPos1,atomPos2, r) 
                        r = (-1.0_8) * r / abs(r)
                        D = tempVec .x. r
                        !D = Dz*r
                        oldEnergy = oldEnergy + ((J* S*S_prime) + (D*(S .x. S_prime)))
                        
                        newEnergy = newEnergy + ((J* S_proposed*S_prime) + (D*(S_proposed .x. S_prime)))
                end do 

                ! Now calculate contribution from next nearest neighbours J coupling 
                if (size(chainMesh%atomShells,2) > 1) then 

                        
                        do i = 1,size(chainMesh%atomShells(atomIndex,2)%NNList)
                                atomIndexTemp = chainMesh%atomShells(atomIndex,2)%NNList(i)

                                call OMP_SET_LOCK(lockArray(atomIndexTemp))
                                        S_prime = chainMesh%atomSpins(atomIndexTemp,:)
                                call OMP_UNSET_LOCK(lockArray(atomIndexTemp))
                                oldEnergy = oldEnergy + J_prime*(S*S_prime)
                                newEnergy = newEnergy + J_prime*(S_proposed*S_prime)

                        end do 
                else 
                       write(stderr,*) "Warning: Not calculating next to nearest neighbour interactions becuase numShells < 2"
                end if   
                oldEnergy = oldEnergy - g*Bohr_magneton*B*s%coords(3) 
                newEnergy = newEnergy - g*Bohr_magneton*B*S_proposed%coords(3) 
                if (calculate_demag) then 
                        oldEnergy = oldEnergy - mu_0*MdotH ! Half not needed due to it being a local change 

                        newEnergy = newEnergy - mu_0*MdotH_proposed
                end if 

                ! Don't need to add magnetic field contributions again as they have already been added 

        end subroutine calculateHeisenbergEnergy


        subroutine calculateTotalHeisenbergEnergyHelper(chainMesh,atomIndex, J, J_prime, Dz, Dz_prime ,B , & 
                        lockArray, energy, demag) ! Helper with correct factors of 1/2 to avoid double counting
                implicit none
                type(ChainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: atomIndex 
                real(kind=8), intent(in) :: J, J_prime, Dz,Dz_prime, B 
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
                real(kind=8), intent(out) ::  energy 
                logical, optional, intent(in) :: demag 

                real(kind = 8) :: x,y,z, tempEnergy, Hx, Hy, Hz, MdotH, MdotH_proposed
                integer :: i, atomIndexTemp, dim_i, nn , nn_index  
                type(vecNd_t) :: S, S_prime, atomPos1, atomPos2, tempVec,r, D, D_prime
                logical :: calculate_demag 

                calculate_demag = .True.
                if (present(demag)) calculate_demag = demag
                !call OMP_SET_LOCK(lockArray(atomIndex))
                energy = 0.0_8
                S = makeVecNdCheck(S, dble(chainMesh%atomSpins(atomIndex,:)))
                x = chainMesh%atoms(atomIndex)%x
                y = chainMesh%atoms(atomIndex)%y 
                z = chainMesh%atoms(atomIndex)%z 
                atomPos1 = makeVecNd([x,y,z]) 
                tempVec = makeVecNdCheck(tempVec, [0.0_8, 0.0_8, Dz])
                Hx = chainMesh%demagnetisation_array(atomIndex,1)
                Hy = chainMesh%demagnetisation_array(atomIndex,2)
                Hz = chainMesh%demagnetisation_array(atomIndex,3)
                MdotH = S%coords(1)* Hx + S%coords(2)* Hy + S%coords(3)* Hz 

                do i = 1,size(chainMesh%atoms(atomIndex)%NeighborList)
                        atomIndexTemp = chainMesh%atoms(atomIndex)%NeighborList(i)
                        if (atomIndexTemp == atomIndex) error stop "Encountered self interaction"
                        call OMP_SET_LOCK(lockArray(atomIndexTemp))
                                S_prime = makeVecNdCheck(S_prime,dble(chainMesh%atomSpins(atomIndexTemp,:)))
                        call OMP_UNSET_LOCK(lockArray(atomIndexTemp))
                        x = chainMesh%atoms(atomIndexTemp)%x
                        y = chainMesh%atoms(atomIndexTemp)%y 
                        z = chainMesh%atoms(atomIndexTemp)%z
                        atomPos2 = makeVecNdCheck(atomPos2, [x,y,z])
                        !r = atomPos1 - atomPos2
                        call distance_points_vec(chainMesh,atomPos1,atomPos2, r) 
                        r = (-1.0_8) * r / abs(r)
                        D = tempVec .x. r
                        !D = Dz*r
                        energy = energy + 0.5_08*((J* S*S_prime) + (D*(S .x. S_prime)))
                        
                end do 
                energy = energy - (g*Bohr_magneton*B*s%coords(3)) 
                if (calculate_demag) then 
                        energy = energy - (0.5_8*mu_0*MdotH) ! Half needed as this is used to calculate total energy 
                                                             ! So there is double counting
                end if 


        end subroutine calculateTotalHeisenbergEnergyHelper




        subroutine Metropolis_mcs(chainMesh, beta,numMCSSweeps, J, J_prime,  Dz, Dz_prime, B, r, lockArray, &
                                demag)
                ! Each thread will randomly select an atom nsteps times and determine whether to flip the spin 
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), intent(in) :: beta, J, J_prime, Dz, Dz_prime, B, r
                integer, intent(in) :: numMCSSweeps
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
                logical, optional, intent(in) :: demag 

                type(random) :: rand 
                integer :: threadID, time, i, atomIndex, numThreads
                real(kind=8) :: rand_num, oldEnergy, newEnergy, P, Z  
                type(vecNd_t) :: S, S_proposed
                integer :: nsteps, MCScounter, demag_update_interval
                real(kind=8), parameter :: pi = 3.14159265358979323846_8
                logical :: calculate_demag 
                integer :: itemp, jtemp, ktemp, N, L, M, stride, iupper, jupper, kupper, cellIndex, atomIndexTemp 
                integer :: atomIncell, counterX, counterY, counterZ, stat
                real(kind=8), dimension(:,:), allocatable :: snapshot
                calculate_demag = .True.
                if (present(demag)) calculate_demag = demag


                demag_update_interval = int((pi * sqrt(5.0_8)) / (2.0_8*r))
                call calculate_demagnetisation_field(chainMesh,chainMesh%demagnetisation_array)
                allocate(snapshot(chainMesh%numAtoms,3),stat=stat)
                if (stat /= 0) error stop "Error allocating spin snapshot"
                !$omp parallel default(private) & 
                !$omp& firstprivate(nsteps,beta, J,J_prime, Dz,Dz_prime, B, r, demag_update_interval, numMCSSweeps) & 
                !$omp&  shared(chainMesh, lockArray, calculate_demag, counterX, counterY, counterZ,snapshot)
                ! parallel book keeping setup
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY
                M = chainMesh%numCellsZ 
                
                stride = 2 ! For now I will decide this statically, later on this should be decided on the fly
                           ! and computed from the interaction range
                iupper = N/stride + modulo(N,stride)
                jupper = L/stride + modulo(L,stride)
                kupper = M/stride + modulo(M,stride)

                numThreads = omp_get_num_threads()
                threadID = omp_get_thread_num()

                call system_clock(time)
                rand = makeRandom(time*threadID + modulo(threadID,time))

                do i = 1,100
                        rand_num = algor_uniform_random(rand) ! Warm up the generator
                end do 
                nsteps = chainMesh%numAtoms

                do MCScounter = 1,numMCSSweeps
                !$omp single 
                        counterX = 1 + nint(algor_uniform_random(rand)*(stride-1))
                        counterY = 1 + nint(algor_uniform_random(rand)*(stride-1))
                        counterZ = 1 + nint(algor_uniform_random(rand)*(stride-1))
                        

                !$omp end single 

                !$omp do 
                do i=1,chainMesh%numAtoms
                        snapshot(i,:) = dble(chainMesh%atomSpins(i,:))
                end do 
                !$omp end do
                !$omp do 
                do i = 1,nsteps
                                !atomIndex = int((algor_uniform_random(rand)*dble(size(chainMesh%atoms)))/&
                                        !(dble(size(chainMesh%atoms))+1) * dble(size(chainMesh%atoms)-1)) + 1 

                                itemp = nint(algor_uniform_random(rand)*(iupper-1))
                                jtemp = nint(algor_uniform_random(rand)*(jupper-1))
                                ktemp = nint(algor_uniform_random(rand)*(kupper-1))

                                itemp = min(itemp*stride + mod(counterX,stride),N-1)
                                jtemp = min(jtemp*stride + mod(counterY,stride),L-1)
                                ktemp = min(ktemp*stride + mod(counterZ,stride),M-1)
                                cellIndex = IndexFromCoordinates(chainMesh,itemp+1,jtemp+1,ktemp+1) ! Select Cell
                                atomInCell = 1 + &
                                        nint(algor_uniform_random(rand)*(chainMesh%chainMeshCells(cellIndex)%NumAtomsPerUnitCell-1))
                                        
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                do atomIndexTemp=1,atomInCell-1
                                        atomIndex = chainMesh%atoms(atomIndex)%nextAtom
                                end do 
                                ! Given atomIndex, propose a new spin then accept/reject based on the energy
                                call OMP_SET_LOCK(lockArray(atomIndex))
                                S = makeVecNdCheck(S,dble(chainMesh%atomSpins(atomIndex,:)))
                                call OMP_UNSET_LOCK(lockArray(atomIndex))
                                S_proposed = S
                                call UniformRandomInSphere(S_proposed,rand,r)
                                S_proposed = S_proposed + S
                                S_proposed = S_proposed / abs(S_proposed)
                                if (any(S_proposed%coords /= S_proposed%coords)) error stop "NaN in MetropolisMCS"
                                if (any(S%coords /= S%coords)) error stop "NaN in MetropolisMCS"
  
                                call calculateHeisenbergEnergy(chainMesh,atomIndex,J,J_prime,Dz,Dz_prime,&
                                        B,lockArray,S_proposed,oldEnergy,newEnergy,chainMesh%demagnetisation_array, calculate_demag)

                                if (newEnergy < oldEnergy) then 
                                        Z = 1.0_8 
                                else 
                                        Z = exp(- beta * ((newEnergy - oldEnergy)))
   

                                end if
                                p = algor_uniform_random(rand)
                                if (z <= 0.6) then 
                                !print *, "Z = ", Z, "New Energy - Old energy ", NewEnergy - OldEnergy, &
                                !                NewEnergy, OldEnergy, "beta = ", beta
                                end if
                                if (Z >= p) then 
                                        call OMP_SET_LOCK(lockArray(atomIndex))
                                        snapshot(atomIndex,:) = S_proposed%coords
                                        call OMP_UNSET_LOCK(lockArray(atomIndex))
                                end if 
                end do
                !$omp end do

                !$omp do 
                do i = 1,chainMesh%numAtoms
                        chainMesh%atomSpins(i,:) = snapshot(i,:)
                end do 
                !$omp end do

                !$omp barrier
                if (mod(MCScounter,demag_update_interval) == 0 .and. calculate_demag) then                         
                        !$omp single
                        call calculate_demagnetisation_field(chainMesh,chainMesh%demagnetisation_array)
                        write(*,*) "Demag maxval, minval = ", &
                                minval(chainMesh%demagnetisation_array), maxval(chainMesh%demagnetisation_array)
                        !$omp end single
                end if 

                end do 
                !$omp end parallel 

        end subroutine Metropolis_mcs


        subroutine MetropolisMixed(chainMesh, beta, nsteps, J, J_prime,  Dz, Dz_prime, B, lockArray, demagnetisation_array)
                ! Each thread will randomly select an atom nsteps times and determine whether to flip the spin 
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), intent(in) :: beta, J, J_prime, Dz, Dz_prime, B 
                integer, intent(in) :: nsteps 
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
                real(kind=8), dimension(:,:), allocatable, intent(inout) :: demagnetisation_array

                type(random) :: rand 
                integer :: threadID, time, i, atomIndex, counter
                real(kind=8) :: rand_num, oldEnergy, newEnergy, P, Z  
                type(vecNd_t) :: S, S_proposed

                if (allocated(demagnetisation_array)) then 
                        if (any(shape(demagnetisation_array) /= [chainMesh%numAtoms,3])) then 
                                deallocate(demagnetisation_array)
                                allocate(demagnetisation_array(chainMesh%numAtoms,3))
                        end if 
                else 
                        allocate(demagnetisation_array(chainMesh%numAtoms,3))
                end if 

                !$omp parallel default(private) firstprivate(nsteps,beta, J,J_prime, Dz,Dz_prime, B) shared(chainMesh, lockArray,&
                !$omp&  demagnetisation_array)
                block 

                threadID = omp_get_thread_num()
                call system_clock(time)
                rand = makeRandom(time*threadID + modulo(threadID,time))
                !rand = makeRandom(123456 + threadID)
                do i = 1,100
                        rand_num = algor_uniform_random(rand) ! Warm up the generator
                end do 
                
                counter = 0
                do i = 1,nsteps
                                atomIndex = int((algor_uniform_random(rand)*dble(size(chainMesh%atoms)))/&
                                        (dble(size(chainMesh%atoms))+1) * dble(size(chainMesh%atoms)-1)) + 1 
                                ! Given atomIndex, propose a new spin then accept/reject based on the energy
                                call OMP_SET_LOCK(lockArray(atomIndex))
                                S = makeVecNdCheck(S,dble(chainMesh%atomSpins(atomIndex,:)))
                                call OMP_UNSET_LOCK(lockArray(atomIndex))
                                S_proposed = S 
                                if (counter == 0) then 
                                        call UniformRandomSpin(S_proposed, rand)
                                else if (counter < 3) then 
                                        call GaussianStep(S_proposed,rand,beta)
                                end if 
                                S_proposed = S_proposed / abs(S_proposed)
                                if (any(S_proposed%coords /= S_proposed%coords)) error stop "NaN in MetropolisMixed"
                                if (any(S%coords /= S%coords)) error stop "NaN in MetropolisMixed"
                                !oldEnergy = calculateHeisenbergEnergy(chainMesh, atomIndex,J,Dz,B,lockArray)
                                !chainMesh%atoms(atomIndex)%AtomParameters = S_proposed%coords 
                                !newEnergy = calculateHeisenbergEnergy(chainMesh, atomIndex,J,Dz,B,lockArray)
                                !chainMesh%atoms(atomIndex)%AtomParameters = S%coords
                                call calculateHeisenbergEnergy(chainMesh,atomIndex,J,J_prime,Dz,Dz_prime,&
                                        B,lockArray,S_proposed,oldEnergy,newEnergy,demagnetisation_array)
                                if (newEnergy < oldEnergy) then 
                                        Z = 1.0_8 
                                else 
                                        Z = exp(- beta * ((newEnergy - oldEnergy)))
   

                                end if
                                p = algor_uniform_random(rand)
                                if (z <= 0.6) then 
                                !print *, "Z = ", Z, "New Energy - Old energy ", NewEnergy - OldEnergy, &
                                !                NewEnergy, OldEnergy, "beta = ", beta
                                end if
                                if (Z >= p) then 
                                        call OMP_SET_LOCK(lockArray(atomIndex))
                                        chainMesh%atomSpins(atomIndex,:) = S_proposed%coords
                                        call OMP_UNSET_LOCK(lockArray(atomIndex))
                                end if 

                        counter = counter + 1 
                        if (counter >= 3) counter = 0
                end do
                end block
                !$omp end parallel 
                

        end subroutine MetropolisMixed

        subroutine Metropolis_demag(chainMesh, beta, nsteps, nsteps_total, J, J_prime,  Dz, Dz_prime, B, lockArray)
                ! Each thread will randomly select an atom nsteps times and determine whether to flip the spin 
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), intent(in) :: beta, J, J_prime, Dz, Dz_prime, B 
                integer, intent(in) :: nsteps, nsteps_total
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
                real(kind=8), dimension(:,:), allocatable, save :: demagnetisation_array

                type(random) :: rand 
                integer :: threadID, time, i, atomIndex, counter
                real(kind=8) :: rand_num, oldEnergy, newEnergy, P, Z  
                type(vecNd_t) :: S, S_proposed

                counter = 0
                if (allocated(demagnetisation_array)) then 
                        if (any(shape(demagnetisation_array) /= [chainMesh%numAtoms,3])) then 
                                deallocate(demagnetisation_array)
                                allocate(demagnetisation_array(chainMesh%numAtoms,3))
                        end if 
                else 
                        allocate(demagnetisation_array(chainMesh%numAtoms,3))
                end if 

                do while (counter < nsteps_total)

                        

                        call calculate_demagnetisation_field(chainMesh,demagnetisation_array)
                        call MetropolisMixed(chainMesh, beta, nsteps, J, J_prime,  Dz, Dz_prime, B, lockArray,demagnetisation_array)
                        counter = counter + nsteps
                        print *, "Completed metropolis step with counter, total = ", counter, nsteps_total
                end do 

        end subroutine Metropolis_demag
        

        subroutine TotalHeisenbergEnergy(chainMesh, J, J_prime, Dz, Dz_prime, B, lockArray, totalEnergy,demag)
        implicit none

        type(ChainMesh_t), intent(inout)            :: chainMesh
        real(kind=8),     intent(in)                :: J, J_prime, Dz, Dz_prime, B
        integer(kind=OMP_LOCK_KIND), intent(inout)  :: lockArray(:)
        real(kind=8),     intent(out)               :: totalEnergy
        logical, intent(in) :: demag

        integer :: atomIndex
        real(kind=8) :: energyTemp

        call calculate_demagnetisation_field(chainMesh,chainMesh%demagnetisation_array)
        ! Initialize
        totalEnergy = 0.0_8
        energyTemp = 0.0_8
        ! create a dummy 3-vector (only used for newEnergy, which we ignore)

        do atomIndex = 1, size(chainMesh%atoms)
                call calculateTotalHeisenbergEnergyHelper( &
                     chainMesh, atomIndex,      &
                     J, J_prime, Dz, Dz_prime, B, &
                     lockArray,                 &
                     energyTemp,                & ! only oldE is used
                     demag                      &
                )
                totalEnergy = totalEnergy + energyTemp
         end do

        end subroutine TotalHeisenbergEnergy

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


subroutine write_demagnetisation_field_to_file(chainMesh,demag_field,filepath)
        implicit none
        type(ChainMesh_t), intent(in) :: chainMesh 
        real(kind=8), dimension(:,:), intent(in) :: demag_field 
        character(len=*),  intent(in) :: filepath

        integer :: iostat, i
        real(kind=8) :: x,y,z, x_component, y_component, z_component
        open(unit=100, status="replace", action="write",file=trim(filepath), iostat = iostat)
        if (iostat /= 0) error stop "Error opening demagnetisation file"
        
        write(100, '(A)') "x,y,z,fx,fy,fz"
        do i = 1,size(chainMesh%atoms)
                x = chainMesh%atoms(i)%x 
                y = chainMesh%atoms(i)%y 
                z = chainMesh%atoms(i)%z
                x_component = demag_field(i,1)
                y_component = demag_field(i,2)
                z_component = demag_field(i,3)
                write(100, '(F8.4,"," F8.4,",", F8.4,",", F8.4, ",", F8.4, ",", F8.4)') x,y,z, x_component, y_component, z_component

        end do 
        close(100)
end subroutine write_demagnetisation_field_to_file

end module EnergyMin 
