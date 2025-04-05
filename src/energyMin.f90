

module EnergyMin
        use Atom 
        use ChainMesh 
        use ChainMeshCell
        use Rand, only: makeRandom, algor_uniform_random, random, Normal
        use cube_partition, only: bounds_type, partition_unit_cube, dp 
        use constants, only: gyromagnetic_ratio, Kb, Bohr_magneton
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

        subroutine calculateHeisenbergEnergy(chainMesh,atomIndex, J, Dz,B, lockArray,S_proposed, oldEnergy, newEnergy) 
                implicit none
                type(ChainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: atomIndex 
                real(kind=8), intent(in) :: J, Dz, B 
                real(kind=8), intent(inout) ::  oldEnergy, newEnergy 
                type(vecNd_t), intent(in) :: S_proposed 
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)

                real(kind = 8) :: Energy, x,y,z 
                integer :: i, atomIndexTemp
                type(vecNd_t) :: S, S_prime, atomPos1, atomPos2, tempVec,r, D 
                !call OMP_SET_LOCK(lockArray(atomIndex))
                oldEnergy = 0.0_8
                newEnergy = 0.0_8
                S = makeVecNdCheck(S, dble(chainMesh%atoms(atomIndex)%AtomParameters))
                x = chainMesh%atoms(atomIndex)%x
                y = chainMesh%atoms(atomIndex)%y 
                z = chainMesh%atoms(atomIndex)%z 
                atomPos1 = makeVecNd([x,y,z]) 
                tempVec = makeVecNdCheck(tempVec, [0.0_8, 0.0_8, Dz])
                do i = 1,size(chainMesh%atoms(atomIndex)%NeighborList)
                        atomIndexTemp = chainMesh%atoms(atomIndex)%NeighborList(i)
                        if (atomIndexTemp == atomIndex) error stop "Encountered self interaction"
                        call OMP_SET_LOCK(lockArray(atomIndexTemp))
                                S_prime = makeVecNdCheck(S_prime,dble(chainMesh%atoms(atomIndexTemp)%AtomParameters))
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
                        oldEnergy = oldEnergy + (J* S*S_prime) + (D*(S .x. S_prime)) + &
                                B*S%coords(3)
                        
                        newEnergy = newEnergy + (J* S_proposed*S_prime) + (D*(S_proposed .x. S_prime)) + &
                                B*S_proposed%coords(3)
                end do 

                !call OMP_UNSET_LOCK(lockArray(atomIndex))
        end subroutine calculateHeisenbergEnergy


        subroutine MetropolisMixed(chainMesh, beta, nsteps, J, Dz, B, lockArray)
                ! Each thread will randomly select an atom nsteps times and determine whether to flip the spin 
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), intent(in) :: beta, J, Dz, B 
                integer, intent(in) :: nsteps 
                integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)

                type(random) :: rand 
                integer :: threadID, time, i, atomIndex, counter
                real(kind=8) :: rand_num, oldEnergy, newEnergy, P, Z  
                type(vecNd_t) :: S, S_proposed
                !$omp parallel default(private) firstprivate(nsteps,beta, J, Dz, B) shared(chainMesh, lockArray)
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
                                S = makeVecNdCheck(S,dble(chainMesh%atoms(atomIndex)%AtomParameters))
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
                                call calculateHeisenbergEnergy(chainMesh,atomIndex,J,Dz,B,lockArray,S_proposed,oldEnergy,newEnergy)
                                if (newEnergy < oldEnergy) then 
                                        Z = 1.0_8 
                                else 
                                        Z = exp(- beta * ((newEnergy - oldEnergy)))
   

                                end if
                                p = algor_uniform_random(rand)
                                !if (z <= 0.6) then 
                                !print *, "Z = ", Z, "New Energy - Old energy ", NewEnergy - OldEnergy, &
                                !                NewEnergy, OldEnergy, "beta = ", beta
                                !end if
                                if (Z >= p) then 
                                        call OMP_SET_LOCK(lockArray(atomIndex))
                                        chainMesh%atoms(atomIndex)%AtomParameters = S_proposed%coords
                                        call OMP_UNSET_LOCK(lockArray(atomIndex))
                                end if 

                        counter = counter + 1 
                        if (counter >= 3) counter = 0
                end do
                end block
                !$omp end parallel 
                

        end subroutine MetropolisMixed

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
