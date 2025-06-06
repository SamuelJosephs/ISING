module llg
        ! This module will provide the tools to work with a chain Mesh, initialise a Skyrmion, evolve it over time using the LLG
        ! equation, and write the output in a series of files that are are suited to animation and analysis in python.
        use chainMesh 
        use vecNd 
        use stereographicProjection, only: NSphereProjection
        use omp_lib
        use constants 
        use reciprocal_space_processes
        implicit none 
        
        abstract interface 
        function H_eff_class(Mesh, atomIndex,lockArray, J, Dz, B, demagnetisation_array) result(E)
            use chainMesh, only: chainMesh_t
            use vecNd, only: vecNd_t
            use OMP_LIB 
            type(chainMesh_t), intent(inout) :: Mesh
            integer, intent(in) :: atomIndex 
            type(vecNd_t) :: E
            integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
            real(kind=8), intent(in) :: J, Dz, B 
            real(kind=8), dimension(:,:), intent(in) :: demagnetisation_array
        end function H_eff_class
        end interface 
        contains  


        subroutine initialise_skyrmion_sp(chainMesh, r, R_s, chi, N)
                
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                type(vecNd_t), intent(in) :: r
                real(kind=8), intent(in) :: R_s, chi  
                integer, intent(in) :: N ! winding number
                integer :: i
                real(kind=8), allocatable :: angle(:)
                type(vecNd_t) :: AtomPos, SkyrmionCentre, m 
                real(kind=8) :: maxTheta, maxPhi, minTheta, minPhi
                ! For each atom at position (x_A, y_A, z_A), calculate the angles from a stereographic projection 
                ! from a sphere at position (r_x, r_y, z_A). This is a 2d spin texture "tornadoised" down the z axis.
                ! Then from this projection set 
                ! x = r sin(theta) cos(psi + chi)
                ! y = r sin(theta) sin(psi + chi)
                ! z = r cos(theta)
                AtomPos = makeVecNd([0.0_08,0.0_08,0.0_08])
                SkyrmionCentre = makeVecNd([r%coords(1), r%coords(2), 0.0_08])
                maxTheta = TINY(maxTheta)
                maxPhi = TINY(maxPhi)
                minTheta = HUGE(minTheta)
                minPhi = HUGE(minPhi)
                do i = 1,size(chainMesh%atoms)
                        AtomPos%coords(1) = chainMesh%atoms(i)%x 
                        AtomPos%coords(2) = chainMesh%atoms(i)%y
                        AtomPos%coords(3) = chainMesh%atoms(i)%z 

                        SkyrmionCentre%coords(3) = chainMesh%atoms(i)%z + R_s 

                        angle = NSphereProjection(AtomPos,SkyrmionCentre,R_s)
                        if (any(angle /= angle)) error stop "NaN encountered"
                        m = makeVecNd([sin(angle(1))*cos(N*(angle(2) + chi)), &
                                       sin(angle(1))*sin(N*(angle(2)+chi)), &
                                       cos(angle(1))])
                        if (any(m%coords /= m%coords)) error stop "NaN encountered in m"
                        chainMesh%atoms(i)%AtomParameters = m%coords
                        if (angle(1) > maxTheta) maxTheta = angle(1)
                        if (angle(2) > maxPhi) maxPhi = angle(2)
                        
                        if (angle(1) < minTheta) minTheta = angle(1)
                        if (angle(2) < minPhi) minPhi = angle(2)
                end do 
                print *, "maxTheta, minTheta =  ", maxTheta, mintheta  
                print *, "maxPhi, minphi = ", maxPhi, minphi
        end subroutine initialise_skyrmion_sp


        subroutine normalise_spin(chainMesh,i, val)
                type(ChainMesh_t), intent(inout) :: chainMesh
                integer, intent(in) :: i 
                real(kind=8), intent(in) :: val
                real(kind=8) :: mag
                mag = sqrt(chainMesh%atoms(i)%AtomParameters(1)**2 + &
                                chainMesh%atoms(i)%AtomParameters(2)**2 + &
                                        chainMesh%atoms(i)%AtomParameters(3)**2)
                if (mag < 1e-8) then 
                        return
                end if 
                chainMesh%atoms(i)%AtomParameters(1) = chainMesh%atoms(i)%AtomParameters(1)/mag * val  
                chainMesh%atoms(i)%AtomParameters(2) = chainMesh%atoms(i)%AtomParameters(2)/mag * val 
                chainMesh%atoms(i)%AtomParameters(3) = chainMesh%atoms(i)%AtomParameters(3)/mag * val 

        end subroutine normalise_spin

        subroutine write_spins_to_file(chainMesh, filename, append, include_header, frame_index)
            type(ChainMesh_t), intent(in) :: chainMesh
            character(len=*), intent(in) :: filename
            logical, intent(in), optional :: append        ! Append to existing file instead of replacing
            logical, intent(in), optional :: include_header ! Include header row
            integer, intent(in), optional :: frame_index   ! For time-evolution animations
            
            integer :: i, fileunit, iostat
            logical :: do_append, do_include_header
            character(len=20) :: frame_str
            do i=1,size(chainMesh%atoms)
                if (any(chainMesh%atoms(i)%AtomParameters /= chainMesh%atoms(i)%AtomParameters)) then 
                        print *, "NaN value encoutered in atom ", i, "with atomParameters:", chainMesh%atoms(i)%AtomParameters
                        error stop "NaN value encountered"
                end if
            end do 
            ! Set defaults for optional parameters
            do_append = .false.
            if (present(append)) do_append = append
            
            do_include_header = .true.
            if (present(include_header)) do_include_header = include_header
            
            ! Open the file
            if (do_append) then
                open(newunit=fileunit, file=filename, status='unknown', position='append', action='write', iostat=iostat)
            else
                open(newunit=fileunit, file=filename, status='replace', action='write', iostat=iostat)
            end if
            
            if (iostat /= 0) then
                print *, "Error opening file: ", filename
                return
            end if
            
            ! Write frame index if provided
            if (present(frame_index)) then
                write(frame_str, '(I0)') frame_index
                write(fileunit, '(A)', iostat=iostat) "# Frame: " // trim(frame_str)
                if (iostat /= 0) then
                    print *, "Error writing frame index to file: ", filename
                    close(fileunit)
                    return
                end if
            end if
            
            ! Write header if required
            if (do_include_header) then
                write(fileunit, '(A)', iostat=iostat) "x,y,z,Sx,Sy,Sz"
                if (iostat /= 0) then
                    print *, "Error writing header to file: ", filename
                    close(fileunit)
                    return
                end if
            end if
            
            ! Write data for each atom
            do i = 1, size(chainMesh%atoms)
                if (any(chainMesh%atoms(i)%AtomParameters /= chainMesh%atoms(i)%AtomParameters)) error stop "Nan Encountered"
                write(fileunit, '(ES22.14, 5(",", ES22.14))', iostat=iostat) &
                    chainMesh%atoms(i)%x, &
                    chainMesh%atoms(i)%y, &
                    chainMesh%atoms(i)%z, &
                    chainMesh%atoms(i)%AtomParameters(1), &
                    chainMesh%atoms(i)%AtomParameters(2), &
                    chainMesh%atoms(i)%AtomParameters(3)
                if (iostat /= 0) then
                    print *, "Error writing data to file: ", filename
                    close(fileunit)
                    return
                end if
            end do
            
            ! Close the file
            close(fileunit)
            
            print *, "Successfully wrote spin data to: ", trim(filename)
        end subroutine write_spins_to_file


subroutine LLGStep(chainMesh, dt, A, B, C, D, H)
    ! This subroutine evolves the spin lattice on the chain mesh using the LLG equation
    ! 
    ! Parameters:
    ! chainMesh - The chain mesh containing the atoms with their spins
    ! dt - Time step for evolution
    ! A, B, C - Exchange interaction parameters for x, y, z components
    ! D - Anisotropy parameter
    ! H - External magnetic field vector (3 components)
    
    implicit none
    type(ChainMesh_t), intent(inout) :: chainMesh
    real(kind=8), intent(in) :: dt, A, B, C, D
    real(kind=8), intent(in) :: H(3)  ! Assuming H is a 3D vector
    
    integer :: i, j, neighbor_idx
    real(kind=8) :: effective_field(3), cross_product(3), dS(3)
    type(vecNd_t) :: S_vec, H_vec, eff_field_vec, dS_vec, cross_vec
    
    ! Create vector for H
    H_vec = makeVecNd(H)
    
    ! Loop over all atoms in the chain mesh
    do i = 1, size(chainMesh%atoms)
        ! Reset effective field for each atom
        effective_field = 0.0d0
        
        ! Get the current spin vector
        S_vec = makeVecNd((/dble(chainMesh%atoms(i)%AtomParameters(1)),  &
                           dble(chainMesh%atoms(i)%AtomParameters(2)),  &
                           dble(chainMesh%atoms(i)%AtomParameters(3))/))
        
        ! Sum contributions from nearest neighbors
        do j = 1, size(chainMesh%atoms(i)%NeighborList)
            neighbor_idx = chainMesh%atoms(i)%NeighborList(j)
            
            ! Skip self-interaction
            if (neighbor_idx == i) then
                cycle
            end if
            
            ! Add x-component contribution
            effective_field(1) = effective_field(1) + A * chainMesh%atoms(neighbor_idx)%AtomParameters(1)
            
            ! Add y-component contribution
            effective_field(2) = effective_field(2) + B * chainMesh%atoms(neighbor_idx)%AtomParameters(2)
            
            ! Add z-component contribution
            effective_field(3) = effective_field(3) + C * chainMesh%atoms(neighbor_idx)%AtomParameters(3)
        end do
        
        ! Add anisotropy term (assuming z-axis anisotropy)
        effective_field(3) = effective_field(3) + 2.0d0 * D * chainMesh%atoms(i)%AtomParameters(3)
        
        ! Add external magnetic field
        effective_field(1) = effective_field(1) + H(1)
        effective_field(2) = effective_field(2) + H(2)
        effective_field(3) = effective_field(3) + H(3)
        
        ! Create vector for effective field
        eff_field_vec = makeVecNd(effective_field)
        
        ! Calculate cross product S x H_eff (first term in LLG equation)
        ! Using the cross product operation from vecNd module
        cross_vec = S_vec .x. eff_field_vec
        
        ! Extract components from cross product
        dS(1) = cross_vec%coords(1)
        dS(2) = cross_vec%coords(2)
        dS(3) = cross_vec%coords(3)
        
        ! Update spin using forward Euler method
        chainMesh%atoms(i)%AtomParameters(1) = chainMesh%atoms(i)%AtomParameters(1) + dt * dS(1)
        chainMesh%atoms(i)%AtomParameters(2) = chainMesh%atoms(i)%AtomParameters(2) + dt * dS(2)
        chainMesh%atoms(i)%AtomParameters(3) = chainMesh%atoms(i)%AtomParameters(3) + dt * dS(3)
        
        ! Normalize the spin to maintain unit length
        call normalise_spin(chainMesh, i, 1.0d0)
    end do
    
end subroutine LLGStep

function H_eff_Heisenberg(Mesh, atomIndex,lockArray, J, Dz, B, demagnetisation_array) result(H_temp)
    use chainMesh, only: chainMesh_t
    use vecNd, only: vecNd_t
    type(chainMesh_t), intent(inout) :: Mesh
    integer, intent(in) :: atomIndex 
    integer(kind=OMP_LOCK_KIND), intent(inout) :: lockArray(:)
    real(kind=8), intent(in) :: J, Dz, B
    real(kind=8), dimension(:,:), intent(in) :: demagnetisation_array
    type(vecNd_t) :: E
    integer :: atomIndexTemp,i, threadNum
    type(vecNd_t) :: H_temp
    type(vecNd_t) :: D, atomPos1, atomPos2,r, tempVec, S_temp
    real(kind=8) :: demag_x, demag_y, demag_z 

    real(kind=8) :: x,y,z

    threadNum = omp_get_thread_num()

    atomIndexTemp = atomIndex 

    H_temp = makeVecNdCheck(H_temp,[0.0_8,0.0_8,0.0_8])
    atomPos2 = makeVecNdCheck(atomPos2,[0.0_8,0.0_8,0.0_8])

    D = makeVecNdCheck(D,[0.0_8,0.0_8,0.0_8])
    x = dble(Mesh%atoms(atomIndexTemp)%x)
    y = dble(Mesh%atoms(atomIndexTemp)%y)
    z = dble(Mesh%atoms(atomIndexTemp)%z)
    atomPos1 = makeVecNdCheck(atomPos1,[x,y,z])
    tempVec = makeVecNdCheck(tempVec,[0.0_8, 0.0_8, Dz])
    demag_x = demagnetisation_array(atomIndex,1)
    demag_y = demagnetisation_array(atomIndex,2)
    demag_z = demagnetisation_array(atomIndex,3)
    do i = 1,size(Mesh%atoms(atomIndex)%NeighborList)
        !call OMP_SET_LOCK(lockArray(Mesh%atoms(atomIndex)%NeighborList(i)))
        atomIndexTemp = Mesh%atoms(atomIndex)%NeighborList(i) !NeighborList is not written to so don't need to guard against race
                                                              ! conditions
        !call OMP_UNSET_LOCK(lockArray(Mesh%atoms(atomIndex)%NeighborList(i)))
        if (atomIndexTemp == atomIndex) cycle
        x = dble(Mesh%atoms(atomIndexTemp)%x)
        y = dble(Mesh%atoms(atomIndexTemp)%y)
        z = dble(Mesh%atoms(atomIndexTemp)%z)
        atomPos2%coords = [x,y,z]
        call distance_points_vec(Mesh,atomPos1,atomPos2, r) 
        r =  (-1.0_8) * r / abs(r)
        D = tempVec .x. r
        call OMP_SET_LOCK(lockArray(atomIndexTemp))
        S_temp = makeVecNdCheck(S_temp,dble(Mesh%atoms(atomIndexTemp)%atomParameters))
        call OMP_UNSET_LOCK(lockArray(atomIndexTemp))
        !H_temp = H_temp - (J*S_temp + (D .x. S_temp))
        H_temp = H_temp + (J*S_temp + (S_temp .x. D))
    end do
    H_temp%coords(3) = H_temp%coords(3) + B
    H_temp%coords = H_temp%coords - 0.5_8 * [demag_x, demag_y, demag_z]
    H_temp = (-1.0_8)*H_temp ! Testing minus signs

end function H_eff_Heisenberg

subroutine HeunStep(chainMesh, numSteps, dt, H_eff_method, lambda,gamma, J, Dz, B)
    implicit none
    type(chainMesh_t), intent(inout) :: chainMesh
    integer, intent(in) :: numSteps
    real(kind=8), intent(in) :: dt
    real(kind=8), intent(in) :: lambda, gamma, J, Dz, B 
    procedure(H_eff_class), pointer, intent(in) :: H_eff_method
    type(vecNd_t) :: S_prime, S_next, S_temp, H, delta_S, delta_S_prime, H_prime, test_temp
    integer :: atomIndex, i, threadNum, counter, stat
    integer(kind=OMP_LOCK_KIND), allocatable :: lockArray(:)
    real(kind=8), allocatable, dimension(:,:) :: demagnetisation_array
    allocate(lockArray(size(chainMesh%atoms)))
    allocate(demagnetisation_array(chainMesh%numAtoms,3),stat=stat)
    if (stat /= 0) error stop "Error allocating demagnetisation_array"
    do i = 1,size(lockArray)
        call OMP_INIT_LOCK(lockArray(i))
    end do 
    counter = 0
    do i = 1, numSteps
    print *, "Entered Heun Step" 
    call calculate_demagnetisation_field(chainMesh,demagnetisation_array)
    print *, "Heun Step calculated demagnetisation array"
    !$omp parallel do shared(chainMesh,lockArray, H_eff_method,demagnetisation_array) &
    !$omp&  default(private) firstprivate(counter, lambda, gamma,dt,J,Dz,B)

        do atomIndex = 1,size(chainMesh%atoms)
            S_prime = makeVecNdCheck(S_prime,[0.0_8,0.0_8,0.0_8])
            S_next = S_prime
            S_temp = S_prime 
            H = S_prime
            delta_S = S_prime
            delta_S_prime = S_prime
            threadNum = OMP_GET_THREAD_NUM()
            ! For each atom, calculate H_eff using H_eff_method. Then calculate S' and S'' before doing updating the spins
            S_temp = makeVecNdCheck(S_temp,dble(chainMesh%atoms(atomIndex)%atomParameters))
            H = H_eff_method(chainMesh,atomIndex,lockArray, J, Dz, B,demagnetisation_array)
            delta_S = (- gamma / (1 + lambda**2))*((S_temp .x. H )+ ((lambda*S_temp) .x. (S_temp .x. H))) 
            if (.not. allocated(S_temp%coords)) print *, "S_temp is not allocated for i = ",i
            if (.not. allocated(S_prime%coords)) print *, "S_prime is not allocated for i = ",i

            S_prime = S_temp + delta_s*dt
            S_prime%coords = S_prime%coords / abs(S_prime)
            !call OMP_SET_LOCK(lockArray(atomIndex))
                chainMesh%atoms(atomIndex)%AtomParameters = S_prime%coords
                
            !call OMP_UNSET_LOCK(lockArray(atomIndex))
            H_prime = H_eff_method(chainMesh,atomIndex,lockArray, J, Dz, B,demagnetisation_array)

            !call OMP_SET_LOCK(lockArray(atomIndex))
                chainMesh%atoms(atomIndex)%AtomParameters = S_temp%coords
                
            !call OMP_UNSET_LOCK(lockArray(atomIndex))
            delta_S_prime = (- gamma / (1 + lambda**2))*((S_prime .x. H_prime) + ((lambda*S_prime) .x. (S_prime .x. H_prime))) 

            S_next = S_temp + 0.5_8*(delta_S + delta_S_prime)*dt
            S_next%coords = S_next%coords / abs(S_next) 
            !call OMP_SET_LOCK(lockArray(atomIndex))
            chainMesh%atoms(atomIndex)%atomParameters = S_next%coords
            !print *, "delta_S = ", delta_S%coords 
            !print *, "delta_S_prime = ", delta_S_prime%coords
            !test_temp = S_next - S_temp 
            !print *, "S_next - S_temp = ", test_temp%coords 
            !call OMP_UNSET_LOCK(lockArray(atomIndex))
            counter = counter + 1
            

        end do
        !$omp end parallel do
    end do

    do i = 1,size(lockArray)
        call OMP_DESTROY_LOCK(lockArray(i))
    end do 
end subroutine HeunStep

end module llg 
