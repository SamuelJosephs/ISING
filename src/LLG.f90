module llg
        ! This module will provide the tools to work with a chain Mesh, initialise a Skyrmion, evolve it over time using the LLG
        ! equation, and write the output in a series of files that are are suited to animation and analysis in python.
        use chainMesh 
        use vecNd 
        implicit none 
        

        contains  



        subroutine initialise_skyrmion(chainMesh, r, R_s)
                
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                type(vecNd_t), intent(in) :: r 
                real(kind=8), intent(in) :: R_s 
                real(kind=8) :: d, dx, dy, widthX, widthY, widthZ, HalfWidthX, HalfwidthY, HalfwidthZ, psi, phi, gaussian_factor
                integer :: i
                type(vecNd_t) :: atomPos
                ! r: Skyrmion centre 
                ! R_s: Skyrmion radius
                if (size(r) /= 3) then 
                        error stop "Skyrmion position must have three components"
                end if 
                HalfwidthX = chainMesh%numCellsX*chainMesh%latticeParameter / 2
                HalfwidthY = chainMesh%numCellsY*chainMesh%latticeParameter / 2
                HalfwidthZ = chainMesh%numCellsZ*chainMesh%latticeParameter / 2
                widthX = chainMesh%numCellsX*chainMesh%latticeParameter
                widthY = chainMesh%numCellsY*chainMesh%latticeParameter
                widthZ = chainMesh%numCellsZ*chainMesh%latticeParameter
                do i = 1,size(chainMesh%atoms)
       
                        atomPos = makeVecNd((/dble(chainMesh%atoms(i)%x),dble(chainMesh%atoms(i)%y),dble(chainMesh%atoms(i)%z)/))
                        d = distance_points(chainMesh,r,atomPos) 
                        dx = abs(r%coords(1) - atomPos%coords(1))
                        dy = abs(r%coords(2) - atomPos%coords(2))
                         
                        if (dx > HalfwidthX) dx = widthX - dx 
                        if (dy > HalfwidthY) dy = widthY - dy

                        psi = 2*atan2(d,R_s)
                        phi = 2*atan2(dy,dx)
                        gaussian_factor = exp(- d**2 / sqrt(R_s/2))
                        if (gaussian_factor > 1e-4) then 
                                chainMesh%atoms(i)%AtomParameters(1) = gaussian_factor*sin(psi)*sin(phi)
                                chainMesh%atoms(i)%AtomParameters(2) = gaussian_factor*sin(psi)*cos(phi)
                                chainMesh%atoms(i)%AtomParameters(3) = gaussian_factor*cos(psi)
                        end if
                        call normalise_spin(chainMesh,i,1.0_8)
                end do 
        end subroutine initialise_skyrmion

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
                write(fileunit, '(F15.8, 5(",", F15.8))', iostat=iostat) &
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

end module llg 
