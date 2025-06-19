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
                        chainMesh%atomSpins(i,:) = m%coords
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
                mag = sqrt(chainMesh%atomSpins(i,1)**2 + &
                                chainMesh%atomSpins(i,2)**2 + &
                                        chainMesh%atomSpins(i,3)**2)
                if (mag < 1e-8) then 
                        return
                end if 
                chainMesh%atomSpins(i,1) = chainMesh%atomSpins(i,1)/mag * val  
                chainMesh%atomSpins(i,2) = chainMesh%atomSpins(i,2)/mag * val 
                chainMesh%atomSpins(i,3) = chainMesh%atomSpins(i,3)/mag * val 

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
            character(len=100) :: command_buff
            do i=1,size(chainMesh%atoms)
                if (any(chainMesh%atomSpins(i,:) /= chainMesh%atomSpins(i,:))) then 
                        print *, "NaN value encoutered in atom ", i, "with atomParameters:", chainMesh%atomSpins(i,:)
                        error stop "NaN value encountered"
                end if
            end do 
            ! Set defaults for optional parameters
            do_append = .false.
            if (present(append)) do_append = append
            
            do_include_header = .true.
            if (present(include_header)) do_include_header = include_header
            
            ! Open the file
         
            call EXECUTE_COMMAND_LINE("mkdir -p " // filename)
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
                if (any(chainMesh%atomSpins(i,:) /= chainMesh%atomSpins(i,:))) error stop "Nan Encountered"
                write(fileunit, '(ES22.14, 5(",", ES22.14))', iostat=iostat) &
                    chainMesh%atoms(i)%x, &
                    chainMesh%atoms(i)%y, &
                    chainMesh%atoms(i)%z, &
                    chainMesh%atomSpins(i,1), &
                    chainMesh%atomSpins(i,2), &
                    chainMesh%atomSpins(i,3)
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

end module llg 
