module llg
        ! This module will provide the tools to work with a chain Mesh, initialise a Skyrmion, evolve it over time using the LLG
        ! equation, and write the output in a series of files that are are suited to animation and analysis in python.
        use chainMesh 
        use vecNd 
        implicit none 
        

        contains  
        
        subroutine initialise_skyrmion(chainMesh, r)
               type(ChainMesh_t), intent(inout) :: chainMesh 
               type(vecNd_t), intent(in) :: r
               real(kind=8) :: theta, phi 
               integer :: i

               if (size(r) /= 3) then 
                       error stop "R must have three spacial components"
               end if 

               ! Loop through every atom and initialise the atoms spins  
               do i = 1,size(chainMesh%atoms) !TODO: This isn't right
                    ! Calculate position relative to skyrmion center
                   rx = abs(chainMesh%atoms(i)%x - r(1))
                   ry = abs(chainMesh%atoms(i)%y - r(2))
                   
                   ! Calculate distance from skyrmion center in the xy-plane
                   rho = sqrt(rx**2 + ry**2)
                   
                   ! Calculate azimuthal angle (in-plane angle)
                   phi = atan2(ry, rx)
                   
                   ! Calculate polar angle (out-of-plane angle)
                   ! Use a profile function for the skyrmion
                   R_s = 5.0  ! Skyrmion radius parameter
                   theta = 3.14159265359 * (1.0 - exp(-rho/R_s))
                   
                   ! Set the spin components in spherical coordinates
                   chainMesh%atoms(i)%AtomParameters(1) = sin(theta) * cos(phi)  ! Sx
                   chainMesh%atoms(i)%AtomParameters(2) = sin(theta) * sin(phi)  ! Sy
                   chainMesh%atoms(i)%AtomParameters(3) = cos(theta)             ! Sz  
               end do 
        end subroutine initialise_skyrmion

end module llg 
