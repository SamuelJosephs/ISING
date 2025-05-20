module reciprocal_space_processes 
        use chainMesh 
        use vecNd
        use constants
        use, intrinsic :: iso_c_binding 
        !include 'fftw3.f03' !This is included in chainMesh, annoying I would like header files with #pragma once
        
        contains
        
        subroutine interpolate_to_fft_array(chainMesh) 
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer :: atom, iCell, jCell, kCell, chainCellIndex, atomIndex, N,L,M , stat, cellIndex, cellIndexTemp  
                integer, dimension(27) :: nearestNeighborCellList
                type(vecNd_t) :: cellCentre, atomPos, d 
                real(kind = 8) :: x,y,z, xWeight, yWeight, zWeight, cellWidth, weight  
                cellWidth = dble(chainMesh%latticeParameter)
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 
                chainMesh%fft_array_x = real(0.0,C_DOUBLE)
                chainMesh%fft_array_y = real(0.0,C_DOUBLE) 
                chainMesh%fft_array_z = real(0.0,C_DOUBLE)
                do CellIndex = 1,chainMesh%numChainMeshCells
                        call getNeighboringCells(chainMesh,CellIndex,nearestNeighborCellList)
                        call coordinatesFromIndex(chainMesh,CellIndeX, iCell, jCell, kCell)
                        iCell = iCell + 1 ! coordinatesFromIndex works for 0 based indexing 
                        jCell = jCell + 1 
                        kCell = kCell + 1 
                        x = chainMesh%chainMeshCells(CellIndex)%centreX
                        y = chainMesh%chainMeshCells(cellIndex)%centreY
                        z = chainMesh%chainMeshCells(cellIndex)%centreZ
                        cellCentre = makeVecNdCheck(cellCentre,[x,y,z])
                        do cellIndexTemp = 1,size(nearestNeighborCellList)
                               atom = chainMesh%chainMeshCells(nearestNeighborCellList(cellIndexTemp))%firstAtomInMeshCell
                               do while (atom /= -1)
                                   x = dble(chainMesh%atoms(atom)%x) 
                                   y = dble(chainMesh%atoms(atom)%y)
                                   z = dble(chainMesh%atoms(atom)%z) 
                                   atomPos = makeVecNdCheck(atomPos,[x,y,z])
                                   call distance_points_vec(chainMesh,cellCentre,atomPos,d) !Takes periodic boundaries into account.
                                   
                                   !Now calculate the cloud in cell interpolation weights 

                                   xWeight = 1.0_8 - (abs(d%coords(1))/cellWidth)
                                   yWeight = 1.0_8 - (abs(d%coords(2))/cellWidth) 
                                   zWeight = 1.0_8 - (abs(d%coords(3))/cellWidth)
                                   weight = xWeight*yWeight*zWeight
                                   chainMesh%fft_array_x(iCell, jCell, kCell) = chainMesh%fft_array_x(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atoms(atom)%AtomParameters(1) 
                                   chainMesh%fft_array_y(iCell, jCell, kCell) = chainMesh%fft_array_y(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atoms(atom)%AtomParameters(2) 
                                   chainMesh%fft_array_z(iCell, jCell, kCell) = chainMesh%fft_array_z(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atoms(atom)%AtomParameters(3) 
 
                                   atom = chainMesh%atoms(atom)%nextAtom 
                               end do 
                        end do 
                        
                end do 
        end subroutine interpolate_to_fft_array 
        
        subroutine interpolate_fft_to_atoms(chainMesh, outputArray)
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), allocatable, dimension(:,:), intent(out) :: outputArray ! i,j in [1,size(chainmesh%atoms)], [1,3]
                integer :: atom, iCell, jCell, kCell, chainCellIndex, atomIndex, N,L,M , stat, cellIndex, cellIndexTemp  
                integer, dimension(27) :: nearestNeighborCellList
                type(vecNd_t) :: cellCentre, atomPos, d 
                real(kind = 8) :: x,y,z, xWeight, yWeight, zWeight, cellWidth, weight  
                cellWidth = dble(chainMesh%latticeParameter)
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 
                allocate(outputArray(size(chainMesh%atoms),3), stat=stat) 
                outputArray(:,:) = 0.0_8
                if (stat /= 0) error stop "failed to allocate output array"
                do CellIndex = 1,chainMesh%numChainMeshCells
                        call getNeighboringCells(chainMesh,CellIndex,nearestNeighborCellList)
                        call coordinatesFromIndex(chainMesh,CellIndeX, iCell, jCell, kCell)
                        iCell = iCell + 1 ! coordinatesFromIndex works with 0 based indexing 
                        jCell = jCell + 1 
                        kCell = kCell + 1 
                        x = chainMesh%chainMeshCells(CellIndex)%centreX
                        y = chainMesh%chainMeshCells(cellIndex)%centreY
                        z = chainMesh%chainMeshCells(cellIndex)%centreZ
                        cellCentre = makeVecNdCheck(cellCentre,[x,y,z])
                        do cellIndexTemp = 1,size(nearestNeighborCellList)
                               atom = chainMesh%chainMeshCells(nearestNeighborCellList(cellIndexTemp))%firstAtomInMeshCell
                               do while (atom /= -1)
                                   x = dble(chainMesh%atoms(atom)%x) 
                                   y = dble(chainMesh%atoms(atom)%y)
                                   z = dble(chainMesh%atoms(atom)%z) 
                                   atomPos = makeVecNdCheck(atomPos,[x,y,z])
                                   call distance_points_vec(chainMesh,cellCentre,atomPos,d) !Takes periodic boundaries into account.
                                   
                                   !Now calculate the cloud in cell interpolation weights 

                                   xWeight = 1.0_8 - (abs(d%coords(1))/cellWidth)
                                   yWeight = 1.0_8 - (abs(d%coords(2))/cellWidth) 
                                   zWeight = 1.0_8 - (abs(d%coords(3))/cellWidth)
                                   weight = xWeight*yWeight*zWeight
                                        
                                   outputArray(atom,1) = outputArray(atom,1) +  weight*chainMesh%fft_array_x(iCell,jCell,kCell)
                                   outputArray(atom,2) = outputArray(atom,2) +  weight*chainMesh%fft_array_y(iCell,jCell,kCell)
                                   outputArray(atom,3) = outputArray(atom,3) +  weight*chainMesh%fft_array_z(iCell,jCell,kCell)
                                   atom = chainMesh%atoms(atom)%nextAtom 
                               end do 
                        end do 
                        
                end do        

        end subroutine interpolate_fft_to_atoms 


        subroutine fft_forward_chainMesh(chainMesh)
                type(chainMesh_t), intent(inout) :: chainMesh 

                call fftw_execute_dft_r2c(chainMesh%forwardPlanX, chainMesh%fft_array_x, chainMesh%fft_c_view_x)
                call fftw_execute_dft_r2c(chainMesh%forwardPlanY, chainMesh%fft_array_y, chainMesh%fft_c_view_y)
                call fftw_execute_dft_r2c(chainMesh%forwardPlanZ, chainMesh%fft_array_z, chainMesh%fft_c_view_z)

        end subroutine fft_forward_chainMesh 
        
        subroutine fft_backwards_chainMesh(chainMesh)
                type(chainMesh_t), intent(inout) :: chainMesh 

                call fftw_execute_dft_c2r(chainMesh%backwardPlanX, chainMesh%fft_c_view_x, chainMesh%fft_array_x)
                call fftw_execute_dft_c2r(chainMesh%backwardPlanY, chainMesh%fft_c_view_y, chainMesh%fft_array_y)
                call fftw_execute_dft_c2r(chainMesh%backwardPlanZ, chainMesh%fft_c_view_z, chainMesh%fft_array_z)

        end subroutine fft_backwards_chainMesh


        subroutine calculate_demagnetisation_field(chainMesh, outputArray)
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), allocatable, dimension(:,:), intent(out) :: outputArray
                integer :: N,L,M, stat, i, j, k , waveIndexX, waveIndexY, waveIndexZ
                complex(kind=C_DOUBLE_COMPLEX) :: kx, ky, kz, displacement_phase
                real(kind=C_DOUBLE) :: scaleFactorX, scaleFactorY, scaleFactorZ, displacement_vector
                complex(kind = C_DOUBLE_COMPLEX) :: Mx, My, Mz, kdotM, k_squared

                integer :: startClock, endClock, clockRate
                real(kind = C_DOUBLE) :: elapsed_time
                call system_clock(startClock, clockRate)



                displacement_vector = chainMesh%latticeParameter / 2.0_8
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numcellsZ 

                if (allocated(outputArray)) then 
                        if (any(shape(outputArray) /= [chainMesh%numAtoms,3])) then 
                                deallocate(outputArray)
                                allocate(outputArray(chainMesh%numAtoms,3))
                        end if 
                end if 

                scaleFactorX = real(2.0_08,C_DOUBLE) * real(3.14159265358979323846_08, C_DOUBLE) / & !2 pi / N is
                                        real(N,C_DOUBLE)                               ! used to
                                                                                       ! calculate k values         
                                        

                scaleFactorY = real(2.0,C_DOUBLE) * real(3.14159265358979323846, C_DOUBLE) / &
                                        real(L,C_DOUBLE)

                scaleFactorZ = real(2.0_08,C_DOUBLE) * real(3.14159265358979323846, C_DOUBLE) / &
                                        real(M,C_DOUBLE)



                call interpolate_to_fft_array(chainMesh)
                call fft_forward_chainMesh(chainMesh)
                ! Now process data using the complex array view into the in place fft array 
                
                do i = 1, N/2 + 1
                        waveIndexX = i-1
                        kx = scaleFactorX*(waveIndexX)
                        do j = 1,L
                                waveIndexY = j-1
                                if (waveIndexY <= L/2)  then 
                                        ky = scaleFactorY*waveIndexY

                                else  
                                        ky = scaleFactorY*(waveIndexY-L)
                                end if 
                                do k = 1,M 
                                    waveIndexZ = k-1
                                    if (waveIndexZ <= M/2) then 
                                        kz = scaleFactorZ*waveIndexZ
                                    else  
                                        kz = scaleFactorZ*(waveIndexZ - M)
                                    end if  
                                    
                                    displacement_phase = exp(-cmplx(0.0_8,1.0_8) * (kx + ky + kz)*displacement_vector)
                                    Mx = chainMesh%fft_c_view_x(i,j,k)
                                    My = chainMesh%fft_c_view_y(i,j,k)
                                    Mz = chainMesh%fft_c_view_z(i,j,k)
                                    ! demagnetisation kernel is given by - (k.m / k^2) k
                                    kdotM = (kx*Mx + ky*My + kz*Mz)
                                    k_squared = kx*kx + ky*ky + kz*kz
                                    if ((i == 1 .and. j == 1 .and. k == 1)) cycle
                                    chainMesh%fft_c_view_x(i,j,k) = - displacement_phase*(kdotM / k_squared) * kx
                                    chainMesh%fft_c_view_y(i,j,k) = - displacement_phase*(kdotM / k_squared) * ky
                                    chainMesh%fft_c_view_z(i,j,k) = - displacement_phase*(kdotM / k_squared) * kz



                                end do 
                        end do 
                end do 
                chainMesh%fft_c_view_x(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                chainMesh%fft_c_view_y(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                chainMesh%fft_c_view_z(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                call fft_backwards_chainMesh(chainMesh) ! Demagnetising field is now in the three fft_arrays in chainMesh.
                chainMesh%fft_array_x = chainMesh%fft_array_x / (N*L*M)
                chainMesh%fft_array_y = chainMesh%fft_array_y / (N*L*M)
                chainMesh%fft_array_z = chainMesh%fft_array_z / (N*L*M)

                call interpolate_fft_to_atoms(chainMesh,outputArray)
                call system_clock(endClock, clockRate)
                elapsed_time = real(endClock - startClock, C_DOUBLE) / real(clockRate, C_DOUBLE)
        end subroutine calculate_demagnetisation_field

        subroutine calculate_magnetisation_gradient(chainMesh, outputArray)
                ! outputArray contains the gradient of M, with output(atomIndex,i,j) returning 
                ! partial_i M_j for the specified atom index
                ! This routine was not optimised to be called in a loop so do so at your own performance peril.
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), allocatable, dimension(:,:,:), intent(out) :: outputArray
                integer :: N,L,M, stat, i, j, k , waveIndexX, waveIndexY, waveIndexZ
                complex(kind=C_DOUBLE_COMPLEX) :: kx, ky, kz, displacement_phase
                real(kind=C_DOUBLE) :: scaleFactorX, scaleFactorY, scaleFactorZ, displacement_vector
                complex(kind = C_DOUBLE_COMPLEX) :: Mx, My, Mz, kdotM, k_squared

                integer :: startClock, endClock, clockRate, d, atomIndex, stat1, stat2, stat3  
                real(kind = C_DOUBLE) :: elapsed_time
                complex(kind = C_DOUBLE_COMPLEX), dimension(3) :: k_array
                real(kind=8), allocatable, dimension(:,:) :: tempArray
                complex(kind=C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:) :: temp_x, temp_y, temp_z
                call system_clock(startClock, clockRate)

                
                displacement_vector = chainMesh%latticeParameter / 2.0_8
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numcellsZ 
                allocate(temp_x(N/2 + 1, L, M), stat=stat1)
                allocate(temp_y(N/2 + 1, L, M), stat=stat2)
                allocate(temp_z(N/2 + 1, L, M), stat=stat3)

                if(any([stat1,stat2,stat3] /= 0)) error stop "Failed to allocate temporary reciprocal space buffers"

                if (allocated(outputArray)) then 
                        if (any(shape(outputArray) /= [chainMesh%numAtoms,3,3])) then 
                                deallocate(outputArray)
                                allocate(outputArray(chainMesh%numAtoms,3,3))
                        end if 
                else 
                        allocate(outputArray(chainMesh%numAtoms,3,3))
                end if 
                allocate(tempArray(chainMesh%numAtoms,3))
                scaleFactorX = real(2.0_08,C_DOUBLE) * real(3.14159265358979323846_08, C_DOUBLE) / & !2 pi / N is
                                        real(N,C_DOUBLE)                               ! used to
                                                                                       ! calculate k values         
                                        

                scaleFactorY = real(2.0,C_DOUBLE) * real(3.14159265358979323846, C_DOUBLE) / &
                                        real(L,C_DOUBLE)

                scaleFactorZ = real(2.0_08,C_DOUBLE) * real(3.14159265358979323846, C_DOUBLE) / &
                                        real(M,C_DOUBLE)


                call interpolate_to_fft_array(chainMesh) ! interpolates the spins to a grid 
                call fft_forward_chainMesh(chainMesh) ! FFT's these internal chainmesh grids for x,y,z 
                
                temp_x = chainMesh%fft_c_view_x
                temp_y = chainMesh%fft_c_view_y
                temp_z = chainMesh%fft_c_view_z
                do d = 1,3 

                        do i = 1, N/2 + 1
                                waveIndexX = i-1
                                kx = scaleFactorX*(waveIndexX)
                                do j = 1,L
                                        waveIndexY = j-1
                                        if (waveIndexY <= L/2)  then 
                                                ky = scaleFactorY*waveIndexY

                                        else  
                                                ky = scaleFactorY*(waveIndexY-L)
                                        end if 
                                        do k = 1,M 
                                            waveIndexZ = k-1
                                            if (waveIndexZ <= M/2) then 
                                                kz = scaleFactorZ*waveIndexZ
                                            else  
                                                kz = scaleFactorZ*(waveIndexZ - M)
                                            end if  
                                            
                                            k_array = [kx,ky,kz]
                                            displacement_phase = exp(-cmplx(0.0_8,1.0_8) * (kx + ky + kz)*displacement_vector)
                                            Mx = temp_x(i,j,k)
                                            My = temp_y(i,j,k)
                                            Mz = temp_z(i,j,k)
                                            ! gradient is given by i k_d M_j where d and j are indices
                                            chainMesh%fft_c_view_x(i,j,k) = cmplx(0.0,1.0,C_DOUBLE_COMPLEX)*&
                                                                        displacement_phase*k_array(d)*Mx
                                            chainMesh%fft_c_view_y(i,j,k) = cmplx(0.0,1.0,C_DOUBLE_COMPLEX)*&
                                                                        displacement_phase*k_array(d)*My 
                                            chainMesh%fft_c_view_z(i,j,k) = cmplx(0.0,1.0,C_DOUBLE_COMPLEX)*&
                                                                        displacement_phase*k_array(d)*Mz



                                        end do 
                                end do 
                        end do 
                        
                        call fft_backwards_chainMesh(chainMesh)
                        chainMesh%fft_array_x = chainMesh%fft_array_x / (N*L*M)
                        chainMesh%fft_array_y = chainMesh%fft_array_y / (N*L*M)
                        chainMesh%fft_array_z = chainMesh%fft_array_z / (N*L*M)
                        call interpolate_fft_to_atoms(chainMesh,tempArray)
                        
                        outputArray(:,:,d) = tempArray(:,:)
                        
                        
                        

                end do 
                

                call system_clock(endClock, clockRate)
                elapsed_time = real(endClock - startClock, C_DOUBLE) / real(clockRate, C_DOUBLE)
        end subroutine calculate_magnetisation_gradient 

        function calculate_winding_number(chainMesh,Z_index) result(winding_number)
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: Z_index
                real(kind=8) :: winding_number 

                real(kind=8), allocatable, dimension(:,:,:) :: grad_array 
                real(kind = 8) :: acc 
                integer :: numCellsX, numCellsY, numCellsZ, iCell, jCell, i,j, atomIndex, chainCellIndex
                type(vecNd_t) :: spin, spin_grad_x, spin_grad_y  
                call calculate_magnetisation_gradient(chainMesh,grad_array)

                ! For now we will just calculate it at Z = numCellsZ / 2 
                numCellsZ = chainMesh%numCellsZ 
                numCellsX = chainMesh%numCellsX 
                numCellsY = chainMesh%numCellsY
                if (Z_index < 1 .or. Z_index > chainMesh%numCellsZ) error stop "Z index out of range"                
                acc = 0.0_08
                do i = 1, numCellsX 
                        do j = 1,numCellsY 
                                chainCellIndex = IndexFromCoordinates(chainMesh,i,j,Z_index)
                                atomIndex = chainMesh%chainMeshCells(chainCellIndex)%firstAtomInMeshCell
                                spin = makeVecNdCheck(spin,dble(chainMesh%atoms(atomIndex)%AtomParameters))
                                spin_grad_x = makeVecNdCheck(spin_grad_x,[grad_array(atomIndex,1,1), & ! x derivative 
                                                                                grad_array(atomIndex,2,1), & 
                                                                                        grad_array(atomIndex,3,1)])
                                spin_grad_y = makeVecNdCheck(spin_grad_y,[grad_array(atomIndex,1,2), & ! y derivative 
                                                                                grad_array(atomIndex,2,2), &
                                                                                        grad_array(atomIndex,3,2)])

                                acc = acc + spin * (spin_grad_x .x. spin_grad_y)
                        end do 
                end do 
                acc = acc * chainMesh%latticeParameter * chainMesh%latticeParameter ! dx*dy 
                acc = acc / (4*pi)
                winding_number = acc
        end function calculate_winding_number

        function arc_winding(s1,s2,s3) result(sigma_area)
                implicit none
                type(vecNd_t), intent(in) :: s1, s2, s3
                real(kind=8) :: sigma_area
                complex(kind=8) :: num

                
                num = 1 + s1*s2 + s2*s3 + s3*s1 + cmplx(0.0_8,s1*(s2 .x. s3))
                if (abs(real(num)) < 1e-15 .and. abs(aimag(num)) < 1e-15) then 
                        sigma_area = 0.0_8 
                        return 
                end if 
                sigma_area = 2*atan2(aimag(num),real(num))
        end function arc_winding

        function calculate_winding_number2(chainMesh,Z_index) result(winding_number)
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: Z_index
                
                integer :: i, j, i_index, j_index, N, L, M, atomIndex, cellIndex 
                integer :: atom2, atom3, atom4
                type(VecNd_t) :: s1, s2, s3, s4
                real(kind=8) :: x1,x2_1,x2_2,x3,y1,y2_1,y2_2,y3,z1,z2_1,z2_2,z3
                real(kind=8) :: sigma1_area1, sigma2_area2, winding_number

                N = chainMesh%numCellsX 
                L = chainMesh%numcellsY 
                M = chainMesh%numCellsZ
                if (Z_index < 1 .or. Z_index > M) error stop "Z_index out of bounds"
                if (.not. allocated(chainMesh%derivativeList)) error stop "Derivative List not initialised"
                if (any(shape(chainMesh%derivativeList) /= [chainMesh%numAtoms,3,2])) &
                                                error stop "DerivativeList not initialised with proper shape"


                winding_number = 0.0_8
                do i = 1, chainMesh%numCellsX
                        do j = 1, chainMesh%numCellsY
                                cellIndex = IndexFromCoordinates(chainMesh,i,j,Z_index)
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                atom2 = chainMesh%derivativeList(atomIndex,1,2)
                                atom4 = chainMesh%derivativeList(atom2,2,2)
                                atom3 = chainMesh%derivativeList(atomIndex,2,2)
                               s1 = dble(chainMesh%atoms(atomIndex)%AtomParameters)
                               s2 = dble(chainMesh%atoms(atom2)%AtomParameters)
                               s3 = dble(chainMesh%atoms(atom3)%AtomParameters)
                               s4 = dble(chainMesh%atoms(atom4)%AtomParameters)

                               s1 = s1 / abs(s1)
                               s2 = s2 / abs(s2)
                               s3 = s3 / abs(s3) 
                               s4 = s4 / abs(s4)
                               sigma1_area1 = arc_winding(s1,s2,s4)
                               sigma2_area2 = arc_winding(s1,s4,s3)
                               
                               winding_number = winding_number + sigma1_area1 + sigma2_area2 
                        end do 
                end do 
                winding_number = winding_number / (4*pi) ! Leaving this at the end hoping the compiler will vectorise it
        end function calculate_winding_number2


        subroutine calculate_winding_number_density(chainMesh,Z_index, density_matrix)
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: Z_index
                real(kind=8), dimension(:,:), allocatable, intent(out) :: density_matrix 
                
                integer :: i, j, i_index, j_index, N, L, M, atomIndex, cellIndex 
                integer :: atom2, atom3, atom4
                type(VecNd_t) :: s1, s2, s3, s4
                real(kind=8) :: x1,x2_1,x2_2,x3,y1,y2_1,y2_2,y3,z1,z2_1,z2_2,z3
                real(kind=8) :: sigma1_area1, sigma2_area2

                N = chainMesh%numCellsX 
                L = chainMesh%numcellsY 
                M = chainMesh%numCellsZ
                if (Z_index < 1 .or. Z_index > M) error stop "Z_index out of bounds"
                if (.not. allocated(chainMesh%derivativeList)) error stop "Derivative List not initialised"
                if (any(shape(chainMesh%derivativeList) /= [chainMesh%numAtoms,3,2])) &
                                                error stop "DerivativeList not initialised with proper shape"
                if (allocated(density_matrix)) then 
                        if (any(shape(density_matrix) /= [N,L])) then 
                                deallocate(density_matrix)
                                allocate(density_matrix(N,L))
                        end if 
                else 
                        allocate(density_matrix(N,L))
                end if
                density_matrix(:,:) = 0.0_8


                do i = 1, chainMesh%numCellsX
                        do j = 1, chainMesh%numCellsY
                                cellIndex = IndexFromCoordinates(chainMesh,i,j,Z_index)
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                atom2 = chainMesh%derivativeList(atomIndex,1,2)
                                atom4 = chainMesh%derivativeList(atom2,2,2)
                                atom3 = chainMesh%derivativeList(atomIndex,2,2)
                               s1 = dble(chainMesh%atoms(atomIndex)%AtomParameters)
                               s2 = dble(chainMesh%atoms(atom2)%AtomParameters)
                               s3 = dble(chainMesh%atoms(atom3)%AtomParameters)
                               s4 = dble(chainMesh%atoms(atom4)%AtomParameters)

                               s1 = s1 / abs(s1)
                               s2 = s2 / abs(s2)
                               s3 = s3 / abs(s3) 
                               s4 = s4 / abs(s4)
                               sigma1_area1 = arc_winding(s1,s2,s4)
                               sigma2_area2 = arc_winding(s1,s4,s3)
                               
                               density_matrix(i,j) = sigma1_area1 + sigma2_area2 
                        end do 
                end do 
                density_matrix = density_matrix / (4*pi) ! Leaving this at the end hoping the compiler will vectorise it
        end subroutine calculate_winding_number_density


        subroutine add_neighbors_to_stack(chainMesh, i,j,stack_array, stack_ptr, visited_array, density_mask, in_stack_array, &
                                                density_matrix)
                implicit none
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: i,j
                integer, dimension(:,:), intent(inout) :: stack_array
                integer, intent(inout) :: stack_ptr 
                logical, dimension(:,:), intent(in) :: visited_array, density_mask 
                logical, dimension(:,:), intent(inout) :: in_stack_array
                real(kind=8), dimension(:,:), intent(in) :: density_matrix
                integer, dimension(2) :: stack_array_shape
                integer :: i_neighbor, j_neighbor, itemp, jtemp, N, L 
                logical :: downhill
                
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY
                stack_array_shape = shape(stack_array)
                  do i_neighbor = -1,1 
                        do j_neighbor = -1,1 
                                itemp = i + i_neighbor
                                jtemp = j + j_neighbor
                                if (itemp < 1) itemp = N ! Periodic boundaries 
                                if (itemp > N) itemp = 1 
                                
                                if (jtemp < 1) jtemp = L 
                                if (jtemp > L) jtemp = 1 
                                if (in_stack_array(itemp,jtemp)) cycle 
                                downhill = density_matrix(i,j) >= density_matrix(itemp,jtemp) - 1e-10_8
                                if ((.not. visited_array(itemp,jtemp)) .and. density_mask(itemp,jtemp) .and. &
                                                                         downhill) then 
                                        if (stack_ptr + 1 > stack_array_shape(1)) error stop "Stack array capacity exceeded"
                                        stack_array(stack_ptr,1) = itemp ! stack_ptr points to where the next entry
                                        ! should be written and equals the existing number of entries + 1.
                                        stack_array(stack_ptr,2) = jtemp 
                                        stack_ptr = stack_ptr + 1
                                        in_stack_array(itemp,jtemp) = .True.
                                        
                                end if 
                        end do 
                  end do
        end subroutine add_neighbors_to_stack

        function calculate_skyrmion_number(chainMesh,Z_index,q_threshold,particle_number, sigma) result(skyrmion_number)
                implicit none
                type(chainMesh_t), intent(inout) :: chainmesh
                integer, intent(in) :: Z_index, particle_number
                real(kind=8), intent(in) :: q_threshold, sigma
                
                real(kind=8), dimension(:,:), allocatable :: density_matrix
                logical, dimension(:,:), allocatable :: density_mask, visited_array, in_stack_array
                integer :: N,L, stat, skyrmion_number
                
                integer, allocatable, dimension(:,:) :: stack_array ! stack_array(stack_ptr, (i,j))
                integer :: stack_ptr
                integer, parameter :: stack_len = 10000

                integer :: i, j, i_neighbor, j_neighbor, itemp, jtemp, candidate_counter, Nx, Ny, &
                                        xIndex, yIndex
                real(kind=8) :: acc
                logical :: is_candidate, is_local_maxima
                real(kind=8) :: upper_threshold
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 

                allocate(density_matrix(N,L), stat=stat)
                if (stat /= 0) error stop "Error allocating topological charge density matrix"

                allocate(density_mask(N,L), stat=stat)
                if (stat/=0) error stop "Error allocating density mask"

                allocate(visited_array(N,L), stat=stat)
                if (stat/=0) error stop "Error allocating visited_array"

                allocate(in_stack_array(N,L),stat=stat)
                if (stat/=0) error stop "Error allocating in_stack_array"
                
                allocate(stack_array(stack_len,2), stat=stat)
                if (stat /= 0) error stop "Error allocating stack array"
                stack_ptr = 1 

                if (abs(q_threshold) > 1.0_8 .or. q_threshold < 0.0_8) error stop "q_threshold must be between 0 and 1"

                call calculate_winding_number_density(chainMesh, Z_index, density_matrix)
                !print *, "Minval / Maxval in winding array = ", minval(abs(density_matrix)) / maxval(abs(density_matrix))
                ! density_mask(:,:) = density_matrix > (q_threshold * maxval(abs(density_matrix)))
                call Gaussian_filter_2d(density_matrix,sigma)
                density_mask = abs(density_matrix) > (q_threshold * maxval(abs(density_matrix)))
                visited_array(:,:) = .False.
                in_stack_Array(:,:) = .False.
                stack_array(:,:) = 0 
                skyrmion_number = 0 
                candidate_counter = 0
                is_candidate = .False.
                Nx = chainMesh%numcellsX 
                Ny = chainMesh%numCellsY
                do i = 1,N 
                        do j = 1,L 
                                is_local_maxima = .True.
                                do itemp = -1,1
                                        do jtemp = -1,1
                                                xIndex = modulo(i + itemp - 1,Nx) + 1
                                                yIndex = modulo(j + jtemp - 1,Ny) + 1
                                                if (itemp == 0 .and. jtemp == 0) cycle
                                                if (density_matrix(i,j) < density_matrix(xIndex,yIndex)) is_local_maxima = .False.
                                        end do 
                                end do 
                                if (density_mask(i,j) .and. (.not. visited_array(i,j)) .and. is_local_maxima) then 
                                  acc = density_matrix(i,j)
                                  visited_array(i,j) = .True.
 
                                  ! Loop through neighbors, if they are above the theshold and have not been visited or put in the
                                  ! stack then add them to the stack

                                
 
                                  stack_ptr = 1
                                  in_stack_array = .False.
                                  
                                  call add_neighbors_to_stack(chainMesh,i,j,stack_array,stack_ptr, &
                                                 visited_array,density_mask, in_stack_array, density_matrix)
                                  
                                  ! Now for each neighbor in the stack, visit it if not already visited and add it's density to the
                                  ! accumulator. Then add it's unvisited neighbors over the density threshold to the stack, repeat
                                  ! until stack_ptr - 1 == 0
                                  do while (stack_ptr - 1 /= 0)
                                        ! pop cell off of the stack 
                                        itemp = stack_array(stack_ptr - 1,1)
                                        jtemp = stack_array(stack_ptr - 1,2)
                                        stack_ptr = stack_ptr - 1
                                        !print *, "itemp, jtemp = ", itemp, jtemp
                                        if (.not. visited_array(itemp,jtemp)) then

                                                acc = acc + density_matrix(itemp,jtemp)
                                                visited_array(itemp,jtemp) = .True.
                                                call add_neighbors_to_stack(chainMesh,itemp,jtemp,& 
                                                        stack_array,stack_ptr,visited_array,density_mask, &
                                                                in_stack_array, density_matrix)
                                        end if 
                                  end do
                                !   print *, "acc, skyrmion_number = ", acc, skyrmion_number, "q_theshold = ", q_threshold, &
                                !                candidate_counter
                                  if (abs(abs(acc) - particle_number) < 0.2) skyrmion_number = skyrmion_number + 1
                                end if 
                                stack_array = 0
                        end do
                end do 
                call write_2d_real_array_to_file(density_matrix, "./density_matrix.csv")
                call write_2d_logical_array_to_file(visited_array,"./visited_array.csv")
                call write_2d_logical_array_to_file(density_mask,"./density_mask.csv")
        end function calculate_skyrmion_number

        subroutine Gaussian_filter_2d(array,sigma)
                real(kind=8), dimension(:,:), intent(inout) :: array 
                real(kind=8), intent(in) :: sigma
                real(kind=8), dimension(:,:), allocatable :: buffer
                real(kind = 8), dimension(:,:), allocatable :: G 
                integer, parameter :: k = 2 ! k functions as a cutoff
                real(kind = 8) :: acc
                integer :: stat, i, Nx, Ny, index_i, index_j, l,m
                acc = 0.0_8 
                allocate(G(2*k + 1,2*k+1),stat = stat)
                if (stat /= 0) error stop "Error allocating G array"

                Nx = size(array,1)
                Ny = size(array,2)

                allocate(buffer(Nx, Ny), stat = stat)
                if (stat /= 0) error stop "Error allocating Buffer"
                do i = 1, 2*k + 1
                        do j = 1, 2*k + 1 
                                G(i,j) = exp(- ((i - k - 1)**2 +(j - k - 1)**2)/(2*sigma**2))
                                acc = acc + G(i,j)
                        end do 
                end do 
                G = G / acc
                buffer = 0.0_8
                do i = 1, Nx
                        do j = 1, Ny 
                                do l = 1, 2*k + 1 
                                        do m = 1, 2*k + 1
                                        ! Going to compute array \conv G 
                                        ! This is given mathematically (for a 1d array) as (array conv G)[i] = sum_n array[i-n] * G[n]
                                        ! We just need to handle periodic boundaries for i - n that are out of bounds 
                                        index_i = i - l 
                                        index_j = j - m 
                                        if (index_i > Nx) index_i = index_i - Nx 
                                        if (index_i < 1) index_i = Nx - abs(index_i)

                                        if (index_j > Ny) index_j = index_j - Ny 
                                        if (index_j < 1) index_j = Ny - abs(index_j)
                                        buffer(i,j) = buffer(i,j) + array(index_i, index_j)*G(l,m)
                                        end do 
                                end do 
                        end do 
                end do 
                array = buffer
        end subroutine Gaussian_filter_2d

        subroutine write_2d_real_array_to_file(array, filename)
                implicit none
                real(kind=8), intent(in), dimension(:,:) :: array 
                character(len=*), intent(in) :: filename 
                
                integer :: N, L, i, j 
                integer, dimension(2) :: shape_array
                character(len=100) :: real_buffer
                character(len=:), allocatable :: buffer
                shape_array = shape(array)
                
                N = shape_array(1)
                L = shape_array(2)
                open(unit=231,file=trim(filename),status="replace",action="write")                 
                allocate(character(len=100) :: buffer)
                real_buffer = ' '
                buffer = ' '


                do i = 1,N 
                        buffer = ' '
                        do j = 1,L
                                real_buffer = ' '
                                write(real_buffer,'(ES12.5)') array(i,j) 
                                buffer = buffer // trim(adjustl(real_buffer))
                                if (j /= L) buffer = buffer // ','
                        end do 
                        write(231,'(A)') buffer
                end do 
                close(231)
        end subroutine write_2d_real_array_to_file


        subroutine write_2d_logical_array_to_file(array, filename)
                implicit none
                logical, intent(in), dimension(:,:) :: array 
                character(len=*), intent(in) :: filename 
                
                integer :: N, L, i, j 
                integer, dimension(2) :: shape_array
                character(len=100) :: real_buffer
                character(len=:), allocatable :: buffer
                shape_array = shape(array)
                
                N = shape_array(1)
                L = shape_array(2)
                
                open(unit=231,file=trim(filename),status="replace",action="write")                 
                allocate(character(len=100) :: buffer)
                real_buffer = ' '
                buffer = ' '


                do i = 1,N 
                        buffer = ' '
                        do j = 1,L 
                                real_buffer = ' '
                                write(real_buffer,'(I5.5)') merge(1,0,array(i,j)) 
                                buffer = buffer // trim(adjustl(real_buffer))
                                if (j /= L) buffer = buffer // ','
                        end do 
                        write(231,'(A)') buffer
                end do 
                close(231)
        end subroutine write_2d_logical_array_to_file




        subroutine compute_skyrmion_distribution(chainMesh,N,winding_array,min_threshold,max_threshold,num_thresholds, Z_index)
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: N ! Will compute the number of particles with winding numbers 1...N
                integer, dimension(:), allocatable, intent(inout) :: winding_array
                real(kind=8), intent(in) :: min_threshold, max_threshold 
                integer,intent(in) :: num_thresholds, Z_index

                integer :: i, j, sigma_index
                real(kind=8) :: total_charge, threshold, winding, sigma
                
                if (max_threshold < min_threshold) error stop "max_threshold must be greater than min_threshold"
                if (allocated(winding_array)) then 
                        if (size(winding_array) /= N) then 
                                deallocate(winding_array)
                                allocate(winding_array(N))
                        end if 

                else 
                        allocate(winding_array(N))
                end if 

                total_charge = calculate_winding_number2(chainMesh, Z_index)
                winding_array = 0.0_8
                do i = 1,num_thresholds
                        do sigma_index = 1,20
                        ! sigma = (dble(sigma_index)/20.0_8) * (1.0_8 - 0.001_8) + 0.001_8 
                        sigma = 1.0_8 - (dble(sigma_index)/20.0_8)*(1.0_8 - 0.001_8)
                        winding = 0.0_8
                        !threshold = (dble(i) / dble(num_thresholds)) * (max_threshold - min_threshold) + min_threshold
                        threshold = max_threshold - (dble(i) / dble(num_thresholds)) * (max_threshold - min_threshold)
                        do j = 1, N 
                                winding_array(j) = calculate_skyrmion_number(chainMesh,Z_index,threshold,j,sigma)
                                winding = winding + j * winding_array(j)
                        end do 
                        !print *, "Winding array = ", winding_array, "threshold = ", threshold, "sigma = ", sigma
                        if (abs(abs(winding) - abs(total_charge)) < 0.1) then 
                                print *, "Solution found at q_threshold = ", threshold, "threshold = ", threshold, "sigma = ", sigma
                                return
                        end if 
                        end do
                end do 
                winding_array = 0 
        end subroutine compute_skyrmion_distribution
end module reciprocal_space_processes 
