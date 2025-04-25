module reciprocal_space_processes 
        use chainMesh 
        use vecNd
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
                if (any(chainMesh%fft_c_view_x /= chainMesh%fft_c_view_x)) then
                        print *, chainMesh%fft_c_view_x 
                        error stop "c_view_x contains NaN's"
                end if 
                if (any(chainMesh%fft_c_view_y /= chainMesh%fft_c_view_y)) then
                        print *, chainMesh%fft_c_view_y 
                        error stop "c_view_y contains NaN's"
                end if 
                if (any(chainMesh%fft_c_view_z /= chainMesh%fft_c_view_z)) then
                        print *, chainMesh%fft_c_view_z 
                        error stop "c_view_z contains NaN's"
                end if 
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
                print *, "Computed demag field in ", elapsed_time, "seconds"
        end subroutine calculate_demagnetisation_field
end module reciprocal_space_processes 
