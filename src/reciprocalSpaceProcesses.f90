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
                outputArray = 0.0_8
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
end module reciprocal_space_processes 
