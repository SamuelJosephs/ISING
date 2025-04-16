module reciprocal_space_processes 
        use chainMesh 
        use, intrinsic :: iso_c_binding 
        !include 'fftw3.f03' !This is included in chainMesh, annoying I would like header files with #pragma once
        
        contains
        
        subroutine compute_chainMesh_gradient(chainMesh,outputArray, plan) 
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), intent(inout), allocatable :: outputArray(:,:) ! (atomIndex, component) 
                type(C_ptr), intent(in) :: plan ! plan for 3d fftw fourier transform
                integer :: atom, iCell, jCell, kCell, chainCellIndex, atomIndex, N,L,M , stat, d  
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 

                if (.not. allocated(chainMesh%derivativeList)) error stop "Derivative List not allocated"
                if (.not. all([N,L,M] == shape(chainMesh%DerivativeList))) error stop "Derivative List is the wrong shape"
                if (stat /= 0) error stop "Allocation Failed"
                do d = 1,3 ! dimension = 1,3 
                        do atom = 1,chainMesh%chainMeshCells(1)%NumAtomsPerUnitCell
                                do iCell = 1,N 
                                        do jCell = 1, L 
                                                do kCell = 1, M 
                                                        chainCellIndex = &
                                                                IndexFromCoordinates(chainMesh,iCell,jCell,kCell)
                                                        atomIndex = chainMesh%chainMeshCells(chainCellIndex)%firstAtomInMeshCell
                                                        do i = 1,atom - 1 ! Traverse linked list to get atom  
                                                                atomIndex = chainMesh%atoms(atomIndex)%nextAtom 
                                                        end do 
                                                        chainMesh%fft_array(iCell,jCell,kCell) = &
                                                                        real(chainMesh%atoms(atomIndex)%AtomParameters(d),kind=c_double) 
                                                end do 
                                        end do
                                end do

                        end do 
                end do 
        end subroutine compute_chainMesh_gradient

end module reciprocal_space_processes 
