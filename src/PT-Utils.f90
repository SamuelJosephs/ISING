module PT_Utils 
        use chainMesh
        use reciprocal_space_processes
        use constants
        implicit none 
        
        public :: indicesFromSlot, SlotFromIndices
        contains 

                subroutine indicesFromSlot(Slot,NJ,ND,NB,Jindex,Dindex,BIndex)
                        integer, intent(in) :: Slot,NJ,ND,NB ! Number of values for each parameter, 1 based indexing
                        integer, intent(out) :: Jindex, Dindex, Bindex
                        
                        integer :: res 
                        ! index = (Jindex-1)*ND*NB + (DIndex - 1)*NB + (B_Index - 1) 

                        Jindex = (Slot-1) / (ND*NB) + 1 
                        res = mod(Slot-1,ND*NB)
                        DIndex = res / NB + 1 
                        Bindex = mod(res,NB) + 1 
                end subroutine indicesFromSlot

                function SlotFromIndices(NJ,ND,NB,JIndex,DIndex,BIndex) result(res)
                        integer, intent(in) :: NJ, ND, NB, JIndex, DIndex, BIndex
                        integer :: res ! Index in array for a given set of indices, 1 based indexing

                        res = (JIndex-1)*ND*NB + (DIndex-1)*NB + BIndex
                        ! Don't need (Bindex - 1) + 1 as the ones cancel
                end function SlotFromIndices


                subroutine chainMesh_statistics(chainMesh,skyrmion_number_middle,winding_number_middle,winding_number_spread,&
                                                                magnetisation)
                        type(ChainMesh_t), intent(inout) :: chainMesh 
                        real(kind=dp), intent(out) :: winding_number_middle, winding_number_spread
                        integer, intent(out) :: skyrmion_number_middle
                        real(kind=dp), allocatable, dimension(:) :: winding_array 
                        type(vecNd_t), intent(out) :: magnetisation
                        integer, allocatable, dimension(:) :: skyrmion_array
                        
                        integer :: i
                        real(kind=dp), parameter :: lower_bound = 1e-5_dp
                        real(kind=dp), parameter :: upper_bound = 0.95_dp
                        integer, parameter :: num_thresholds = 60
                        type(vecNd_t) :: temp
                        allocate(winding_array(chainMesh%numCellsZ))

                        do i = 1,chainmesh%numCellsZ 
                                winding_array(i) = calculate_winding_number2(chainMesh,i)
                        end do 

                        winding_number_spread = abs(maxval(winding_array) - minval(winding_array))
                        call compute_skyrmion_distribution(chainMesh,1,skyrmion_array, &
                                                        lower_bound, upper_bound,num_thresholds,chainMesh%numCellsZ/2)
                        skyrmion_number_middle = skyrmion_array(1)
                        winding_number_middle = winding_array(chainMesh%numCellsZ/2)

                        magnetisation = [0.0_dp, 0.0_dp, 0.0_dp]
                        
                        do i = 1,chainMesh%numAtoms 
                               temp = chainMesh%atomSpins(i,:)
                               magnetisation = magnetisation + temp
                        end do 
                        magnetisation = magnetisation / dble(chainMesh%numAtoms)
                end subroutine chainMesh_statistics 
end module PT_Utils 
