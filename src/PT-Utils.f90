module PT_Utils 
        implicit none 
        
        public :: indicesFromSlot, SlotFromIndices
        contains 

                subroutine indicesFromSlot(Slot,NJ,ND,NB,Jindex,Dindex,BIndex)
                        integer, intent(in) :: Slot,NJ,ND,NB ! Number of values for each parameter, 1 based indexing
                        integer, intent(out) :: Jindex, Dindex, Bindex
                        
                        integer :: res 
                        ! index = (Jindex-1)*ND*NB + (DIndex - 1)*NB + (B_Index - 1) 

                        Jindex = Slot / (ND*NB) + 1 
                        res = mod(Slot,ND*NB)
                        DIndex = res / NB + 1 
                        Bindex = mod(res,NB) + 1 
                end subroutine indicesFromSlot

                function SlotFromIndices(NJ,ND,NB,JIndex,DIndex,BIndex) result(res)
                        integer, intent(in) :: NJ, ND, NB, JIndex, DIndex, BIndex
                        integer :: res ! Index in array for a given set of indices, 1 based indexing

                        res = (JIndex-1)*ND*NB + (DIndex-1)*NB + BIndex
                        ! Don't need (Bindex - 1) + 1 as the ones cancel
                end function SlotFromIndices

end module PT_Utils 
