module algo 


        interface quicksort 
                module procedure quicksort_dp 
        end interface quicksort 

        public :: quicksort 

contains 
        recursive subroutine quicksort_dp(array)
                use, intrinsic :: iso_fortran_env, dp=>real64 
                implicit none 

                real(kind=dp), dimension(:), intent(inout) :: array 

                integer :: lesserPartitionIndex, GreaterPartitionIndex, &
                        pivotIndex, pivot  
                real(kind=dp) :: temp

                if (size(array) == 1) return 
                if (size(array) == 2) then 
                        if (array(1) <= array(2)) return 
                end if 
                ! Choose pivot to be in the middle of the array 
                pivotIndex = size(array)/2 
                
                lesserPartitionIndex = 1
                GreaterPartitionIndex = size(array) 

                pivot = array(pivotIndex)

                do while (lesserPartitionIndex < GreaterPartitionIndex)
                        
                        if (array(lesserPartitionIndex) > pivot) then 
                                temp = array(lesserPartitionIndex)
                                array(lesserPartitionIndex) = pivot 
                                array(pivotIndex) = temp 
                        else if (array(GreaterPartitionIndex) < pivot) then 
                                temp = array(GreaterPartitionIndex)
                                array(GreaterPartitionIndex) = pivot 
                                array(pivotIndex) = temp 
                        end if 
                
                        lesserPartitionIndex = lesserPartitionIndex + 1
                        GreaterPartitionIndex = GreaterPartitionIndex - 1
                end do 

                print *, "pivot = ", pivot 
                print *, "Lower partition = ", array(1:GreaterPartitionIndex)
                print *, "Upper Partition = ", array(GreaterPartitionIndex:size(array))
                
                call quicksort_dp(array(1:GreaterPartitionIndex))
                call quicksort_dp(array(GreaterPartitionIndex:size(array))) 
        end subroutine quicksort_dp 
end module algo 
