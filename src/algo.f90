module algo 


        interface quicksort 
                module procedure quicksort_dp 
        end interface quicksort 

        public :: quicksort 

contains 
        recursive subroutine quicksort_dp(array,integer_companion) 
                ! Integer Companion exists so that it's entries are kept in line with the sorted array entries 
                use, intrinsic :: iso_fortran_env, dp=>real64 
                implicit none 

                real(kind=dp), dimension(:), intent(inout) :: array 
                integer, dimension(:), optional, intent(inout) :: integer_companion
                integer :: lesserPartitionIndex, GreaterPartitionIndex, &
                        pivotIndex, pivot  
                real(kind=dp) :: temp

                if (present(integer_companion)) then 
                        if (size(integer_companion) /= size(array)) error stop "Error: &
                                integer_companion must be the same size as array"
                end if 
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
                                if (present(integer_companion)) then 
                                        temp = integer_companion(lesserPartitionIndex)
                                        integer_companion(lesserPartitionIndex) = integer_companion(pivotIndex) 
                                        integer_companion(pivotIndex) = temp 
                                end if 
                                pivotIndex = lesserPartitionIndex 
                        end if 
                        if (array(GreaterPartitionIndex) < pivot) then 
                                temp = array(GreaterPartitionIndex)
                                array(GreaterPartitionIndex) = pivot 
                                array(pivotIndex) = temp 
                                if (present(integer_companion)) then 
                                        temp = integer_companion(GreaterPartitionIndex)
                                        integer_companion(GreaterPartitionIndex) = integer_companion(pivotIndex) 
                                        integer_companion(pivotIndex) = temp 
                                end if 
                                pivotIndex = GreaterPartitionIndex 
                        end if 
                
                        lesserPartitionIndex = lesserPartitionIndex + 1
                        GreaterPartitionIndex = GreaterPartitionIndex - 1
                end do 

                if (.not. present(integer_companion)) then  
                        call quicksort_dp(array(1:GreaterPartitionIndex))
                        call quicksort_dp(array(GreaterPartitionIndex + 1:size(array))) 
                else 
                        call quicksort_dp(array(1:GreaterPartitionIndex),&
                                integer_companion=integer_companion(1:GreaterPartitionIndex))
                        call quicksort_dp(array(GreaterPartitionIndex + 1:size(array)),&
                                integer_companion=integer_companion(GreaterPartitionIndex + 1:size(integer_companion))) 
                end if 
        end subroutine quicksort_dp 
end module algo 
