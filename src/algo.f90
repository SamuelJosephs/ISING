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
                        pivotIndex  
                real(kind=dp) :: temp, pivot
                integer :: temp_int 
                
                if (size(array) == 1) return 


                print *, "entering quicksort with array :", array
                ! Choose pivot 

                pivotIndex = size(array) / 2
                pivot = array(pivotIndex)

                ! Now partition into elements below element and above element 

                lesserPartitionIndex = 1
                GreaterPartitionIndex = size(array)

                do while (lesserPartitionIndex < GreaterPartitionIndex) 
                        ! Move lesser index until we find element that shouldn't be in the lesser partition 
                        ! Same for the right hand partition
                        if (array(lesserPartitionIndex) < pivot) lesserPartitionIndex = lesserPartitionIndex + 1 
                        if (array(GreaterPartitionIndex) > pivot) GreaterPartitionIndex = GreaterPartitionIndex - 1 

                        ! Swap if suitable 

                        if ((array(lesserPartitionIndex) >= pivot) .and. (array(GreaterPartitionIndex) <= pivot)) then 
                                if (lesserPartitionIndex == pivotIndex) then 
                                        pivotIndex = GreaterPartitionIndex
                                else if (GreaterPartitionIndex == pivotIndex) then 
                                        pivotIndex = lesserPartitionIndex 
                                end if 
                                call swap_dp_array_elements(array,lesserPartitionIndex,GreaterPartitionIndex)
                                if (present(integer_companion)) call swap_integer_array_elements(integer_companion,&
                                                                        lesserPartitionIndex, GreaterPartitionIndex)
                        end if 

                end do 
                if (pivotIndex == size(array)) pivotIndex = pivotIndex - 1
                
                if (.not. present(integer_companion)) then 
                        call quicksort_dp(array(1:pivotIndex))
                        call quicksort_dp(array(pivotIndex + 1:size(array)))
                else 
                        call quicksort_dp(array(1:pivotIndex),integer_companion=integer_companion(1:pivotIndex))
                        call quicksort_dp(&
                                array(pivotIndex + 1:size(array)),&
                                integer_companion=integer_companion(pivotIndex + 1:size(array)))
                end if 
        end subroutine quicksort_dp 

        subroutine swap_dp_array_elements(array,i,j)
                use iso_fortran_env, only: dp=>real64 
                implicit none 
                real(kind=dp), dimension(:), intent(inout) :: array 
                integer, intent(in) :: i,j
                real(kind=dp) :: temp 

                temp = array(i) 
                array(i) = array(j) 
                array(j) = temp 
                return 
        end subroutine swap_dp_array_elements 

        subroutine swap_integer_array_elements(array,i,j)
                implicit none 
                integer, dimension(:), intent(inout) :: array 
                integer, intent(in) :: i,j
                integer :: temp 

                temp = array(i) 
                array(i) = array(j) 
                array(j) = temp 
                return 
        end subroutine swap_integer_array_elements 
end module algo 
