module algo 


        interface quicksort 
                module procedure quicksort_dp 
        end interface quicksort 

        interface mergesort 
                module procedure mergesort_dp 
        end interface mergesort 

        public :: quicksort, mergesort 

contains 
        recursive subroutine quicksort_dp(array,integer_companion) 
                ! WARNING: Seems very tempermental and doesn't work with duplicates, only use the mergesort implementation please
                ! future self.

                ! Integer Companion exists so that it's entries are kept in line with the sorted array entries 
                use, intrinsic :: iso_fortran_env, dp=>real64 
                implicit none 

                real(kind=dp), dimension(:), intent(inout) :: array 
                integer, dimension(:), optional, intent(inout) :: integer_companion
                integer :: lesserPartitionIndex, GreaterPartitionIndex, &
                        pivotIndex  
                real(kind=dp) :: temp, pivot
                integer :: temp_int, i, j 
                
                if (size(array) <= 1) return 


                !print *, "entering quicksort with array :", array
                ! Choose pivot 

                pivotIndex = size(array) / 2
                pivot = array(pivotIndex)

                ! Now partition into elements below element and above element 

                lesserPartitionIndex = 1
                GreaterPartitionIndex = size(array)

                !print *, "Lesser Parition Index before loop = ", lesserPartitionIndex

                !print *, "Greater Parition Index before loop = ", GreaterPartitionIndex
                print *, "Pivot Index, Pivot value before loop = ", pivotIndex, pivot
                do while (lesserPartitionIndex < GreaterPartitionIndex) 
                        ! Move lesser index until we find element that shouldn't be in the lesser partition 
                        ! Same for the right hand partition
                        !print *, "Entering loop with lesserIndex, GreaterIndex = ", lesserPartitionIndex, GreaterPartitionIndex
                        ! Swap if suitable 

                        do while (array(lesserPartitionIndex) < pivot)
                                lesserPartitionIndex = lesserPartitionIndex + 1
                        end do 

                        do while (array(GreaterPartitionIndex) > pivot)
                                !print *, "GreaterPartitionIndex = ", GreaterPartitionIndex, GreaterPartitionIndex - 1
                                GreaterPartitionIndex = GreaterPartitionIndex - 1
                        end do  

                        if (lesserPartitionIndex >= GreaterPartitionIndex) exit 
                        !print *, "Array before swap : ", array
                        call swap_dp_array_elements(array,lesserPartitionIndex,GreaterPartitionIndex)
                        if (present(integer_companion)) call swap_integer_array_elements(integer_companion,lesserPartitionIndex,&
                                                                GreaterPartitionIndex)

                        if (array(lesserPartitionIndex) == array(GreaterPartitionIndex)) then 
                                lesserPartitionIndex = lesserPartitionIndex + 1 
                        end if 
                        !print *, "Array after  swap : ", array
                        !if (array(lesserPartitionIndex) < pivot) lesserPartitionIndex = lesserPartitionIndex + 1 
                        !if (array(GreaterPartitionIndex) > pivot) GreaterPartitionIndex = GreaterPartitionIndex - 1 


                        !print *, "pivotIndex, pivotValue = ", pivotIndex, pivot
                        !print *, "lesserPartitionIndex, value = ", lesserPartitionIndex, array(lesserPartitionIndex)
                        !print *, "GreaterPartitionIndex, value = ", GreaterPartitionIndex, array(GreaterPartitionIndex)
                end do 
                
                !if (GreaterPartitionIndex >= size(array)) GreaterPartitionIndex = size(array) - 1 
               
                !print *, "About to recurse with lowerParititonIndex, GreaterPartitionIndex = ", lesserPartitionIndex, &
                !        GreaterPartitionIndex
                if (.not. present(integer_companion)) then 
                        call quicksort_dp(array(1:GreaterPartitionIndex))
                        call quicksort_dp(array(GreaterPartitionIndex + 1:size(array)))
                else 
                        call quicksort_dp(array(1:GreaterPartitionIndex),&
                                integer_companion=integer_companion(1:GreaterPartitionIndex))
                        call quicksort_dp(&
                                array(GreaterPartitionIndex + 1:size(array)),&
                                integer_companion=integer_companion(GreaterPartitionIndex + 1:size(array)))
                end if 
        end subroutine quicksort_dp 

        subroutine mergesort_dp(array)
                use iso_fortran_env, only: dp=>real64 
                implicit none 

                real(kind=dp), dimension(:), intent(out), target :: array 
                
                real(kind=dp), dimension(:), allocatable :: scratchSpace  
                real(kind=dp), dimension(:), pointer :: A,B
                logical :: sorted 
                integer :: stat, NBins, i, strideIndex, NBinsTemp, counter  
                integer, parameter :: stride = 2   
                integer, dimension(:,:), allocatable :: Bins, BinsTemp  

                sorted = .False.
                allocate(scratchSpace(size(array)),stat=stat )
                if (stat /= 0) error stop "Failed to allocate scratchSpace" 
                
                NBins = size(array) 
                
                allocate(bins(NBins,2),stat=stat)
                if (stat /= 0) error stop "failed to allocate bins array"
                bins = -1 
                
                allocate(binsTemp(NBins,2),stat=stat)
                if (stat /= 0) error stop "failed to allocate binsTemp array"
                binsTemp = -1

                do i = 1,size(array)
                        Bins(i,1) = i ! Initially each element is it's own bin  
                        Bins(i,2) = i
                end do 

                do while (Nbins > 1)
                        NBinsTemp = NBins 
                        counter = 1
                        do i = 1,NBins/2 + mod(Nbins,2) 

                                strideIndex = (i-1)*stride + 1 
                                if ((mod(NBins,2) /= 0) .and. (i == NBins/2 + mod(Nbins,2))) then 
                                        ! Just need to record unswapped bins into binsTemp in the correct position 
                                        binsTemp(counter,1) = bins(NBins,1)
                                        binsTemp(counter,2) = bins(Nbins,2)
                                        counter = counter + 1
                                        cycle   

                                end if   
                                A => array(Bins(strideIndex,1):Bins(strideIndex,2))
                                B => array(Bins(strideIndex + 1,1):Bins(strideIndex + 1,2))
                                call merge_dp(A,B,scratchSpace)
                                array(Bins(strideIndex,1):Bins(strideIndex + 1,2)) = scratchSpace(1:size(A) + size(B))
                                NBinsTemp = NBinsTemp - 1 ! Each Merge reduces the number of subarrays by 1
                                binsTemp(counter,1) =bins(strideIndex,1) ! Start is the start 
                                binsTemp(counter,2) = bins(strideIndex + 1,2)! End of the new merged array is the end of the RHS
                                                                             ! unmerged array
                                counter = counter + 1

                        end do 
                        BinsTemp(counter:size(Bins,1),:) = Bins(counter:size(Bins,1),:)
                        Bins = BinsTemp 
                        NBins = NBinsTemp 

                end do 

        end subroutine mergesort_dp 

        recursive subroutine merge_dp(A,B,scratchSpace,scratchIndex)
                use iso_fortran_env, only: dp=>real64 
                implicit none 
                real(kind=dp), dimension(:), intent(inout) :: A,B,scratchSpace 
                integer, optional, intent(in) :: scratchIndex
                integer :: SI 

                if (size(scratchSpace) < size(A) + size(B)) error stop "Error: Insufficient scratchspace"

                if (present(scratchIndex)) then 
                        SI = scratchIndex
                else 
                        SI = 1
                end if 
                
                
                if (size(A) == 0) then 
                        scratchSpace(SI:SI + size(B)) = B(:)
                        return 
                else if (size(B) == 0) then 
                        scratchSpace(SI:SI + size(A)) = A(:)
                        return 
                end if 
                
                if (A(1) < B(1)) then 
                        scratchSpace(SI) = A(1)
                        call merge_dp(A(2:size(A)),B,scratchSpace,SI+1)
                else 
                        scratchSpace(SI) = B(1)
                        call merge_dp(A,B(2:size(B)),scratchSpace,SI+1)
                end if 

        end subroutine merge_dp 
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
