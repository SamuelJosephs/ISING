
program tests 
        use algo, only: quicksort, mergesort  
        use iso_fortran_env, only: dp=>real64 
        implicit none 

        real(kind=8), dimension(8) :: array 
        integer, dimension(8) :: qs_companion_array


        array = (/9.0,1.0,5.0,7.0,1.0,4.0,2.0,1.0/)
        qs_companion_array = (/1,2,3,4,5,6,7,8/)
        !call quicksort(array,integer_companion=qs_companion_array) 
        !print *, "Sorted Array = ", array
        !print *, "sorted integer companion = ", qs_companion_array
        print *, "Array Before Sorting = ", array 
        call mergesort(array, integer_companion=qs_companion_array) 
        print *, "Array After Sorting = ", array
        print *, "Companion Array After Sorting = ", qs_companion_array
end program tests 
