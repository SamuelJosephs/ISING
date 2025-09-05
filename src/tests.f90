
program tests 
        use algo, only: quicksort 
        use iso_fortran_env, only: dp=>real64 
        implicit none 

        real(kind=8), dimension(3) :: array 
        integer, dimension(3) :: qs_companion_array


        array = (/9,3,6/)
        qs_companion_array = (/1,2,3/)
        call quicksort(array,integer_companion=qs_companion_array) 
        print *, "Sorted Array = ", array
        print *, "sorted integer companion = ", qs_companion_array
end program tests 
