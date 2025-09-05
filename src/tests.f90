
program tests 
        use algo, only: quicksort 
        use iso_fortran_env, only: dp=>real64 
        implicit none 

        real(kind=dp), dimension(3) :: array 


        array = (/9,3,6/)
        call quicksort(array) 
        print *, "Sorted Array = ", array
end program tests 
