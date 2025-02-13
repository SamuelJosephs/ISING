module vecNd  

        implicit none 
        type vecNd_t  
                real, allocatable:: coords(:) 
        end type vecNd_t   
        

        interface operator(+) 
                module procedure add_Nvec
        end interface 

        interface operator(-)
                module procedure sub_Nvec
        end interface 

        interface abs 
                module procedure abs_nvec
        end interface 

        interface operator(==)
                module procedure vecNd_eq
        end interface
        contains 
        function makeVecNd(input) result(res)
                real, intent(in) :: input(:) 
                type(vecNd_t):: res 

                allocate(res%coords(size(input))) 
                res%coords = input 
                
        end function makeVecNd 

        function add_Nvec(vec1, vec2) result(res) 
                type(vecNd_t), intent(in) :: vec1, vec2 
                type(vecNd_t) :: res 
                integer :: i 
                if ((.not. allocated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then 
                        error stop "Error: Adding unitialised vectors."

                else if (size(vec1%coords) /= size(vec2%coords)) then 
                        error stop "Error: Adding vectors with different dimensions."
                end if 


                allocate(res%coords(size(vec1%coords)))
                
                do i = 1,size(res%coords)
                        res%coords(i) = vec1%coords(i) + vec2%coords(i)
                end do 
                

        end function add_NVec

        function sub_Nvec(vec1, vec2) result(res) 
                type(vecNd_t), intent(in) :: vec1, vec2 
                type(vecNd_t) :: res 
                integer :: i 
                if ((.not. allocated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then 
                        error stop "Error: Subtracting unitialised vectors."

                else if (size(vec1%coords) /= size(vec2%coords)) then 
                        error stop "Error: Subtracting vectors with different dimensions."
                end if 


                allocate(res%coords(size(vec1%coords)))
                
                do i = 1,size(res%coords)
                        res%coords(i) = vec1%coords(i) - vec2%coords(i)
                end do 
                

        end function sub_NVec

        function abs_nvec(input) result(res)
                type(vecNd_t), intent(in) :: input 
                real :: res 
                integer :: i 
                res = 0.0
                
                do i=1,size(input%coords)
                        res = res + input%coords(i)**2 
                end do 
                res = sqrt(res)
        end function abs_nvec

        function vecNd_eq(vec1,vec2) result(res)
                type(vecNd_t), intent(in) :: vec1, vec2 
                logical :: res 
                integer :: i 
                res = .True.
                if ((.not. allocated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then 
                        error stop "Error: equating unitialised vectors."

                else if (size(vec1%coords) /= size(vec2%coords)) then 
                        error stop "Error: equating vectors with different dimensions."
                end if                 

                do i = 1, size(vec1%coords) 
                        if (abs(vec1%coords(i) - vec2%coords(i)) > 1e-8) then 
                                res = .False.
                                return 
                        end if 
                end do 
        
        end function vecNd_eq
end module vecNd 
