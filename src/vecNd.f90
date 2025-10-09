
module vecNd  
        implicit none 
        
        type vecNd_t  
                real(kind=8), allocatable:: coords(:) 
        end type vecNd_t   
        
        ! Secondary type for shared views
        type vecNd_view_t
                real(kind=8), pointer :: coords(:) => null()
        end type vecNd_view_t
        
        interface operator(+) 
                module procedure add_Nvec
        end interface 

        interface operator(-)
                module procedure sub_Nvec
        end interface 

        interface operator(*)
                module procedure vecNd_scalar_product
                module procedure vecNd_scalar_multiplication
                module procedure vecNd_scalar_multiplication2

        end interface 

        interface assignment(=)
                module procedure vecNDINIT
                module procedure VECNdINIT_array
        end interface
        
        interface operator(/)
                module procedure vecNd_scalar_division
        end interface

        interface operator(.x.)
                module procedure cross_product
                module procedure cross_product_view
                module procedure cross_product_mixed1
                module procedure cross_product_mixed2
        end interface

        interface abs 
                module procedure abs_nvec
                module procedure abs_nvec_array
                module procedure abs_nvec_view
        end interface 

        interface operator(==)
                module procedure vecNd_eq
                module procedure vecNd_eq_array
                module procedure vecNd_view_eq
        end interface

        interface size
                module procedure vecNdSize
                module procedure vecNdViewSize

        end interface size 

        ! New interface to create a view
        interface view
                module procedure createView
                module procedure createViewFromArray
        end interface view

        contains 

        function vecSTP(A,B,C) result(res)
                type(vecNd_t), intent(in) :: A,B,C
                real(kind=dp) :: res

                res = A*(B.x.C)
        
        end function vecSTP

        function vecNdSize(input) result(res)
                type(vecNd_t), intent(in) :: input
                integer :: res 

                res = size(input%coords)
        end function vecNdSize 
        
        function vecNdViewSize(input) result(res)
                type(vecNd_view_t), intent(in) :: input
                integer :: res 

                res = size(input%coords)
        end function vecNdViewSize
        function normalSize(input) result(res)
                real(kind=8) :: input(:)
                integer :: res 

                res = size(input)

        end function normalSize
        function makeVecNd(input) result(res)
                real(kind=8), intent(in) :: input(:) 
                type(vecNd_t):: res 
                
                allocate(res%coords(size(input))) 
                res%coords = input 
        end function makeVecNd

        function makeVecNdCheck(self, input) result(res) ! Checks for allocation and avoids any unneccessary allocations
                type(vecNd_t), intent(inout) :: self
                type(vecNd_t) :: res
                real(kind = 8), intent(in) :: input(:)

                if (allocated(self%coords)) then
                        if (size(self%coords) == size(input)) then
                                call move_alloc(self%coords,res%coords)
                                res%coords = input
                                return
                        else 
                                deallocate(self%coords)
                                allocate(res%coords(size(input)))
                                res%coords = input
                                return
                        end if 
                end if 
                allocate(res%coords(size(input)))
                res%coords = input
        end function makeVecNdCheck
        
        ! Function to create a view of a vector without copying data
        function createView(src) result(res)
                type(vecNd_t), target, intent(in) :: src
                type(vecNd_view_t) :: res
                
                ! Check if source vector is allocated
                if (.not. allocated(src%coords)) then
                        error stop "Error: Creating view of uninitialized vector."
                end if
                
                ! Point to the source data without copying
                res%coords => src%coords
        end function createView
        
        ! Function to create a view from a stack-allocated array
        function createViewFromArray(array) result(res)
                real(kind=8), target, intent(in) :: array(:)
                type(vecNd_view_t) :: res
                
                ! Point directly to the stack-allocated array
                res%coords => array
        end function createViewFromArray
        
        ! Function to convert a view back to a standard vector (with copying)
        function vecFromView(v) result(res)
                type(vecNd_view_t), intent(in) :: v
                type(vecNd_t) :: res
                
                if (.not. associated(v%coords)) then
                        error stop "Error: Converting uninitialized view."
                end if
                
                allocate(res%coords(size(v%coords)))
                res%coords = v%coords
        end function vecFromView

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
        
        ! Cross product function (primarily for 3D vectors)
        function cross_product(vec1, vec2) result(res)
                type(vecNd_t), intent(in) :: vec1, vec2
                type(vecNd_t) :: res
                integer :: dim1, dim2
                
                dim1 = size(vec1%coords)
                dim2 = size(vec2%coords)
                
                ! Check if dimensions are appropriate for cross product
                if ((.not. allocated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then
                        error stop "Error: Cross product with uninitialized vectors."
                else if (dim1 /= dim2) then
                        error stop "Error: Cross product with vectors of different dimensions."
                else if (dim1 /= 3) then
                        error stop "Error: Cross product is only implemented for 3D vectors."
                end if
                
                ! Allocate result
                allocate(res%coords(3))
                
                ! Calculate cross product for 3D vectors
                res%coords(1) = vec1%coords(2) * vec2%coords(3) - vec1%coords(3) * vec2%coords(2)
                res%coords(2) = vec1%coords(3) * vec2%coords(1) - vec1%coords(1) * vec2%coords(3)
                res%coords(3) = vec1%coords(1) * vec2%coords(2) - vec1%coords(2) * vec2%coords(1)
        end function cross_product

        function abs_nvec(input) result(res)
                type(vecNd_t), intent(in) :: input 
                real(kind=8) :: res 
                integer :: i 
                res = 0.0
                
                do i=1,size(input%coords)
                        res = res + input%coords(i)**2 
                end do 
                res = sqrt(res)
        end function abs_nvec
        
        function abs_nvec_view(input) result(res)
                type(vecNd_view_t), intent(in) :: input 
                real(kind=8) :: res 
                integer :: i 
                res = 0.0
                
                do i=1,size(input%coords)
                        res = res + input%coords(i)**2 
                end do 
                res = sqrt(res)
        end function abs_nvec_view

        function abs_nvec_array(input) result(output_array)
                type(vecNd_t), intent(in) :: input(:) 
                real(kind = 8), allocatable :: output_array(:)
                integer :: i 
                allocate(output_array(size(input)))
                
                do i = 1, size(input)
                        output_array(i) = abs(input(i))
                end do
        end function abs_nvec_array

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
        
        function vecNd_view_eq(vec1,vec2) result(res)
                type(vecNd_view_t), intent(in) :: vec1
                type(vecNd_t), intent(in) :: vec2 
                logical :: res 
                integer :: i 
                res = .True.
                if ((.not. associated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then 
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
        
        end function vecNd_view_eq
        
        ! Cross product for view type
        function cross_product_view(vec1, vec2) result(res)
                type(vecNd_view_t), intent(in) :: vec1
                type(vecNd_view_t), intent(in) :: vec2
                type(vecNd_t) :: res
                integer :: dim1, dim2
                
                dim1 = size(vec1%coords)
                dim2 = size(vec2%coords)
                
                ! Check if dimensions are appropriate for cross product
                if ((.not. associated(vec1%coords)) .or. (.not. associated(vec2%coords))) then
                        error stop "Error: Cross product with uninitialized vector views."
                else if (dim1 /= dim2) then
                        error stop "Error: Cross product with vector views of different dimensions."
                else if (dim1 /= 3) then
                        error stop "Error: Cross product is only implemented for 3D vectors."
                end if
                
                ! Allocate result
                allocate(res%coords(3))
                
                ! Calculate cross product for 3D vectors
                res%coords(1) = vec1%coords(2) * vec2%coords(3) - vec1%coords(3) * vec2%coords(2)
                res%coords(2) = vec1%coords(3) * vec2%coords(1) - vec1%coords(1) * vec2%coords(3)
                res%coords(3) = vec1%coords(1) * vec2%coords(2) - vec1%coords(2) * vec2%coords(1)
        end function cross_product_view
        
        ! Cross product mixing vecNd_t and vecNd_view_t
        function cross_product_mixed1(vec1, vec2) result(res)
                type(vecNd_t), intent(in) :: vec1
                type(vecNd_view_t), intent(in) :: vec2
                type(vecNd_t) :: res
                integer :: dim1, dim2
                
                dim1 = size(vec1%coords)
                dim2 = size(vec2%coords)
                
                ! Check if dimensions are appropriate for cross product
                if ((.not. allocated(vec1%coords)) .or. (.not. associated(vec2%coords))) then
                        error stop "Error: Cross product with uninitialized vectors."
                else if (dim1 /= dim2) then
                        error stop "Error: Cross product with vectors of different dimensions."
                else if (dim1 /= 3) then
                        error stop "Error: Cross product is only implemented for 3D vectors."
                end if
                
                ! Allocate result
                allocate(res%coords(3))
                
                ! Calculate cross product for 3D vectors
                res%coords(1) = vec1%coords(2) * vec2%coords(3) - vec1%coords(3) * vec2%coords(2)
                res%coords(2) = vec1%coords(3) * vec2%coords(1) - vec1%coords(1) * vec2%coords(3)
                res%coords(3) = vec1%coords(1) * vec2%coords(2) - vec1%coords(2) * vec2%coords(1)
        end function cross_product_mixed1
        
        ! Cross product mixing vecNd_view_t and vecNd_t
        function cross_product_mixed2(vec1, vec2) result(res)
                type(vecNd_view_t), intent(in) :: vec1
                type(vecNd_t), intent(in) :: vec2
                type(vecNd_t) :: res
                integer :: dim1, dim2
                
                dim1 = size(vec1%coords)
                dim2 = size(vec2%coords)
                
                ! Check if dimensions are appropriate for cross product
                if ((.not. associated(vec1%coords)) .or. (.not. allocated(vec2%coords))) then
                        error stop "Error: Cross product with uninitialized vectors."
                else if (dim1 /= dim2) then
                        error stop "Error: Cross product with vectors of different dimensions."
                else if (dim1 /= 3) then
                        error stop "Error: Cross product is only implemented for 3D vectors."
                end if
                
                ! Allocate result
                allocate(res%coords(3))
                
                ! Calculate cross product for 3D vectors
                res%coords(1) = vec1%coords(2) * vec2%coords(3) - vec1%coords(3) * vec2%coords(2)
                res%coords(2) = vec1%coords(3) * vec2%coords(1) - vec1%coords(1) * vec2%coords(3)
                res%coords(3) = vec1%coords(1) * vec2%coords(2) - vec1%coords(2) * vec2%coords(1)
        end function cross_product_mixed2

        function vecNd_eq_array(input_array, vec) result(output_Array) 
                type(vecNd_t), intent(in) :: input_array(:) 
                type(vecNd_t), intent(in) :: vec 
                logical, allocatable :: output_array(:)
                integer :: i 
                allocate(output_array(size(input_array)))
                
                do i = 1, size(input_array)
                        output_array(i) = input_array(i) == vec 
                end do                
        end function vecNd_eq_array 

        function vecNd_scalar_product(vec1,vec2) result(res)
                type(vecNd_t), intent(in) :: vec1, vec2 
                real(kind=8) :: res 
                integer :: i 
                res = 0.0_8
                do i = 1,size(vec1)
                        res = res + vec1%coords(i) * vec2%coords(i)
                end do 
                
        end function vecNd_scalar_product

        function vecNd_scalar_division(vec1,a) result(res) 
                type(vecNd_t), intent(in) :: vec1 
                real(kind=8), intent(in) :: a 
                type(vecNd_t) :: res 
                integer :: i 
                res = vec1  
                do i = 1,size(vec1)
                        res%coords(i) = res%coords(i) / a 
                end do 
                
        end function vecNd_scalar_division

        function vecNd_scalar_multiplication(vec1,a) result(res) 
                type(vecNd_t), intent(in) :: vec1 
                real(kind=8), intent(in) :: a 
                type(vecNd_t) :: res 
                integer :: i 
                res = vec1  
                do i = 1,size(vec1)
                        res%coords(i) = res%coords(i) * a 
                end do 
                
        end function vecNd_scalar_multiplication

        function vecNd_scalar_multiplication2(a,vec1) result(res) 
                type(vecNd_t), intent(in) :: vec1 
                real(kind=8), intent(in) :: a 
                type(vecNd_t) :: res 
                integer :: i 
                res = vec1  
                do i = 1,size(vec1)
                        res%coords(i) = res%coords(i) * a 
                end do 
                
        end function vecNd_scalar_multiplication2

        subroutine vecNdINIT(lhs,rhs) ! input1 = input2
                type(vecNd_t), intent(inout) :: lhs
                type(vecNd_t), intent(in) :: rhs 
                if (.not. allocated(rhs%coords)) error stop "RHS is not allocated"               
                if (allocated(lhs%coords)) then 
                        if (size(lhs) == size(rhs)) then 
                                lhs%coords = rhs%coords
                                return 
                        else 
                                deallocate(lhs%coords)
                                allocate(lhs%coords(size(rhs)))
                                lhs%coords = rhs%coords
                                return 
                        end if

                end if
                allocate(lhs%coords(size(rhs)))

                lhs%coords = rhs%coords
                return
                
        end subroutine vecNDINIT

        subroutine vecNdINIT_array(lhs,rhs) ! input1 = input2
                type(vecNd_t), intent(inout) :: lhs
                real(kind=8), intent(in) :: rhs(:) 
                if (allocated(lhs%coords)) then 
                        if (size(lhs%coords) == size(rhs)) then 
                                lhs%coords = rhs 
                                return 
                        else 
                                deallocate(lhs%coords)
                                allocate(lhs%coords(size(rhs)))
                                lhs%coords = rhs 
                                return 
                        end if

                end if
                allocate(lhs%coords(size(rhs)))

                lhs%coords = rhs
                return
                
        end subroutine vecNdINIT_array 
end module vecNd
