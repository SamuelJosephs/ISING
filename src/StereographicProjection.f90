module stereographicProjection
        use vecNd
        implicit none

        contains 

                function planeToSphere(coord, centre, radius,free) result(angles)
                        type(vecNd_t), intent(in) :: coord, centre ! Coordinate on the plane in question and the centre of the
                                                                   ! sphere. 
                        real(kind=8), intent(in) :: radius
                        logical, intent(in), optional :: free 
                        type(vecNd_t) :: gradient, temp, point1, point2, targetPoint, top
                        real(kind=8) :: b,c,t1,t2, denominator  
                        real(kind=8), allocatable, save :: angles
                        integer :: i, j 
                        if (size(coord) /= size(centre)) error stop "input and central coordinates must have the same dimension."
                        if (size(coord) == 0) error stop "Empty input coordinates."
                        if (.not. allocated(angles)) allocate(angles(size(coord)))
                        angles = 0.0_8
                        top = centre
                        top%coords(size(top)) = top%coords(size(coord)) + radius 
                        gradient = (top - coord)
                        gradient = gradient / abs(gradient)  
                        temp = coord - centre 

                        b = gradient*temp
                        if (b**2 - 4*c < 0.0_8) error stop "negative square root"
                        c = temp*temp - radius**2

                        t1 = (-b + sqrt(b**2 - 4*c)) / 2 
                        t2 = (-b - sqrt(b**2 - 4*c)) / 2
                
                        point1 = gradient*t1 + coord
                        point2 = gradient*t2 + coord 

                        if (point1 == top) then 
                                targetPoint = point2 
                        else 
                                targetPoint = point1 
                        end if

                        
                        do i = 1,size(angles)
                                denominator = radius 
                                do j = 1,i
                                        denominator = denominator * sin(angles(j))   
                                end do 
                                angles(i) = acos(targetPoint(i) / denominator)
                        end do 
        
                        
                        if (present(free)) then 
                                if (free == .true.) deallocate(angles)
                        end if
                end function planeToSphere 


end module stereographicProjection 
