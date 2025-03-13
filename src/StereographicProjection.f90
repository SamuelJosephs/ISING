module stereographicProjection
        use vecNd
        implicit none

        public :: NSphereProjection
        contains 

                function planeToSphere(coord, centre, radius,free) result(angles)
                        type(vecNd_t), intent(in) :: coord, centre ! Coordinate on the plane in question and the centre of the
                                                                   ! sphere. 
                        real(kind=8), intent(in) :: radius
                        logical, intent(in), optional :: free 
                        type(vecNd_t) :: gradient, temp, point1, point2, targetPoint, top
                        real(kind=8) :: b,c,t1,t2, denominator  
                        real(kind=8), allocatable :: angles(:)
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
                                angles(i) = acos(targetPoint%coords(i) / denominator)
                        end do 
        
                        
                        if (present(free)) then 
                                if (free) deallocate(angles)
                        end if
                end function planeToSphere 

                function NSphereProjection(coord, centre, radius) result(angles)
                        type(vecNd_t), intent(in) :: coord, centre ! Coordinate on the plane in question and the centre of the
                                                                   ! sphere. 
                        real(kind=8), intent(in) :: radius

                        real(kind=8), allocatable :: angles(:) 
                        type(vecNd_t) :: point1 ,point2, targetPoint, top
                        real(kind=8) :: a, b, c, denominator, qOpt1,qOpt2, discriminant, k 
                        integer :: i, j 
                        if (size(coord) /= size(centre)) error stop "input and central coordinates must have the same dimension."
                        if (size(coord) == 0) error stop "Empty input coordinates."
                        if (.not. allocated(angles)) allocate(angles(size(coord)))
                        allocate(targetPoint%coords(size(coord)+1))
                        angles = 0.0_8
                        top = centre
                        top%coords(size(top)) = top%coords(size(coord)) + radius 
                        

                        k = 0.0
                        do i=1,size(coord)
                                k = k + coord%coords(i)**2 
                        end do 
                        k = k / (4*radius**2)

                        a = k + 1 
                        b = 4*radius*k 
                        c = 4*k - 1 
                        discriminant = b**2 - 4*a*c 
                        if (discriminant < 0.0_8) error stop "Encountered negative square root."
                        qOpt1 = (-b + sqrt(discriminant))/(2*a)
                        qOpt2 = (-b - sqrt(discriminant))/(2*a)
                        if (abs(qOpt1 - 2*radius) < 1e-4) then 
                                targetPoint%coords(size(targetPoint)) = qOpt2 
                        else if (abs(qOpt2 - 2*radius) < 1e-4) then 
                                targetPoint%coords(size(targetPoint)) = qOpt1 
                        else
                                error stop "Unable to correctly compute projection."
                        end if

                        do i = 1, size(targetPoint)-1
                                targetPoint%coords(i) = (coord%coords(i)/(2*radius)) * (2*radius -  targetPoint%coords(size(targetPoint)))
                        end do 

                        do i = 1,size(angles)
                                denominator = radius 
                                do j = 1,i
                                        denominator = denominator * sin(angles(j))   
                                end do 
                                angles(i) = acos(targetPoint%coords(i) / denominator)
                        end do 
        
                end function NSphereProjection 


                

end module stereographicProjection 
