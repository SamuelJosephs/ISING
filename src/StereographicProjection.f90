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

                function NSphereProjection(coord, centre, radius) result(angles) ! This was a pain to figure out
                        type(vecNd_t), intent(in) :: coord, centre ! Coordinate on the plane in question and the centre of the
                                                                   ! sphere. 
                        real(kind=8), intent(in) :: radius

                        type(vecNd_t) :: q
                        real(kind=8), allocatable :: angles(:) 
                        type(vecNd_t) :: point1 ,point2, targetPoint, top
                        real(kind=8) :: a, b, c, denominator, qOpt1,qOpt2, discriminant, k,r 
                        integer :: i, j, N
                        if (size(coord) /= size(centre)) error stop "input and central coordinates must have the same dimension."
                        if (size(coord) == 0) error stop "Empty input coordinates."
                        if (.not. allocated(angles)) allocate(angles(size(coord)-1))
                        allocate(targetPoint%coords(size(coord)))
                        angles = 0.0_8
                        ! targetPoint will be a N+1 dimensional point, q will be an N dimensional point. 
                        q = coord - centre 
                        q%coords(size(q)) = 0.0_8
                        top = q 
                        do i = 1, size(top) - 1 
                                top%coords(i) = 0.0_8 
                        end do 
                        top%coords(size(top)) = 2*radius 

                        k = 0.0_8 
                        do i = 1, size(q) - 1 
                            k = k + q%coords(i)**2 
                        end do 
                        k = k / (4*(radius**2))

                        a = k + 1 
                        b = -(4*k + 2)*radius 
                        c = 4*radius**2 * k
                        discriminant = b**2 - 4*a*c 
                        
                        qOpt1 = (-b + sqrt(b**2 - 4*a*c))/(2*a) ! This yields 2*R which passes the sanity check!
                        qOpt2 = (-b - sqrt(b**2 - 4*a*c))/(2*a) ! This yields the coordinate we actually want.
                        targetPoint = coord
                        targetPoint%coords = 0.0_8
                        do i = 1, size(targetPoint) - 1 
                                targetPoint%coords(i) = (q%coords(i)/(2*radius))*(2*radius - qOpt2)
                        end do 
                        targetPoint%coords(size(targetPoint)) = qOpt2 - radius !q - r_C to centre the circle at the origin
                        N = size(targetPoint)
                        do i = 1,size(angles) ! for a 2 sphere in the range [1,2]
                                r = 0.0_8
                                if (i == 1) then 
                                        angles(N-1) = atan2(targetPoint%coords(2),targetPoint%coords(1))
                                        if (angles(N-1) < 0) angles(N-1) = angles(N-1) + 2*3.14159265358979323846 
                                        cycle 
                                end if 

                                do j = 1,i 
                                        r = r + targetPoint%coords(j)**2 
                                end do 
                                r = sqrt(r)
                                angles(N-i) = atan2(r,targetPoint%coords(i+1)) 
                        end do 
                        if (any(angles /= angles)) angles = 0.0_8 !test for NaN's
                end function NSphereProjection 


                

end module stereographicProjection 
