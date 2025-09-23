module io
        
        integer, save :: io_NJ, io_ND, io_NB

        public :: io_NJ, io_ND, io_NB

        contains 


                subroutine io_parsefile(filename)
                        character(len=*), intent(in) :: filename
                        
                        integer :: fileunit, stat, counter
                        character(len=:), allocatable :: lineBuffer, tagBuffer
                        open(newunit=fileunit, file=filename, action="read", iostat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to open input file"

                        ! Read file line by line:

                        stat = 1 
                        counter = 1
                        do while (stat /= 0)
                               read(fileunit,*, iostat=stat) lineBuffer
                               ! Now parse each line 
                               do i = 1,len(lineBuffer)
                                   if (isWhiteSpace(lineBuffer(i:i))) then
                                           cycle
                                   else 
                                   ! When we encounter a non whitespace space character place it into          
                                        counter = i
                                        do while (.not. isWhiteSpace(lineBuffer(counter:counter)))
                                                if (counter > len(lineBuffer)) exit ! Check if we have reached the end of the line 
                                                if (lineBuffer(counter:counter) == ':') exit
                                                counter = counter + 1

                                        end do  

                                        tagBuffer = lineBuffer(i:counter)
                                        print *, "TagBuffer  = ", tagBuffer
                                   end if  
                               end do 
                        end do 
                        close(fileunit)
                end subroutine io_parsefile

                function isWhiteSpace(a) result(ret)
                        character, intent(in) :: a
                        logical :: ret
                        
                        if (a == achar(32) .or. a == achar(9)) then
                                ret = .True.
                        else 
                                ret = .False.
                        end if 
                 end function isWhiteSpace
                        

end module io
