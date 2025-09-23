module io
        
        integer, save :: io_NJ, io_ND, io_NB

        public :: io_NJ, io_ND, io_NB
        public :: io_parsefile
        contains 


                subroutine io_parsefile(filename)
                        implicit none 
                        character(len=*), intent(in) :: filename
                        
                        integer :: fileunit, stat, counter, temp 
                        character(len=:), allocatable :: lineBuffer, tagBuffer, valBuffer
                        logical :: foundSeperator
                        open(newunit=fileunit, file=filename, action="read", iostat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to open input file"

                        
                        allocate(character(len=1024) :: lineBuffer,stat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to allocate array"
                        ! Read file line by line:

                        counter = 1
                        do 
                               read(fileunit,'(A)', iostat=stat) lineBuffer
                               print *, "Read in line : ", lineBuffer
                               if (stat /= 0) exit
                               ! Now parse each line 
                               call io_Token(lineBuffer,counter,temp,tagBuffer)
                               counter = temp 
                               print *, "Counter after io_Token = ", counter
                               if (counter < len(lineBuffer) - 1) then 
                                        ! Find : 
                                        do 
                                                if (counter > len(lineBuffer)) then 
                                                        foundSeperator = .False.
                                                        exit 
                                                else if (lineBuffer(counter:counter) == ':') then 
                                                        foundSeperator = .True.
                                                        exit 
                                                else 
                                                        counter = counter + 1 
                                                        cycle 
                                                end if 
                                                
                                        end do 
                                        if (foundSeperator .and. (counter < len(lineBuffer))) then 
                                                call io_Token(lineBuffer,counter,temp,valBuffer)
                                                counter = temp
                                        end if 
                               end if 

                        print *, "Tag = ", trim(tagBuffer), "val =", trim(valBuffer)
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
                        
                 subroutine io_Token(string,startIndex,finishPointer,res) 
                         implicit none 
                         character(len=*), intent(in) :: string 
                         integer, intent(in) :: startIndex
                         integer, intent(out) :: finishPointer 
                         character(len=:), allocatable, intent(out) :: res 

                         integer :: i, stat, counter   
                         ! From startIndex return the first token, that is a string of consecutive non whitespace characters 

                         !allocate(character(len=1024) :: res, stat=stat)
                         
                         if (startIndex > len(string)) then 
                                 finishPointer = startIndex 
                                 return 
                         end if 
                         i = startIndex 

                         do 

                                if (i > len(string)) exit 
                                if (isWhiteSpace(string(i:i))) then 
                                        i = i + 1 
                                        cycle 
                                end if 
                                
                                ! Now that we know that string(i:i) is our first non whitespace character start appending to res 
                                counter = i
                                do while(.not. isWhiteSpace(string(counter:counter)))
                                        counter = counter + 1 
                                        if (counter > len(string)) exit
                                end do 
                                res = string(i:counter)
                                exit 
                         end do 
                         finishPointer = counter 
                         return 
                 end subroutine io_Token 

end module io
