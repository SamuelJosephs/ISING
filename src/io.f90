module io
        
        integer, save :: io_NJ, io_ND, io_NB

        public :: io_NJ, io_ND, io_NB
        public :: io_parsefile

        type stringWrapper
                character(len=:), allocatable :: string
        end type stringWrapper 
        contains 

                subroutine io_parsefile(filename)
                        implicit none 
                        character(len=*), intent(in) :: filename 

                        integer :: fileunit, stat, i, temp, startPos, endPos, j
                        character(len=:), allocatable :: TokenBuffer, lineBuffer
                        type(stringWrapper), allocatable, dimension(:) :: TokenArray
                        type(stringWrapper) :: tempStringWrapper

                        open(file=filename, newunit=fileunit,action="read",iostat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to open file"
                        allocate(character(len=1024) :: TokenBuffer, stat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to allocate TokenBuffer"
                        allocate(character(len=1024) :: lineBuffer, stat=stat)
                        if (stat /= 0) error stop "Error: io_parsefile failed to allocate lineBuffer"
                        
                        outermost_do: do 
                                
                              i = 1
                              read(fileunit,'(A)',iostat=stat) lineBuffer
                              if (stat /= 0) exit outermost_do  
                              print *, "Read in line ", lineBuffer
                              
                               
                              call skipWhiteSpace(lineBuffer,i,temp)
                              i = temp
                              if (i == -1) cycle outermost_do 
                              print *, "First non whitespace Character Found = ", lineBuffer(i:i)
                              
                              startPos = i
                              call nextWhiteSpacePosition(lineBuffer,i,endPos)
                              if (endPos == -1) then 
                                      endPos = len(lineBuffer)
                              else 
                                      endPos = endPos - 1
                              end if 

                                
                              ! Now Find Token
                              tempStringWrapper%string = lineBuffer(startPos:endPos)
                              print *, "TokenBuffer = ", TokenBuffer, " len(TokenBuffer) = ", len(TokenBuffer)
                              TokenArray = [TokenArray, tempStringWrapper]
                              
                              do j = 1,size(TokenArray)
                                  print *, "Collected Token: ", TokenArray(j)%string
                              end do 

                        end do outermost_do 


                        close(fileunit)
                end subroutine io_parsefile 

                subroutine nextWhiteSpacePosition(string,position,FirstNonWhiteSpacePosition)
                        character(len=*), intent(in) :: string 
                        integer, intent(in) :: position ! Starting Index to begin search 
                        integer, intent(out) :: FirstNonWhiteSpacePosition
                        
                        integer :: i
                        
                        if (position > len(string)) then 
                                FirstNonWhiteSpacePosition = -1
                                return 
                        end if 

                        i = position 
                        do while (.not. isWhiteSpace(string(i:i)))
                                i = i + 1
                                if (i > len(string)) then 
                                        FirstNonWhiteSpacePosition = -1
                                        return 
                                end if 
                        end do 

                        FirstNonWhiteSpacePosition = i 
                end subroutine nextWhiteSpacePosition
                subroutine skipWhiteSpace(string,position,FirstNonWhiteSpacePosition)
                        character(len=*), intent(in) :: string 
                        integer, intent(in) :: position ! Starting Index to begin search 
                        integer, intent(out) :: FirstNonWhiteSpacePosition
                        
                        integer :: i
                        
                        if (position > len(string)) then 
                                FirstNonWhiteSpacePosition = -1
                                return 
                        end if 

                        i = position 
                        do while (isWhiteSpace(string(i:i)))
                                i = i + 1
                                if (i > len(string)) then 
                                        FirstNonWhiteSpacePosition = -1
                                        return 
                                end if 
                        end do 

                        FirstNonWhiteSpacePosition = i 
                end subroutine skipWhiteSpace
                function isWhiteSpace(a) result(ret)
                        character, intent(in) :: a
                        logical :: ret
                        
                        if (a == achar(32) .or. a == achar(9)) then
                                ret = .True.
                        else 
                                ret = .False.
                        end if 
                end function isWhiteSpace 

                function isChar(a) result(res)
                        character, intent(in) :: a
                        integer :: a_char 
                        logical :: res 
                        
                        a_char = iachar(a)

                        if (( a_char >= iachar('a')) .and. (a_char <= iachar('z'))) then 
                              res = .True.
                        else if ((a_char >= iachar('A')) .and. (a_char <= iachar('Z'))) then
                              res = .True.
                        else 
                              res = .False.
                        end if    

                        return
                end function isChar
end module io
