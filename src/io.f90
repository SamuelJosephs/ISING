module io
        use iso_fortran_env, only: dp=>real64
        integer, save :: io_NJ, io_ND, io_NB
        real(kind=dp) :: io_JMIN, io_JMAX, io_DMIN, io_DMAX, io_BMIN, io_BMAX

        character, parameter, dimension(2) :: SeperatorArray = (/'=', ':'/)
        character, parameter, dimension(2) :: OperatorArray  = (/'=', ':'/)
        logical, parameter, dimension(2) :: IsBinaryOperator = (/.True., .True./)

        public :: io_NJ, io_ND, io_NB, io_JMIN, io_JMAX, io_DMIN, io_DMAX
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
                              print *, "Read in line: ", trim(lineBuffer)
                              
                              call TokensInLine(lineBuffer,TokenArray)

                              call ParseTokensInLine(TokenArray)

                              
                        end do outermost_do 

                        print *,"io_NJ = ", io_NJ
                        print *,"io_ND = ", io_ND
                        print *,"io_NB = ", io_NB
                        print *,"io_JMin = ", io_JMIN
                        print *,"io_JMax = ", io_JMAX
                        print *, "io_DMin = ", io_DMIN
                        print *, "io_DMax = ", io_DMAX
                        print *, "io_BMin = ", io_BMIN
                        print *, "io_BMax = ", io_BMAX

                        close(fileunit)
                end subroutine io_parsefile 

                subroutine ParseTokensInLine(TokenArray)
                        type(stringWrapper), dimension(:), intent(in) :: TokenArray
                        
                        integer :: i, j, operatorIndex
                        character :: operatorChar
                        outer_loop: do i = 1,size(TokenArray)
                                ! Test to see if it is in operatorArray 
                                do j = 1,size(operatorArray)
                                        if (trim(adjustl(TokenArray(i)%string)) == operatorArray(j)) then 
                                               
                                               if (isBinaryOperator(j)) then 
                                                       operatorIndex = i
                                                       call ParseBinaryOperator(TokenArray,OperatorIndex)
                                                       cycle outer_loop
                                               end if 
                                        end if 
                                end do 
                        end do outer_loop
                end subroutine ParseTokensInLine
                
                subroutine ParseBinaryOperator(TokenArray,OperatorIndex)
                        type(stringWrapper), dimension(:), intent(in) :: TokenArray
                        integer, intent(in) :: OperatorIndex

                        character(len=:), allocatable :: VarArray, valArray
                        integer :: stat
                        if (TokenArray(OperatorIndex)%string == '=') then 
                                if (OperatorIndex - 1 < 1) error stop "Error: Invalid input string encountered"
                                if (OperatorIndex + 1 > size(TokenArray)) error stop "Error: Invalid input string encountered"
                                
                                varArray = TokenArray(OperatorIndex - 1)%string
                                valArray = TokenArray(OperatorIndex + 1)%string 

                                call to_upper(varArray)
                                if (varArray == "NJ") then 
                                        read(valArray,*,iostat=stat) io_NJ
                                else if (varArray == "ND") then 
                                        read(valArray,*,iostat=stat) io_ND
                                else if (varArray == "NB") then 
                                        read(valArray,*,iostat=stat) io_NB
                                else if (varArray == "JMIN") then 
                                        read(valArray,*,iostat=stat) io_JMIN
                                else if (varArray == "JMAX") then 
                                        read(valArray,*,iostat=stat) io_JMAX
                                else if (varArray == "DMIN") then 
                                        read(valArray,*,iostat=stat) io_DMIN
                                else if (varArray == "DMAX") then 
                                        read(valArray,*,iostat=stat) io_DMAX
                                else if (varArray == "BMIN") then 
                                        read(valArray,*,iostat=stat) io_BMIN
                                else if (varArray == "BMAX") then 
                                        read(valArray,*,iostat=stat) io_BMAX
                                else 
                                        error stop "Error: Unrecognised LHS Value for the = operator"
                                end if 
                                
                                if (stat /= 0) error stop "Error: Failed to parse input parameter"


                        end if 
                end subroutine ParseBinaryOperator
                subroutine TokensInLine(lineBuffer, TokenArray)
                        implicit none 

                        character(len=*), intent(in) :: lineBuffer
                        type(stringWrapper), allocatable, dimension(:), intent(inout) :: TokenArray
                        integer :: i, stat, temp
                        type(stringWrapper) :: TokenBuffer
                        if (allocated(TokenArray)) deallocate(TokenArray)
                        allocate(TokenArray(0), stat=stat)
                        if (stat /= 0) error stop "Error: TokensInLine Failed to allocate TokenArray"
                       
                        i = 1

                        do while (i < len(lineBuffer))
                                
                                call skipWhiteSpace(lineBuffer,i,temp)
                                if (temp == -1) return 
                                i = temp
                                if (any(lineBuffer(i:i) == seperatorArray)) then 
                                        TokenBuffer%string = lineBuffer(i:i)
                                        TokenArray = [TokenArray,TokenBuffer]
                                        if (i == len(lineBuffer)) then
                                                return 
                                        else 
                                                i = i + 1
                                                cycle 
                                        end if 
                                end if 
                                call endOfTokenPosition(lineBuffer,i,temp)
                                if (temp == -1) then 
                                       temp = len(lineBuffer)
                                else
                                      temp = temp - 1
                                end if 

                                TokenBuffer%string = lineBuffer(i:temp)

                                TokenArray = [TokenArray, TokenBuffer]
                                
                                i = i + len(TokenBuffer%string)
                        end do  

                        
                end subroutine TokensInLine

                subroutine endOfTokenPosition(string,position,FirstWhiteSpacePosition)
                        implicit none
                        character(len=*), intent(in) :: string 
                        integer, intent(in) :: position ! Starting Index to begin search 
                        integer, intent(out) :: FirstWhiteSpacePosition
                        
                        integer :: i
                        if (position > len(string)) then 
                                FirstWhiteSpacePosition = -1
                                return 
                        end if 

                        i = position 
                        do while ((.not. isWhiteSpace(string(i:i))) .and. &
                                (.not. (any(string(i:i) == SeperatorArray)))&
                                )
                                i = i + 1
                                if (i > len(string)) then 
                                        FirstWhiteSpacePosition = -1
                                        return 
                                end if 
                        end do 

                        FirstWhiteSpacePosition = i 
                end subroutine endOfTokenPosition

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

                subroutine to_upper(s)
                        character(len=*), intent(inout) :: s
                        integer :: i, ich

                        do i = 1, len_trim(s)
                                ich = iachar(s(i:i))
                                if (ich >= iachar('a') .and. ich <= iachar('z')) then
                                        s(i:i) = achar(ich - 32)   ! shift lowercase to uppercase
                                end if
                        end do
                end subroutine to_upper
end module io
