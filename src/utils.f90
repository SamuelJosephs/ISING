module utils

      public :: utils_to_upper
      contains
              subroutine utils_to_upper(s)
                      character(len=*), intent(inout) :: s
                      integer :: i, ich
                              
                      do i = 1, len_trim(s)
                              ich = iachar(s(i:i)) 
                              if (ich >= iachar('a') .and. ich <= iachar('z')) then
                                      s(i:i) = achar(ich - 32)   ! shift lowercase to uppercase
                              end if
                      end do  
              end subroutine utils_to_upper



end module utils 
