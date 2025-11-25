
!!!!!!!!!!! Module Instructions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FFT planning and array book keeping should be kept seperate from actual FFT's. All the book keeping is encapsulated by the type
! fft_object, on which all FFT routines should act.





module fft
use iso_c_binding
include 'fftw3.f03' 

type fft_object
        ! At the moment these objects can be initialised with create_plan_nd routines, only inplace are currently supported due to
        ! time constraints.
         
        type(C_ptr) :: plan_forward, plan_backward
        type(C_ptr) :: fft_array_real_base_ptr ! Points to the first element of the real space array 
        type(C_ptr) :: fft_array_recip_base_ptr ! Points to the first element of the reciprocal space array 
        integer :: real_buffer_len ! Length of the array buffer, proper C programming in Fortran I'm sure there is no risk of stack
                                   ! smashing myself again.
        integer :: recip_buffer_len                            ! These are type agnositc, stores either the number of real or
                                                               ! complex entries in the respective buffer the C pointers address.
        integer, dimension(:), allocatable :: RealBufferShape  
        integer, dimension(:), allocatable :: RecipBufferShape ! (shape of the reciprocal space array, can differ from that of the
                                                               ! real space shape for an in place transform), size(BufferShape) is
                                                               ! the dimension of the FFT for the object, the entries give the
                                                               ! shape.
        logical :: in_place ! Needed for cleanup

end type fft_object

public :: fft_2d_r2c, fft_2d, fft_object, create_plan_2d_r2c, create_plan_2d, c2f_indexing_3d, c2f_indexing_2d, fft_alloc_real 

contains 

        subroutine create_plan_2d_r2c(fft_obj,Nx,Ny,outputBufferReal,outputBufferComplex, inplace)
                ! This routine is used to initialise an fft_object such that an in place transformation can be done in outputBuffer
                implicit none

                type(fft_object), intent(inout) :: fft_obj
                integer, intent(in) :: Nx, Ny
                real(kind=c_double), intent(out), pointer, dimension(:,:) :: outputBufferReal ! This is just a reference to the memory
                complex(kind=c_double_complex), intent(out), pointer, dimension(:,:) :: outputBufferComplex
                logical, intent(in), optional :: inplace

                logical :: myInPlace
                integer(kind=c_size_t) :: Nx_c, Ny_c
                integer :: stat 
                

                myInPlace = .False.
                if (present(inplace)) myInPlace = inplace
                fft_obj%in_place = myInPlace
                Nx_c = int(Nx,c_size_t)
                Ny_c = int(Ny,c_size_t)

                ! We only need Nx/2 +1 complex numbers in the transformed, +1 for the 0 frequency, then N up to the Nyquist
                ! frequency.
                ! complex space due to symmetry above the Nyquist frequency FFTW reuses these values for lower memory
                ! overhead and faster FFT's. 

                ! Periodicity requires that in reciprocal space every N entries are identical. Thus f^twiddle(N-k) = f(k)^*
                ! For frequencies abovet the Nyquist frequency this maps to  
                ! Note: This is only true for real valued funcitons, where f = f^*
                fft_obj%fft_array_real_base_ptr = fftw_alloc_real((2*(Nx_c/2 + 1))*Ny_c) 
                if (myInPlace) then 
                        fft_obj%fft_array_recip_base_ptr = fft_obj%fft_array_real_base_ptr
                else 
                        fft_obj%fft_array_recip_base_ptr = fftw_alloc_real((2*(Nx_c/2 + 1))*Ny_c)  
                end if 
                                                                                               
                call C_F_POINTER(fft_obj%fft_array_real_base_ptr,outputBufferReal,[Nx,Ny])

                call C_F_POINTER(fft_obj%fft_array_recip_base_ptr,outputBufferComplex,[Nx/2 + 1,Ny])
                 ! Now outputBuffer will behave like an ordinary 2d array, 
                outputBufferReal(:,:) = 0.0_c_double ! Zero the array to avoid uninitialised data later.

                fft_obj%real_buffer_len = product(shape(outputBufferReal))
                fft_obj%recip_buffer_len = product(shape(outputBufferComplex))

                allocate(fft_obj%RealBufferShape(2),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate real buffer shape"

                allocate(fft_obj%RecipBufferShape(2),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate reciprocal buffer shape"

                fft_obj%RealBufferShape = shape(outputBufferReal)
                fft_obj%RecipBufferShape = shape(outputBufferComplex)

                ! Now we need to initialise the plans 

                ! We have to transpose the shapes in the actual C calls.
                fft_obj%plan_forward = fftw_plan_dft_r2c_2d(Ny,Nx,outputBufferReal,outputBufferComplex,FFTW_ESTIMATE)
                fft_obj%plan_backward = fftw_plan_dft_c2r_2d(Ny,Nx,outputBufferComplex,outputBufferReal,FFTW_ESTIMATE)
                

        end subroutine create_plan_2d_r2c

        subroutine create_plan_2d(fft_obj, Nx, Ny,inplace)
                use iso_c_binding
                implicit none 
                type(fft_object), intent(out) :: fft_obj 
                integer, intent(in) :: Nx, Ny
                logical, intent(in), optional :: inPlace

                logical :: myInPlace
                integer :: stat
                complex(c_double_complex), dimension(:,:), pointer :: inBuff, outBuff
                

                myInPlace = .False.
                if (present(inplace)) myInPlace = inplace 

                fft_obj%in_place = myInPlace

                fft_obj%fft_array_real_base_ptr = fftw_alloc_complex(int(Nx*Ny,C_SIZE_T))
                if (myInPlace) then 
                        fft_obj%fft_array_recip_base_ptr = fft_obj%fft_array_real_base_ptr        
                else 
                        fft_obj%fft_array_recip_base_ptr = fftw_alloc_complex(int(Nx*Ny,C_SIZE_T))
                end if 
                fft_obj%real_buffer_len = Nx*Ny
                fft_obj%recip_buffer_len = Nx*Ny

                allocate(fft_obj%RealBufferShape(2),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate fft_obj%RealBufferShape"
 
                allocate(fft_obj%RecipBufferShape(2),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate fft_obj%RealBufferShape"               


               
                fft_obj%RealBufferShape = (/Nx, Ny/)
                fft_obj%RecipBufferShape = (/Nx, Ny/)

                call C_F_POINTER(fft_obj%fft_array_real_base_ptr,inBuff,fft_obj%RealBufferShape)
                call C_F_POINTER(fft_obj%fft_array_recip_base_ptr,outBuff,fft_obj%RecipBufferShape)
                 
                fft_obj%plan_forward = fftw_plan_dft_2d(Ny,Nx,inBuff,outBuff,FFTW_FORWARD,FFTW_ESTIMATE)
                fft_obj%plan_backward = fftw_plan_dft_2d(Ny,Nx,outBuff,inBuff,FFTW_BACKWARD,FFTW_ESTIMATE)
        end subroutine create_plan_2d

        subroutine fft_2d_r2c(fft_obj,direction)
                ! This function assumes that the real and complex buffers have been initialised with create_plan_2d_inplace, not
                ! doing this could very easily cause a book keeping error!!!!!!
                use utils, only: utils_to_upper
                implicit none 
                type(fft_object), intent(in) :: fft_obj
                character(len=1) :: direction

                real(kind=c_double), dimension(:,:), pointer :: RealBuffer
                complex(kind=c_double_complex), dimension(:,:), pointer :: ComplexBuffer
                character(len=1) :: tempChar
               
                if (.not. fft_obj%in_place) error stop "Error: fft_2d_r2c can only be used for in place transformations"
                tempChar = direction 
                call utils_to_upper(tempChar)

                call c_f_pointer(fft_obj%fft_array_real_base_ptr,RealBuffer,fft_obj%RealBufferShape)
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,ComplexBuffer,fft_obj%RecipBufferShape)
                if (tempChar == 'F') then 
                        call fftw_execute_dft_r2c(fft_obj%plan_forward,RealBuffer,ComplexBuffer)
                else if (tempChar == 'B') then 
                        call fftw_execute_dft_c2r(fft_obj%plan_backward,ComplexBuffer,RealBuffer)
                else 
                        error stop "Error: Invalid option passed to fft_2d_r2c"
                end if 
                
        end subroutine fft_2d_r2c

        subroutine fft_2d(fft_obj,direction)
                ! This function assumes that the real and complex buffers have been initialised with create_plan_2d_inplace, not
                ! doing this could very easily cause a book keeping error!!!!!!
                use utils, only: utils_to_upper
                implicit none 
                type(fft_object), intent(in) :: fft_obj
                character(len=1) :: direction

               
                complex(kind=c_double_complex), dimension(:,:), pointer :: realBuffer, recipBuffer
                character(len=1) :: tempChar
               
                if (fft_obj%in_place) error stop "Error: fft_2d cannot be used for in place transformations"

                tempChar = direction 
                call utils_to_upper(tempChar)

                call c_f_pointer(fft_obj%fft_array_real_base_ptr,realBuffer,fft_obj%RealBufferShape)
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,recipBuffer,fft_obj%RecipBufferShape)
                if (tempChar == 'F') then 
                        call fftw_execute_dft(fft_obj%plan_forward,realBuffer,recipBuffer)
                else if (tempChar == 'B') then 
                        call fftw_execute_dft(fft_obj%plan_backward,recipBuffer,realBuffer)
                else 
                        error stop "Error: Invalid option passed to fft_2d"
                end if 
                
        end subroutine fft_2d

        subroutine fft_alloc_real(ptr,numElems)
                ! allocate with alignment for SIMD correctness, ptr should be a pointer to the first element of the array
                use iso_fortran_env, only: dp => real64
                use iso_c_binding
                implicit none 
                real(kind=c_double), pointer, intent(inout) :: ptr 
                integer, intent(in) :: numElems
                type(C_ptr) :: tmp

                tmp = fftw_alloc_real(int(numElems,c_size_t)) ! fftw_alloc_real takes integers of type c_size_t
                call c_f_pointer(tmp,ptr)
        end subroutine fft_alloc_real


        subroutine fft_alloc_cmplx(ptr,numElems)
                ! allocate with alignment for SIMD correctness, ptr should be a pointer to the first element of the array 
                use iso_fortran_env, only: dp => real64
                use iso_c_binding
                implicit none 
                complex(kind=c_double_complex), pointer, intent(inout) :: ptr 
                integer, intent(in) :: numElems
                type(C_ptr) :: tmp

                tmp = fftw_alloc_complex(int(numElems,c_size_t)) ! fftw_alloc_complex takes integers of type c_size_t
                call c_f_pointer(tmp,fft_alloc_ptr)
        end subroutine fft_alloc_cmplx
        ! Subroutines for converting between Fortran and C array indexes

        subroutine c2f_indexing_2d(i_c,j_c,i_f,j_f)
                implicit none
                integer, intent(in) :: i_c,j_c
                integer, intent(out) :: i_f, j_f

                i_f = j_c
                j_f = i_f ! Transpose 
        end subroutine c2f_indexing_2d

         subroutine c2f_indexing_3d(i_c,j_c,k_c,i_f,j_f,k_f)
                implicit none 
                integer, intent(in) :: i_c,j_c,k_c
                integer, intent(out) :: i_f, j_f, k_f

                i_f = k_c
                j_f = j_f ! Transpose 
                k_f = i_c
         end subroutine c2f_indexing_3d
        
         

end module fft
