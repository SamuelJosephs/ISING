
module fft
use iso_c_binding
include 'fftw3.f03' 

type fft_object
        type(C_ptr) :: plan_forward, plan_backward
        type(C_ptr) :: fft_array_real_base_ptr ! Points to the first element of the real space array 
        type(C_ptr) :: fft_array_recip_base_ptr ! Points to the first element of the reciprocal space array 
        integer :: real_buffer_len ! Length of the array buffer, proper C programming in Fortran I'm sure there is no risk of stack
                                   ! smashing myself again.
        integer :: recip_buffer_len
        integer, dimension(:), allocatable :: RealBufferShape  ! (shape of the real space array)
        integer, dimension(:), allocatable :: RecipBufferShape ! (shape of the reciprocal space array, can differ from that of the
                                                               ! real space shape for an in place transform)

end type fft_object

public :: fft_2d_r2c, fft_object, create_plan_2d_inplace, c2f_indexing_3d, c2f_indexing_2d 

contains 

        subroutine create_plan_2d_inplace(fft_obj,Nx,Ny,outputBufferReal,outputBufferComplex)
                ! This routine is used to initialise an fft_object such that an in place transformation can be done in outputBuffer
                implicit none

                type(fft_object), intent(inout) :: fft_obj
                integer, intent(in) :: Nx, Ny
                real(kind=c_double), intent(out), pointer, dimension(:,:) :: outputBufferReal ! This is just a reference to the memory
                complex(kind=c_double_complex), intent(out), pointer, dimension(:,:) :: outputBufferComplex

                integer(kind=c_size_t) :: Nx_c, Ny_c
                integer :: stat 
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
                fft_obj%fft_array_recip_base_ptr = fft_obj%fft_array_real_base_ptr
                                                                                               
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
                


        end subroutine create_plan_2d_inplace


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
               
                tempChar = direction 
                call utils_to_upper(tempChar)

                call c_f_pointer(fft_obj%fft_array_real_base_ptr,RealBuffer,fft_obj%RealBufferShape)
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,ComplexBuffer,fft_obj%RecipBufferShape)
                if (tempChar == 'F') then 
                        call fftw_execute_dft_r2c(fft_obj%plan_forward,RealBuffer,ComplexBuffer)
                else if (tempChar == 'B') then 
                        call fftw_execute_dft_c2r(fft_obj%plan_backward,ComplexBuffer,RealBuffer)
                else 
                        error stop "Error: Invalid option passed to fft_2d_realBuff"
                end if 
                
        end subroutine fft_2d_r2c


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
