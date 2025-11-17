
module fft
use iso_c_binding
include 'fftw3.f03' 

type fft_object
        type(C_ptr) :: plan_forward, plan_backward
        real(kind=c_double), pointer :: fft_array 
        type(C_ptr) :: fft_array_base_ptr ! Points to the first element of the array 
        integer :: buffer_len ! Length of the array buffer, proper C programming in Fortran I'm sure there is no risk of stack
                              ! smashing myself again.


end type fft_object


contains 

        subroutine create_plan_2d_inplace(fft_obj,Nx,Ny,outputBufferReal,outputBufferComplex)
                ! This routine is used to initialise an fft_object such that an in place transformation can be done in outputBuffer
                implicit none

                type(fft_object), intent(inout) :: fft_obj
                integer, intent(in) :: Nx, Ny
                real(kind=c_double), intent(out), pointer, dimension(:,:) :: outputBufferReal ! This is just a reference to the memory
                complex(kind=c_double_complex), intent(out), pointer, dimension(:,:) :: outputBufferComplex

                integer(kind=c_size_t) :: Nx_c, Ny_c
                
                Nx_c = int(Nx,c_size_t)
                Ny_c = int(Ny,c_size_t)

                ! We only need Nx/2 complex numbers in the transformed complex space, then add 1 and multiply by 2 to make sure
                ! there is enough space for the transformed buffer
                fft_obj%fft_array_base_ptr = fftw_alloc_real((2*(Nx_c/2 + 1))*Ny_c) 
                                                                                               
                call C_F_POINTER(fft_obj%fft_array_base_ptr,outputBufferReal,[(2*(Nx/2 + 1)),Ny])

                call C_F_POINTER(fft_obj%fft_array_base_ptr,outputBufferComplex,[Nx,Ny])
                 ! Now outputBuffer will behave like an ordinary 2d array, 
                outputBufferReal(:,:) = 0.0_c_double ! Zero the array to avoid uninitialised data later.

                fft_obj%buffer_len = product(shape(outputBufferReal))

                ! Now we need to initialise the plans 

                fft_obj%plan_forward = fftw_plan_dft_r2c_2d(Ny,Nx,outputBufferReal,outputBufferComplex,FFTW_ESTIMATE)
                fft_obj%plan_backward = fftw_plan_dft_c2r_2d(Ny,Nx,outputBufferComplex,outputBufferReal,FFTW_ESTIMATE)
                


        end subroutine create_plan_2d_inplace


end module fft
