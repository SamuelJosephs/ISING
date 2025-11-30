
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
        integer, dimension(:), allocatable :: num_elems_without_padding_real  ! This array will store the number of elements in the
                                                                              ! real buffer that are meaningful. Only this one
                                                                              ! should be used for normalising transforms. 
        integer, dimension(:), allocatable :: num_elems_without_padding_recip ! This array will store the number of elements in the
                                                                              ! recip buffer that are meaningful. 
        logical :: in_place ! Needed for cleanup
        logical :: is_r2c ! Easier book keeping, if .True. then the real space array is of type real(kind=c_double), else the real
                          ! space array is of type complex(kind=c_double_complex). The reciprocal space array is always of type
                          ! complex(kind=c_double_complex).


end type fft_object

public ::  fft_2d, fft_object, create_plan_2d_r2c_many, create_plan_2d_many, c2f_indexing_3d, c2f_indexing_2d, fft_alloc_real 

contains 

        subroutine create_plan_2d_r2c_many(fft_obj, Nx, Ny, vdim, inplace, zero_arrays)
                ! use this routine for vectors 
                ! The resulting arrays will have shape (Nx,Ny,vdim), the (:,:,i) subarrays will be contiguous.
                implicit none 
                type(fft_object), intent(inout) :: fft_obj
                integer, intent(in) :: Nx, Ny, vdim ! vdim is the dimension of the vector field we wish to transform. 
                logical, optional, intent(in) :: inplace
                logical, optional, intent(in) :: zero_arrays
                ! For a three dimensional vector (like spin) rank = 3 

                integer :: numX_real, numY_real, numX_recip, numY_recip
                logical :: myInPlace, my_zero_arrays
                integer, dimension(2) :: n_real, n_recip, real_embed, recip_embed
                integer :: real_stride, recip_stride, real_dist, recip_dist
                real(kind=c_double), dimension(:), pointer :: real_ptr
                complex(kind=c_double_complex), dimension(:), pointer :: recip_ptr
                integer :: stat


                my_zero_arrays = .True.
                if (present(zero_arrays)) my_zero_arrays = zero_arrays

                myInPlace = .False.
                if (present(inplace)) myinplace = inplace 

                numX_real = 2*(Nx/2 + 1) ! Amount of real entries required for there to be enough space to take advantage of the r2c
                                         ! symmetry.
                numY_real = Ny
                
                numX_recip = Nx/2 + 1
                numY_recip = Ny

                ! Now do the allocations 

                fft_obj%fft_array_real_base_ptr = fftw_alloc_real(int(numX_real*numY_real*vdim,c_size_t)) ! Multiply by vdim to
                                                                                                           ! allocate enough memory for all of the batched FFT's.
                if (.not. c_associated(fft_obj%fft_array_real_base_ptr)) error stop "Error:&
                                & failed to allocate real base ptr"
                if (myInPlace) then
                        fft_obj%fft_array_recip_base_ptr = fft_obj%fft_array_real_base_ptr
                else 
                        fft_obj%fft_array_recip_base_ptr = fftw_alloc_complex(int(numX_recip*numY_recip*vdim,c_size_t))
                        if (.not. c_associated(fft_obj%fft_array_recip_base_ptr)) error stop "Error:&
                                & Failed to allocate recip base ptr"
                end if 



                ! Now do the Forward plan, we will need the TRANSPOSED shape,
                n_real = (/numY_real, numX_real/)
                n_recip = (/numY_recip, numX_recip/)

                real_embed = n_real
                recip_embed = n_recip ! No subarrays, so the shape of the sub arrays are that of the entire input 
                ! take in an array of shape (Nx,Ny,rank), the stride will be Nx*Ny
                real_stride = 1
                recip_stride = 1 ! Arrays are contiguous, we could have stride = rank if we packed the array as (rank, Nx, Ny) instead.

                real_dist = numX_real*numY_real ! The size of each subarray we wish to transform 
                recip_dist = numX_recip*numY_recip 
                allocate(fft_obj%RealBufferShape(3),stat = stat)
                if (stat /= 0 ) error stop "Failed to allocate RealBufferShape"

                allocate(fft_obj%RecipBufferShape(3),stat = stat)
                if (stat /= 0 ) error stop "Failed to allocate RecipBufferShape"

                ! These are for array usage in Fortran so column major is fine.
                fft_obj%RealBufferShape = (/numX_real, numY_real, vdim/)
                fft_obj%RecipBufferShape = (/numX_recip, numY_recip, vdim/)

                allocate(fft_obj%num_elems_without_padding_real(2),stat=stat)
                if (stat /= 0) error stop "Failed to allocate array"

                allocate(fft_obj%num_elems_without_padding_recip(2),stat=stat)
                if (stat /= 0) error stop "Failed to allocate array"

                fft_obj%num_elems_without_padding_real = (/ Nx, Ny /) ! This will be used for array normalisation.
                fft_obj%num_elems_without_padding_recip = (/ Nx/2 + 1, Ny /) ! Only have up to the Nyquist frequency in the x

                fft_obj%real_buffer_len = product(fft_obj%RealBufferShape)
                fft_obj%recip_buffer_len = product(fft_obj%RecipBufferShape)
                fft_obj%in_place = myInPlace
                fft_obj%is_r2c = .True. ! This is a r2c plan creation helper                                                                            
                ! Bind the in_forward and out_forward fortran pointers to the actual memory addresses,
                ! don't need shapes they are made up, in fact they are so made up we can flatten them.
                call c_f_pointer(fft_obj%fft_array_real_base_ptr,real_ptr,(/fft_obj%real_buffer_len/))
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,recip_ptr,(/fft_obj%recip_buffer_len/))

                if (my_zero_arrays) then 
                        real_ptr(:) = 0.0_c_double 
                        if (.not. myInPlace) recip_ptr(:) = cmplx(0.0,0.0,c_double_complex)
                end if 

                !TODO: Check everything here is the correct way round.
                ! 2 dimensional trasnform so n = 2
                ! In the Forward transform:
                ! n is the shape of the transform so that will be the row major in_forward (real space dimension)
                ! howmany: We are doing one transform for each dimension of the input vector, which we have called vdim
                ! in: The input array, which for the real space data we call in_forward 
                ! inembed: No subarrays so this is just the (row major) shape of each whole transform 
                ! istride: We define the arrays of each transform to be contiguous so istride = 1.
                ! idist: The size of the stride for each input trasnform (How far in memory do we need to jump to reach the base of
                ! address of the next transform, as each trasnform is contiguous we will only need the products of the shapes of
                ! each transform
                !
                ! backward_embed: again no padding so just the shape of the reciprical space transform 
                fft_obj%plan_forward = fftw_plan_many_dft_r2c(2,n_real,vdim,real_ptr,real_embed,real_stride,real_dist,&
                        recip_ptr,recip_embed,recip_stride,recip_dist,FFTW_ESTIMATE)
                ! For the output trasnform n = n_real (not n_recip) because fftw knows to cut the last row major length in half
                ! (+1). TLDR: FFTW always takes the full real size as the size of the transform.
                fft_obj%plan_backward = fftw_plan_many_dft_c2r(2,n_real,vdim,recip_ptr,recip_embed,recip_stride,recip_dist,&
                        real_ptr,real_embed,real_stride,real_dist,FFTW_ESTIMATE)

                
                if (.not. c_associated(fft_obj%plan_forward)) error stop "Error: Failed to create forward plan"
                if (.not. c_associated(fft_obj%plan_backward)) error stop "Error: Failed to create backward plan"
        end subroutine create_plan_2d_r2c_many


        subroutine create_plan_2d_many(fft_obj, Nx, Ny, vdim, inplace,zero_arrays)
                ! use this routine for vectors 
                implicit none 
                type(fft_object), intent(inout) :: fft_obj
                integer, intent(in) :: Nx, Ny, vdim ! vdim is the dimension of the vector field we wish to transform. 
                logical, optional, intent(in) :: inplace
                logical, optional, intent(in) :: zero_arrays
                ! For a three dimensional vector (like spin) rank = 3 

                integer :: numX_real, numY_real, numX_recip, numY_recip
                logical :: myInPlace, my_zero_arrays
                integer, dimension(2) :: n_real, n_recip, real_embed, recip_embed
                integer :: real_stride, recip_stride, real_dist, recip_dist
                complex(kind=c_double_complex), dimension(:), pointer :: real_ptr, recip_ptr
                integer :: stat
           

                my_zero_arrays = .True.
                if (present(zero_arrays)) my_zero_arrays = zero_arrays

                myInPlace = .False.
                if (present(inplace)) myinplace = inplace 

                numX_real = Nx
                numY_real = Ny
                
                numX_recip = Nx
                numY_recip = Ny

                ! Now do the allocations 

                fft_obj%fft_array_real_base_ptr = fftw_alloc_complex(int(numX_real*numY_real*vdim,c_size_t)) ! Multiply by vdim to
                                                                                                             ! allocate enough memory for all of the batched FFT's.
                if (.not. c_associated(fft_obj%fft_array_real_base_ptr)) error stop "Error: Failed to allocate&
                                & real base ptr"
                if (myInPlace) then
                        fft_obj%fft_array_recip_base_ptr = fft_obj%fft_array_real_base_ptr
                else 
                        fft_obj%fft_array_recip_base_ptr = fftw_alloc_complex(int(numX_recip*numY_recip*vdim,c_size_t))
                        if (.not. c_associated(fft_obj%fft_array_recip_base_ptr)) error stop "Error: Failed to allocate&
                                & recip base ptr"
                end if 



                ! Now do the Forward plan, we will need the TRANSPOSED shape,
                n_real = (/numY_real, numX_real/)
                n_recip = (/numY_recip, numX_recip/)

                real_embed = n_real
                recip_embed = n_recip ! No subarrays, so the shape of the sub arrays are that of the entire input 
                ! take in an array of shape (Nx,Ny,rank), the stride will be Nx*Ny
                real_stride = 1
                recip_stride = 1 ! Arrays are contiguous, we could have stride = rank if we packed the array as (rank, Nx, Ny) instead.

                real_dist = numX_real*numY_real ! The size of each subarray we wish to transform 
                recip_dist = numX_recip*numY_recip 
                allocate(fft_obj%RealBufferShape(3),stat = stat)
                if (stat /= 0 ) error stop "Failed to allocate RealBufferShape"

                allocate(fft_obj%RecipBufferShape(3),stat = stat)
                if (stat /= 0 ) error stop "Failed to allocate RecipBufferShape"

                ! These are for array usage in Fortran so column major is fine.
                fft_obj%RealBufferShape = (/numX_real, numY_real, vdim/)
                fft_obj%RecipBufferShape = (/numX_recip, numY_recip, vdim/)

                allocate(fft_obj%num_elems_without_padding_real(2),stat=stat)
                if (stat /= 0) error stop "Failed to allocate array"

                allocate(fft_obj%num_elems_without_padding_recip(2),stat=stat)
                if (stat /= 0) error stop "Failed to allocate array"

                fft_obj%num_elems_without_padding_real = (/ Nx, Ny /) ! This will be used for array normalisation.
                fft_obj%num_elems_without_padding_recip = (/ Nx, Ny /) 
                ! Bind the in_forward and out_forward fortran pointers to the actual memory addresses,
                ! don't need shapes they are made up, in fact they are so made up we can flatten them 
                call c_f_pointer(fft_obj%fft_array_real_base_ptr,real_ptr,(/fft_obj%RealBufferShape/))
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,recip_ptr,(/fft_obj%RecipBufferShape/))

                if (my_zero_arrays) then
                        real_ptr(:) = cmplx(0.0,0.0,c_double_complex)
                        ! If in place then we don't need to zero the same memory twice.
                        if (.not. myInPlace) recip_ptr(:) = cmplx(0.0,0.0,c_double_complex)
                end if 
                fft_obj%real_buffer_len = product(fft_obj%RealBufferShape)
                fft_obj%recip_buffer_len = product(fft_obj%RecipBufferShape)
                fft_obj%in_place = myInPlace
                fft_obj%is_r2c = .False. 
                !TODO: Check everything here is the correct way round.
                ! 2 dimensional trasnform so n = 2
                ! In the Forward transform:
                ! n is the shape of the transform so that will be the row major in_forward (real space dimension)
                ! howmany: We are doing one transform for each dimension of the input vector, which we have called vdim
                ! in: The input array, which for the real space data we call in_forward 
                ! inembed: No subarrays so this is just the (row major) shape of each whole transform 
                ! istride: We define the arrays of each transform to be contiguous so istride = 1.
                ! idist: The size of the stride for each input trasnform (How far in memory do we need to jump to reach the base of
                ! address of the next transform, as each trasnform is contiguous we will only need the products of the shapes of
                ! each transform
                !
                ! backward_embed: again no padding so just the shape of the reciprical space transform 
                fft_obj%plan_forward = fftw_plan_many_dft(2,n_real,vdim,real_ptr,real_embed,real_stride,real_dist,&
                        recip_ptr,recip_embed,recip_stride,recip_dist,FFTW_FORWARD,FFTW_ESTIMATE)
                ! For the output trasnform n = n_real (not n_recip) because fftw knows to cut the last row major length in half
                ! (+1). TLDR: FFTW always takes the full real size as the size of the transform.
                fft_obj%plan_backward = fftw_plan_many_dft(2,n_real,vdim,real_ptr,real_embed,real_stride,real_dist,&
                        recip_ptr,recip_embed,recip_stride,recip_dist,FFTW_BACKWARD,FFTW_ESTIMATE)
                
                if (.not. c_associated(fft_obj%plan_forward)) error stop "Error: Failed to create forward plan"
                if (.not. c_associated(fft_obj%plan_backward)) error stop "Error: Failed to create backward plan"

        end subroutine create_plan_2d_many


        subroutine fft_2d(fft_obj,direction)
                type(fft_object), intent(inout) :: fft_obj
                character(len=1) :: direction 


                if (fft_obj%is_r2c) then 
                        call fft_2d_r2c_wrapped(fft_obj,direction)
                else 
                        call fft_2d_wrapped(fft_obj,direction)
                endif 


        end subroutine fft_2d

        subroutine fft_2d_r2c_wrapped(fft_obj,direction)
                ! This function assumes that the real and complex buffers have been initialised with create_plan_2d_inplace, not
                ! doing this could very easily cause a book keeping error!!!!!!
                use utils, only: utils_to_upper
                implicit none 
                type(fft_object), intent(inout) :: fft_obj
                character(len=1) :: direction

                real(kind=c_double), dimension(:,:,:), pointer :: RealBuffer ! (Nx, Ny, vdim)
                complex(kind=c_double_complex), dimension(:,:,:), pointer :: ComplexBuffer
                character(len=1) :: tempChar
               
                tempChar = direction 
                call utils_to_upper(tempChar)

                call c_f_pointer(fft_obj%fft_array_real_base_ptr,RealBuffer,fft_obj%RealBufferShape)
                call c_f_pointer(fft_obj%fft_array_recip_base_ptr,ComplexBuffer,fft_obj%RecipBufferShape)
                if (tempChar == 'F') then 
                        call fftw_execute_dft_r2c(fft_obj%plan_forward,RealBuffer,ComplexBuffer)
                else if (tempChar == 'B') then 
                        call fftw_execute_dft_c2r(fft_obj%plan_backward,ComplexBuffer,RealBuffer)
                        RealBuffer = RealBuffer / product(fft_obj%num_elems_without_padding_real)
                else 
                        error stop "Error: Invalid option passed to fft_2d_r2c"
                end if 


                
                
        end subroutine fft_2d_r2c_wrapped

        subroutine fft_2d_wrapped(fft_obj,direction)
                ! This function assumes that the real and complex buffers have been initialised with create_plan_2d_inplace, not
                ! doing this could very easily cause a book keeping error!!!!!!
                use utils, only: utils_to_upper
                implicit none 
                type(fft_object), intent(inout) :: fft_obj
                character(len=1) :: direction

               
                complex(kind=c_double_complex), dimension(:,:,:), pointer :: realBuffer, recipBuffer
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
                        realBuffer = realBuffer / product(fft_obj%num_elems_without_padding_real)
                else 
                        error stop "Error: Invalid option passed to fft_2d"
                end if 
                
        end subroutine fft_2d_wrapped

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
                call c_f_pointer(tmp,ptr)
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
