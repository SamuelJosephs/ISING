module reciprocal_space_processes 
        use chainMesh 
        use vecNd
        use constants
        use, intrinsic :: iso_c_binding 
        use iso_fortran_env, only: dp=>real64
        !include 'fftw3.f03' !This is included in chainMesh, annoying I would like header files with #pragma once
        

        contains
       
        subroutine spectral_interpolate_to_fine(chainMesh)
                ! This routine interpolates chainMesh%fft_obj_std to chainMesh%fft_obj_fine.
                ! This routine assumes that the standard grid has the relevant data already on it.
                use fft, only: fft_object, fft_2d
                use iso_fortran_env, only: dp => real64
                use iso_c_binding 
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh
                ! We need to write the output of the standard grid fft to the fine grid, so will need handles for both 
                complex(kind=c_double_complex), dimension(:,:,:), pointer :: std_recip, fine_recip
                integer :: i,j ! Loop variables
                integer :: scl_fctr_x, scl_fctr_y ! scaling factors, these are equal the number of channels to pad with 0. 
                integer :: std_num_elems_x, std_num_elems_y, fine_num_elems_x, fine_num_elems_y
                ! Associate the pointers
                call c_f_pointer(chainmesh%fft_obj_std%fft_array_recip_base_ptr,&
                                 std_recip,&
                                 chainMesh%fft_obj_std%RecipBufferShape) ! For a 2d transform this is a size 3 array (Nx,Ny,vdim).
                                                                         ! For an N dimensiona interpolation I would have to treat
                                                                         ! the arrays as flattened and manually do the indexing, for
                                                                         ! now this is easier.

                call c_f_pointer(chainmesh%fft_obj_fine%fft_array_recip_base_ptr,&
                                 fine_recip,&
                                 chainMesh%fft_obj_fine%RecipBufferShape)
                std_num_elems_x = chainMesh%fft_obj_std%num_elems_without_padding_recip(1)
                std_num_elems_y = chainMesh%fft_obj_std%num_elems_without_padding_recip(2)

                fine_num_elems_x = chainMesh%fft_obj_fine%num_elems_without_padding_recip(1)
                fine_num_elems_y = chainMesh%fft_obj_fine%num_elems_without_padding_recip(2)

                
                scl_fctr_x = fine_num_elems_x / std_num_elems_x
                scl_fctr_y = fine_num_elems_y / std_num_elems_y
                ! Test that these are exact multiples 
                if (modulo(fine_num_elems_x,std_num_elems_x) /= 0) error stop "Error: Inexact multiples in spectral interpolation &
                                                                                & along x"
                if (modulo(fine_num_elems_y,std_num_elems_y) /= 0) error stop "Error: Inexact multiples in spectral interpolation &
                                                                                & along y"
                if (scl_fctr_x <= 1) error stop "Error: Interpolation scale factor along x must be two or greater"
                if (scl_fctr_y <= 1) error stop "Error: Interpolation scale factor along y must be two or greater"

                ! Now perform the forward transform 
                call fft_2d(chainMesh%fft_obj_std,"F")
                fine_recip(:,:,:) = cmplx(0.0,0.0,c_double_complex) ! Zero the buffer
                do j = 1, chainMesh%fft_obj_std%num_elems_without_padding_recip(2)
                        do i = 1, chainMesh%fft_obj_std%num_elems_without_padding_recip(1)
                                   fine_recip((i-1)*scl_fctr_x + 1,&
                                               (j-1)*scl_fctr_y + 1,:)  = std_recip(i,j,:)     
                        end do 
                end do 
                ! Do the backwards transform 
                call fft_2d(chainMesh%fft_obj_fine,"B")
        end subroutine spectral_interpolate_to_fine

             



        subroutine interpolate_to_fft_array(chainMesh) 
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer :: atom, iCell, jCell, kCell, chainCellIndex, atomIndex, N,L,M , stat, cellIndex, cellIndexTemp  
                integer, dimension(27) :: nearestNeighborCellList
                type(vecNd_t) :: cellCentre, atomPos, d 
                real(kind = 8) :: x,y,z, xWeight, yWeight, zWeight, weight  
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 
                chainMesh%fft_array_x = real(0.0,C_DOUBLE)
                chainMesh%fft_array_y = real(0.0,C_DOUBLE) 
                chainMesh%fft_array_z = real(0.0,C_DOUBLE)
                do CellIndex = 1,chainMesh%numChainMeshCells
                        call getNeighboringCells(chainMesh,CellIndex,nearestNeighborCellList)
                        call coordinatesFromIndex(chainMesh,CellIndeX, iCell, jCell, kCell)
                        ! iCell, jCell, and kCell corrospond to the a,b,c lattice vector coefficients. The interpolation weight
                        ! function will need to be max(1 - a_coeff/a,0.0)
                        iCell = iCell + 1 ! coordinatesFromIndex works for 0 based indexing 
                        jCell = jCell + 1 
                        kCell = kCell + 1 
                        x = chainMesh%chainMeshCells(CellIndex)%centreX
                        y = chainMesh%chainMeshCells(cellIndex)%centreY
                        z = chainMesh%chainMeshCells(cellIndex)%centreZ
                        cellCentre = makeVecNdCheck(cellCentre,[x,y,z])
                        do cellIndexTemp = 1,size(nearestNeighborCellList)
                               atom = chainMesh%chainMeshCells(nearestNeighborCellList(cellIndexTemp))%firstAtomInMeshCell
                               do while (atom /= -1)
                                   x = dble(chainMesh%atoms(atom)%x) 
                                   y = dble(chainMesh%atoms(atom)%y)
                                   z = dble(chainMesh%atoms(atom)%z) 
                                   atomPos = makeVecNdCheck(atomPos,[x,y,z])
                                   call distance_points_vec_bravais(chainMesh,cellCentre,atomPos,d) !Takes periodic boundaries into
                                                                                                    !account and returns the distance in the Bravais lattice basis 
                                   
                                   !Now calculate the cloud in cell interpolation weights 

                                   xWeight = max(1.0_8 - (abs(d%coords(1))),0.0_8)
                                   yWeight = max(1.0_8 - (abs(d%coords(2))),0.0_8) 
                                   zWeight = max(1.0_8 - (abs(d%coords(3))),0.0_8)
                                   weight = xWeight*yWeight*zWeight
                                   chainMesh%fft_array_x(iCell, jCell, kCell) = chainMesh%fft_array_x(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atomSpins(atom,1) 
                                   chainMesh%fft_array_y(iCell, jCell, kCell) = chainMesh%fft_array_y(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atomSpins(atom,2) 
                                   chainMesh%fft_array_z(iCell, jCell, kCell) = chainMesh%fft_array_z(iCell,jCell,kCell) + &
                                                   weight*chainMesh%atomSpins(atom,3) 
 
                                   atom = chainMesh%atoms(atom)%nextAtom 
                               end do 
                        end do 
                        
                end do 
        end subroutine interpolate_to_fft_array 
        
        subroutine interpolate_fft_to_atoms(chainMesh, outputArray)
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), allocatable, dimension(:,:), intent(out) :: outputArray ! i,j in [1,size(chainmesh%atoms)], [1,3]
                integer :: atom, iCell, jCell, kCell, chainCellIndex, atomIndex, N,L,M , stat, cellIndex, cellIndexTemp  
                integer, dimension(27) :: nearestNeighborCellList
                type(vecNd_t) :: cellCentre, atomPos, d 
                real(kind = 8) :: x,y,z, xWeight, yWeight, zWeight,  weight  
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 
                allocate(outputArray(size(chainMesh%atoms),3), stat=stat) 
                outputArray(:,:) = 0.0_8
                if (stat /= 0) error stop "failed to allocate output array"
                do CellIndex = 1,chainMesh%numChainMeshCells
                        call getNeighboringCells(chainMesh,CellIndex,nearestNeighborCellList)
                        call coordinatesFromIndex(chainMesh,CellIndeX, iCell, jCell, kCell)
                        iCell = iCell + 1 ! coordinatesFromIndex works with 0 based indexing 
                        jCell = jCell + 1 
                        kCell = kCell + 1 
                        x = chainMesh%chainMeshCells(CellIndex)%centreX
                        y = chainMesh%chainMeshCells(cellIndex)%centreY
                        z = chainMesh%chainMeshCells(cellIndex)%centreZ
                        cellCentre = makeVecNdCheck(cellCentre,[x,y,z])
                        do cellIndexTemp = 1,size(nearestNeighborCellList)
                               atom = chainMesh%chainMeshCells(nearestNeighborCellList(cellIndexTemp))%firstAtomInMeshCell
                               do while (atom /= -1)
                                   x = dble(chainMesh%atoms(atom)%x) 
                                   y = dble(chainMesh%atoms(atom)%y)
                                   z = dble(chainMesh%atoms(atom)%z) 
                                   atomPos = makeVecNdCheck(atomPos,[x,y,z])
                                   call distance_points_vec_bravais(chainMesh,cellCentre,atomPos,d) !Takes periodic boundaries into account.
                                   
                                   !Now calculate the cloud in cell interpolation weights 

                                   xWeight = max(1.0_8 - (abs(d%coords(1))),0.0_8)
                                   yWeight = max(1.0_8 - (abs(d%coords(2))),0.0_8) 
                                   zWeight = max(1.0_8 - (abs(d%coords(3))),0.0_8)
                                   weight = xWeight*yWeight*zWeight
                                        
                                   outputArray(atom,1) = outputArray(atom,1) +  weight*chainMesh%fft_array_x(iCell,jCell,kCell)
                                   outputArray(atom,2) = outputArray(atom,2) +  weight*chainMesh%fft_array_y(iCell,jCell,kCell)
                                   outputArray(atom,3) = outputArray(atom,3) +  weight*chainMesh%fft_array_z(iCell,jCell,kCell)
                                   atom = chainMesh%atoms(atom)%nextAtom 
                               end do 
                        end do 
                        
                end do        

        end subroutine interpolate_fft_to_atoms 


        subroutine fft_forward_chainMesh(chainMesh)
                type(chainMesh_t), intent(inout) :: chainMesh 

                call fftw_execute_dft_r2c(chainMesh%forwardPlanX, chainMesh%fft_array_x, chainMesh%fft_c_view_x)
                call fftw_execute_dft_r2c(chainMesh%forwardPlanY, chainMesh%fft_array_y, chainMesh%fft_c_view_y)
                call fftw_execute_dft_r2c(chainMesh%forwardPlanZ, chainMesh%fft_array_z, chainMesh%fft_c_view_z)

        end subroutine fft_forward_chainMesh 
        
        subroutine fft_backwards_chainMesh(chainMesh)
                type(chainMesh_t), intent(inout) :: chainMesh 

                call fftw_execute_dft_c2r(chainMesh%backwardPlanX, chainMesh%fft_c_view_x, chainMesh%fft_array_x)
                call fftw_execute_dft_c2r(chainMesh%backwardPlanY, chainMesh%fft_c_view_y, chainMesh%fft_array_y)
                call fftw_execute_dft_c2r(chainMesh%backwardPlanZ, chainMesh%fft_c_view_z, chainMesh%fft_array_z)

        end subroutine fft_backwards_chainMesh


        subroutine calculate_demagnetisation_field(chainMesh, outputArray)
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=8), allocatable, dimension(:,:), intent(out) :: outputArray
                integer :: N,L,M, stat, i, j, k , waveIndexX, waveIndexY, waveIndexZ
                real(kind=C_DOUBLE) :: scaleFactorX, scaleFactorY, scaleFactorZ
                complex(kind = C_DOUBLE_COMPLEX) :: Mx, My, Mz, kdotM, k_squared

                integer :: startClock, endClock, clockRate
                real(kind = C_DOUBLE) :: elapsed_time, sincX, sincY, sincZ, tmp
                type(vecNd_t) :: kx, ky, kz, wg,k_vec
                real(kind=8) :: wg_squared
                call system_clock(startClock, clockRate)



                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 
                M = chainMesh%numcellsZ 
                kx = makeVecNdCheck(kx,[0.0_8,0.0_8,0.0_8])
                ky = makeVecNdCheck(kx,[0.0_8,0.0_8,0.0_8])
                kz = makeVecNdCheck(kx,[0.0_8,0.0_8,0.0_8])
                wg = makeVecNdCheck(wg,[pi/chainMesh%a, pi/chainMesh%b, pi/chainMesh%c])
                wg_squared = wg*wg
                if (allocated(outputArray)) then 
                        if (any(shape(outputArray) /= [chainMesh%numAtoms,3])) then 
                                deallocate(outputArray)
                                allocate(outputArray(chainMesh%numAtoms,3))
                        end if 
                end if 

                call interpolate_to_fft_array(chainMesh)
                call fft_forward_chainMesh(chainMesh)
                ! Now process data using the complex array view into the in place fft array 
                
                do i = 1, N/2 + 1
                        waveIndexX = i-1
                        kx = 2*pi*(dble(waveIndexX)/dble(N))*chainMesh%ar_vec
                        do j = 1,L
                                waveIndexY = j-1
                                if (waveIndexY > L/2) waveIndexY = waveIndexY - L 
                                ky = 2*pi*(dble(waveIndexY)/dble(L))*chainMesh%br_vec
                                do k = 1,M 
                                    waveIndexZ = k-1
                                    if (waveIndexZ > M/2) waveIndexZ = waveIndexZ - M 
                                    kz = 2*pi*(dble(waveIndexZ)/dble(M))*chainMesh%cr_vec
                                    Mx = chainMesh%fft_c_view_x(i,j,k)
                                    My = chainMesh%fft_c_view_y(i,j,k)
                                    Mz = chainMesh%fft_c_view_z(i,j,k)

                                    ! demagnetisation kernel is given by - (k.m / k^2) k
                                    k_vec = kx + ky + kz
                                    kdotM = k_vec%coords(1)*Mx + k_vec%coords(2)*My + k_vec%coords(3)*Mz
                                    k_squared = k_vec*k_vec
                                    if ((i == 1 .and. j == 1 .and. k == 1)) cycle
                                    chainMesh%fft_c_view_x(i,j,k) = - (kdotM / k_squared) * k_vec%coords(1)
                                    chainMesh%fft_c_view_y(i,j,k) = - (kdotM / k_squared) * k_vec%coords(2) 
                                    chainMesh%fft_c_view_z(i,j,k) = - (kdotM / k_squared) * k_vec%coords(3)

                                    if ((abs(k_vec%coords(1)) > pi/chainMesh%a) .or. (abs(k_vec%coords(2)) > pi/chainMesh%b) &
                                                        .or. (abs(k_vec%coords(3)) > pi/chainMesh%c)) then 
                                            chainMesh%fft_c_view_x(i,j,k) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX) ! Impose wg cutoff to prevent aliasing
                                            chainMesh%fft_c_view_y(i,j,k) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX) 
                                            chainMesh%fft_c_view_z(i,j,k) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                                    end if 

                                end do 
                        end do 
                end do 
                chainMesh%fft_c_view_x(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                chainMesh%fft_c_view_y(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                chainMesh%fft_c_view_z(1,1,1) = cmplx(0.0,0.0,C_DOUBLE_COMPLEX)
                call fft_backwards_chainMesh(chainMesh) ! Demagnetising field is now in the three fft_arrays in chainMesh.
                chainMesh%fft_array_x = chainMesh%fft_array_x / (N*L*M)
                chainMesh%fft_array_y = chainMesh%fft_array_y / (N*L*M)
                chainMesh%fft_array_z = chainMesh%fft_array_z / (N*L*M)

                call interpolate_fft_to_atoms(chainMesh,outputArray)
                call system_clock(endClock, clockRate)
                elapsed_time = real(endClock - startClock, C_DOUBLE) / real(clockRate, C_DOUBLE)
        end subroutine calculate_demagnetisation_field


        function arc_winding(s1,s2,s3) result(sigma_area)
                ! taken from https://doi.org/10.1016/0550-3213(81)90568-X
                use constants, only: pi
                implicit none
                type(vecNd_t), intent(in) :: s1, s2, s3
                real(kind=8) :: sigma_area
                complex(kind=8) :: num
                real(kind=dp) :: a, b, c, interiorA, interiorB, interiorC, s, norm
                ! This is the numerically unstable approach, there are issues when s1*s2 + s2*s3 + s3*s1 \approx -1
                num = 1 + s1*s2 + s2*s3 + s3*s1 + cmplx(0.0_8,s1*(s2 .x. s3))
                ! normalisation needed for unitary definition in paper.
                norm = sqrt(2.0_dp*(1.0_dp + s1*s2)*(1.0_dp + s2*s3)*(1.0_dp + s3*s1))
                num = num / norm 
                sigma_area = 2*atan2(aimag(num),real(num))
                ! The numerically stable approach computes the interior angles from the input vectors 
                
                ! We assume that s1, s2, and s3 are unit vectors 
                
                !a = vec_angle(s2,s3)
                !b = vec_angle(s1,s3)
                !c = vec_angle(s1,s2)
                !s = 0.5_dp*(a + b + c) ! semi perimeter

                ! vecSTP(A,B,C) = A\cdot (B\times C), it is the scalar triple product.
                
                
               
                ! Now using L'Hulliers theorem
                !sigma_area = (tan(s/2.0_dp)*tan((s-a)/2.0_dp)*tan((s-b)/2.0_dp)*tan((s-c)/2.0_dp))
                !sigma_area = max(sigma_area,0.0_dp) ! Don't want it to be negative due to the square root and floating point noise
                                                    ! can apparently make it slightly negative.
                !sigma_area = sign(1.0_dp,vecSTP(s1,s2,s3))*4.0_dp*atan(sqrt(sigma_area))
                ! vec(STP) gives the sign of the area

        end function arc_winding

        function calculate_winding_number2(chainMesh,Z_index) result(winding_number)
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: Z_index
                
                integer :: i, j, i_index, j_index, N, L, M, atomIndex, cellIndex 
                integer :: atom2, atom3, atom4
                type(VecNd_t) :: s1, s2, s3, s4
                real(kind=8) :: x1,x2_1,x2_2,x3,y1,y2_1,y2_2,y3,z1,z2_1,z2_2,z3
                real(kind=8) :: sigma1_area1, sigma2_area2, winding_number

                N = chainMesh%numCellsX 
                L = chainMesh%numcellsY 
                M = chainMesh%numCellsZ
                if (Z_index < 1 .or. Z_index > M) error stop "Z_index out of bounds"
                if (.not. allocated(chainMesh%derivativeList)) error stop "Derivative List not initialised"
                if (any(shape(chainMesh%derivativeList) /= [chainMesh%numAtoms,3,2])) &
                                                error stop "DerivativeList not initialised with proper shape"


                winding_number = 0.0_8
                do i = 1, chainMesh%numCellsX
                        do j = 1, chainMesh%numCellsY
                                cellIndex = IndexFromCoordinates(chainMesh,i,j,Z_index)
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                atom2 = chainMesh%derivativeList(atomIndex,1,2)
                                atom4 = chainMesh%derivativeList(atom2,2,2)
                                atom3 = chainMesh%derivativeList(atomIndex,2,2)
                               s1 = dble(chainMesh%atomSpins(atomIndex,:))
                               s2 = dble(chainMesh%atomSpins(atom2,:))
                               s3 = dble(chainMesh%atomSpins(atom3,:))
                               s4 = dble(chainMesh%atomSpins(atom4,:))

                               s1 = s1 / abs(s1)
                               s2 = s2 / abs(s2)
                               s3 = s3 / abs(s3) 
                               s4 = s4 / abs(s4)
                               sigma1_area1 = arc_winding(s1,s2,s4)
                               sigma2_area2 = arc_winding(s1,s4,s3)
                               
                               winding_number = winding_number + sigma1_area1 + sigma2_area2 
                        end do 
                end do 
                winding_number = winding_number / (4*pi) ! Leaving this at the end hoping the compiler will vectorise it
        end function calculate_winding_number2

        subroutine calculate_winding_density_interpolated(chainMesh,density)
                use iso_fortran_env, only: dp => real64
                use iso_c_binding
                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                real(kind=dp), dimension(:,:), allocatable, intent(inout) :: density ! The output winding density 
                real(kind=c_double), dimension(:,:,:), pointer :: fine_real ! The fine and standard grids work with real to complex
                                                                            ! transforms, so real is the appropriate type here (not
                                                                            ! complex). For 2d transforms they have the shape
                                                                            ! (Nx,Ny,vdim).
                integer :: i, j, itemp, jtemp
                integer :: numX, numY
                integer :: stat
                type(VecNd_t) :: s1, s2, s3, s4
                real(kind=dp) :: sigma_area_1, sigma_area_2
                ! First we need the data on the standard grid 

                call map_spins_to_std_grid(chainMesh)

                ! Then we need to interpolate from the standard grid on to the fine grid.

                call spectral_interpolate_to_fine(chainMesh)

                ! In order to work with the interpolated data we will need a handle to work with it in Fortran land. 
                call c_f_pointer(chainMesh%fft_obj_fine%fft_array_real_base_ptr,&
                                 fine_real,&
                                 chainMesh%fft_obj_fine%RealBufferShape)

                numX = chainMesh%fft_obj_fine%num_elems_without_padding_real(1)
                numY = chainMesh%fft_obj_fine%num_elems_without_padding_real(2)
                
                ! Make sure density is allocated with the correct shape.
                if (allocated(density)) then 
                        if (any(shape(density) /= chainMesh%fft_obj_fine%num_elems_without_padding_real)) then 
                                deallocate(density)
                                allocate(density(numX, numY),stat=stat)
                                if (stat /= 0) error stop "Error: Failed to allocate interpolated density array"
                        end if 
                else 
                        allocate(density(numX, numY),stat=stat)
                        if (stat /= 0) error stop "Error: Failed to allocate interpolated density array"
                end if 
                density(:,:) = 0.0_dp
                do j = 1,numY
                        jtemp = j + 1
                        if (jtemp > numY) jtemp = 1 ! Periodicity
                        do i = 1,numX
                                itemp = i + 1
                                if (itemp > numX) itemp = 1 ! periodicity 

                                s1 = fine_real(i,j,:)
                                s2 = fine_real(itemp,j,:)
                                s3 = fine_real(itemp,jtemp,:)
                                s4 = fine_real(i,jtemp,:)

                                s1 = s1 / abs(s1)
                                s2 = s2 / abs(s2)
                                s3 = s3 / abs(s3)
                                s4 = s4 / abs(s4) 
                                sigma_area_1 = arc_winding(s1,s2,s3) ! Anti-Clockwise
                                sigma_area_2 = arc_winding(s1,s3,s4)
                                
                                density(i,j) = sigma_area_1 + sigma_area_2

                        end do 
                end do 

                density = density / (4*pi)
        end subroutine calculate_winding_density_interpolated

        subroutine calculate_winding_number_density(chainMesh,Z_index, density_matrix)
                use iso_fortran_env, only: dp => real64
                use iso_c_binding, only: c_double, c_double_complex 
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: Z_index
                real(kind=dp), dimension(:,:), allocatable, intent(out) :: density_matrix 
                
                integer :: i, j, i_index, j_index, N, L, M, atomIndex, cellIndex, comp 
                integer :: atom2, atom3, atom4
                type(VecNd_t) :: s1, s2, s3, s4, s
                real(kind=dp) :: x1,x2_1,x2_2,x3,y1,y2_1,y2_2,y3,z1,z2_1,z2_2,z3
                real(kind=dp) :: sigma1_area1, sigma2_area2
                integer :: sclx, scly
                


                sclx = chainMesh%sclx
                scly = chainMesh%scly

                N = chainMesh%numCellsX 
                L = chainMesh%numcellsY 
                M = chainMesh%numCellsZ
                if (Z_index < 1 .or. Z_index > M) error stop "Z_index out of bounds"
                if (.not. allocated(chainMesh%derivativeList)) error stop "Derivative List not initialised"
                if (any(shape(chainMesh%derivativeList) /= [chainMesh%numAtoms,3,2])) &
                                                error stop "DerivativeList not initialised with proper shape"
                if (allocated(density_matrix)) then 
                        if (any(shape(density_matrix) /= [N,L])) then 
                                deallocate(density_matrix)
                                allocate(density_matrix(N,L))
                        end if 
                else 
                        allocate(density_matrix(N,L))
                end if
                density_matrix(:,:) = 0.0_8

                ! Need to set up the grid that the winding density will be calculated from 

                do i = 1, chainMesh%numCellsX
                        do j = 1, chainMesh%numCellsY
                                cellIndex = IndexFromCoordinates(chainMesh,i,j,Z_index)
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                atom2 = chainMesh%derivativeList(atomIndex,1,2)
                                atom4 = chainMesh%derivativeList(atom2,2,2)
                                atom3 = chainMesh%derivativeList(atomIndex,2,2)
                               s1 = dble(chainMesh%atomSpins(atomIndex,:))
                               s2 = dble(chainMesh%atomSpins(atom2,:))
                               s3 = dble(chainMesh%atomSpins(atom3,:))
                               s4 = dble(chainMesh%atomSpins(atom4,:))

                               s1 = s1 / abs(s1)
                               s2 = s2 / abs(s2)
                               s3 = s3 / abs(s3) 
                               s4 = s4 / abs(s4)
                               sigma1_area1 = arc_winding(s1,s2,s4)
                               sigma2_area2 = arc_winding(s1,s4,s3)
                               
                               density_matrix(i,j) = sigma1_area1 + sigma2_area2 
                        end do 
                end do 
                density_matrix = density_matrix / (4*pi) ! Leaving this at the end hoping the compiler will vectorise it
        end subroutine calculate_winding_number_density

        subroutine write_winding_number_density(chainMesh, filepath)
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh 
                character(len=*), intent(in) :: filepath 

                real(kind=8), dimension(:,:), allocatable :: density_matrix 
                integer :: i, unit, stat, j, k 
                type(vecNd_t) :: pos
                character(len=400) :: string_buff
                character(len=:), allocatable :: outputString
                character(len=1), parameter :: newline = NEW_LINE('a')

                

                string_buff(1:len(string_buff)) = ' '

                open(newunit=unit, file=filepath, status="replace",action="write",iostat=stat)
                if (stat/=0) error stop "Error opening file to write density information" 

                outputString = "x,y,z,Winding_Density" // NEW_LINE('a')
                do i=1,chainMesh%numcellsZ 
                        call calculate_winding_number_density(chainMesh,i,density_matrix)
                        do j = 1,chainMesh%numCellsX 
                                do k = 1,chainMesh%numCellsY 
                                        pos = (dble(j)*chainMesh%a_vec) + (dble(k)*chainMesh%b_vec) + (dble(i)*chainMesh%c_vec)
                                        string_buff(1:len(string_buff)) = ' '

                                        write(string_buff,'(F0.8, A, F0.8, A, F0.8, A, F0.8, A1)') &
                pos%coords(1), ",", pos%coords(2), ",", pos%coords(3),",", density_matrix(j,k), newline ! Thank you gfortran parser
                                        outputString = outputString // trim(adjustl(string_buff)) // NEW_LINE('a')
                                end do 
                        end do 
                end do 
                write(unit,'(A)') outputString 
                
                close(unit)
        end subroutine write_winding_number_density
        subroutine add_neighbors_to_stack(chainMesh, i,j,stack_array, stack_ptr, visited_array, density_mask, in_stack_array, &
                                                density_matrix)
                implicit none
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: i,j
                integer, dimension(:,:), intent(inout) :: stack_array
                integer, intent(inout) :: stack_ptr 
                logical, dimension(:,:), intent(in) :: visited_array, density_mask 
                logical, dimension(:,:), intent(inout) :: in_stack_array
                real(kind=8), dimension(:,:), intent(in) :: density_matrix
                integer, dimension(2) :: stack_array_shape
                integer :: i_neighbor, j_neighbor, itemp, jtemp, N, L 
                logical :: downhill
                
                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY
                stack_array_shape = shape(stack_array)
                  do i_neighbor = -1,1 
                        do j_neighbor = -1,1 
                                itemp = i + i_neighbor
                                jtemp = j + j_neighbor
                                if (itemp < 1) itemp = N ! Periodic boundaries 
                                if (itemp > N) itemp = 1 
                                
                                if (jtemp < 1) jtemp = L 
                                if (jtemp > L) jtemp = 1 
                                if (in_stack_array(itemp,jtemp)) cycle 
                                downhill = abs(density_matrix(i,j)) >= abs(density_matrix(itemp,jtemp)) - 1e-10_8
                                if ((.not. visited_array(itemp,jtemp)) .and. density_mask(itemp,jtemp) .and. &
                                                                         downhill) then 
                                        if (stack_ptr + 1 > stack_array_shape(1)) error stop "Stack array capacity exceeded"
                                        stack_array(stack_ptr,1) = itemp ! stack_ptr points to where the next entry
                                        ! should be written and equals the existing number of entries + 1.
                                        stack_array(stack_ptr,2) = jtemp 
                                        stack_ptr = stack_ptr + 1
                                        in_stack_array(itemp,jtemp) = .True.
                                        
                                end if 
                        end do 
                  end do
        end subroutine add_neighbors_to_stack




        function calculate_skyrmion_number(chainMesh,Z_index,q_threshold,particle_number, sigma,&
                        sclx, scly) result(skyrmion_number)
                use iso_fortran_env, only: dp=>real64
                implicit none
                type(chainMesh_t), intent(inout) :: chainmesh
                integer, intent(in) :: Z_index, particle_number
                real(kind=8), intent(in) :: q_threshold, sigma
                integer, intent(in), optional :: sclx, scly
                real(kind=8), dimension(:,:), allocatable :: density_matrix
                logical, dimension(:,:), allocatable :: density_mask, visited_array, in_stack_array
                integer :: N,L, stat, skyrmion_number
                
                integer, allocatable, dimension(:,:) :: stack_array ! stack_array(stack_ptr, (i,j))
                integer :: stack_ptr
                integer, parameter :: stack_len = 10000

                integer :: i, j, i_neighbor, j_neighbor, itemp, jtemp, candidate_counter, Nx, Ny, &
                                        xIndex, yIndex
                real(kind=dp) :: acc
                logical :: is_candidate, is_local_maxima
                real(kind=dp) :: upper_threshold

                real(kind=dp), allocatable, dimension(:) :: accs
                integer, allocatable, dimension(:) :: coordsi 
                integer, allocatable, dimension(:) :: coordsj
                logical, allocatable, dimension(:) :: merged
                integer, allocatable, dimension(:,:) :: regionID
                integer :: regionCounter, ID1, ID2
                integer :: my_sclx, my_scly
                real(kind=dp), parameter :: skThreshold = 0.05_dp ! 95% of a skyrmion must be accounted for 


                my_sclx = 1
                my_scly = 1

                if (present(sclx)) my_sclx = sclx 
                if (present(scly)) my_scly = scly 

                N = chainMesh%numCellsX 
                L = chainMesh%numCellsY 

                allocate(density_matrix(N,L), stat=stat)
                if (stat /= 0) error stop "Error allocating topological charge density matrix"

                allocate(density_mask(N,L), stat=stat)
                if (stat/=0) error stop "Error allocating density mask"

                allocate(visited_array(N,L), stat=stat)
                if (stat/=0) error stop "Error allocating visited_array"

                allocate(in_stack_array(N,L),stat=stat)
                if (stat/=0) error stop "Error allocating in_stack_array"
                
                allocate(stack_array(stack_len,2), stat=stat)
                if (stat /= 0) error stop "Error allocating stack array"

                allocate(accs(0))
                allocate(merged(0))
                allocate(coordsi(0))
                allocate(coordsj(0))
                allocate(regionID(N,L))
                regionID = 0
                regionCounter = 0

                stack_ptr = 1 

                if (abs(q_threshold) > 1.0_8 .or. q_threshold < 0.0_8) error stop "q_threshold must be between 0 and 1"

                call calculate_winding_number_density(chainMesh, Z_index, density_matrix)
                !print *, "Minval / Maxval in winding array = ", minval(abs(density_matrix)) / maxval(abs(density_matrix))
                ! density_mask(:,:) = density_matrix > (q_threshold * maxval(abs(density_matrix)))
                call Gaussian_filter_2d(density_matrix,sigma)
                density_mask = abs(density_matrix) > (q_threshold * maxval(abs(density_matrix)))
                visited_array(:,:) = .False.
                in_stack_Array(:,:) = .False.
                stack_array(:,:) = 0 
                skyrmion_number = 0 
                candidate_counter = 0
                is_candidate = .False.
                Nx = chainMesh%numcellsX 
                Ny = chainMesh%numCellsY
                do i = 1,N 
                        do j = 1,L 
                                is_local_maxima = .True.
                                do itemp = -1,1
                                        do jtemp = -1,1
                                                xIndex = modulo(i + itemp - 1,Nx) + 1
                                                yIndex = modulo(j + jtemp - 1,Ny) + 1
                                                if (itemp == 0 .and. jtemp == 0) cycle
                                                if (abs(density_matrix(i,j)) < &
                                                      abs(density_matrix(xIndex,yIndex))) is_local_maxima = .False.
                                        end do 
                                end do 
                                if (density_mask(i,j) .and. (.not. visited_array(i,j)) .and. is_local_maxima) then 
                                acc = density_matrix(i,j)
                                visited_array(i,j) = .True.

                                coordsi = [coordsi,i]
                                coordsj = [coordsj,j]
                                accs = [accs,density_matrix(i,j)]
                                merged = [merged,.False.]  
                                regioncounter = regionCounter + 1
                                regionID(i,j) = regionCounter

 
                                  ! Loop through neighbors, if they are above the theshold and have not been visited or put in the
                                  ! stack then add them to the stack

                                
 
                                  stack_ptr = 1
                                  in_stack_array = .False.
                                  
                                  call add_neighbors_to_stack(chainMesh,i,j,stack_array,stack_ptr, &
                                                 visited_array,density_mask, in_stack_array, density_matrix)
                                  
                                  ! Now for each neighbor in the stack, visit it if not already visited and add it's density to the
                                  ! accumulator. Then add it's unvisited neighbors over the density threshold to the stack, repeat
                                  ! until stack_ptr - 1 == 0
                                  do while (stack_ptr - 1 /= 0)
                                        ! pop cell off of the stack 
                                        itemp = stack_array(stack_ptr - 1,1)
                                        jtemp = stack_array(stack_ptr - 1,2)
                                        stack_ptr = stack_ptr - 1
                                        !print *, "itemp, jtemp = ", itemp, jtemp
                                        if (.not. visited_array(itemp,jtemp)) then

                                                acc = acc + density_matrix(itemp,jtemp)
                                                accs(regionCounter) = accs(regionCounter) + density_matrix(itemp,jtemp)
                                                
                                                visited_array(itemp,jtemp) = .True.
                                                call add_neighbors_to_stack(chainMesh,itemp,jtemp,& 
                                                        stack_array,stack_ptr,visited_array,density_mask, &
                                                                in_stack_array, density_matrix)
                                                regionID(itemp,jtemp) = regionCounter
                                        end if 
                                  end do

                                end if 
                                stack_array = 0
                        end do
                end do 

                skyrmion_number = 0 
                ! Loop to merge regions
                stack_array = 0
                in_stack_Array = .False.
                do i = 1,N 
                        do j = 1,L 
                                if (regionID(i,j) == 0) cycle 
                        
                                ID1 = regionID(i,j)
                                if (merged(ID1)) cycle 

                                acc = accs(ID1) 

                                if (abs(abs(acc) - particle_number) < skThreshold) then 
                                        skyrmion_number = skyrmion_number + nint(sign(1.0_dp,acc)) 
                                        cycle 
                                end if 
                                do itemp = -1,1
                                        do jtemp = -1,1
                                                xIndex = modulo(i + itemp - 1,Nx) + 1
                                                yIndex = modulo(j + jtemp - 1,Ny) + 1
                                                if (itemp == 0 .and. jtemp == 0) cycle
                                                ID2 = regionID(xIndex,yIndex)

                                                if (ID2 == 0) cycle
                                                if (merged(ID2)) cycle 
                                                if (ID2 == 0) cycle

                                                if (abs(abs(accs(ID1) + accs(ID2)) - particle_number) < skThreshold) then 
                                                        !skyrmion_number = skyrmion_number + nint(sign(1.0_dp,accs(ID1) + accs(ID2)))
                                                        ! For now remove mergin logic it is probably far too lenient.
                                                        merged(ID1) = .True.
                                                        merged(ID2) = .True.
                                                end if 
                                        end do 
                                end do 

                        end do 
                end do 
                !call write_2d_real_array_to_file(density_matrix, "./density_matrix.csv")
                !call write_2d_logical_array_to_file(visited_array,"./visited_array.csv")
                !call write_2d_logical_array_to_file(density_mask,"./density_mask.csv")
        end function calculate_skyrmion_number

        subroutine Gaussian_filter_2d(array,sigma)
                real(kind=8), dimension(:,:), intent(inout) :: array 
                real(kind=8), intent(in) :: sigma
                real(kind=8), dimension(:,:), allocatable :: buffer
                real(kind = 8), dimension(:,:), allocatable :: G 
                integer, parameter :: k = 2 ! k functions as a cutoff
                real(kind = 8) :: acc
                integer :: stat, i, Nx, Ny, index_i, index_j, l,m
                acc = 0.0_8 
                allocate(G(2*k + 1,2*k+1),stat = stat)
                if (stat /= 0) error stop "Error allocating G array"

                Nx = size(array,1)
                Ny = size(array,2)

                allocate(buffer(Nx, Ny), stat = stat)
                if (stat /= 0) error stop "Error allocating Buffer"
                do i = 1, 2*k + 1
                        do j = 1, 2*k + 1 
                                G(i,j) = exp(- ((i - k - 1)**2 +(j - k - 1)**2)/(2*sigma**2))
                                acc = acc + G(i,j)
                        end do 
                end do 
                G = G / acc
                buffer = 0.0_8
                do i = 1, Nx
                        do j = 1, Ny 
                                do l = 1, 2*k + 1 
                                        do m = 1, 2*k + 1
                                        ! Going to compute array \conv G 
                                        ! This is given mathematically (for a 1d array) as (array conv G)[i] = sum_n array[i-n] * G[n]
                                        ! We just need to handle periodic boundaries for i - n that are out of bounds 
                                        index_i = i - l 
                                        index_j = j - m 
                                        if (index_i > Nx) index_i = index_i - Nx 
                                        if (index_i < 1) index_i = Nx - abs(index_i)

                                        if (index_j > Ny) index_j = index_j - Ny 
                                        if (index_j < 1) index_j = Ny - abs(index_j)
                                        buffer(i,j) = buffer(i,j) + array(index_i, index_j)*G(l,m)
                                        end do 
                                end do 
                        end do 
                end do 
                array = buffer
        end subroutine Gaussian_filter_2d

        subroutine write_2d_real_array_to_file(array, filename)
                implicit none
                real(kind=8), intent(in), dimension(:,:) :: array 
                character(len=*), intent(in) :: filename 
                
                integer :: N, L, i, j 
                integer, dimension(2) :: shape_array
                character(len=100) :: real_buffer
                character(len=:), allocatable :: buffer
                shape_array = shape(array)
                
                N = shape_array(1)
                L = shape_array(2)
                open(unit=231,file=trim(filename),status="replace",action="write")                 
                allocate(character(len=100) :: buffer)
                real_buffer = ' '
                buffer = ' '


                do i = 1,N 
                        buffer = ' '
                        do j = 1,L
                                real_buffer = ' '
                                write(real_buffer,'(ES12.5)') array(i,j) 
                                buffer = buffer // trim(adjustl(real_buffer))
                                if (j /= L) buffer = buffer // ','
                        end do 
                        write(231,'(A)') buffer
                end do 
                close(231)
        end subroutine write_2d_real_array_to_file


        subroutine write_2d_logical_array_to_file(array, filename)
                implicit none
                logical, intent(in), dimension(:,:) :: array 
                character(len=*), intent(in) :: filename 
                
                integer :: N, L, i, j 
                integer, dimension(2) :: shape_array
                character(len=100) :: real_buffer
                character(len=:), allocatable :: buffer
                shape_array = shape(array)
                
                N = shape_array(1)
                L = shape_array(2)
                
                open(unit=231,file=trim(filename),status="replace",action="write")                 
                allocate(character(len=100) :: buffer)
                real_buffer = ' '
                buffer = ' '


                do i = 1,N 
                        buffer = ' '
                        do j = 1,L 
                                real_buffer = ' '
                                write(real_buffer,'(I5.5)') merge(1,0,array(i,j)) 
                                buffer = buffer // trim(adjustl(real_buffer))
                                if (j /= L) buffer = buffer // ','
                        end do 
                        write(231,'(A)') buffer
                end do 
                close(231)
        end subroutine write_2d_logical_array_to_file




        subroutine compute_skyrmion_distribution(chainMesh,N,winding_array,min_threshold,max_threshold,num_thresholds, Z_index)
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: N ! Will compute the number of particles with winding numbers 1...N
                integer, dimension(:), allocatable, intent(inout) :: winding_array
                real(kind=8), intent(in) :: min_threshold, max_threshold 
                integer,intent(in) :: num_thresholds, Z_index
                integer, dimension(:), allocatable :: closest_skyrmion_number
                integer :: i, j, sigma_index, current_winding_number,stat
                real(kind=8) :: total_charge, threshold, winding, sigma
                
                if (max_threshold < min_threshold) error stop "max_threshold must be greater than min_threshold"
                if (allocated(winding_array)) then 
                        if (size(winding_array) /= N) then 
                                deallocate(winding_array)
                                allocate(winding_array(N),stat=stat)
                        end if 

                else 
                        allocate(winding_array(N),stat=stat)
                end if 
                allocate(closest_skyrmion_number(N),stat=stat)
                if (stat /= 0) error stop "Error allocating skyrmion distribution arrays"

                total_charge = calculate_winding_number2(chainMesh, Z_index)
                winding_array = 0.0_8
                closest_skyrmion_number = 0
                do i = 1,num_thresholds
                        do sigma_index = 1,20
                        ! sigma = (dble(sigma_index)/20.0_8) * (1.0_8 - 0.001_8) + 0.001_8 
                        sigma = 1.0_8 - (dble(sigma_index)/20.0_8)*(1.0_8 - 0.00001_8)
                        winding = 0.0_8
                        current_winding_number = 0
                        !threshold = (dble(i) / dble(num_thresholds)) * (max_threshold - min_threshold) + min_threshold
                        threshold = max_threshold - (dble(i) / dble(num_thresholds)) * (max_threshold - min_threshold)
                        do j = 1, N 
                                winding_array(j) = calculate_skyrmion_number(chainMesh,Z_index,threshold,j,sigma)
                                winding = winding + j * winding_array(j)
                                
                                current_winding_number = current_winding_number + j*closest_skyrmion_number(j)
                        end do 

                        if (abs(nint(winding) - nint(total_charge)) < &
                                abs(current_winding_number - nint(total_charge))) then
                                closest_skyrmion_number = winding_array
                        end if 
                        if (abs(winding - total_charge) < 0.1) then 
                                !print *, "Solution found at q_threshold = ", threshold, "threshold = ", threshold, "sigma = ", sigma
                                return
                        end if 
                        end do
                end do 
                winding_array(:) = closest_skyrmion_number(:) 
        end subroutine compute_skyrmion_distribution
end module reciprocal_space_processes 
