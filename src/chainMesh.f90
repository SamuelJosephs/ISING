module ChainMesh
        use Atom 
        use ChainMeshCell 
        use vecNd 
        use omp_lib
        use constants
        use, intrinsic :: iso_c_binding
        use fft, only: fft_object
        include 'fftw3.f03' ! Mistakes were made, for now the technical debt will have to grow.

        type NNContainer_t
                integer, allocatable, dimension(:) :: NNList ! List of nearest Neighbors 
                real(kind=8) :: distance ! At a specified distance 
        end type NNContainer_t 


        type ChainMesh_t 
                integer :: numAtoms, numChainMeshCells
                type(Atom_t), allocatable :: atoms(:) ! Each atom is in a chain mesh cell and points to the next atom within that cell 
                real(kind=8), allocatable, dimension(:,:) :: atomSpins
                type(ChainMeshCell_t), allocatable :: chainMeshCells(:)
                integer :: numCellsX, numCellsY, numCellsZ  ! NumCellsX/Y/Z will be for the a,b,c Bravais Lattice Vectors
                integer, allocatable, dimension(:,:,:) :: derivativeList !(i,j,k) : i=atomInex, j=dim (1,2,3), k = lower,higher (1,2)
                type(C_ptr) :: forwardPlanX, forwardPlanY, forwardPlanZ, backwardPlanX, backwardPlanY, backwardPlanZ
                real(kind=c_double), pointer :: fft_array_x(:,:,:), fft_array_y(:,:,:), fft_array_z(:,:,:)
                complex(kind=c_double_complex), pointer :: fft_c_view_x(:,:,:), fft_c_view_y(:,:,:), fft_c_view_z(:,:,:)
                type(C_ptr) :: fft_array_ptr
                real(kind=8) :: a, b, c, Bravais_ab, Bravais_bc, Bravais_ca
                type(vecNd_t) :: a_vec, b_vec, c_vec, ar_vec, br_vec, cr_vec ! ar stands for a reciprocal  
                real(kind=8), allocatable, dimension(:,:) :: demagnetisation_array
                
                type(NNContainer_t), allocatable, dimension(:) :: NNNList ! N'th nearest neighbor list OBSOLETE
                type(NNContainer_t), allocatable, dimension(:,:) ::  atomShells ! (atomIndex, shellIndex)


                integer :: sclx, scly ! The scale factors by which the density matrix is being upsampled.
                type(fft_object) :: fft_obj_std, fft_obj_fine ! fine grid for interpolation, standard grid to  
                                                                     
        end type ChainMesh_t 

        interface IndexFromCoordinates
                module procedure IndexFromCoordinatesInteger
        end interface IndexFromCoordinates
        contains


        subroutine getNeighboringCells(chainMesh,cellIndex,neighborCellList)
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(out) :: neighborCellList(27)
                integer, intent(in) :: cellIndex 
                integer :: counter, i,j,k, icoord, jcoord, kcoord, N, tempInt, a,b,c  
                !print *, "DEBIG: Entered getNeighboringCells"
                !print *, "N = ",N
                counter = 1
                a = chainMesh%numCellsX 
                b = chainMesh%numCellsY 
                c = chainMesh%numCellsZ
                !print *, "DEBUG: get Neighboring Cells"
                call coordinatesFromIndex(chainMesh,cellIndex,icoord,jcoord,kcoord)
                 do i = icoord - 1, icoord + 1 
                            do j = jcoord - 1, jcoord + 1 
                                do k = kcoord - 1, kcoord + 1 
                                        ! if i,j,k == N or -1 then wrap around 
                                        itemp = i 
                                        jtemp = j 
                                        ktemp = k 
                                        if (itemp == a) then 
                                                itemp = 0 
                                        else if (itemp == -1) then 
                                                itemp = a-1 
                                        end if 
                                        if (jtemp == b) then 
                                                jtemp = 0 
                                        else if (jtemp == -1) then 
                                                jtemp = b-1 
                                        end if 
                                        if (ktemp == c) then 
                                                ktemp = 0 
                                        else if (ktemp == -1) then 
                                                ktemp = c-1  
                                        end if 
                                        !print *, "cellIndex = ",cellIndex," itemp,jtemp,ktemp = ",itemp,jtemp,ktemp

                                        tempIndex = (itemp)*b*c + (jtemp)*c + (ktemp) + 1
                                        neighborCellList(counter) = tempIndex
                                
                                        counter = counter + 1
                                        
                                end do 
                            end do 
                        end do 

        end subroutine getNeighboringCells

        function distance(chainMesh, atomIndex1, atomIndex2) result(d)
                type(ChainMesh_t), intent(in) :: chainMesh
                integer, intent(in) :: atomIndex1, atomIndex2 
                integer :: ind1,ind2 ! ChainCell Coordinates of atoms 1 and 2 
                integer :: i1,j1,k1,i2,j2,k2 ! Chain Cell i,j,k coordinates of 
                type(Atom_t) :: atom1, atom2 
                real :: d 
                real(kind=8) :: a_coeff, b_coeff, c_coeff, widthX, widthY, widthZ
                real(kind=8), parameter :: eps = 1e-8_8
                
                widthX = dble(chainMesh%numCellsX) ! widthX is the maximum width along the a axis 
                widthY = dble(chainMesh%numCellsY) ! widthY is the maximum width along the b axis 
                widthZ = dble(chainMesh%numCellsZ) ! widthZ is the maximum width aling the c axis

                atom1 = chainMesh%atoms(atomIndex1)
                atom2 = chainMesh%atoms(atomIndex2)
                dx = atom2%x - atom1%x
                dy = atom2%y - atom1%y
                dz = atom2%z - atom1%z
                
                ! Take inner product with respective dual (reciprocal) vectors to extract coefficients in the Bravais basis. 
                a_coeff = dx*chainMesh%ar_vec%coords(1) + dy*chainMesh%ar_vec%coords(2) + dz*chainMesh%ar_vec%coords(3)
                b_coeff = dx*chainMesh%br_vec%coords(1) + dy*chainMesh%br_vec%coords(2) + dz*chainMesh%br_vec%coords(3)
                c_coeff = dx*chainMesh%cr_vec%coords(1) + dy*chainMesh%cr_vec%coords(2) + dz*chainMesh%cr_vec%coords(3)
                
                a_coeff = a_coeff - widthX*nint(a_coeff/widthX)
                b_coeff = b_coeff - widthY*nint(b_coeff/widthY)
                c_coeff = c_coeff - widthZ*nint(c_coeff/widthZ)
                ! Now reconstruct the cartesian distance from these vectors 
                dx = a_coeff*chainMesh%a_vec%coords(1) + b_coeff*chainMesh%b_vec%coords(1) + c_coeff*chainMesh%c_vec%coords(1)
                dy = a_coeff*chainMesh%a_vec%coords(2) + b_coeff*chainMesh%b_vec%coords(2) + c_coeff*chainMesh%c_vec%coords(2)
                dz = a_coeff*chainMesh%a_vec%coords(3) + b_coeff*chainMesh%b_vec%coords(3) + c_coeff*chainMesh%c_vec%coords(3)

                d = sqrt(dx**2 + dy**2 + dz**2)
        end function distance 

        function distance_points(chainMesh, point1, point2) result(d)
            implicit none
            type(ChainMesh_t), intent(in) :: chainMesh
            type(vecNd_t), intent(in) :: point1, point2
            real(kind=8) :: d, dx, dy, dz, widthX, widthY, widthZ
            integer :: a, b, c
            real(kind=8) :: a_coeff, b_coeff, c_coeff
            ! Get system dimensions

                
            widthX = dble(chainMesh%numCellsX) ! widthX is the maximum width along the a axis 
            widthY = dble(chainMesh%numCellsY) ! widthY is the maximum width along the b axis 
            widthZ = dble(chainMesh%numCellsZ) ! widthZ is the maximum width aling the c axis

            dx = point2%coords(1) - point1%coords(1)
            dy = point2%coords(2) - point1%coords(2)
            dz = point2%coords(3) - point1%coords(3)
                
                ! Take inner product with respective dual (reciprocal) vectors to extract coefficients in the Bravais basis. 
            a_coeff = dx*chainMesh%ar_vec%coords(1) + dy*chainMesh%ar_vec%coords(2) + dz*chainMesh%ar_vec%coords(3)
            b_coeff = dx*chainMesh%br_vec%coords(1) + dy*chainMesh%br_vec%coords(2) + dz*chainMesh%br_vec%coords(3)
            c_coeff = dx*chainMesh%cr_vec%coords(1) + dy*chainMesh%cr_vec%coords(2) + dz*chainMesh%cr_vec%coords(3)
                
            a_coeff = a_coeff - widthX*nint(a_coeff/widthX)
            b_coeff = b_coeff - widthY*nint(b_coeff/widthY)
            c_coeff = c_coeff - widthZ*nint(c_coeff/widthZ)
                ! Now reconstruct the cartesian distance from these vectors 
            dx = a_coeff*chainMesh%a_vec%coords(1) + b_coeff*chainMesh%b_vec%coords(1) + c_coeff*chainMesh%c_vec%coords(1)
            dy = a_coeff*chainMesh%a_vec%coords(2) + b_coeff*chainMesh%b_vec%coords(2) + c_coeff*chainMesh%c_vec%coords(2)
            dz = a_coeff*chainMesh%a_vec%coords(3) + b_coeff*chainMesh%b_vec%coords(3) + c_coeff*chainMesh%c_vec%coords(3)

            d = sqrt(dx**2 + dy**2 + dz**2)

        end function distance_points

        subroutine distance_points_vec(chainMesh, point1, point2,d)
            implicit none
            type(ChainMesh_t), intent(in) :: chainMesh
            type(vecNd_t), intent(in) :: point1, point2
            type(vecNd_t), intent(inout) :: d 

            real(kind=8) :: dx, dy, dz, widthX, widthY, widthZ
            real(kind=8) :: a_coeff, b_coeff, c_coeff
            integer :: a, b, c

                
            widthX = dble(chainMesh%numCellsX) ! widthX is the maximum width along the a axis 
            widthY = dble(chainMesh%numCellsY) ! widthY is the maximum width along the b axis 
            widthZ = dble(chainMesh%numCellsZ) ! widthZ is the maximum width aling the c axis

            dx = point2%coords(1) - point1%coords(1)
            dy = point2%coords(2) - point1%coords(2)
            dz = point2%coords(3) - point1%coords(3)
                
                ! Take inner product with respective dual (reciprocal) vectors to extract coefficients in the Bravais basis. 
            a_coeff = dx*chainMesh%ar_vec%coords(1) + dy*chainMesh%ar_vec%coords(2) + dz*chainMesh%ar_vec%coords(3)
            b_coeff = dx*chainMesh%br_vec%coords(1) + dy*chainMesh%br_vec%coords(2) + dz*chainMesh%br_vec%coords(3)
            c_coeff = dx*chainMesh%cr_vec%coords(1) + dy*chainMesh%cr_vec%coords(2) + dz*chainMesh%cr_vec%coords(3)
                
            a_coeff = a_coeff - widthX*nint(a_coeff/widthX)
            b_coeff = b_coeff - widthY*nint(b_coeff/widthY)
            c_coeff = c_coeff - widthZ*nint(c_coeff/widthZ)
                ! Now reconstruct the cartesian distance from these vectors 
            dx = a_coeff*chainMesh%a_vec%coords(1) + b_coeff*chainMesh%b_vec%coords(1) + c_coeff*chainMesh%c_vec%coords(1)
            dy = a_coeff*chainMesh%a_vec%coords(2) + b_coeff*chainMesh%b_vec%coords(2) + c_coeff*chainMesh%c_vec%coords(2)
            dz = a_coeff*chainMesh%a_vec%coords(3) + b_coeff*chainMesh%b_vec%coords(3) + c_coeff*chainMesh%c_vec%coords(3)

            d = makeVecNdCheck(d,[dx,dy,dz]) 


        end subroutine distance_points_vec


        subroutine distance_points_vec_bravais(chainMesh, point1, point2,d)
            implicit none
            type(ChainMesh_t), intent(in) :: chainMesh
            type(vecNd_t), intent(in) :: point1, point2
            type(vecNd_t), intent(inout) :: d 

            real(kind=8) :: dx, dy, dz, widthX, widthY, widthZ
            real(kind=8) :: a_coeff, b_coeff, c_coeff
            integer :: a, b, c
            ! Returns the distance in units of the Bravais lattice vectors.
                
            widthX = dble(chainMesh%numCellsX) ! widthX is the maximum width along the a axis 
            widthY = dble(chainMesh%numCellsY) ! widthY is the maximum width along the b axis 
            widthZ = dble(chainMesh%numCellsZ) ! widthZ is the maximum width aling the c axis

            dx = point2%coords(1) - point1%coords(1)
            dy = point2%coords(2) - point1%coords(2)
            dz = point2%coords(3) - point1%coords(3)
                
                ! Take inner product with respective dual (reciprocal) vectors to extract coefficients in the Bravais basis. 
            a_coeff = dx*chainMesh%ar_vec%coords(1) + dy*chainMesh%ar_vec%coords(2) + dz*chainMesh%ar_vec%coords(3)
            b_coeff = dx*chainMesh%br_vec%coords(1) + dy*chainMesh%br_vec%coords(2) + dz*chainMesh%br_vec%coords(3)
            c_coeff = dx*chainMesh%cr_vec%coords(1) + dy*chainMesh%cr_vec%coords(2) + dz*chainMesh%cr_vec%coords(3)
                
            a_coeff = a_coeff - widthX*nint(a_coeff/widthX)
            b_coeff = b_coeff - widthY*nint(b_coeff/widthY)
            c_coeff = c_coeff - widthZ*nint(c_coeff/widthZ)
                ! Now reconstruct the cartesian distance from these vectors 

            d = makeVecNdCheck(d,[a_coeff,b_coeff,c_coeff]) 


        end subroutine distance_points_vec_bravais 



       subroutine AssignAtomNearestNeighbhors(chainMesh,AtomIndex, AtomCellIndex, NeighborCellList,atomLockArray)
                implicit none
                type(chainMesh_t), intent(inout), target :: chainMesh
                integer, intent(in) :: AtomIndex 
                integer, intent(in) :: AtomCellIndex, NeighborCellList(27)
                integer(kind=OMP_LOCK_KIND), intent(inout) :: atomLockArray(:) 
                type(Atom_t), pointer :: atoms(:)
                type(chainMeshCell_t), pointer :: chainMeshCells(:)
                integer :: cellIndex, NeighborListIter,NumCells, atomIndexTemp  
                type(Atom_t) :: tempAtom,Atom  
                real :: dist, nearestNeighborDistance
                integer, allocatable :: tempList(:)
                
                !print *, "Neighbor Cell list = ", NeighborCellList
                !print *, "DEBUG: AssignAtomNearestNeighbhors"
                NumCells = chainMesh%numChainMeshCells 
                atoms => chainMesh%atoms
                chainMeshCells => chainMesh%chainMeshCells 
                !print *, "Assigned Pointers"
                nearestNeighborDistance = HUGE(0.0)
                ! Get nearest neighbor distance 

                Atom = atoms(AtomIndex)
                ! This loop is looping through all the atoms in the cell and it's neighbors to compute the nearest neighbor distance 
                do cellIndex = 1,size(neighborCellList)
                    ! First atom in the i'th cell in neighbhorCellList
                    atomIndexTemp = chainMeshCells(neighborCellList(cellIndex))%firstAtomInMeshCell
                    do while (atomIndexTemp .ne. -1) 
                        
                        if (atomIndexTemp == AtomIndex) then ! atom skips itself 
                                atomIndexTemp = atoms(atomIndexTemp)%nextAtom
                        
                                cycle
                        end if
                        tempAtom = atoms(atomIndexTemp)
                         dist = distance(chainMesh,AtomIndex,AtomIndexTemp)
                        if (dist .lt. nearestNeighborDistance) then 
                                nearestNeighborDistance = dist 
                        end if 
                        atomIndexTemp = atoms(atomIndexTemp)%nextAtom
                    end do 
                end do
                ! Now actually compute the nearest neighbhor list 
                do cellIndex = 1,size(neighborCellList)
                    ! First atom in the i'th cell in neighbhorCellList
                    atomIndexTemp = chainMeshCells(neighborCellList(cellIndex))%firstAtomInMeshCell
                    do while (atomIndexTemp .ne. -1)
                        if (atomIndexTemp == AtomIndex) then 
                                atomIndexTemp = atoms(atomIndexTemp)%nextAtom ! We don't want to consider the atom we are assigning
                                                                              ! nearest neighbors to 
                               cycle 
                        else if (any(atoms(AtomIndex)%NeighborList == atomIndexTemp)) then 
                                atomIndexTemp = atoms(atomIndexTemp)%nextAtom ! Don't want to consider an atom if it is already in
                                                                              ! nearest neighbors list 
                               cycle 
                       end if
                        tempAtom = atoms(atomIndexTemp)

                        dist = distance(chainMesh,AtomIndex, atomIndexTemp)
                        if (abs(dist - nearestNeighborDistance) < chainMesh%a/100.0_8) then ! Tolerance 
                                call omp_set_lock(atomLockArray(atomIndex))

                                allocate(tempList((size(atoms(AtomIndex)%NeighborList) + 1)))
                                tempList(1:size(atoms(AtomIndex)%NeighborList)) = atoms(AtomIndex)%NeighborList 
                                tempList(size(atoms(AtomIndex)%NeighborList)+1) = atomIndextemp   
                                deallocate(atoms(AtomIndex)%NeighborList)
                                allocate(atoms(AtomIndex)%NeighborList(size(tempList)))
                                atoms(AtomIndex)%NeighborList = tempList
                                !print *, "New Neighbor List: ", atoms(atomIndexTemp)%NeighborList
                                deallocate(tempList)

                                call omp_unset_lock(atomLockArray(atomIndex))
                               !atoms(atomIndexTemp)%NeighborList = [atoms(atomIndexTemp)%NeighborList, &
                                !        atomIndexTemp]
                        end if 
                        atomIndexTemp = atoms(atomIndexTemp)%nextAtom
                    end do 
                end do 
        end subroutine AssignAtomNearestNeighbhors

        subroutine assignNearestNeighbors(chainMesh)
                implicit none
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer :: N, i, j ,k , numchainMeshCells, AtomIdent, icoord,jcoord,kcoord,idex, counter 
                integer :: neighborCellList(27), tempIndex, itemp, jtemp, ktemp, idexTemp   
                integer(kind=OMP_LOCK_KIND), allocatable :: atomLockArray(:) 
                allocate(atomLockArray(size(chainMesh%atoms)))
                do i=1,size(atomLockArray)
                    call omp_init_lock(atomLockArray(i))
                end do 
                !print *, "Entered Assign Nearest Neighbhors"
                !print *, "Assigned N"
                numChainMeshCells = chainMesh%numChainMeshCells
                !print *, "Assigned numChainMeshCells = ", numChainMeshCells
                ! For each chain mesh cell, for each atom in that cell, compute nearest neighbors and     
                ! check neighbor chain mesh cells for nearest neighbor atoms, taking account the
                ! periodicity of the system
                !$omp parallel do shared(chainMesh,numChainMeshCells,atomLockArray) private(idex,AtomIdent,NeighborCellList) &
                !$omp& default(firstprivate)
                do idex = 1,numChainMeshCells
                        !print *, "DEBUG: idex = ",idex
                        AtomIdent = chainMesh%chainMeshCells(idex)%firstAtomInMeshCell
                        ! Compute nearest neighbors for this atom 
                        ! Handle neighboring cells
                        !print *, "DEBUG: about to call getNeighboringCells"
                        call getNeighboringCells(chainMesh,idex,neighborCellList) 
                        !print *, "DEBIG: Succsessfully called getNeighboringCells"
                        ! For each atom in cell Idex, test get it's nearest neighbor list based on 
                        !it's cell and all of the neighboring cells 
                        ! For each atom determine the nearest neighbor list 
                        
                        do while (AtomIdent .ne. -1 )
                                !$omp critical 
                                call AssignAtomNearestNeighbhors(chainMesh,AtomIdent,idex,neighborCellList,atomLockArray) 
                                !$omp end critical
                                AtomIdent = chainMesh%atoms(AtomIdent)%NextAtom
                        end do
                        
                        !print *, "DEBUG: idex = ", idex
                end do 
                !$omp end parallel do 
                !print *, "DEBUG: succsessfully assigned nearest neighbors"
                do i=1,size(atomLockArray)
                    call omp_destroy_lock(atomLockArray(i))
                end do        
        end subroutine assignNearestNeighbors


        subroutine initNNNList(chainMesh,N) ! This routine is suspect and will be replaced with initAtomShells
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: N 

                integer :: i, atomIndex, cellIndex, numAtomsInCell, cellIndexTemp, atomIndexTemp
                integer :: distanceCounter, j
                integer, allocatable, dimension(:) :: distanceLoc
                real(kind=8) :: currentDistance, distance_val, prevDistance
                type(NNContainer_t), allocatable, dimension(:,:) :: NNNArray
                integer, allocatable :: NNeighboringCells(:)
                real(kind=8), allocatable, dimension(:) :: distanceArray
                if (allocated(chainMesh%NNNList)) then 
                        deallocate(chainMesh%NNNList)
                end if 
                
                allocate(chainMesh%NNNList(0))
                
                allocate(distanceArray(N))
                distanceArray(:) = 0.0_8

                allocate(NNNArray(chainMesh%numAtoms,N))
                do i = 1,N 
                        do j = 1, chainMesh%numAtoms
                                allocate(NNNArray(j,i)%NNList(N))
                        end do 
                end do 
                do cellIndex = 1,chainMesh%numChainMeshCells 
                                numAtomsIncell = chainMesh%chainMeshCells(cellIndex)%NumAtomsPerUnitCell
                                ! Compute N'th neighboring cell indexes.
                                call NNearestCells(chainMesh,cellIndex,N,NNeighboringCells)
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell


                                do while (atomIndex /= -1)  

                                        prevDistance = -1.0_8
                                        do distanceCounter = 1,N  

                                                currentDistance = HUGE(currentDistance)
                                                ! For each atom in Cell 
                                                ! For each potential neighbor cell
                                                do i = 1,size(NNeighboringCells)
                                                        atomIndexTemp = chainMesh%chainMeshCells(&
                                                                NNeighboringCells(i))%firstAtomInMeshCell
                                                        ! For each potential neighbor atom in neighbouring cells  
                                                        do while (atomIndexTemp /= -1)
                                                                if (atomIndexTemp == atomIndex) cycle
                                                                distance_val = distance(chainMesh,atomIndex,atomIndexTemp)  

                                                                if (distance_val < currentDistance .and. &
                                                                         distance_val  > prevDistance) currentDistance = &
                                                                                        distance_val 
                                                                atomIndexTemp = chainMesh%atoms(atomIndexTemp)%nextAtom
                                                        end do 
                                                end do 
                                                prevDistance = currentDistance
                                                distanceArray(distanceCounter) = currentDistance 
                                        end do 
                                        ! Now that we have the ranks of distances for the current atom we can do another pass
                                        ! through it's neighbors and assign them into bins by their distance 
                                         
                                        do i = 1,size(NNeighboringCells)
                                                atomIndexTemp = chainMesh%chainMeshCells(&
                                                        NNeighboringCells(i))%firstAtomInMeshCell
                                                ! For each potential neighbor atom in neighbouring cells  
                                                do while (atomIndexTemp /= -1)
                                                        if (atomIndexTemp == atomIndex) cycle
                                                        distance_val = distance(chainMesh,atomIndex,atomIndexTemp)  
                                                        
                                                        distanceLoc = findloc(abs(distanceArray - distance_val) < 1e-5,.True.)
                                                        if (size(distanceLoc) < 1) error stop "Error: Cannot find distance entry &
                                                                & in array"
                                                        NNNArray(atomIndex,distanceLoc(1))%NNList = &
                                                                [NNNArray(atomIndex,distanceLoc(1))%NNList,atomIndexTemp]
                                                        atomIndexTemp = chainMesh%atoms(atomIndexTemp)%nextAtom
                                                end do 
                                        end do 

                                        do i = 1,N 
                                                NNNArray(atomIndex,i)%distance = distanceArray(i)
                                        end do 
                                        atomIndex = chainMesh%atoms(atomIndex)%nextAtom
                                end do 
                end do 
        end subroutine initNNNList 

        subroutine initAtomShells(chainMesh,N,NNeighbourCells,distance_threshold) ! Initialise N shells around each atom 
                use iso_fortran_env, only: dp=>real64
                use algo, only: mergesort  

                implicit none 
                type(chainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: N, NNeighbourCells 
                real(kind=dp), intent(in) :: distance_threshold

                integer :: i, j, atomIndex, atomIndexTemp,  cellIndex, NCellIndex, stat, &
                        NCAtomIndex
                real(kind=dp), allocatable, dimension(:) :: distanceArray 
                integer, allocatable, dimension(:) :: distanceArrayAtomIndices
                integer :: distanceArrayIndex, shellIndex
                real(kind=dp) :: dist, prevDist 

                integer, dimension(:), allocatable :: NCellList

                ! Initialise chainMesh%atomShells 

                if (allocated(chainMesh%atomShells)) then 
                    do i = 1, chainMesh%numAtoms
                        do j = 1, N
                                deallocate(chainMesh%atomShells(i,j)%NNList) ! Initial allocation
                        end do
                    end do
                    deallocate(chainMesh%atomShells)
                end if 

                allocate(chainMesh%atomShells(chainMesh%numAtoms,N),stat=stat)
                if (stat /=0 ) error stop "Error: Failed to allocated atomShells array"
                do i = 1, chainMesh%numAtoms 
                    do j = 1, N
                        allocate(chainMesh%atomShells(i,j)%NNList(0)) ! Initial allocation 
                        chainMesh%atomShells(i,j)%distance = 0.0_dp 
                    end do 
                end do 

                allocate(distanceArray((2*NNeighbourCells+1)**3*chainMesh%chainMeshCells(1)%numAtomsPerUnitCell),stat=stat)
                if (stat /= 0) error stop "Error: failed to allocate distanceArray"
                allocate(distanceArrayAtomIndices((2*NNeighbourCells+1)**3*chainMesh%chainMeshCells(1)%numAtomsPerUnitCell),&
                        stat=stat)
                if (stat /= 0) error stop "Error: failed to allocate distanceArray"
                               
                ! Now start iterating through each atom, each atoms neighbors in N neibouring cells, and recording the distances in
                ! an array 
                
                do cellIndex = 1, chainMesh%numchainMeshCells
                        call NNearestCells(chainMesh,cellIndex,NNeighbourCells,NCellList)
                        ! At this stage the atoms should have been assigned to each cell 
                        atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell

                        do while (atomIndex /= -1)
                                distanceArrayIndex = 1
                                distanceArray(:) = HUGE(distanceArray(1)) 
                                distanceArrayAtomIndices(:) = -1                               
                                do i = 1, size(NCellList)
                                        NCellIndex = NCellList(i) 
                                        NCAtomIndex = chainMesh%chainMeshCells(NCellIndex)%firstAtomInMeshCell 
                                        do while (NCAtomIndex /= -1)
                                                if (NCAtomIndex == atomIndex) then 
                                                        NCAtomIndex = chainMesh%atoms(NCAtomIndex)%nextAtom 
                                                        cycle 
                                                end if 
                                                dist = distance(chainMesh,atomIndex,NCAtomIndex) 
                                                distanceArray(distanceArrayIndex) = dist 
                                                distanceArrayAtomIndices(distanceArrayIndex) = NCAtomIndex 
                                                distanceArrayIndex = distanceArrayIndex + 1
                                                NCAtomIndex = chainMesh%atoms(NCAtomIndex)%nextAtom 
                                        end do 
                                end do         
                                ! Now for the atom given by atomIndex, we have a set of distances, we can sort these and compute
                                ! "shells" 
                                ! Naively there should be no excess elements but we will put in a check just in case 

                                call mergesort(distanceArray,integer_companion=distanceArrayAtomIndices)

                                shellIndex = 1
                                chainMesh%atomShells(atomIndex,shellIndex)%NNList = [chainMesh%atomShells(atomIndex,&
                                        shellIndex)%NNList,&
                                        distanceArrayAtomIndices(1)] ! Add the closest atom to the first shell 
                                do j = 2,size(distanceArray)
                                        dist=distanceArray(j)
                                        prevDist=distanceArray(j - 1)
                                        if (abs(dist - prevDist) > distance_threshold) shellIndex = shellIndex + 1 
                                        if (shellIndex > N) exit
                                         chainMesh%atomShells(atomIndex,shellIndex)%NNList = &
                                                    [chainMesh%atomShells(atomIndex,shellIndex)%NNList,&
                                                                distanceArrayAtomIndices(j)]
                                end do 
                                
                                atomIndex = chainMesh%atoms(atomIndex)%nextAtom

                        end do 
                                
                end do 
        end subroutine initAtomShells 

        subroutine NNearestCells(chainMesh,cellIndex,N,list) ! Computed the cell indices of the N nearest cells and put them in list 
                implicit none
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: cellIndex, N 
                integer, allocatable, dimension(:), intent(inout) :: list 

                integer :: aCoord,bCoord,cCoord, atemp,btemp,ctemp
                integer :: i,j,k,cellIndexOut, counter, stat
                if (allocated(list)) then 
                        if (size(list) /= (2*N+1)**3) then 
                                deallocate(list)
                                allocate(list((2*N+1)**3), stat=stat)
                                if (stat /= 0) error stop "Error: Failed to allocate list array"
                        end if 
                else 
                        allocate(list((2*N+1)**3),stat=stat)
                        if (stat /= 0) error stop "Error: Failed to allocate list array"
                end if 
                
                call coordinatesFromIndex(chainMesh,cellIndex,aCoord,bCoord,cCoord)

                aCoord = aCoord  ! 1 - numCellsX 
                bCoord = bCoord 
                cCoord = cCoord 
                counter = 1
                do i = -N,N 
                        do j = -N,N 
                                do k = -N,N 
                                        atemp = aCoord + i 
                                        btemp = bCoord + j 
                                        ctemp = cCoord + k 

                                        atemp = modulo(atemp,chainMesh%numCellsX) + 1 
                                        btemp = modulo(btemp,chainMesh%numCellsY) + 1 
                                        ctemp = modulo(ctemp,chainMesh%numCellsZ) + 1 
                                        
                                        cellIndexOut = IndexFromCoordinates(chainmesh, atemp, btemp, ctemp)
                                        list(counter) = cellIndexOut 
                                        counter = counter + 1
                                end do 
                        end do 
                end do 
        end subroutine NNearestCells  

        subroutine enumerateChainMeshCells(chainMesh)
                type(chainMesh_t), intent(inout), target :: chainMesh 
                type(ChainMeshCell_t), pointer :: chainMeshCells(:)
                type(Atom_t), pointer :: atoms(:)
                integer :: atomIdent, i, atomCounter, idex ,N 
                chainMeshCells => chainMesh%chainMeshCells
                atoms => chainMesh%atoms
                atomCounter = 0 
                do i = 1,size(chainMeshCells)
                        print *, "###################################"
                        atomIdent = chainMeshCells(i)%firstAtomInMeshCell
                        print *, "MeshCell ",i," contains atoms:"
                        do while (atomIdent .ne. -1)
                                print *, atomIdent,&
                                       " With coordinates: ", "(",atoms(atomIdent)%x,",",&
                                      atoms(atomIdent)%y,",",atoms(atomIdent)%z,")"
                                print *, "tmpx, tmpy, tmpz = ",atoms(atomIdent)%tmpx,atoms(atomIdent)%tmpy,atoms(atomIDent)%tmpz
                                print *, "Atom ", atomIdent," has nearest Neighbors: ", atoms(atomIdent)%NeighborList
                                print *, "Atom ", atomIdent," has ",size(atoms(atomIdent)%NeighborList)," nearest Neighbors: " 
                                atomIdent = atoms(atomIdent)%nextAtom
                                atomCounter = atomCounter + 1 
                        end do 
                end do
                print *, "Mesh cells have a combined ",atomCounter," Atoms"
        end subroutine enumerateChainMeshCells

        subroutine addAtomToChainCell(CellIndex,AtomIndex,chainMesh)
                integer, intent(in) :: CellIndex, AtomIndex
                type(ChainMesh_t), intent(inout) :: chainMesh
                type(ChainMeshCell_t) :: tempChainCell
                integer :: temp 
                temp = chainMesh%chainMeshCells(CellIndex)%firstAtomInMeshCell
                chainMesh%chainMeshCells(CellIndex)%firstAtomInMeshCell=AtomIndex 
                chainMesh%atoms(AtomIndex)%NextAtom=temp
        end subroutine addAtomToChainCell 

        function IndexFromCoordinatesInteger(chainMesh,i,j,k) result(res)
                implicit none
                type(ChainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: i,j,k 
                integer :: res 
                integer :: Amesh, Bmesh, Cmesh 
                Amesh = chainMesh%numCellsX
                Bmesh = chainMesh%numCellsY 
                Cmesh = chainMesh%numCellsZ 
                
                ! i,j, and k will be from loops with 1 based indexing, e.g. do i = 1,size(array)
                res = (i-1)*Bmesh*Cmesh + (j-1)*Cmesh + (k-1) + 1 ! +1 to work with fortran 1 based array indexing
        end function IndexFromCoordinatesInteger

        subroutine coordinatesFromIndex(chainMesh,Ind,i,j,k)
                implicit none
                type(ChainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: Ind
                integer, intent(out) :: i,j,k
                integer :: remainderj, IndTemp, reconstructedIndex, a, b, c, tempInt 
                a = chainMesh%numCellsX 
                b = chainMesh%numCellsY 
                c = chainMesh%numCellsZ 
                tempInt = b*c 
                IndTemp = Ind-1
                i = IndTemp/tempInt  
                remainderj = mod(IndTemp,tempInt)
                j = remainderj / c 
                k = mod(remainderj,c)
                reconstructedIndex = i*b*c + j*b + k + 1 ! Needed for 1 based indexing in Fortran 
        end subroutine coordinatesFromIndex

        subroutine create_chainMesh_plan(chainMesh)
                type(chainmesh_t), intent(inout) :: chainMesh 
                type(C_ptr) :: plan 
                
                integer :: N,L,M, status, plannerThreads
                
                N = chainMesh%numCellsX
                L = chainMesh%numCellsY 
                M = chainMesh%numCellsZ 
                

                call fftw_plan_with_nthreads(omp_get_max_threads())
                chainMesh%forwardPlanX = fftw_plan_dft_r2c_3d(M,L,N,chainMesh%fft_array_x,chainMesh%fft_c_view_x,FFTW_ESTIMATE)
                chainMesh%backwardPlanX = fftw_plan_dft_c2r_3d(M,L,N,chainMesh%fft_c_view_x,chainMesh%fft_array_x,FFTW_ESTIMATE)


                chainMesh%forwardPlanY = fftw_plan_dft_r2c_3d(M,L,N,chainMesh%fft_array_y,chainMesh%fft_c_view_y,FFTW_ESTIMATE)
                chainMesh%backwardPlanY = fftw_plan_dft_c2r_3d(M,L,N,chainMesh%fft_c_view_y,chainMesh%fft_array_y,FFTW_ESTIMATE)


                chainMesh%forwardPlanZ = fftw_plan_dft_r2c_3d(M,L,N,chainMesh%fft_array_z,chainMesh%fft_c_view_z,FFTW_ESTIMATE)
                chainMesh%backwardPlanZ = fftw_plan_dft_c2r_3d(M,L,N,chainMesh%fft_c_view_z,chainMesh%fft_array_z,FFTW_ESTIMATE)
        end subroutine create_chainMesh_plan

        function makeChainMesh(numCellsX, numCellsY, numCellsZ, & 
                        AtomsInUnitCell, &
                        a,b,c,ab_deg,bc_deg,ca_deg,distance_threshold,numShells,debug,&
                        sclx,scly) result(chainMesh)
                use algo, only: mergesort
                use fft, only: fft_object, create_plan_2d_many, create_plan_2d_r2c_many
                use iso_fortran_env, only: dp => real64
                implicit none 
                integer, intent(in) ::  numCellsX, numCellsY, numCellsZ 
                type(Atom_t), intent(in) :: AtomsInUnitCell(:)
                real(kind=8), intent(in) :: a,b,c,ab_deg,bc_deg,ca_deg
                real(kind=8), intent(in), optional :: distance_threshold ! For Neighbour atom shells 
                ! Need to create all the atoms, create all the ChainMeshCells, then allocate all of the atoms to a chain mesh cell
                ! default = 0.1*a where a is the first Bravais lattice dimension length
                integer, intent(in), optional :: numShells ! the number of shells we wish to compute, default = 2
                logical, intent(in), optional :: debug
                integer, intent(in), optional :: sclx, scly
                

                integer :: my_sclx, my_scly
                type(ChainMesh_t), target :: chainMesh 
                integer :: numChainMeshCells, padX, stat, numAtomsPerUnitCell
                integer :: numAtoms, stride 
                integer :: i,j,k, icoord, jcoord, kcoord, my_numShells
                type(ChainMeshCell_t) :: tempChainMeshCell
                type(Atom_t) :: tempAtoms(size(AtomsInUnitCell))
                real :: domainWidth
                type(C_ptr) :: temp_c_ptr
                type(vecNd_t) :: a_vec, b_vec, c_vec, ar_vec, br_vec, cr_vec, tmp_vec, &
                                        pos
                real(kind=8) :: cx, cy, cz, ab, bc,ca,&
                       my_threshold 
                integer, allocatable, dimension(:) :: intBuffer
                logical :: my_debug 

                my_sclx = 1
                my_scly = 1
                if (present(sclx)) my_sclx = sclx
                if (present(scly)) my_scly = scly

                chainMesh%sclx = my_sclx
                chainMesh%scly = my_scly
                ! Initialise fft_objects
                call create_plan_2d_r2c_many(chainMesh%fft_obj_std,numCellsX,numCellsY,3,inplace=.False.) 
                call create_plan_2d_r2c_many(chainMesh%fft_obj_fine,&
                                             my_sclx*numCellsX,&
                                             my_scly*numCellsY,&
                                             3,inplace = .False.)

                my_debug = .False.
                if (present(debug)) my_debug = debug 
                my_threshold = 0.1_dp*a 
                my_numShells = 2
                if (present(distance_threshold)) my_threshold = distance_threshold   
                if (present(numShells)) my_numShells = numShells  
                if (my_numShells < 1) error stop "Error in makeChainMesh: numShells must be greater than or equal to 1"
                numAtomsPerUnitCell = size(AtomsInUnitCell)
                ab = ab_deg*(pi/180.0_dp)
                bc = bc_deg*(pi/180.0_dp)
                ca = ca_deg*(pi/180.0_dp)
                chainMesh%a = a 
                chainMesh%b = b 
                chainMesh%c = c 
                chainMesh%Bravais_ab = ab 
                chainMesh%Bravais_bc = bc 
                chainMesh%Bravais_ca = ca 
                a_vec = makeVecNd([a,0.0_8,0.0_8])
                b_vec = makeVecNd([b*cos(ab),b*sin(ab),0.0_8])
                cx = c*cos(ca)
                cy = (c*b*cos(bc)-cx*b_vec%coords(1))
                cz = sqrt(c**2 - cx**2 - cy**2)
                c_vec = makeVecNd([cx,cy,cz])
                chainMesh%a_vec = a_vec
                chainMesh%b_vec = b_vec
                chainMesh%c_vec = c_vec 
                

                print *, "a_vec = ", chainMesh%a_vec%coords
                print *, "b_vec = ", chainMesh%b_vec%coords
                print *, "c_vec = ", chainMesh%c_vec%coords
                print *, "numCellsX = ", numCellsX, "numCellsY = ", numCellsY, "numCellsZ = ", numCellsZ
                chainMesh%a_vec%coords = merge(chainMesh%a_vec%coords, 0.0_8,abs(chainMesh%a_vec%coords) > 1e-5_8 )
                chainMesh%b_vec%coords = merge(chainMesh%b_vec%coords, 0.0_8,abs(chainMesh%b_vec%coords) > 1e-5_8 )
                chainMesh%c_vec%coords = merge(chainMesh%c_vec%coords, 0.0_8,abs(chainMesh%c_vec%coords) > 1e-5_8 )
                !chainMesh%ar_vec%coords = merge(ar_vec%coords, 0.0_8,abs(ar_vec%coords) > 1e-5 )
                !chainMesh%br_vec%coords = merge(br_vec%coords, 0.0_8,abs(br_vec%coords) > 1e-5 )
                !chainMesh%cr_vec%coords = merge(cr_vec%coords, 0.0_8,abs(cr_vec%coords) > 1e-5 )
                a_vec = chainMesh%a_vec
                b_vec = chainMesh%b_vec 
                c_vec = chainMesh%c_vec 

                chainMesh%ar_vec = (b_vec .x. c_vec)/(a_vec*(b_vec .x. c_vec)) ! reciprocal lattice vectors, form a dual basis to  
                chainMesh%br_vec = (c_vec .x. a_vec)/(b_vec*(c_vec .x. a_vec)) ! the Bravais lattice vectors.
                chainMesh%cr_vec = (a_vec .x. b_vec)/(c_vec*(a_vec .x. b_vec))

                chainMesh%numCellsX = numCellsX 
                chainMesh%numCellsY = numCellsY 
                chainMesh%numCellsZ = numCellsZ
                padX = 2*(numCellsX / 2 + 1)
                !allocate(chainMesh%fft_array(numCellsX, numCellsY, numCellsZ),stat=stat)
                chainMesh%fft_array_ptr = fftw_alloc_real(int(2*(numCellsX/2 + 1)*numCellsY*numCellsZ,C_SIZE_T))
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_array_x,[padX,numCellsY,numCellsZ])
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_c_view_x,[numCellsX/2 + 1,numCellsY,numCellsZ])

                chainMesh%fft_array_ptr = fftw_alloc_real(int(2*(numCellsX/2 + 1)*numCellsY*numCellsZ,C_SIZE_T))
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_array_y,[padX,numCellsY,numCellsZ])
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_c_view_y,[numCellsX/2+1,numCellsY,numCellsZ])

                chainMesh%fft_array_ptr = fftw_alloc_real(int(2*(numCellsX/2 + 1)*numCellsY*numCellsZ,C_SIZE_T))
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_array_z,[padX,numCellsY,numCellsZ])
                call c_f_pointer(chainMesh%fft_array_ptr,chainMesh%fft_c_view_z,[numCellsX/2+1,numCellsY,numCellsZ])

                call create_chainMesh_plan(chainMesh)
                numChainMeshCells = numCellsX*numCellsY*numCellsZ 
                numAtoms = numAtomsPerUnitCell*numChainMeshCells
                chainMesh%numAtoms = numAtoms
                chainMesh%numChainMeshCells = numChainMeshCells
                allocate(chainMesh%chainMeshCells(numChainMeshCells))
                allocate(chainMesh%atoms(numAtoms))
                stride = size(AtomsInUnitCell)
                ! initialise atoms and chain Mesh Cells 
                do j = 1,numChainMeshCells ! do it for each unit cell 
                        tempAtoms = AtomsInUnitCell
                        tempChainMeshCell%NumAtomsPerUnitCell = NumAtomsPerUnitCell
                        tempChainMeshCell%firstAtomInMeshCell = -1 ! Initialise to garbage value so that we don't have to deal with
                        ! it
                        chainMesh%chainMeshCells(j) = tempChainMeshCell ! Just so that it isn't initialised to garbage  
             
                        
                        call coordinatesFromIndex(chainMesh,j,icoord,jcoord,kcoord) ! icoord, jcoord, kcoord integer
                        tmp_vec = chainMesh%a_vec*(dble(icoord)+0.5_8) + chainMesh%b_vec*(dble(jcoord)+0.5_8) & 
                                                + chainMesh%c_vec*(dble(kcoord) + 0.5_8)

                        chainMesh%chainMeshCells(j)%centreX = tmp_vec%coords(1)
                        chainMesh%chainMeshCells(j)%centreY = tmp_vec%coords(2)
                        chainMesh%chainMeshCells(j)%centreZ = tmp_vec%coords(3)
                        tmp_vec = chainMesh%a_vec*(dble(icoord)) + chainMesh%b_vec*(dble(jcoord)) & 
                                                + chainMesh%c_vec*(dble(kcoord))
                        do i=1,size(AtomsInUnitCell)
                                pos = dble([AtomsInUnitCell(i)%x, AtomsInUnitCell(i)%y, AtomsInUnitCell(i)%z])
                                pos = pos%coords(1)*chainMesh%a_vec + pos%coords(2)*chainMesh%b_vec &
                                        + pos%coords(3)*chainMesh%c_vec
                                tempAtoms(i)%x = pos%coords(1) + tmp_vec%coords(1) 
                                tempAtoms(i)%y = pos%coords(2) + tmp_vec%coords(2)
                                tempAtoms(i)%z = pos%coords(3) + tmp_vec%coords(3) 
                                ! Assign Atoms to Unit Cells 
                                chainMesh%atoms((j-1)*stride + i) = tempAtoms(i) ! This line confuses me, I don't know what I was
                                                                                 ! thinking when I wrote it.
                                call addAtomToChainCell(j,(j-1)*stride + i,chainMesh)
                        end do
               end do

               call assignNearestNeighbors(chainMesh) ! TODO: This should be depreciated in favour of initAtomShells 

               allocate(chainMesh%demagnetisation_array(chainMesh%numAtoms,3),stat=stat)
               if (stat /= 0) error stop "Error: Failed to allocate demagnetisation_array"
               chainMesh%demagnetisation_array = 0.0_8
               allocate(chainMesh%atomSpins(chainMesh%numAtoms,3),stat=stat)
               if (stat /= 0) error stop "Error: Failed to allocated atomSpins array"
               call DerivativeList(chainMesh,chainMesh%derivativeList)        
               call initAtomShells(chainMesh,my_numShells,my_numShells,my_threshold) 

               ! DEBUG: Check that initAtomShells first shell is identical to assignNearestNeighbors calculation 
               if (my_debug) then 
                       do i = 1,chainMesh%numAtoms
                                do j = 1,my_numShells
                                        call mergesort(chainMesh%atomShells(i,j)%NNList)         
                                end do 
                                call mergesort(chainMesh%atoms(i)%NeighborList)
                                if (all(chainMesh%atoms(i)%NeighborList == chainMesh%atomShells(i,1)%NNList)) then 
                                        print *, "Atom ", i ," Has correct Neighbours"
                                else 
                                        print *, "ERROR: atom ", i, " has inconsistent neighbours initAtoms gives :",&
                                                chainMesh%atomShells(i,1)%NNList, " The old method gives :",&
                                                         chainMesh%atoms(i)%NeighborList
                                end if 
                        
                       end do 
               end if 
        end function makeChainMesh 


        subroutine deallocateChainMesh(chainMesh)
                type(ChainMesh_t), intent(inout) :: chainMesh
                integer :: i,j 

                do i=1,size(chainMesh%atoms)
                        if (allocated(chainMesh%atoms(i)%NeighborList)) then 
                                deallocate(chainMesh%atoms(i)%NeighborList)
                        end if 
                end do
                if (allocated(chainMesh%atoms)) then 
                        deallocate(chainMesh%atoms)
                end if
                if (allocated(chainMesh%chainMeshCells)) then 
                        deallocate(chainMesh%chainMeshCells)
                end if 
        end subroutine deallocateChainMesh

        subroutine DerivativeList(chainMesh,array)
                implicit none
                type(chainMesh_t), intent(inout) :: chainMesh
                integer, allocatable, intent(inout) :: array(:,:,:) ! Output Array 
                
                integer :: N, d, numChainMeshCells, cellIndex, atomIndex, lowerAtom, HigherAtom   
                if (allocated(array)) deallocate(array)
                if (.not. allocated(chainMesh%atoms)) error stop "ChainMesh not properly initialised"
                N = size(chainMesh%atoms)
                allocate(array(N,3,2)) ! (atomIndex, dim, (low index, high index))
                numChainMeshCells = chainMesh%numChainMeshCells

                do d = 1,3 ! For each dimension  
                        do cellIndex = 1,numChainMeshCells
                                atomIndex = chainMesh%chainMeshCells(cellIndex)%firstAtomInMeshCell
                                do while (atomIndex /= -1)
                                        call computeAdjacentAtoms(chainMesh,d,cellIndex,atomIndex,lowerAtom,HigherAtom)
                                        array(AtomIndex,d,1) = lowerAtom 
                                        array(AtomIndex,d,2) = HigherAtom
                                        atomIndex = chainMesh%atoms(atomIndex)%nextAtom
                                end do 
                        end do 
                end do 
        end subroutine DerivativeList 

        subroutine computeAdjacentAtoms(chainMesh,d,cellIndex,atomIndex,lowerAtom,HigherAtom)
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(in) :: d, cellIndex, atomIndex 
                integer, intent(out) :: lowerAtom, HigherAtom 
                
                integer, dimension(27) :: neighborCellList
                integer :: i, j, cellIndexTemp, atomIndexTemp, currentMinIndex, currentMaxIndex  
                type(vecNd_t) :: distance, atomPos1, atomPos2 
                real(kind=8) :: x,y,z, currentMin, currentMax
                logical :: candidate 
                currentMinIndex = -1
                currentMaxIndex = -1
                currentMin = HUGE(currentMin) ! This is a confusing naming scheme, but we need to find the closest atoms with
                                              !coordinate in dimension d that is greater and less than the current atoms.
                currentMax = HUGE(currentMax)
                call getNeighboringCells(chainMesh,cellIndex,neighborCellList) 
                x = chainMesh%atoms(atomIndex)%x 
                y = chainMesh%atoms(atomIndex)%y 
                z = chainMesh%atoms(atomIndex)%z 
                atomPos1 = makeVecNdCheck(atomPos1, [x,y,z])
                do i = 1,size(neighborCellList)
                        cellIndexTemp = neighborCellList(i)
                        atomIndexTemp = chainMesh%chainMeshCells(cellIndexTemp)%firstAtomInMeshCell

                        if (cellIndexTemp == cellIndex) cycle ! Make sure that we only compute adjacent atoms outside of our current
                                                              ! cell. 

                        do while (atomIndexTemp /= -1)
                                if (atomIndexTemp == atomIndex) then 
                                        atomIndexTemp = chainMesh%atoms(atomIndexTemp)%nextAtom
                                        cycle 
                                end if
                                x = chainMesh%atoms(atomIndexTemp)%x 
                                y = chainMesh%atoms(atomIndexTemp)%y 
                                z = chainMesh%atoms(atomIndexTemp)%z 
                                atomPos2 = makeVecNdCheck(atomPos2, [x,y,z])
                                call distance_points_vec_bravais(chainMesh,atomPos1, atomPos2, distance) 
                                ! d contains the coefficients in the expansion d = a_coeff*a + b_coeff*b + c_coeff*c, we have to
                                ! define planes relative to the Bravais lattice.
                                candidate = .True.
                                do j = 1,3
                                        if (j == d) cycle 
                                        if (abs(distance%coords(j)) > 0.0001) then 
                                                candidate = .False.
                                                exit 
                                        end if 
                                end do 

                                if (candidate) then 

                                        !print *, "Candidate found! , distance between atom is given by: ", distance%coords, "d = ",d
                                        if (distance%coords(d) < 0.0_8) then ! Figure out we put it into the bin for the closest
                                                                             !atom with position less than atomIndex's
                                                if (abs(distance%coords(d)) <= currentMin) then 
                                                        currentMin = abs(distance%coords(d))
                                                        currentMinIndex = atomIndexTemp 
                                                end if 

                                        end if 

                                        if (distance%coords(d) > 0.0_8) then ! Figure out if we put it into the bin for closest atom
                                                                            !with position greater than atomIndex's
                                                if (abs(distance%coords(d)) <= currentMax) then 
                                                        currentMax = abs(distance%coords(d))
                                                        currentMaxIndex = atomIndexTemp 
                                                end if 
                                        end if 
                                end if 
                                atomIndexTemp = chainMesh%atoms(atomIndexTemp)%nextAtom
                                
                        end do 
                end do 
                if ((currentMinIndex == -1) .or. (currentMaxIndex == -1)) error stop "Unable to locate atoms"
                lowerAtom = currentMinIndex 
                HigherAtom = currentMAxIndex 
        end subroutine computeAdjacentAtoms


end module ChainMesh 

