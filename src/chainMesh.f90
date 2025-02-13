

module ChainMesh
        use Atom 
        use ChainMeshCell 
        use vecNd 
        use omp_lib
type ChainMesh_t 
        integer :: numAtoms, numChainMeshCells, numChainMeshCellsPerSide
        type(Atom_t), allocatable :: atoms(:) ! Each atom is in a chain mesh cell and points to the next atom within that cell 
        type(ChainMeshCell_t), allocatable :: chainMeshCells(:)
        real :: domainWidth, latticeParameter
end type ChainMesh_t 

        contains
        subroutine getNeighboringCells(chainMesh,cellIndex,neighborCellList)
                type(chainMesh_t), intent(in) :: chainMesh 
                integer, intent(out) :: neighborCellList(27)
                integer, intent(in) :: cellIndex 
                integer :: counter, i,j,k, icoord, jcoord, kcoord, N, tempInt  
                !print *, "DEBIG: Entered getNeighboringCells"
                N = chainMesh%numChainMeshCellsPerSide
                !print *, "N = ",N
                counter = 1
                !print *, "DEBUG: get Neighboring Cells"
                call coordinatesFromIndex(cellIndex,N,icoord,jcoord,kcoord)
                 do i = icoord - 1, icoord + 1 
                            do j = jcoord - 1, jcoord + 1 
                                do k = kcoord - 1, kcoord + 1 
                                        ! if i,j,k == N or -1 then wrap around 
                                        itemp = i 
                                        jtemp = j 
                                        ktemp = k 
                                        if (itemp == N) then 
                                                itemp = 0 
                                        else if (itemp == -1) then 
                                                itemp = N-1
                                        end if 
                                        if (jtemp == N) then 
                                                jtemp = 0 
                                        else if (jtemp == -1) then 
                                                jtemp = N-1
                                        end if 
                                        if (ktemp == N) then 
                                                ktemp = 0 
                                        else if (ktemp == -1) then 
                                                ktemp = N-1 
                                        end if 
                                        !print *, "cellIndex = ",cellIndex," itemp,jtemp,ktemp = ",itemp,jtemp,ktemp

                                        tempIndex = (itemp)*N**2 + (jtemp)*N + (ktemp) + 1
                                        neighborCellList(counter) = tempIndex
                                
                                        counter = counter + 1
                                        
                                end do 
                            end do 
                        end do 

        end subroutine getNeighboringCells
function getNeighborCells(cellIndex, N) result(neighborCells)
    integer, intent(in) :: cellIndex, N  ! N is cells per side
    integer :: neighborCells(27)
    integer :: i, j, k, counter, icoord, jcoord, kcoord
    integer :: itemp, jtemp, ktemp, tempIndex
    
    ! Get i,j,k coordinates from cell index
    icoord = (cellIndex - 1)/(N*N)
    jcoord = mod((cellIndex - 1)/N, N)
    kcoord = mod((cellIndex - 1), N)
    
    counter = 1
    do i = icoord - 1, icoord + 1 
        do j = jcoord - 1, jcoord + 1 
            do k = kcoord - 1, kcoord + 1 
                ! Handle periodic boundaries
                itemp = i
                jtemp = j
                ktemp = k
                
                if (itemp == N) then 
                    itemp = 0 
                else if (itemp == -1) then 
                    itemp = N-1
                end if 
                
                if (jtemp == N) then 
                    jtemp = 0 
                else if (jtemp == -1) then 
                    jtemp = N-1
                end if 
                
                if (ktemp == N) then 
                    ktemp = 0 
                else if (ktemp == -1) then 
                    ktemp = N-1 
                end if 
                
                tempIndex = itemp*N**2 + jtemp*N + ktemp + 1
                neighborCells(counter) = tempIndex
                counter = counter + 1
            end do 
        end do 
    end do 
        ! Sort neighbor array to prevent deadlocks
        do i = 1, size(neighborCells)-1
            do j = i+1, size(neighborCells)
                if (neighborCells(i) >neighborCells(j)) then
                    tempint =neighborCells(i)
                    neighborCells(i) = neighborCells(j)
                    neighborCells(j) = tempint
                end if
            end do
        end do
end function getNeighborCells
        function distance(chainMesh, atomIndex1, atomIndex2) result(d)
                type(ChainMesh_t), intent(in) :: chainMesh
                integer, intent(in) :: atomIndex1, atomIndex2 
                integer :: ind1,ind2 ! ChainCell Coordinates of atoms 1 and 2 
                integer :: i1,j1,k1,i2,j2,k2 ! Chain Cell i,j,k coordinates of 
                integer :: N 
                type(Atom_t) :: atom1, atom2 
                real :: d, domainWidth 
                domainwidth = chainMesh%domainWidth
                atom1 = chainMesh%atoms(atomIndex1)
                atom2 = chainMesh%atoms(atomIndex2)
                N = chainMesh%numChainMeshCellsPerSide
                dx = abs(atom1%x - atom2%x)
                dy = abs(atom1%y - atom2%y) 
                dz = abs(atom1%z - atom2%z)

                if (dx > domainWidth/2) dx = domainWidth - dx 
                if (dy > domainWidth/2) dy = domainWidth - dy 
                if (dz > domainWidth/2) dz = domainWidth - dz 
                
                d = sqrt(dx**2 + dy**2 + dz**2)
        end function distance 
        function distance_between_lattice_vectors(chainMesh, atomIndex1, atomIndex2) result(res)
                type(ChainMesh_t), intent(in) :: chainMesh
                integer, intent(in) :: atomIndex1, atomIndex2 
                type(vecNd_t) :: res 
                type(Atom_t) :: atom1, atom2 
                real :: domainWidth 
                implicit none 
                domainwidth = chainMesh%domainWidth
                atom1 = chainMesh%atoms(atomIndex1)
                atom2 = chainMesh%atoms(atomIndex2)
                N = chainMesh%numChainMeshCellsPerSide
                dx = abs(atom1%x - atom2%x)
                dy = abs(atom1%y - atom2%y) 
                dz = abs(atom1%z - atom2%z)

                if (dx > domainWidth/2) dx = domainWidth - dx 
                if (dy > domainWidth/2) dy = domainWidth - dy 
                if (dz > domainWidth/2) dz = domainWidth - dz 
                
                res = makeVecNd((/dx,dy,dz/))

        end function distance_between_lattice_vectors
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
                        
                        if (atomIndexTemp == AtomIndex) then 
                                atomIndexTemp = atoms(atomIndexTemp)%nextAtom
                        
                                cycle
                        end if
                        tempAtom = atoms(atomIndexTemp)
                         dist = distance(chainMesh,AtomIndex,AtomIndexTemp)
                        if ( dist .lt. nearestNeighborDistance) then 
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
                        if (abs(dist - nearestNeighborDistance) < chainMesh%chainMeshCells(1)%latticePArameter/10) then ! Tolerance 
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
                N = chainMesh%numChainMeshCellsPerSide
                !print *, "Assigned N"
                numChainMeshCells = chainMesh%numChainMeshCells
                !print *, "Assigned numChainMeshCells = ", numChainMeshCells
                ! For each chain mesh cell, for each atom in that cell, compute nearest neighbors and     
                ! check neighbor chain mesh cells for nearest neighbor atoms, taking account the
                ! periodicity of the system
                !$omp parallel do shared(chainMesh,numChainMeshCells,atomLockArray) private(idex,AtomIdent,NeighborCellList)
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
                                !print *, "DEBUG: Looping with AtomIdent = ", AtomIdent
                                call AssignAtomNearestNeighbhors(chainMesh,AtomIdent,idex,neighborCellList,atomLockArray) 
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

        subroutine enumerateChainMeshCells(chainMesh)
                type(chainMesh_t), intent(inout), target :: chainMesh 
                type(ChainMeshCell_t), pointer :: chainMeshCells(:)
                type(Atom_t), pointer :: atoms(:)
                integer :: atomIdent, i, atomCounter, idex ,N 
                N = chainMesh%numChainMeshCellsPerSide
                W = chainMesh%domainWidth 
                chainMeshCells => chainMesh%chainMeshCells
                atoms => chainMesh%atoms
                atomCounter = 0 
                do i = 1,size(chainMeshCells)
                        print *, "###################################"
                        atomIdent = chainMeshCells(i)%firstAtomInMeshCell
                        print *, "MeshCell ",i," contains atoms:"
                        do while (atomIdent .ne. -1)
                                idex = IndexFromCoordinates(atoms(atomIdent)%x,atoms(atomIdent)%y,atoms(atomIdent)%z,N,W)
                                print *, atomIdent,&
                                       " With coordinates: ", "(",atoms(atomIdent)%x,",",&
                                      atoms(atomIdent)%y,",",atoms(atomIdent)%z,")"
                                print *, "tmpx, tmpy, tmpz = ",atoms(atomIdent)%tmpx,atoms(atomIdent)%tmpy,atoms(atomIDent)%tmpz
                                print *, "Corresponding to Index: ", idex
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
        function IndexFromCoordinates(Ain,Bin,Cin,N,Win) result(res)
                real, intent(in) :: Ain,Bin,Cin, Win 
                real(kind=8) :: A,B,C,W, latticeParam  
                real(kind=8), parameter :: eps = 1e-8
                integer :: res 
                integer :: i,j,k
                latticeParam = W/dble(N)
                A = abs(dble(Ain)) 
                B = abs(dble(Bin)) 
                C = abs(dble(Cin))
                W = dble(Win)
                !print *, "A/W = ", A/W, "B/W = ", B/W, "C/W = ", C/W
                
                i = floor((A/W) * real(N))  
                j = floor((B/W) * real(N)) 
                k = floor((C/W) * real(N))

                !if ((latticeParam - (A - dble(i)*latticeParam)) < eps .and. (A - dble(i)*latticeParam > 0.0) .and. i < N-1) then
                 !       i = i + 1 
                !end if
                !if ((latticeParam - (B - dble(i)*latticeParam)) < eps .and. (B - dble(i)*latticeParam > 0.0) .and. j < N-1) then
                 !       j = j + 1 
                !end if                       
                !if ((latticeParam - (C - dble(i)*latticeParam)) < eps .and. (C - dble(i)*latticeParam > 0.0) .and. k < N-1) then
                 !       k = k + 1 
                !end if            !    i = nint((((dble(N)*A)/W))) 
                 !j = nint((((dble(N)*B)/W))) 
          !      k = nint((((dble(N)*C)/W)))
                !print *, "i,j,k = ", i,j,k 

                res = i*N**2 + j*N + k + 1 ! +1 to work with fortran 1 based array indexing
                !print *, "W = ", W 
                !print *, "A/W, B/W, C/W = ", A/W, B/W, C/W 
                !print *, "Assigning coordinates: ", A,B,C, " to ijk = ", i,j,k
        end function IndexFromCoordinates

        subroutine coordinatesFromIndex(Ind,N,i,j,k)
                implicit none
                integer, intent(in) :: Ind,N 
                integer, intent(out) :: i,j,k
                integer :: remainderj, IndTemp, reconstructedIndex
                IndTemp = Ind-1
                i = int(floor(real(IndTemp)/real(N**2)))
                remainderj = mod(IndTemp,N**2)
                j = int(real(remainderj)/real(N))
                k = mod(remainderj,N)
                reconstructedIndex = i*N**2 + j*N + k + 1 ! Needed for 1 based indexing in Fortran 
                !!print *, "Index ", Ind, " Corresponds to (",i,",",j,",",k,") & 
                       ! reconstructed Index: ",reconstructedIndex 
        end subroutine coordinatesFromIndex
        subroutine AssignAtomsToUnitCells(chainMesh,domainWidth,N)
                type(ChainMesh_t), intent(inout) :: chainMesh
                integer :: ChainMeshCellsLen
                integer :: atomsLen
                integer, intent(in) :: N 
                real, intent(in) :: domainWidth 
                integer :: i,j, chainCellIndex
                real :: tmpx, tmpy, tmpz, eps 
                eps = chainMesh%latticeParameter / 20.0
                ChainMeshCellsLen = size(chainMesh%chainMeshCells)
                atomsLen = size(chainMesh%atoms)
                
                do i=1,atomsLen 
                        tmpx = chainMesh%atoms(i)%x
                        tmpy = chainMesh%atoms(i)%y
                        tmpz = chainMesh%atoms(i)%z
                        if (modulo(chainMesh%atoms(i)%x,chainMesh%latticeParameter) < eps .and. i > 0 ) then 
                                tmpx = tmpx + eps  
                        end if 
                        if (modulo(chainMesh%atoms(i)%y,chainMesh%latticeParameter) < eps .and. j > 0) then 
                                tmpy = tmpy + eps 
                        end if 
                        if (modulo(chainMesh%atoms(i)%z,chainMesh%latticeParameter) < eps .and. k > 0) then 
                                tmpz = tmpz + eps 
                        end if 
                        chainMesh%atoms(i)%tmpx=tmpx
                        chainMesh%atoms(i)%tmpy=tmpy
                        chainMesh%atoms(i)%tmpz=tmpz
                        chainCellIndex = IndexFromCoordinates(tmpx,tmpy,&
                                tmpz,N,domainWidth)
                        call addAtomToChainCell(chainCellIndex,i,chainMesh)
                !TODO: Complete this subroutine
                end do 

        end subroutine AssignAtomsToUnitCells

        function makeChainMesh(numAtomsPerUnitCell, numUnitCellsPerSideLength, & 
                        latticeParameter,   AtomsInUnitCell ) result(chainMesh)
                implicit none 
                integer, intent(in) :: numAtomsPerUnitCell, numUnitCellsPerSideLength
                real, intent(in) :: latticeParameter
                type(Atom_t), intent(in) :: AtomsInUnitCell(numAtomsPerUnitCell)
                ! Need to create all the atoms, create all the ChainMeshCells, then allocate all of the atoms to a chain mesh cell 
                type(ChainMesh_t) :: chainMesh 
                integer :: numChainMeshCells 
                integer :: numAtoms, stride 
                integer :: i,j,k, icoord, jcoord, kcoord
                type(ChainMeshCell_t) :: tempChainMeshCell
                type(Atom_t) :: tempAtoms(size(AtomsInUnitCell))
                real :: domainWidth
                chainMesh%latticeParameter = latticeParameter
                chainMesh%numChainMeshCellsPerSide = numUnitCellsPerSideLength
                numChainMeshCells = numUnitCellsPerSideLength**3
                domainWidth = latticeParameter*numUnitCellsPerSideLength
                chainMesh%domainWidth = domainWidth
                numAtoms = numAtomsPerUnitCell*numChainMeshCells
                numChainMeshCells = numUnitCellsPerSideLength**3
                chainMesh%numChainMeshCells = numChainMeshCells
                allocate(chainMesh%chainMeshCells(numChainMeshCells))
                allocate(chainMesh%atoms(numAtoms))
                stride = size(AtomsInUnitCell)
                ! initialise atoms and chain Mesh Cells 
                do j = 1,numChainMeshCells ! do it for each unit cell 
                        tempAtoms = AtomsInUnitCell
                        tempChainMeshCell%latticeParameter = latticeParameter
                        tempChainMeshCell%NumAtomsPerUnitCell = NumAtomsPerUnitCell
                        tempChainMeshCell%firstAtomInMeshCell = -1 ! Initialise to garbage value so that we don't have to deal with
                        ! it
                        chainMesh%chainMeshCells(j) = tempChainMeshCell ! Just so that it isn't initialised to garbage  
             

                        call coordinatesFromIndex(j,numUnitCellsPerSideLength,icoord,jcoord,kcoord) ! icoord, jcoord, kcoord integer
                        do i=1,size(AtomsInUnitCell)
                                !chainMesh%atoms(j,j+atoms) = AtomsInUnitCell ! Assuming they are all initialised properly
                                ! coordinates for the chain mesh cell 
                                tempAtoms(i)%x = AtomsInUnitCell(i)%x + real(icoord)*latticeParameter 
                                tempAtoms(i)%y = AtomsInUnitCell(i)%y + real(jcoord)*latticeParameter
                                tempAtoms(i)%z = AtomsInUnitCell(i)%z + real(kcoord)*latticeParameter 
                                ! Assign Atoms to Unit Cells 
                                chainMesh%atoms((j-1)*stride + i) = tempAtoms(i)
                                call addAtomToChainCell(j,(j-1)*stride + i,chainMesh)
                        end do
                        !!print *, "j debug = ", j, " stride = ",stride
                        !chainMesh%atoms((j-1)*stride + 1:(j-1)*stride+size(AtomsInUnitCell)) = tempAtoms
                        !Suspicious! ^^
               end do
                ! Todo: map atoms into their respective chain mesh cells
                !call AssignAtomsToUnitCells(chainMesh,domainWidth,numUnitCellsPerSideLength) 
        end function makeChainMesh 

function H(chainMesh,sigma) result(Energy)
    type(chainMesh_t), intent(in), target :: chainMesh 
    real, intent(in) :: sigma 
    real(kind=8) :: Energy
    type(Atom_t), pointer :: atoms(:)
    integer :: i, j, neighbor_idx
    
    atoms => chainMesh%atoms 
    Energy = 0.0d0
    
    do i = 1, size(atoms)
        do j = 1, size(atoms(i)%NeighborList)
            ! Get the neighbor's index first
            neighbor_idx = atoms(i)%NeighborList(j)
            ! Then use it to access the neighbor's parameters
            Energy = Energy + dble(sigma) * dble(atoms(i)%AtomParameters(1)) * &
                    dble(atoms(neighbor_idx)%AtomParameters(1))
        end do 
    end do 
end function H



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


        function compute_unique_lattice_vectors(chainMesh, outputArray) result(distanceArray)
                type(ChainMesh_t), intent(in), target :: chainMesh 
                type(Atom_t), pointer :: tempAtom
                integer :: i 
                type(vecNd_t) :: G
                type(vecNd_t), allocatable :: distanceArray(:)
                
                allocate(distanceArray(1))
                distanceArray(1) = makeVecNd((/0.0,0.0,0.0/))
                tempAtom => chainMesh%atoms(1) 
                
                do i = 1, size(tempAtom%NeighborList)
                                ! distance_between_lattice_vectors correctly handles the periodicity of the lattice.
                                G = distance_between_lattice_vectors(chainMesh, 1,chainMesh%atoms(tempAtom%NeighborList(i)))
                                if (any(distanceArray == G)) then 
                                        cycle 
                                end if 
                                distanceArray = [distanceArray, G] 
                end do 
        end function compute_unique_lattice_vectors

end module ChainMesh 

