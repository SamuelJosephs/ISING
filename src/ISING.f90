
module Atom 
        use omp_lib
        type Atom_t 
                real :: x,y,z, tmpx, tmpy, tmpz
                integer :: NumAtomParameters
                real , allocatable :: AtomParameters(:) ! e.g. spin
                integer :: nextAtom
                integer, allocatable :: NeighborList(:)

        end type
contains
        function makeAtom(x,y,z,atomParameters,NumAtomParameters,  NextAtom) result(res)
                real, intent(in) :: x,y,z
                integer, intent(in) :: NumAtomParameters
                real, intent(inout), target :: AtomParameters(NumAtomParameters)
                type(Atom_t) :: res
               res%x = x
               res%y = y
               res%z = z
               res%NumAtomParameters = NumAtomParameters
               res%NextAtom = -1
               allocate(res%atomParameters(size(atomParameters)))
               res%atomParameters = atomParameters

               allocate(res%NeighborList(0)) 
               return
        end function makeAtom

        function AtomDeepCopy(atom) result(res)
        type(Atom_t), intent(inout), pointer :: atom 
        type(Atom_t) :: res 
        res%nextAtom = atom%nextAtom 
        res%x = atom%x 
        res%y = atom%y
        res%z = atom%z
        res%numAtomParameters = atom%numAtomParameters 
        allocate(res%atomParameters(size(atom%atomPArameters)))
        res%atomParameters = atom%atomParameters 
        allocate(res%NeighborList((size(atom%NeighborList))))
        res%neighborList = atom%neighborList
        end function atomDeepCopy 
end module Atom 

module ChainMeshCell
        type ChainMeshCell_t 
                real :: latticeParameter
                integer :: NumAtomsPerUnitCell
                integer :: firstAtomInMeshCell 
        end type ChainMeshCell_t 
end module ChainMeshCell 

module ChainMesh
        use Atom 
        use ChainMeshCell 
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
end module ChainMesh 



module EnergyMin
        use Atom 
        use ChainMesh 
        use ChainMeshCell
        use Rand, only: makeRandom, algor_uniform_random, random
        use cube_partition, only: bounds_type, partition_unit_cube, dp 
        use omp_lib 
        implicit none 
        contains 

function AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex) result(Energy)
    type(chainMesh_t), intent(inout), target :: chainMesh 
    real(kind=8), intent(in) :: sigma
    type(Atom_t), pointer, intent(inout) :: atom
    integer, intent(in) :: threadID, atomIndex
    real(kind=8) :: Energy
    type(Atom_t), pointer :: atoms(:)
    integer :: j, neighbor_idx
    if (.not. associated(atom)) then 
        error stop "atom is not a valid pointer, critical error: aborting"
    end if

    atoms => chainMesh%atoms 
    Energy = 0.0d0
    
    do j = 1, size(atom%NeighborList)
        ! Get the neighbor's index first
        neighbor_idx = atom%NeighborList(j)
        ! Then use it to access the neighbor's parameters

        Energy = Energy - dble(sigma) * dble(atom%AtomParameters(1)) * &
                dble(atoms(neighbor_idx)%AtomParameters(1))
    end do 
!     print *, "Thread ", threadID, "Calculating energy of atomIndex ", atomIndex, " Energy = ", Energy

end function AtomEnergy 

                subroutine coordsFromThreadIndex(Ind,b,c, i,j,k) ! factored Ind into a*b*c. 
                        integer, intent(in) :: Ind,b,c 
                        integer, intent(out) :: i,j,k
                        integer :: remainder 
                        real(kind=dp), dimension(3) :: min_point, max_point
                       ! Ind = i (bc) + j(c) + k. This is just the formula for stride. This can be solved for i,j,k using the
                        ! following: 
                        i = Ind/(b*c)
                        remainder = modulo(Ind,b*c)
                        j = remainder / c 
                        k = modulo(remainder,c)
                       ! i,j,k range from (0,a-1), (0,b-1), and (0,c-1) respectively
                end subroutine coordsFromThreadIndex 

               ! Do the metropolis algorithm for the chain mesh 
                subroutine Metropolis(chainMesh, sigma, beta, nsteps)
                        implicit none
                        type(chainMesh_t), intent(inout), target :: chainMesh
                        real(kind=8), intent(in) :: sigma 
                        real(kind=8), intent(in) :: beta !1 / k_B T 
                        integer, intent(in) :: nsteps
                        integer :: threadID, numThreads, binsPerSide, a,b,c 
                        type(bounds_type), allocatable :: thread_bounds(:)
                        integer :: i, imin, imax, jmin, jmax, kmin, kmax 
                        type(bounds_type) :: thread_bound
                        logical :: critical_error, success
                        type(random) :: random
                        real(kind=dp) :: randomNumber, randi, randj, randk
                        real(kind=dp), dimension(3) :: min_point, max_point
                        integer :: N, randAtom, itemp, atomIndex, cell_index, num_atoms_in_cell, first_atom_index
                        integer :: tempint
                        type(Atom_t), pointer :: atom, tempAtomPointer 
                        real(kind=8) :: Energy, NewEnergy, Z, prob  
                        integer(kind=OMP_LOCK_KIND), allocatable :: lockArray(:), cell_lock_array(:)
                        real, allocatable :: atomParameterArray(:)
                        integer :: clock_count, clock_rate, attempt, j 
                        type(Atom_t), target :: tempAtom 
                        integer, dimension(27) :: cell_neighbor_array, sorted_cells 
                        call system_clock(clock_count,clock_rate)
                        allocate(lockArray(chainMesh%numAtoms))
                        allocate(cell_lock_array(chainMesh%numChainMeshCells))
                        allocate(atomParameterArray(size(chainMesh%atoms(1)%AtomParameters)))
                        do i=1,size(lockArray)
                            call omp_init_lock(lockArray(i))
                        end do
                        do i=1,size(cell_lock_array)
                           
                            call omp_init_lock(cell_lock_array(i))
                        end do 
                        ! Algorithm:
                        ! Each thread will get it's thread ID and figure out how many threads there are 
                        ! We will factor the number of threads into a,b,c such that a*b*c = num threads 
                        ! From this we will assign a range of chain mesh cells to each thread 
                        ! Each thread will do Metropolis on each of it's chain mesh cells 
                        critical_error = .false.
                        print *, "Entering parralell section"
                        !$omp parallel &
                        !$omp shared(chainMesh, critical_error, lockArray, beta, nsteps, sigma,cell_lock_array) &
                        !$omp private(threadID, numThreads, a, b, c, binsPerSide, thread_bounds, & 
                        !$omp         random, N, randAtom, itemp, i, imin, imax, jmin, jmax, &
                        !$omp         kmin, kmax, atomIndex, cell_index, Energy, NewEnergy, Z, prob, &
                        !$omp         randi, randj, randk, randomNumber, min_point, max_point, tempint, atom, tempAtom, &
                        !$omp         tempAtomPointer, atomParameterArray,cell_neighbor_array, sorted_cells, attempt,j,success) & 
                        !$omp default(private)
                       block 
                                !binsPerSide = chainMesh%numChainMeshCellsPerSide
                                N = chainMesh%numChainMeshCellsPerSide
                                threadID = omp_get_thread_num()
                                numThreads = omp_get_num_threads()
                                ! Compute which range of x,y,z for atoms within which the thread will randomly sample
                                thread_bounds = partition_unit_cube(numThreads)
                                ! if (threadID == 0) then !testing that it works 
                                !         print *, "num threads: ", numThreads 
                                !         do i=1,size(thread_bounds)
                                !                 print *, "i = ", i 
                                !                 print *, "lower bound: ", thread_bounds(i)%min_point
                                !                 print *, "upper bound: ", thread_bounds(i)%max_point 
                                                
                                !         end do
                                ! end if
                                thread_bound = thread_bounds(threadID + 1) ! +1 due to Fortran indexing
                                min_point = thread_bound%min_point
                                max_point = thread_bound%max_point 
                                imin = int(min_point(1)*N)
                                imax = int(max_point(1)*N)
                                !if (imax == N) print *, "imax == N"
                                !if (imin == 0) print *, "imin == 0"
                                jmin = int(min_point(2)*N)
                                jmax = int(max_point(2)*N)
                               !if (jmax == N) print *, "jmax == N"
                               !if (jmin == 0) print *, "jmin == 0"
                                kmin = int(min_point(3)*N)
                                kmax = int(max_point(3)*N)
                               !if (kmax == N) print *, "kmax == N" 
                               !if (kmin == 0) print *, "kmin == 0"
                                ! Check for sensible values
                                if (imin == imax) then 
                                        print *, "imin == imax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                if (jmin == jmax) then 
                                        print *, "jmin == jmax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                if (kmin == kmax) then 
                                        print *, "kmin == kmax, please use less threads or more mesh cells"
                                        !$omp atomic write 
                                        critical_error = .true.
                                        !$omp end atomic 
                                end if 
                                ! Read only so no race conditions within this barrier 
                                if (critical_error) then 
                                        error stop "Critical Error detected, ending execution on all threads"
                                end if 

                                ! Avoid overlaps between threads 
                                if (imin /= 0 ) then 
                                        imin = imin + 1
                                end if 
                                if (jmin /= 0 ) then 
                                        jmin = jmin + 1
                                end if 
                                if (kmin /= 0) then 
                                        kmin = kmin + 1 
                                end if 
                                ! Now need to select a mesh cell i,j,k coordinate within the threads assigned bounds, then select an
                                ! atom within it 
                                ! To do this we will need the lobotomised CASTEP random number generator as the standard library one
                                ! is probably not thread safe (better safe than sorry).
                                ! This is a bitshift random number generator which returns numbers in the range [0,1], we will also need to warm up
                                ! the generator otherwise threads will make very similar choices for the first few steps.
                                ! We will use the threadID as the seed 
                                !print *, "Thread ", threadID, " is about to make a random number"
                                random=makeRandom(threadID + clock_count)
                                !print *, "Thread ", threadID, " has made a  random number"

                                !Warm up the generator 
                                do i=1,30
                                    randomNumber = algor_uniform_random(random)
                                end do
                                !print *, "Thread ", threadID, " has warmed up the random number generator"

                                do itemp=1,nsteps
                                        !$omp barrier
                                        !Randomly select i,j, and k 
                                        !randi = imin + floor(algor_uniform_random(random)*real((imax - imin)))
                                        !randj = jmin + floor(algor_uniform_random(random)*real((jmax - jmin)))
                                        !randk = kmin + floor(algor_uniform_random(random)*real((kmax - kmin)))
                                        !if (int(randi) == N - 1) print *, "randi == N - 1"
                                        !if (int(randi) == 0) print *, "randi == 0"
                                        !Now randomly select atom in this chain mesh cell
                                        !claude: 
                                       ! cell_index = nint(randi*N**2) + nint(randj*N) + nint(randk) + 1
                                        !cell_index = int(randi)*N**2 + int(randj)*N + int(randk) + 1 TODO: This is a temporary
                                        !exploration using an alternative sampling algorithm 
                                        randi = floor(algor_uniform_random(random) * N) 
                                        randj = floor(algor_uniform_random(random) * N)
                                        randk = floor(algor_uniform_random(random) * N)
                                        cell_index = int(randi)*N*N + int(randj)*N + int(randk) + 1
                                        !print *, "Thread ", threadID, " Is trying to acquire the lock of cell number ", cell_index

                                        call omp_set_lock(cell_lock_array(cell_index))
                                        !print *, "Thread ", threadID, "Has acquired the lock of cell: ", cell_index 
                                        !print *, "cell_index = ", cell_index, " out of ", N**3 
                                        !print *, "Thread ", threadID, "is setting it's locks: ", cell_neighbor_array
        
                                        num_atoms_in_cell = chainMesh%chainMeshCells(cell_index)%NumAtomsPerUnitCell
                                        randAtom = floor(algor_uniform_random(random) * real(num_atoms_in_cell))
                                        !randAtom = floor(algor_uniform_random(random)*real(chainMesh%chainMeshCells(nint(randi*N**2) + &
                                                !nint(randj*N) +  nint(randk) + 1 )%NumAtomsPerUnitCell))                                               ! randomly 
                                        !claude: 

                                        !cell_index = int(randi*N**2) + int(randj*N) + int(randk) + 1
!                                        cell_index = int(randi)*N**2 + int(randj)*N + int(randk) + 1 
                                        first_atom_index = chainMesh%chainMeshCells(cell_index)%firstAtomInMeshCell
                                        atom => chainMesh%atoms(first_atom_index)
                                        !print *, "Thread ", threadID, "has initialised atom"
                                        !atom = chainMesh%atoms(chainMesh%chainMeshCells(int(randi*N**2) + int(randj*N) + int(randk) + 1 &
                                        !)%firstAtomInMeshCell)
                                        !atomIndex = chainMesh%chainMeshCells(int(randi*N**2) + int(randj*N) + int(randk) + 1 &
                                        !)%firstAtomInMeshCell  
                                        atomIndex = first_atom_index
                                        if (atom%nextAtom == -1) then !DEBUG  
                                                print *, "Atom has next atom -1, critical error: aborting"
                                                error stop "Critical Error"
                                        end if
                                        !TODO: Do metropolis step then put in loop
                                        ! Check that all of the first elements of the chainMeshCell's first atoms are initialised:
                                        do i=1,size(chainMesh%chainMeshCells) !DEBUG 
                                            if (chainMesh%chainMeshCells(1)%firstAtomInMeshCell == -1) then 
                                                    print *, "chainMeshCell ", i, " is not properly initialised"
                                                    error stop "Critical error, aborting"
                                            end if 
                                        end do 
                                        do i = 1,randAtom ! Don't need + 1 as if randAtom=0 the od loop won't run, if randAtom=1 it will run
                                                          !once
                                            !print *, "i = ", i, "atom%nextatom = ", atom%nextAtom, " randAtom = ", randAtom
                                            atomIndex = atom%NextAtom !This should never be -1, if it is we would get an out of
                                                                      ! Out of bounds memory accsess.
                                            
                                            atom => chainMesh%atoms(atom%nextAtom)
                                            
                                        end do
                                        ! print *, "Thread ", threadID , " has selected atom with index", atomIndex
                                        !print*, "Thread ", threadID, "Has randomly selected atom ", randAtom, " In it's unit cell: ", &
                                        !int(randi), int(randj), int(randk)
                                        !TODO: Finish metropolis on these randomly selected atoms, maybe make atom a pointer such that we
                                        !can change it's value within the array
                                        if (.not. associated(atom)) then 
                                                print *, "Thread ", threadID, " has a not valid atom pointer with index ", atomIndex
                                                error stop "Threadm atom not a valid pointer"
                                        end if 
                                        !print *, "Thread ", threadID, " is passing atomIndex ", atomIndex
                                        Energy = AtomEnergy(chainMesh, sigma, atom, threadID, atomIndex)
                                        ! print *, "Thread ", threadID, " is about to index and write to atomParameters of atom "
                                        ! print *, atomIndex

                                        ! atom%atomParameters(1) = - atom%atomParameters(1) ! TODO:Atom is a pointer so should make a copy of atomparameters(1)
                                        atomParameterArray = - atom%atomParameters !This is not being used but should be
                                        
                                        ! print *, "Thread ", threadID, " has indexed atomParameters"
                                        ! In order to avoid race conditions while limiting the use of locks, a copy of the atom
                                        ! being considered must be made. Annoyingly this 
                                        !tempAtom = makeAtom(atom%x,atom%y,atom%z,atomParameterArray,atom%numAtomParameters,atom%nextAtom)
                                        tempAtom = atomDeepCopy(atom)  
                                        tempAtomPointer => tempAtom
                                        tempAtom%atomParameters = -atom%atomParameters
                                        NewEnergy = AtomEnergy(chainMesh, sigma, tempAtomPointer,threadID,atomIndex)
                                        ! tempAtom should be automatically deallocated at the end of the scope
                                        !print *, "Thread ", threadID, " has obtained a probability"
                                        
                                        if (NewEnergy < Energy) then 
                                                Z = 1.0 
                                        else 
                                                Z = exp(- beta * ((NewEnergy - Energy)/sigma))
                                        end if 
                                        prob = algor_uniform_random(random)
                                        ! print *, "Thread ", threadID, " has obtained a probability"
                                        ! tempint = atomIndex
                                        ! print *, "Thread ", threadID, " has selected atomIndex to be ",tempint
                                        ! Barrier here to prevent race conditions from threads calculating energies (Reading) and
                                        ! updating the mesh (writing)
                                        !print *, "Debug: thread ", threadID, " Z = ", Z, " prob = ", prob, " atomIndex = ",atomIndex, "New Energy: ", NewEnergy, " Old Energy ", Energy
                                        if (Z >= prob) then 
                                                !call omp_set_lock(lockArray(atomIndex))
                                                do i = 1,size(atom%atomParameters)
                                                !     print *, "Thread ", threadID, " is about to write to atom"
                                                    chainMesh%atoms(atomIndex)%atomParameters(i) = - atom%atomParameters(i) 
                                                    !print *, "Thread ", threadID, " Writing to mesh with value: ", chainMesh%atoms(atomIndex)%atomParameters(i)
                                                end do 
                                                !call omp_unset_lock(lockArray(atomIndex))
                                                !print *, "Thread ", threadID, "Has unset the lock"
                                                ! print *, "Thread ", threadID, " has updated the mesh for atom index", atomIndex

                                        end if 
                                        call omp_unset_lock(cell_lock_array(cell_index))
                                        !print *, "Thread ", threadID, " has released the lock of cell number ", cell_index
                                        !$omp barrier
 
                                        !print *, "######################################################"
                                end do 

                      end block 

                     !$omp end parallel 
                     do i = 1,size(lockArray)
                        call omp_destroy_lock(lockArray(i))
                     end do 

                     do i=1,size(cell_lock_array)
                        call omp_destroy_lock(cell_lock_array(i))
                     end do 
            end subroutine Metropolis 
subroutine WriteMagnetization(chainMesh, filename)
    implicit none
    type(ChainMesh_t), intent(in) :: chainMesh
    character(len=*), intent(in) :: filename
    real(kind = 8) :: total_mag, mag_per_atom
    integer :: i, unit_num
    
    ! Calculate average magnetization
    total_mag = 0.0
    do i = 1, size(chainMesh%atoms)
        total_mag = total_mag + chainMesh%atoms(i)%AtomParameters(1)
    end do
    mag_per_atom = total_mag / real(size(chainMesh%atoms))
    
    ! Open file in append mode
    open(newunit=unit_num, file=trim(filename), position='append', action='write')
    
    ! Write magnetization per atom to new line
    write(unit_num, *) mag_per_atom
    
    close(unit_num)
    
end subroutine WriteMagnetization
end module EnergyMin 
program main 
        use Atom
        use chainMeshCell 
        use ChainMesh
        use EnergyMin
        implicit none
        integer :: numCells, numsteps, i 
        character(len=30) :: arg
        type(ChainMesh_t) :: testMesh
        type(Atom_t) :: AtomsInUnitCell(2)
        real :: AtomParam1(1), AtomParam2(1), a, energy 
        real(kind=8) :: idble, T, Kb,beta, sigma, upperT, lowerT  
        !Kb = 8.62e-5 ! ev/ J 
        Kb = 1.38e-23_dp 
        !Kb = 8.617e-5_dp
        a = 2.8
        numsteps = 60
        upperT = 1500 
        lowerT = 500  
        if (command_argument_count() /= 1) then 
                print *, "Please specify the number of unit cells per side length"
                stop
        end if
        call get_command_argument(1,arg)
        read(arg,*) numCells 
        AtomParam1 = (/1.0/)
        AtomParam2 = (/1.0/)
        AtomsInUnitCell(1) = makeAtom(0.0,0.0,0.0,AtomParam1,1,-1) 
        AtomsInUnitCell(2) = makeAtom(a/2,a/2,a/2,AtomParam2,1,-1)
        testMesh = makeChainMesh(2,numCells,a,AtomsInUnitCell)
        
        call assignNearestNeighbors(testMesh)
        !print *, "Succsess!"
        call enumerateChainMeshCells(testMesh) 
        !!print *, "Num atoms = ",size(testMesh%atoms)

        energy = H(testMesh,0.04565376186997809)
        print *, "Energy: ", energy 
        print *, "Num Atoms = ", size(testMesh%atoms)
        print *, "Energy / atom= ", energy / float(size(testMesh%atoms))
        !sigma = 0.04565376186997809_dp*1.6e-19_dp  
        !sigma = 0.04565376186997809_dp*1.6e-19_dp  
        sigma = 0.01565376186997809_dp*1.6e-19_dp   
        !call WriteMagnetization(testMesh,"Magnetisation3.dat") 
        do i=1,numsteps  
                idble = dble(i)
                idble = (idble / numsteps)
                T = idble * (upperT - lowerT) + lowerT ! Scan T in range(500,900)
                !beta = 1.0_dp / (Kb*T)
                beta = sigma / (Kb * T) !reduced 
                print *, "Reduced beta = ", beta
                call Metropolis(testMesh,sigma,beta,1000000)
                print *, "Writing Magnetisation to file:", i , "/" , numsteps
                call WriteMagnetization(testMesh,"Magnetisation4.dat")
        end do 
        !call deallocateChainMesh(testMesh)
end program main
