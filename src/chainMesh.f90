

module ChainMesh
        use Atom 
        use ChainMeshCell 
        use vecNd 
        use omp_lib
type ChainMesh_t 
        integer :: numAtoms, numChainMeshCells
        type(Atom_t), allocatable :: atoms(:) ! Each atom is in a chain mesh cell and points to the next atom within that cell 
        type(ChainMeshCell_t), allocatable :: chainMeshCells(:)
        real ::  latticeParameter
        integer :: numCellsX, numCellsY, numCellsZ   
        integer, allocatable, dimension(:,:,:) :: derivativeList !(i,j,k) : i=atomInex, j=dim (1,2,3), k = lower,higher (1,2)
end type ChainMesh_t 

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
                integer :: a,b,c  
                type(Atom_t) :: atom1, atom2 
                real :: d, widthX, widthY, widthZ  
                
                a = chainMesh%numCellsX 
                b = chainMesh%numCellsY 
                c = chainMesh%numCellsZ 
                widthX = a*chainMesh%latticeParameter
                widthY = b*chainMesh%latticeParameter 
                widthZ = c*chainMesh%latticeParameter

                atom1 = chainMesh%atoms(atomIndex1)
                atom2 = chainMesh%atoms(atomIndex2)
                dx = abs(atom1%x - atom2%x)
                dy = abs(atom1%y - atom2%y) 
                dz = abs(atom1%z - atom2%z)

                if (dx > widthX/2) dx = widthX - dx 
                if (dy > widthY/2) dy = widthY - dy 
                if (dz > widthZ/2) dz = widthZ - dz 
                
                d = sqrt(dx**2 + dy**2 + dz**2)
        end function distance 

        function distance_points(chainMesh, point1, point2) result(d)
            type(ChainMesh_t), intent(in) :: chainMesh
            type(vecNd_t), intent(in) :: point1, point2
            real(kind=8) :: d, dx, dy, dz, widthX, widthY, widthZ
            integer :: a, b, c

            ! Get system dimensions
            a = chainMesh%numCellsX
            b = chainMesh%numCellsY
            c = chainMesh%numCellsZ
            
            ! Calculate total width in each dimension
            widthX = a * chainMesh%latticeParameter
            widthY = b * chainMesh%latticeParameter
            widthZ = c * chainMesh%latticeParameter
            
            ! Calculate coordinate differences
            dx = abs(point2%coords(1) - point1%coords(1))
            dy = abs(point2%coords(2) - point1%coords(2))
            dz = abs(point2%coords(3) - point1%coords(3))
            
            ! Apply periodic boundary conditions to find the minimum distance
            if (dx > widthX/2) dx = widthX - dx
            if (dy > widthY/2) dy = widthY - dy
            if (dz > widthZ/2) dz = widthZ - dz
            
            ! Calculate Euclidean distance
            d = sqrt(dx**2 + dy**2 + dz**2)
        end function distance_points

        subroutine distance_points_vec(chainMesh, point1, point2,d)
            type(ChainMesh_t), intent(in) :: chainMesh
            type(vecNd_t), intent(in) :: point1, point2
            type(vecNd_t), intent(inout) :: d 

            real(kind=8) :: dx, dy, dz, widthX, widthY, widthZ
            integer :: a, b, c
            ! Get system dimensions
            a = chainMesh%numCellsX
            b = chainMesh%numCellsY
            c = chainMesh%numCellsZ
            
            ! Calculate total width in each dimension
            widthX = a * chainMesh%latticeParameter
            widthY = b * chainMesh%latticeParameter
            widthZ = c * chainMesh%latticeParameter
            
            ! Calculate coordinate differences
            dx = (point2%coords(1) - point1%coords(1)) ! Warning! May need to do point 1 - point 2 
            dy = (point2%coords(2) - point1%coords(2))
            dz = (point2%coords(3) - point1%coords(3))
            
            ! Apply periodic boundary conditions to find the minimum distance
            !if (dx > widthX/2) dx = widthX - dx
            !if (dy > widthY/2) dy = widthY - dy
            !if (dz > widthZ/2) dz = widthZ - dz


            dx = dx - widthX*nint(dx/widthX)
            dy = dy - widthY*nint(dy/widthY)
            dz = dz - widthZ*nint(dz/widthZ)
            
            
            ! Calculate Euclidean distance
            d = makeVecNdCheck(d,[dx,dy,dz])
        end subroutine distance_points_vec

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
                                !print *, "DEBUG: Looping with AtomIdent = ", AtomIdent
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
                                idex = IndexFromCoordinates(chainMesh,atoms(atomIdent)%x,atoms(atomIdent)%y,atoms(atomIdent)%z)
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
        function IndexFromCoordinates(chainMesh,Ain,Bin,Cin) result(res)
                type(ChainMesh_t), intent(in) :: chainMesh 
                real, intent(in) :: Ain,Bin,Cin 
                real(kind=8) :: A,B,C,W, latticeParam  
                real(kind=8), parameter :: eps = 1e-8
                integer :: res 
                integer :: i,j,k, Amesh, Bmesh, Cmesh 
                Amesh = chainMesh%numCellsX
                Bmesh = chainMesh%numCellsY 
                Cmesh = chainMesh%numCellsZ 
                
                latticeParam = chainMesh%latticeParameter
                i = int(Ain / latticeParam)
                j = int(Bin / latticeParam)
                k = int(Cin / latticeParam)
                res = i*Bmesh*Cmesh + j*Cmesh + k + 1 ! +1 to work with fortran 1 based array indexing
        end function IndexFromCoordinates

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
                        chainCellIndex = IndexFromCoordinates(chainMesh,tmpx,tmpy,&
                                tmpz)
                        call addAtomToChainCell(chainCellIndex,i,chainMesh)
                !TODO: Complete this subroutine
                end do 

        end subroutine AssignAtomsToUnitCells

        function makeChainMesh(numAtomsPerUnitCell, numCellsX, numCellsY, numCellsZ, & 
                        latticeParameter,   AtomsInUnitCell ) result(chainMesh)
                implicit none 
                integer, intent(in) :: numAtomsPerUnitCell, numCellsX, numCellsY, numCellsZ 
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

                chainMesh%numCellsX = numCellsX 
                chainMesh%numCellsY = numCellsY 
                chainMesh%numCellsZ = numCellsZ 

                numChainMeshCells = numCellsX*numCellsY*numCellsZ 
                numAtoms = numAtomsPerUnitCell*numChainMeshCells
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
             

                        call coordinatesFromIndex(chainMesh,j,icoord,jcoord,kcoord) ! icoord, jcoord, kcoord integer
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
                        
                        do while (atomIndexTemp /= -1)
                                if (atomIndexTemp == atomIndex) then 
                                        atomIndexTemp = chainMesh%atoms(atomIndexTemp)%nextAtom
                                        cycle 
                                end if
                                x = chainMesh%atoms(atomIndexTemp)%x 
                                y = chainMesh%atoms(atomIndexTemp)%y 
                                z = chainMesh%atoms(atomIndexTemp)%z 
                                atomPos2 = makeVecNdCheck(atomPos2, [x,y,z])
                                call distance_points_vec(chainMesh,atomPos1, atomPos2, distance)
                                
                                candidate = .True.
                                do j = 1,3
                                        if (j == d) cycle 
                                        if (abs(distance%coords(j)) > chainMesh%latticeParameter / 20.0_8) then 
                                                candidate = .False.
                                                exit 
                                        end if 
                                end do 
                                if (candidate) then 
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

        subroutine calculatePartialDerivative(chainMesh,atomIndex,d,outputVec)
                type(ChainMesh_t), intent(inout) :: chainMesh 
                integer, intent(in) :: atomIndex, d , outputVec

                if (.not. allocated(chainMesh%derivativeList)) then 
                        print *, "Warning: derivativeList not allocated, allocating and computing"
                        call DerivativeList(chainMesh,chainMesh%derivativeList)
                end if 
                
        end subroutine calculatePartialDerivative
end module ChainMesh 

