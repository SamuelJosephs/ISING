! This file will contain an MPI-openMP program for parallel tempering of the Heisenberg model as outlined in
! https://doi.org/10.1016/j.cor.2025.107000


program PT 
        use mpi_f08
        use PT_Utils
        use chainMesh
        use Atom
        use iso_fortran_env, only: error_unit, real64
        implicit none 
        integer :: MPI_ierr, MPI_rank, MPI_num_procs 
        integer, parameter :: dp = real64
        integer, parameter :: NJ = 2 ! The number of J, D, and B values to perform the parallel tempering at.  
        integer, parameter :: ND = 2 ! Hard code these for now but eventually they should be taken as input
        integer, parameter :: NB = 1 
        real(kind=dp), parameter :: JMin = 0.0_dp
        real(kind=dp), parameter :: JMax = 1.5_dp 
        real(kind=dp), parameter :: DMin = 0.0_dp 
        real(kind=dp), parameter :: DMax = 1.5_dp
        real(kind=dp), parameter :: BMin = 1.0_dp 
        real(kind=dp), parameter :: Bmax = 1.0_dp

        ! Set up constants for the lattice, for now they will be hardcoded but eventually they should be taken as input.
        type(Atom_t), dimension(2) :: atomsInUnitCell
        real, dimension(3), parameter :: atomParams = (/1.0, 0.0, 0.0/)
        integer, parameter :: numCellsX = 40 
        integer, parameter :: numCellsY = 40
        integer, parameter :: numCellsZ = 6
        real(kind=dp), parameter :: a_bravais = 2.8
        real(kind=dp), parameter :: b_bravais = 2.8
        real(kind=dp), parameter :: c_bravais = 2.8 
        real(kind=dp), parameter :: ab = 90
        real(kind=dp), parameter :: bc = 90 
        real(kind=dp), parameter :: ca = 90


        integer :: NumSlots, BasePtr, TopPtr, NumParams 
        integer :: stat, i, j, JIndex, DIndex, BIndex
        integer, allocatable, dimension(:) :: ParamIndexArray ! Contains the global index of each parameter set
        real(kind=dp), allocatable, dimension(:,:) :: ParamArray ! (paramIndex, (J,D,B))
        ! For each Parameter set we need a buffer for the spins 
        type(ChainMesh_t), allocatable, dimension(:,:) :: meshBuffer ! (paramIndex, Temp index) 
        
        integer, allocatable, dimension(:,:) :: TemperatureMeshArray ! (paramIndex, Temp index)
        real(kind=dp), allocatable, dimension(:) :: TemperatureArray ! (Temp)
        real(kind=dp), parameter :: TMax = 2.0_dp 
        real(kind=dp), parameter :: TMin = 0.0000001_dp 
        integer, parameter :: numTemps = 10

        call MPI_Init(MPI_ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD,MPI_rank)
        call MPI_Comm_Size(MPI_COMM_WORLD,MPI_num_procs)

        
        AtomsInUnitCell(1) = makeAtom(0.0, 0.0, 0.0, atomParams, -1) 
        AtomsInUnitCell(2) = makeAtom(0.5, 0.5, 0.5, atomParams, -1)
        
        

        NumSlots = NJ*ND*NB
        NumParams = NumSlots / MPI_num_procs
        if (NumSlots < MPI_num_procs .or. MPI_num_procs < 1) then 
                if (MPI_rank == 0) then 
                        write(error_unit,'(A)') "ERROR: Number of Parameter Slots must be greater than the number of MPI Ranks, either &
                        increase the number of slots or decrease the number of MPI Ranks for the given input."  
                end if 
                call MPI_Abort(MPI_COMM_WORLD,1,MPI_ierr)
        end if 
        if (MPI_Rank == MPI_num_procs - 1) NumParams = NumParams + mod(NumSlots,MPI_num_procs)
        
        BasePtr = MPI_rank*(NumSlots/MPI_num_procs)
        TopPtr = BasePtr + (NumParams-1)
 
        allocate(ParamIndexArray(NumParams),stat=stat)
        if (stat /= 0) error stop "Error: Falure to allocate ParamIndexArray"

        do i = 1,NumParams
                ParamIndexArray(i) = BasePtr + i
        end do 
        
        allocate(ParamArray(NumParams,3),stat=stat) 
        do i = 1,NumParams
                call indicesFromSlot(ParamIndexArray(i),NJ,ND,NB,JIndex,DIndex,BIndex)
                ParamArray(i,1) = dble(JIndex)/dble(NJ)*(Jmax - Jmin) + Jmin 
                ParamArray(i,2) = dble(DIndex)/dble(ND)*(Dmax - Dmin) + Dmin
                ParamArray(i,3) = dble(BIndex)/dble(NB)*(Bmax - Bmin) + Bmin
        end do 
        if (stat /= 0) error stop "Error allocating ParamArray"

        allocate(meshBuffer(NumParams,numTemps),stat=stat)
        if (stat /= 0) error stop "Error: Failed to allocate mesh buffer"
        
        do i = 1,NumParams
                do j = 1,numTemps 
                        meshBuffer(i,j) = makeChainMesh(2, numCellsX, numCellsY, numCellsZ, AtomsInUnitCell,&
                                                a_bravais,b_bravais,c_bravais,ab,bc,ca)
                end do 
        end do 

        ! Each temperature is being simulated on a specific mesh, we will keep track of this in the TemperatureMesh array 
        ! We will keep track of what temperature each index in the TemperatureMesh corresponds to in TemperatureArray
        
        allocate(TemperatureMeshArray(NumPArams,numTemps),stat=stat)
        if (stat /= 0) error stop "Error: Failed to allocate TemperatureMesh array" 

        allocate(TemperatureArray(numTemps),stat=stat)
        if (stat/=0) error stop "Error: Failed to allocate TemperatureArray"

        do i = 1, numParams 
                do j = 1, numTemps 
                        TemperatureMeshArray(i,j) = j ! Start out with a simple mapping, these will be swapped in parallel tempering later
                        TemperatureArray(j) = ((dble(j)/dble(numTemps)) * (TMax - TMin)) + TMin 
                end do  
        end do 
        print *, "All okay from rank", MPI_rank
        call MPI_Finalize() 
        

end program PT 
