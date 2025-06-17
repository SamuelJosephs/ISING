! This file will contain an MPI-openMP program for parallel tempering of the Heisenberg model as outlined in
! https://doi.org/10.1016/j.cor.2025.107000


program PT 
        use mpi_f08
        use PT_Utils
        use Rand

        use chainMesh
        use Atom
        use EnergyMin
        use constants
        use iso_fortran_env, only: error_unit
        use omp_lib
        implicit none 
        integer :: MPI_ierr, MPI_rank, MPI_num_procs 
        integer, parameter :: NJ = 4 ! The number of J, D, and B values to perform the parallel tempering at.  
        integer, parameter :: ND = 5 ! Hard code these for now but eventually they should be taken as input
        integer, parameter :: NB = 1 
        real(kind=dp), parameter :: JMin = 0.0_dp
        real(kind=dp), parameter :: JMax = 2.5_dp 
        real(kind=dp), parameter :: DMin = 0.0_dp 
        real(kind=dp), parameter :: DMax = 2.5_dp
        real(kind=dp), parameter :: BMin = 1.0_dp 
        real(kind=dp), parameter :: Bmax = 1.0_dp

        real(kind=dp), parameter :: TMax = 2.50_dp 
        real(kind=dp), parameter :: TMin = 0.01_dp 
        integer, parameter :: numTemps = 5
        ! Set up constants for the lattice, for now they will be hardcoded but eventually they should be taken as input.
        type(Atom_t), dimension(2) :: atomsInUnitCell
        real, dimension(3), parameter :: atomParams = (/1.0, 0.0, 0.0/)
        integer, parameter :: numCellsX = 30 
        integer, parameter :: numCellsY = 30
        integer, parameter :: numCellsZ = 6
        real(kind=dp), parameter :: a_bravais = 2.8
        real(kind=dp), parameter :: b_bravais = 2.8
        real(kind=dp), parameter :: c_bravais = 2.8 
        real(kind=dp), parameter :: ab = 90
        real(kind=dp), parameter :: bc = 90 
        real(kind=dp), parameter :: ca = 90

        integer, parameter :: numSwaps = numTemps ! Number of swaps to do per iteration
        integer, parameter :: numIterations = 20 
        integer, parameter :: numMCSSweepsPerSwap = 300
         
        integer :: NumSlots, BasePtr, TopPtr, NumParams, Iteration, meshIndex, swapIndex
        integer :: meshIndex1, meshIndex2 ! intuitive naming requires more variables than are strictly needed
        integer :: stat, i, j, JIndex, DIndex, BIndex
        integer, allocatable, dimension(:) :: ParamIndexArray ! Contains the global index of each parameter set
        real(kind=dp), allocatable, dimension(:,:) :: ParamArray ! (paramIndex, (J,D,B))
        ! For each Parameter set we need a buffer for the spins 
        type(ChainMesh_t), allocatable, dimension(:,:) :: meshBuffer ! (paramIndex, Temp index) 
        
        integer, allocatable, dimension(:,:) :: TemperatureMeshArray ! (paramIndex, Temp index)
        real(kind=dp), allocatable, dimension(:) :: TemperatureArray ! (Temp)

        real(kind=dp) :: beta, J_H,D_H,B_H ! underscore for Heisenberg

        integer(kind=OMP_LOCK_KIND), allocatable, dimension(:) :: lockArray
        type(random) :: rand_gen
        real(kind=dp) :: u, E1, E2, beta1, beta2, Delta, z, p, temp
        integer :: Index1, Index2, tempInt

        integer :: fileunit
        type(MPI_FILE) :: mpi_file_handle 
        type(MPI_INFO) :: mpi_info_handle 
        type(MPI_STATUS) :: mpi_status_handle
        integer :: mpi_file_ierr
        character(len=:), allocatable :: output_string
        character(len=500) :: string_buff
        integer :: mpi_mode
        integer :: skyrmion_number_middle 
        real(kind=dp) :: winding_number_middle, winding_number_spread
        type(vecNd_t) :: magnetisation
        call MPI_Init(MPI_ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD,MPI_rank)
        call MPI_Comm_Size(MPI_COMM_WORLD,MPI_num_procs)

        
        AtomsInUnitCell(1) = makeAtom(0.0, 0.0, 0.0, atomParams, -1) 
        AtomsInUnitCell(2) = makeAtom(0.5, 0.5, 0.5, atomParams, -1)
        
        

        NumSlots = NJ*ND*NB
        NumParams = NumSlots / MPI_num_procs
        if (NumSlots < MPI_num_procs .or. MPI_num_procs < 1) then 
                if (MPI_rank == 0) then 
                        write(error_unit,'(A)') "ERROR: Number of Parameter Slots must be greater than the number &
                                                of MPI Ranks, either increase the number of slots or decrease the number &
                                                        of MPI Ranks for the given input."  
                end if 
                call MPI_Abort(MPI_COMM_WORLD,1,MPI_ierr)
        end if 
        if (MPI_Rank == MPI_num_procs - 1) NumParams = NumParams + mod(NumSlots,MPI_num_procs)
        
        BasePtr = MPI_rank*(NumSlots/MPI_num_procs)
        TopPtr = BasePtr + (NumParams-1)
 
        allocate(ParamIndexArray(NumParams),stat=stat)
        if (stat /= 0) error stop "Error: Falure to allocate ParamIndexArray"

        do i = 1,NumParams
                ParamIndexArray(i) = BasePtr + i ! The global index of each swap
        end do 
        
        allocate(ParamArray(NumParams,3),stat=stat) 
        do i = 1,NumParams
                call indicesFromSlot(ParamIndexArray(i),NJ,ND,NB,JIndex,DIndex,BIndex)
                print *, "MPI_rank ", MPI_rank, "Has Jindex,DIndex,BIndex = ", JIndex,DIndex,BIndex
                ParamArray(i,1) = dble(JIndex)/dble(NJ)*(Jmax - Jmin) + Jmin 
                ParamArray(i,2) = dble(DIndex)/dble(ND)*(Dmax - Dmin) + Dmin
                ParamArray(i,3) = dble(BIndex)/dble(NB)*(Bmax - Bmin) + Bmin
        end do 
        print *, "DEBUG, MPI_rank ", MPI_rank, "Has paramArray: ", ParamArray(1,:)
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
                        TemperatureArray(j) = TMin*(TMax/TMin)**(dble(j-1)/dble(numTemps-1))
                end do  
        end do 

        allocate(lockArray(meshBuffer(1,1)%numAtoms),stat=stat) ! All meshes have the same number of atoms
        if (stat /= 0) error stop "Error: Failed to allocate lock array"

        do i = 1,size(lockArray)
                call OMP_INIT_LOCK(lockArray(i))
        end do 

        ! Initialse and warm up random number generator
        rand_gen = makeRandom(MPI_rank)
        do i = 1,50
                u = algor_uniform_random(rand_gen)
        end do 

        do Iteration = 1,numIterations
                ! First need to perform the specified number of MCS sweeps on each chain mesh on this rank 
                
                do i = 1,numParams 
                                do j = 1,numTemps 
                                        beta = 1.0_dp / TemperatureArray(j)
                                        meshIndex = TemperatureMeshArray(i,j) ! The j'th temperature is being computed at this index
                                        J_H = ParamArray(i,1)
                                        D_H = ParamArray(i,2)
                                        B_H = ParamArray(i,3)
                                        call Metropolis_mcs(meshBuffer(i,meshIndex),beta,numMCSSweepsPerSwap,&
                                                J_H,0.0_8,D_H,0.0_8,B_H,0.2_8, lockArray,demag=.False.)  
                                end do 
                end do 

                ! After the MCS updates we must attempt to swap adjacent temperatures 
                do i = 1,numParams 
                        do swapIndex = 1,numSwaps 
                                Index1 = nint(algor_uniform_random(rand_gen)*(numtemps-1)) + 1
                                if (Index1 == 1) then 
                                        Index2 = 2
                                else if (Index1 == numTemps) then
                                        Index2 = numtemps - 1
                                else 
                                        Index2 = Index1 + 1
                                end if 

                                J_H = ParamArray(i,1)
                                D_H = ParamArray(i,2)
                                B_H = ParamArray(i,3)
                                meshIndex1 = TemperatureMeshArray(i,Index1) ! Need to find the meshes at the temperatures indexed by 
                                meshIndex2 = TemperatureMeshArray(i,Index2) ! Index1 and Index2
                                call totalHeisenbergEnergy(meshBuffer(i,meshIndex1),J_H,0.0_dp,D_H,0.0_dp,B_H,lockArray,E1)
                                call totalHeisenbergEnergy(meshBuffer(i,meshIndex2),J_H,0.0_dp,D_H,0.0_dp,B_H,lockArray,E2)

                                beta1 =  1.0_dp / TemperatureArray(Index1)
                                beta2 =  1.0_dp / TemperatureArray(Index2)

                                Delta = (beta2 - beta1)*(E1 - E2)
                                u = algor_uniform_random(rand_gen)
                                if (u <= 0.0_dp) u = 1e-10_dp
                                u = log(u)

                                if (u < Delta) then 
                                        print *, "MPI_rank", MPI_rank, "Has accepted a swap between temperatures:", index1, index2
                                        tempInt = TemperatureMeshArray(i,Index1)
                                        TemperatureMeshArray(i,index1) = TemperatureMeshArray(i,Index2) 
                                        TemperatureMeshArray(i,index2) = tempInt
                                        
                                end if 

                                print *, "TempuratureMeshArray from MPI Rank:", MPI_rank, " = ", TemperatureMeshArray(i,:)
                        end do  
                end do 
        end do

        ! Now need to collect statistics from each slot and write them to a file

        
        if (MPI_rank == 0) then 
                open(newunit=fileunit,status="replace",action="write",file="output.csv") ! Just to wipe the file
                ! write csv header
                write(fileunit, '(A)') "J,D,B,T,winding_number_middle,skyrmion_number_middle,winding_number_spread,&
                                        mx,my,mz"
                close(unit=fileunit)

        end if 

        call MPI_barrier(MPI_COMM_WORLD)


        ! TODO calculate statistics and write them to a file

        print *, "All okay from rank", MPI_rank
        mpi_info_handle = MPI_INFO_NULL 
        mpi_mode = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        mpi_mode = ior(MPI_MODE_APPEND,mpi_mode)
        call MPI_FILE_OPEN(MPI_COMM_WORLD,"output.csv",mpi_mode,mpi_info_handle, &
                        mpi_file_handle,mpi_file_ierr)
        if (mpi_file_ierr /= 0) error stop "Error: Failed to open output file"
        write(string_buff,'(I8.1)') MPI_rank  
        output_string = "Hello from MPI_rank " // trim(adjustl(string_buff)) // new_line('a') 

        string_buff = " "
        output_string = ""
        do i = 1,NumParams
                do j = 1,numTemps
                        string_buff = " "
                        meshIndex = TemperatureMeshArray(i,j)
                        call chainMesh_statistics(meshBuffer(i,meshIndex),skyrmion_number_middle,winding_number_middle,&
                                winding_number_spread, magnetisation)
                        temp = TemperatureArray(j)
                        ! write model parameters J, D, B, T
                        write(string_buff,'((F8.4,",",F8.4,","F8.4,",",F8.4,","))') ParamArray(i,1), &
                                        ParamArray(i,2), ParamArray(i,3),temp
                        output_string = output_string // trim(adjustl(string_buff))
                        string_buff = " "
                        ! write model statistics 
                        write(string_buff,'((F8.4,",",I8.4,",",F8.2,",",F8.2,",",F8.2,",",F8.2))') winding_number_middle, &
                                skyrmion_number_middle,winding_number_spread,magnetisation%coords(1),&
                                magnetisation%coords(2), magnetisation%coords(3)
                        output_string = output_string // string_buff // new_line('a')

                end do 
        end do 
        call MPI_FILE_WRITE_SHARED(mpi_file_handle,output_string,len(output_string),MPI_CHARACTER,mpi_status_handle,mpi_file_ierr)
        ! Cleanup 
        do i = 1,size(lockArray)
                call OMP_DESTROY_LOCK(lockArray(i))
        end do 
        call MPI_Finalize() 
        

end program PT 
