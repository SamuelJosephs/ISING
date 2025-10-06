! This file will contain an MPI-openMP program for parallel tempering of the Heisenberg model as outlined in
! https://doi.org/10.1016/j.cor.2025.107000


program PT 
        use mpi_f08
        use PT_Utils
        use Rand
        use LLG
        use chainMesh
        use Atom
        use EnergyMin
        use constants
        use io
        use iso_fortran_env, only: error_unit, dp=>real64
        use omp_lib
        implicit none 
        integer :: MPI_ierr, MPI_rank, MPI_num_procs 
        integer :: NJ = 10 ! The number of J, D, and B values to perform the parallel tempering at.  
        integer :: ND = 10 ! Hard code these for now but eventually they should be taken as input
        integer :: NB = 1 
        real(kind=dp) :: JMin = -2.5_dp
        real(kind=dp) :: JMax = 2.5_dp 
        real(kind=dp) :: DMin = -2.5_dp 
        real(kind=dp) :: DMax = 2.5_dp
        real(kind=dp) :: BMin = 1.5_dp 
        real(kind=dp) :: BMax = 1.5_dp
        logical :: demag = .False.
        real(kind=dp), parameter :: TMax = 2000_dp 
        real(kind=dp), parameter :: TMin = 100_dp 
        integer, parameter :: numTemps = 10
        ! Set up constants for the lattice, for now they will be hardcoded but eventually they should be taken as input.
        type(Atom_t), dimension(2) :: atomsInUnitCell
        real, dimension(3), parameter :: atomParams = (/1.0, 0.0, 0.0/)
        integer, parameter :: numCellsX = 30 
        integer, parameter :: numCellsY = 30
        integer, parameter :: numCellsZ = 6
        real(kind=dp):: a_bravais = 2.8
        real(kind=dp) :: b_bravais = 2.8
        real(kind=dp) :: c_bravais = 2.8 
        real(kind=dp) :: ab = 90
        real(kind=dp) :: bc = 90 
        real(kind=dp) :: ca = 90

        integer, parameter :: numSwaps = numTemps ! Number of swaps to do per iteration
        integer, parameter :: numIterations = 100 
        integer, parameter :: numMCSSweepsPerSwap = 1000
         
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
        character(len=:), allocatable :: output_string, filename_string, filename_string_spins,&
                filename_string_density
        character(len=500) :: string_buff
        integer :: mpi_mode
        integer :: skyrmion_number_middle 
        real(kind=dp) :: winding_number_middle, winding_number_spread
        type(vecNd_t) :: magnetisation

        integer :: numArgs, fileHandle 
        character(len=400) :: outputPath, outputPathSpins, outputPathDensity 
        call MPI_Init(MPI_ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD,MPI_rank)
        call MPI_Comm_Size(MPI_COMM_WORLD,MPI_num_procs)

        numArgs = command_argument_count()

        !call io_parsefile("testInput.txt")

        
        if (numArgs /= 10) then 
                write(error_unit,'(A)') "Incorrect number of command line arguments, &
                        input: Jmin Jmax NJ Dmin Dmax ND Bmin Bmax NB outputFilePath"
                call MPI_abort(MPI_COMM_WORLD,1)
        end if 

        ! Not elegant but a fancy parsing system is not a good use of time 
        string_buff = " " 
        call GET_COMMAND_ARGUMENT(1,string_buff)
        read(string_buff,*) JMin

        string_buff = " "
        call GET_COMMAND_ARGUMENT(2,string_buff)
        read(string_buff,*) JMax 

        string_buff = " "
        call GET_COMMAND_ARGUMENT(3,string_buff)
        read(string_buff,*) NJ 

        string_buff = " " 
        call GET_COMMAND_ARGUMENT(4,string_buff)
        read(string_buff,*) DMin 

        string_buff = " " 
        call GET_COMMAND_ARGUMENT(5,string_buff)
        read(string_buff,*) DMax 

        string_buff = " " 
        call GET_COMMAND_ARGUMENT(6,string_buff)
        read(string_buff,*) ND
        

        string_buff = " " 
        call GET_COMMAND_ARGUMENT(7,string_buff)
        read(string_buff,*) BMin

        string_buff = " " 
        call GET_COMMAND_ARGUMENT(8,string_buff)
        read(string_buff,*) BMax
 
        string_buff = " " 
        call GET_COMMAND_ARGUMENT(9,string_buff)
        read(string_buff,*) NB
        
        outputPath = " " 
        call GET_COMMAND_ARGUMENT(10,outputPath)
        
        
        AtomsInUnitCell(1) = makeAtom(0.0,  0.0,  0.0, -1)
        AtomsInUnitCell(2) = makeAtom(0.5,  0.5,  0.5, -1)
        
        a_bravais = 1.0_dp
        b_bravais = 1.0_dp
        c_bravais = 1.0_dp
        
        ab = 90.0_dp
        bc = 90.0_dp
        ca = 90.0_dp 

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
                ParamArray(i,1) = dble(JIndex - 1)/dble(NJ)*(Jmax - Jmin) + Jmin 
                ParamArray(i,2) = dble(DIndex - 1)/dble(ND)*(Dmax - Dmin) + Dmin
                ParamArray(i,3) = dble(BIndex - 1)/dble(NB)*(Bmax - Bmin) + Bmin
        end do 
        print *, "DEBUG, MPI_rank ", MPI_rank, "Has paramArray: ", ParamArray(1,:)
        if (stat /= 0) error stop "Error allocating ParamArray"

        allocate(meshBuffer(NumParams,numTemps),stat=stat)
        if (stat /= 0) error stop "Error: Failed to allocate mesh buffer"
        
        do i = 1,NumParams
                do j = 1,numTemps 
                        meshBuffer(i,j) = makeChainMesh(numCellsX, numCellsY, numCellsZ, AtomsInUnitCell,&
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
                print *, "MPI_rank", MPI_rank, "is starting iteration", Iteration, "out of", numIterations 
                do i = 1,numParams 
                                do j = 1,numTemps 
                                        beta = 1.0_dp / (kb*TemperatureArray(j))
                                        meshIndex = TemperatureMeshArray(i,j) ! The j'th temperature is being computed at this index
                                        J_H = ParamArray(i,1)
                                        D_H = ParamArray(i,2)
                                        B_H = ParamArray(i,3)
                                        call Metropolis_mcs(meshBuffer(i,meshIndex),beta,numMCSSweepsPerSwap,&
                                                J_H,0.0_8,D_H,0.0_8,B_H,0.2_8, lockArray,demag=demag)  
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
                                call totalHeisenbergEnergy(meshBuffer(i,meshIndex1),J_H,0.0_dp,D_H,0.0_dp,B_H,lockArray,E1,&
                                        demag=demag)
                                call totalHeisenbergEnergy(meshBuffer(i,meshIndex2),J_H,0.0_dp,D_H,0.0_dp,B_H,lockArray,E2,&
                                        demag=demag)

                                beta1 =  1.0_dp / TemperatureArray(Index1)
                                beta2 =  1.0_dp / TemperatureArray(Index2)

                                Delta = (beta2 - beta1)*(E1 - E2)
                                u = algor_uniform_random(rand_gen)
                                if (u <= 0.0_dp) u = 1e-10_dp
                                u = log(u)

                                if (u < Delta) then 
                                        tempInt = TemperatureMeshArray(i,Index1)
                                        TemperatureMeshArray(i,index1) = TemperatureMeshArray(i,Index2) 
                                        TemperatureMeshArray(i,index2) = tempInt
                                        
                                end if 

                        end do  
                end do 
        end do

        ! Finally relax each domain
        do i = 1,numParams 
                do j = 1,numTemps 
                
                beta = 1.0_dp / TemperatureArray(j)
                meshIndex = TemperatureMeshArray(i,j) ! The j'th temperature is being computed at this index
                J_H = ParamArray(i,1)
                D_H = ParamArray(i,2)
                B_H = ParamArray(i,3)
                call Metropolis_mcs(meshBuffer(i,meshIndex),beta,1000,&
                        J_H,0.0_8,D_H,0.0_8,B_H,0.2_8, lockArray,demag=demag)  
                end do 
        end do 
        ! Now need to collect statistics from each slot and write them to a file

        print *, "MPI_rank: ", MPI_rank, "Has reached the barrier"
        call MPI_barrier(MPI_COMM_WORLD)


        ! TODO calculate statistics and write them to a 
        print *, "All okay from rank", MPI_rank

        string_buff = " "
        output_string = ""
        outputPathSpins = ' '
        outputPathSpins = trim(adjustl(outputPath)) // "_spins"
        outputPathDensity = trim(adjustl(outputPath)) // "_density"
        if (MPI_RANK == 0) call EXECUTE_COMMAND_LINE("mkdir -p " // trim(adjustl(outputPathSpins)))
                
        if (MPI_RANK == 0) call EXECUTE_COMMAND_LINE("mkdir -p " // trim(adjustl(outputPathDensity)))

        do i = 1,NumParams
                do j = 1,numTemps
                        string_buff = " "
                        meshIndex = TemperatureMeshArray(i,j)
                        call chainMesh_statistics(meshBuffer(i,meshIndex),skyrmion_number_middle,winding_number_middle,&
                                winding_number_spread, magnetisation)
                        temp = TemperatureArray(j)
                        ! write model parameters J, D, B, T
                        output_string = "J,D,B,T,winding_number_middle, &
                                skyrmion_number_middle,winding_number_spread,sx,&
                                sy, sz" // new_line('a')
                        write(string_buff,'((F0.4,",",F0.4,","F0.4,",",F0.10,","))') ParamArray(i,1), &
                                        ParamArray(i,2), ParamArray(i,3),temp
                        output_string = output_string // trim(adjustl(string_buff))
                        string_buff = " "
                        ! write model statistics 
                        write(string_buff,'((F0.4,",",I0.4,",",F0.2,",",F0.2,",",F0.2,",",F0.2))') winding_number_middle, &
                                skyrmion_number_middle,winding_number_spread,magnetisation%coords(1),&
                                magnetisation%coords(2), magnetisation%coords(3)
                        output_string = output_string // string_buff // new_line('a')

                        print *, "J,D,B,T, skyrmion_number_middle, winding_number_middle, winding_number_spread = ", &
                        ParamArray(i,1), ParamArray(i,2), ParamArray(i,3),temp, skyrmion_number_middle, winding_number_middle,&
                                                                 winding_number_spread
                        string_buff = " "
                        write(string_buff,'((F0.4,"_",F0.4,"_",F0.4,"_",e0.4))') ParamArray(i,1), ParamArray(i,2), & 
                                ParamArray(i,3), temp
                        filename_string = ' '
                        filename_string = trim(adjustl(outputPath)) // "/data_"& 
                                // trim(adjustl(string_buff)) // ".csv"
                        filename_string_spins = ' '
                        filename_string_spins = trim(adjustl(outputPathSpins)) // "/"&
                                // "spins_" // trim(adjustl(string_buff)) // ".csv"
                        print *, "MPI_RANK: ", MPI_rank, "Has filename_string: ", filename_string

                        filename_string_density = trim(adjustl(outputPathDensity)) // "/" // "density_" &
                                // trim(adjustl(string_buff)) // ".csv"

                        call write_spins_to_file(meshBuffer(i,meshIndex),&
                                filename_string_spins)

                        call write_winding_number_density(meshBuffer(i,meshIndex),filename_string_density)
                        open(newunit=fileHandle, status="replace", action="write",file=filename_string)

                        write(fileHandle, '(A)') output_string

                        close(fileHandle)
                end do 
        end do 
        ! Cleanup 
        do i = 1,size(lockArray)
                call OMP_DESTROY_LOCK(lockArray(i))
        end do 
        call MPI_Finalize() 
        

end program PT 
