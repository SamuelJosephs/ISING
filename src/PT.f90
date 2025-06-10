! This file will contain an MPI-openMP program for parallel tempering of the Heisenberg model as outlined in
! https://doi.org/10.1016/j.cor.2025.107000


program PT 
        use mpi_f08
        use PT_Utils
        use iso_fortran_env, only: error_unit, real64
        implicit none 
        integer :: MPI_ierr, MPI_rank, MPI_num_procs 
        integer, parameter :: dp = real64
        integer, parameter :: NJ = 2 ! The number of J, D, and B values to perform the parallel tempering at.  
        integer, parameter :: ND = 2 ! Hard code these for now but eventually they should be taken as input
        integer, parameter :: NB = 1 
        real(kind=dp), parameter :: MinJ = 0.0_dp
        real(kind=dp), parameter :: MaxJ = 1.5_dp 
        real(kind=dp), parameter :: MinD = 0.0_dp 
        real(kind=dp), parameter :: MaxD = 1.5_dp

        integer :: NumSlots, BasePtr, TopPtr, NumParams 
        integer :: stat, i, JIndex, DIndex, BIndex
        integer, allocatable, dimension(:) :: ParamIndexArray
        call MPI_Init(MPI_ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD,MPI_rank)
        call MPI_Comm_Size(MPI_COMM_WORLD,MPI_num_procs)


        


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
        if (stat /= 0) error stop "Error: Falure to allocate"

        do i = 1,NumParams
                ParamIndexArray(i) = BasePtr + i
        end do 
        do i = 1,Numparams 
                call indicesFromSlot(ParamIndexArray(i),NJ,ND,NB,JIndex,DIndex,BIndex)
                print *, "Rank, ", MPI_rank, "Has J,D,B indices: ", JIndex, DIndex, BIndex
        end do 
        call MPI_Finalize()
        

end program PT 
