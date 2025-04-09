program main 
        use Atom
        use chainMeshCell 
        use ChainMesh
        use vecNd
        use LLG
        use omp_lib
        use EnergyMin
        use constants, only: Kb, gyromagnetic_ratio, bohr_magneton
        implicit none
        integer :: numCellsX, numCellsY, numCellsZ, i, skyrmion_type, frame, num_frames, numMetropolisSteps
        type(ChainMesh_t) :: testMesh
        type(Atom_t) :: AtomsInUnitCell(2)
        type(vecNd_t) :: skyrmion_center
        real(kind=8), allocatable :: center_coords(:), skyrmion_radius
        real :: AtomParam1(3), AtomParam2(3), latticeParam
        character(len=100) :: output_dir, output_filename, frame_filename
        logical :: z_localized
        procedure(H_eff_class), pointer :: p
        integer :: argc, counter 
        character(len=50) :: arg         
        ! LLG evolution parameters
        real(kind=8) :: dt, total_time
        real(kind=8) :: H_field(3)

        ! Metropolis parameters 
        real(kind=8) :: betaMin, betaMax, beta, J, J_prime,Dz, Dz_prime, B, tempReal, T, Tmax, Tmin  
        integer(kind=OMP_LOCK_KIND), allocatable :: lockArray(:) 
        integer :: numBetaSteps 
        character(len=90) :: filepath_output
        
        
        argc = command_argument_count() 
        if (argc /= 6) error stop "Must have three command line arguments: J J_prime Dz Dz_prime B outputFile"
        call get_command_argument(1,arg)
        read(arg,*) J 
        call get_command_argument(2,arg)
        read(arg,*) J_prime 
        call get_command_argument(3,arg)
        read(arg,*) Dz 
        call get_command_argument(4,arg)
        read(arg,*) Dz_prime 
        call get_command_argument(5,arg)
        read(arg,*) B 
        call get_command_argument(6,arg)
        filepath_output = arg

        print *, "Comand line arguments: "
        print *, "J:",J 
        print *, "J':", J_prime 
        print *, "Dz:", Dz 
        print *, "Dz':", Dz_prime 
        print *, "B:", B 
        print *, "outDir:", filepath_output

        ! Initialize parameters
        latticeParam = 2.8
        ! Create 3D spin parameters for atoms (initialize all spins pointing up)
        AtomParam1 = (/0.0, 0.0, 1.0/)
        AtomParam2 = (/0.0, 0.0, 1.0/)
        
        ! Create atoms in unit cell
        AtomsInUnitCell(1) = makeAtom(0.0, 0.0, 0.0, AtomParam1, 3, -1) 
        AtomsInUnitCell(2) = makeAtom(latticeParam/2, latticeParam/2, latticeParam/2, AtomParam2, 3, -1)
        numCellsX = 80
        numCellsY = 80
        numCellsZ = 6
        ! Create the chain mesh
        testMesh = makeChainMesh(2, numCellsX, numCellsY, numCellsZ, latticeParam, AtomsInUnitCell)
        
        ! Assign nearest neighbors
        call assignNearestNeighbors(testMesh)
        call DerivativeList(testMesh,testMesh%derivativeList)       
        ! Define skyrmion center position (middle of the mesh)
        allocate(center_coords(3))
        center_coords(1) = numCellsX * latticeParam / 2.0d0
        center_coords(2) = numCellsY * latticeParam / 2.0d0
        center_coords(3) = numCellsZ * latticeParam / 2.0d0
        skyrmion_center = makeVecNd(center_coords)
        skyrmion_radius = 1.0_8*testMesh%latticeParameter 
        
        ! Initialize the skyrmion
        !call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,3.12_8/2.0_08,1)
        !call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,0.0_08,1)
        do i = 1, size(testMesh%atoms)
                if (any(testMesh%atoms(i)%AtomParameters /= testMesh%atoms(i)%AtomParameters)) then 
                        print *, "atom ", i, " has atom parameters ", testMesh%atoms(i)%AtomParameters
                        error stop "NaN encountered"
                end if 
        end do 
        ! Set up output directory
        !output_dir = "skyrmion_evolution"
        output_dir = filepath_output

        call system('mkdir -p ' // trim(output_dir))
        
        ! Write initial configuration
        !write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", 0, ".csv"
        !call write_spins_to_file(testMesh, frame_filename)
        
        H_field = (/0.0d0, 0.0d0, 0.1d0/)  ! External magnetic field in z-direction
        
        ! Time evolution parameters
        ! dt = 0.00000000000002_8
        dt = 1e-3
        total_time = 30.0d0
        num_frames = 20
        numMetropolisSteps = 9200
        numBetaSteps = 100
        
        ! Main evolution loop
        p => H_eff_Heisenberg
                
        allocate(lockArray(size(testMesh%atoms)))
        do i = 1,size(lockArray)
                call OMP_INIT_LOCK(lockArray(i))
        end do 
        counter = 1 
        do i = 0,numBetaSteps
                Tmax = 200.0_8 
                !Tmin = 0.1*(0.76*8*J)/(3*Kb)
                Tmin = 0.001_8
                T = Tmax - (Tmax - Tmin)*(dble(i)/dble(numBetaSteps)) 
                beta = 1.0_8 / (T)
        
                call MetropolisMixed(testMesh,beta,numMetropolisSteps,J,J_prime,Dz,Dz_prime,B, lockArray)
                if (mod(i,10) == 0) then 
                        write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", counter-1, ".csv"
                        call write_spins_to_file(testMesh, frame_filename)
                        print *, "Completed metropolis run at beta = ", beta 
                        
                        counter = counter + 1
                end if 
                
        end do 
        print *, "Completed Mixed Metropolis, beggining Heun evolution"
         do frame = 1, num_frames
             ! Calculate how many LLG steps to perform between frames
            
            
            call HeunStep(testMesh,10,dt,p,0.8_08,gyromagnetic_ratio, J, Dz, B)
            
                            
             ! Write current state to file
             print *, "Frame = ", frame, "Counter = ", counter
             write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", frame + counter - 2, ".csv"
             call write_spins_to_file(testMesh, frame_filename)
            
             print *, "Completed frame", frame, "of", num_frames
         end do
        
        ! Write information for Python visualization script
        output_filename = trim(output_dir) // "/info.txt"
        open(unit=10, file=output_filename, status='replace')
        !write(10, *) num_frames + 1  ! Total number of frames (including initial frame)
        write(10, *) counter + num_frames - 1  ! Total number of frames (including initial frame)
        write(10, *) numCellsX, latticeParam     ! Mesh parameters
        close(10)
        
        print *, "Skyrmion evolution completed. Output files in:", trim(output_dir)
        print *, "Run your Python visualization script to create the animation."
        
        ! Clean up memory
        deallocate(center_coords)
        call deallocateChainMesh(testMesh)
                
        do i = 1,size(lockArray)
                call OMP_DESTROY_LOCK(lockArray(i))
        end do 
end program main
