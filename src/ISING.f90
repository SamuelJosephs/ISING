program main 
        use Atom
        use chainMeshCell 
        use ChainMesh
        use vecNd
        use LLG
        use omp_lib
        implicit none
        integer :: numCells, i, skyrmion_type, frame, num_frames
        character(len=30) :: arg
        type(ChainMesh_t) :: testMesh
        type(Atom_t) :: AtomsInUnitCell(2)
        type(vecNd_t) :: skyrmion_center
        real(kind=8), allocatable :: center_coords(:), skyrmion_radius
        real :: AtomParam1(3), AtomParam2(3), latticeParam
        character(len=100) :: output_dir, output_filename, frame_filename
        logical :: z_localized
        procedure(H_eff_class), pointer :: p
        
        ! LLG evolution parameters
        real(kind=8) :: dt, total_time, A, B, C, D
        real(kind=8) :: H_field(3)
        
        ! Initialize parameters
        latticeParam = 2.8
        ! Create 3D spin parameters for atoms (initialize all spins pointing up)
        AtomParam1 = (/0.0, 0.0, 1.0/)
        AtomParam2 = (/0.0, 0.0, -1.0/)
        
        ! Create atoms in unit cell
        AtomsInUnitCell(1) = makeAtom(0.0, 0.0, 0.0, AtomParam1, 3, -1) 
        AtomsInUnitCell(2) = makeAtom(latticeParam/2, latticeParam/2, latticeParam/2, AtomParam2, 3, -1)
        numCells = 15
        ! Create the chain mesh
        testMesh = makeChainMesh(2, numCells, numCells, numCells, latticeParam, AtomsInUnitCell)
        
        ! Assign nearest neighbors
        call assignNearestNeighbors(testMesh)
       
        ! Define skyrmion center position (middle of the mesh)
        allocate(center_coords(3))
        center_coords(1) = numCells * latticeParam / 2.0d0
        center_coords(2) = numCells * latticeParam / 2.0d0
        center_coords(3) = numCells * latticeParam / 2.0d0
        skyrmion_center = makeVecNd(center_coords)
        skyrmion_radius = 1*testMesh%latticeParameter 
        
        ! Initialize the skyrmion
        call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,3.12_8/2.0_08,1)
        !call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,0.0_8,2)
        do i = 1, size(testMesh%atoms)
                if (any(testMesh%atoms(i)%AtomParameters /= testMesh%atoms(i)%AtomParameters)) then 
                        print *, "atom ", i, " has atom parameters ", testMesh%atoms(i)%AtomParameters
                        error stop "NaN encountered"
                end if 
        end do 
        ! Set up output directory
        output_dir = "skyrmion_evolution"
        call system('mkdir -p ' // trim(output_dir))
        
        ! Write initial configuration
        write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", 0, ".csv"
        call write_spins_to_file(testMesh, frame_filename)
        
        ! Set LLG parameters
        A = 1.0d0  ! Exchange parameter for x
        B = 1.0d0  ! Exchange parameter for y
        C = 1.0d0  ! Exchange parameter for z
        D = 0.1d0  ! Anisotropy parameter
        H_field = (/0.0d0, 0.0d0, 0.1d0/)  ! External magnetic field in z-direction
        
        ! Time evolution parameters
        dt = 0.01d0
        total_time = 30.0d0
        num_frames = 30
        
        ! Main evolution loop
        print *, "Starting skyrmion evolution..."
        p => H_eff_Heisenberg
        do frame = 1, num_frames
            ! Calculate how many LLG steps to perform between frames
            do i = 1, int(total_time / dt / num_frames)
                ! Evolve the system using LLG equation
                ! call LLGStep(testMesh, dt, A, B, C, D, H_field)
                call HeunStep(testMesh,4,dt,p,1.0_8,1.0_8)
            end do
            
            ! Write current state to file
            write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", frame, ".csv"
            call write_spins_to_file(testMesh, frame_filename)
            
            print *, "Completed frame", frame, "of", num_frames
        end do
        
        ! Write information for Python visualization script
        output_filename = trim(output_dir) // "/info.txt"
        open(unit=10, file=output_filename, status='replace')
        write(10, *) num_frames + 1  ! Total number of frames (including initial frame)
        write(10, *) numCells, latticeParam     ! Mesh parameters
        close(10)
        
        print *, "Skyrmion evolution completed. Output files in:", trim(output_dir)
        print *, "Run your Python visualization script to create the animation."
        
        ! Clean up memory
        deallocate(center_coords)
        call deallocateChainMesh(testMesh)
end program main
