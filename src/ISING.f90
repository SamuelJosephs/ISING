program main 
        use Atom
        use chainMeshCell 
        use ChainMesh
        use vecNd
        use LLG
        use omp_lib
        use EnergyMin
        use reciprocal_space_processes
        use constants, only: Kb, gyromagnetic_ratio, bohr_magneton
        use iso_fortran_env, only : output_unit
        implicit none
        integer :: numCellsX, numCellsY, numCellsZ, i, skyrmion_type, frame, num_frames, numMetropolisStepsTotal, numMetropolisSteps
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
        real(kind=8) :: totalEnergy1, totalEnergy2, a_bravais,b_bravais,c_bravais, ab, bc, ca
        ! LLG evolution parameters
        real(kind=8) :: dt, total_time
        real(kind=8) :: H_field(3)

        ! Metropolis parameters 
        real(kind=8) :: betaMin, betaMax, beta, J, J_prime,Dz, Dz_prime, B, tempReal, T, Tmax, Tmin
        real(kind=8) :: winding_number_middle, winding_number_bottom, winding_number_top, winding_number_2
        real(kind=8), dimension(3) :: winding_number_array
        integer(kind=OMP_LOCK_KIND), allocatable :: lockArray(:) 
        integer :: numBetaSteps, fftw_status
        character(len=90) :: filepath_output

        real(kind=8), allocatable, dimension(:,:,:) :: test_grad_array
        
        integer ::  skyrmion_index, num_thresholds, skyrmion_number
        integer, allocatable, dimension(:) :: winding_array
        real(kind=8) :: lower_bound, upper_bound
        fftw_status = fftw_init_threads()
        if (fftw_status == 0) error stop "Error initialising fftw threads"
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
        AtomsInUnitCell(2) = makeAtom(0.5, 0.5, 0.5, AtomParam2, 3, -1)
        numCellsX = 40
        numCellsY = 40
        numCellsZ = 6
        a_bravais = latticeParam
        b_bravais = latticeParam
        c_bravais = latticeParam
        ab = 90 
        bc = 90 
        ca = 90 
        ! Create the chain mesh
        testMesh = makeChainMesh(2, numCellsX, numCellsY, numCellsZ, AtomsInUnitCell,&
                                a_bravais,b_bravais,c_bravais,ab,bc,ca)
        print *, "Attempting to allocate ", testMesh%numAtoms, "atoms"
        ! Assign nearest neighbors
        do i = 1,size(testMesh%atoms) 
                print *, "Num nearest Neighbours = ", size(testMesh%atoms(i)%NeighborList)
        end do 

        ! Define skyrmion center position (middle of the mesh)
        allocate(center_coords(3))
        center_coords(1) = numCellsX * latticeParam / 2.0d0
        center_coords(2) = numCellsY * latticeParam / 2.0d0
        center_coords(3) = numCellsZ * latticeParam / 2.0d0
        skyrmion_center = makeVecNd(center_coords)
        skyrmion_radius = 1.0_8*testMesh%a 
        
        ! Initialize the skyrmion
        call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,3.12_8/2.0_08,1)
        winding_number_middle = calculate_winding_number2(testMesh,testMesh%numCellsZ / 2)
        print *, "Test winding number = ", winding_number_middle 
        print *, "*************************************"
        !call initialise_skyrmion_sp(testMesh, skyrmion_center, skyrmion_radius,0.0_08,1)

        
        do i = 1, size(testMesh%atoms)
                if (any(testMesh%atomSpins(i,:) /= testMesh%atomSpins(i,:))) then 
                        print *, "atom ", i, " has atom parameters ", testMesh%atomSpins(i,:)
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
        dt = 1e-17_08
        total_time = 30.0d0
        num_frames = 0
        numMetropolisStepsTotal = 220000
        numMetropolisSteps = 2000
        numBetaSteps = 20
        
        lower_bound = 0.000001
        upper_bound = 0.95
        num_thresholds = 60
        ! Main evolution loop
                
        allocate(lockArray(size(testMesh%atoms)))
        do i = 1,size(lockArray)
                call OMP_INIT_LOCK(lockArray(i))
        end do 
        counter = 1 
        do i = 0,numBetaSteps
                Tmax = 5.0_8 
                !Tmin = 0.1*(0.76*8*J)/(3*Kb)
                Tmin = 0.0000001
                T = Tmax - (Tmax - Tmin)*(dble(i)/dble(numBetaSteps)) 
                beta = 1.0_8 / (T)
                call TotalHeisenbergEnergy(testMesh,J,J_prime,Dz,Dz_prime,B,lockArray,totalEnergy1)
                call Metropolis_mcs(testMesh,beta,numMetropolisSteps,&
                                                J,J_prime,Dz,Dz_prime,B,0.2_8, lockArray,demag=.True.)
                call TotalHeisenbergEnergy(testMesh,J,J_prime,Dz,Dz_prime,B,lockArray,totalEnergy2)
                winding_number_middle = calculate_winding_number2(testMesh,testMesh%numCellsZ/2)
                winding_number_bottom = calculate_winding_number2(testMesh,1)
                winding_number_top = calculate_winding_number2(testMesh,testMesh%numCellsZ)
                winding_number_array = [winding_number_bottom, winding_number_middle, winding_number_top]

                print *, "Winding Numbers = ", winding_number_array 
                print *, "Range of winding numbers = ", maxval(winding_number_array) - minval(winding_number_array)

                print *, "Delta E = ", totalEnergy2 - totalEnergy1, "T = ", T, "oldEnergy, newEnergy = ", totalEnergy1, totalEnergy2
                !call compute_skyrmion_distribution(testMesh,3,winding_array,lower_bound,upper_bound,&
                !                num_thresholds,testmesh%numCellsZ / 2)
                ! skyrmion_number = calculate_skyrmion_number(testMesh,testMesh%numCellsZ/2,0.3_8,1,0.0_8)      
                ! print *, "Skyrmion number = ", skyrmion_number
                if (mod(i,10) == 0) then 
                        if (i /= 0) then
                                call compute_skyrmion_distribution(testMesh,3,winding_array,lower_bound,upper_bound,&
                                       num_thresholds,testmesh%numCellsZ / 2)                 
                                call write_winding_number_density(testMesh,"./density_array_3d.csv")
                                ! skyrmion_number = calculate_skyrmion_number(testMesh,testMesh%numCellsZ/2,0.3_8,1,0.0_8)      
                                ! print *, "Skyrmion number = ", skyrmion_number
                                print *, "Skyrmion Distribution = ", winding_array
                                call write_demagnetisation_field_to_file(testMesh,&
                                                        testMesh%demagnetisation_array,"./demagnetisation_field.csv")
                        end if
                        write(frame_filename, '(A,A,I5.5,A)') trim(output_dir), "/frame_", counter-1, ".csv"
                        call write_spins_to_file(testMesh, frame_filename)
                        print *, "Completed metropolis run at beta = ", beta 
                        
                        counter = counter + 1
                        
                end if 
                print *, "*************************************"
                call flush(output_unit)
                
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
        if (any(winding_array /= 0)) then 
                print *, "Skyrmions Found!"
        end if 
 
end program main
