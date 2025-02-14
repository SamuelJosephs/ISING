program main 
        use Atom
        use chainMeshCell 
        use ChainMesh
        use EnergyMin
        use RGFlow, only: compute_correlation_function, write_correlations_to_file 
        implicit none
        integer :: numCells, numsteps, i, iostat 
        character(len=30) :: arg
        type(ChainMesh_t) :: testMesh
        type(Atom_t) :: AtomsInUnitCell(2)
        real :: AtomParam1(1), AtomParam2(1), a, energy 
        real(kind=8) :: idble, T, Kb,beta, sigma, upperT, lowerT 
        real(kind=8), allocatable :: correlation(:), positions(:)
        real(kind=8) :: max_distance 
        character(len=20) :: str_buff

        !Kb = 8.62e-5 ! ev/ J 
        Kb = 1.38e-23_dp 
        !Kb = 8.617e-5_dp
        a = 2.8
        numsteps = 120
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

        max_distance = dble(testMesh%domainWidth / 2)
        
        call assignNearestNeighbors(testMesh)
        !print *, "Succsess!"
        !call enumerateChainMeshCells(testMesh) 
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
                call Metropolis2(testMesh,sigma,beta,10000000)
                print *, "Finished metropolis2 step"
                !print *, "Writing Magnetisation to file:", i , "/" , numsteps
                call WriteMagnetization(testMesh,"Magnetisation4.dat")
                call compute_correlation_function(testMesh, max_distance, correlation, positions)
                write(str_buff, "(I0)") i 
                !call write_array_to_file(correlation,trim("./correlation_results/correlation_") // trim(str_buff), iostat)
                call write_correlations_to_file(correlation,positions,trim("./correlation_results/correlation_") // trim(str_buff))
        end do 

        !call deallocateChainMesh(testMesh)
end program main
