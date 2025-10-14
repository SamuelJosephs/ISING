module unit_cells
        implicit none 

        public :: hyper_kagome_unit_cell
        contains 

        subroutine bcc_unit_cell(a,b,c,alpha,beta,gamma,atomsInUnitcell)

                use iso_fortran_env, only : dp=> real64
                use atom 
                implicit none 
                real(kind=dp), intent(inout) :: a,b,c,alpha,beta,gamma
                type(atom_t), allocatable, dimension(:), intent(inout) :: atomsInUnitCell 

                integer :: stat

                if (allocated(atomsInUnitCell)) deallocate(atomsInUnitCell)
                allocate(atomsInUnitCell(2),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate atoms in unit cell array"

                a = 2.0_dp
                b = 2.0_dp
                c = 2.0_dp
                alpha = 90.0_dp
                beta = 90.0_dp
                gamma = 90.0_dp
                atomsInUnitCell(1) = makeAtom(0.0,0.0,0.0,-1)
                atomsInUnitCell(2) = makeAtom(0.5,0.5,0.5,-1)

        end subroutine bcc_unit_cell

        subroutine hyper_kagome_unit_cell(a,b,c,alpha,beta,gamma,atomsInUnitcell)
                ! Taken from DOI: https://doi.org/10.1103/PhysRevB.110.054428

                use iso_fortran_env, only : dp=> real64
                use atom 
                implicit none 
                real(kind=dp), intent(out) :: a,b,c,alpha,beta,gamma
                type(atom_t), allocatable, dimension(:), intent(inout) :: atomsInUnitCell 

                integer :: stat

                if (allocated(atomsInUnitCell)) deallocate(atomsInUnitCell)
                allocate(atomsInUnitCell(12),stat=stat)
                if (stat /= 0) error stop "Error: Failed to allocate atoms in unit cell array"

                a = 6.0_dp
                b = 6.0_dp
                c = 6.0_dp
                alpha = 90.0_dp
                beta = 90.0_dp
                gamma = 90.0_dp
                ! r1 = 1/4 (-2,0,2)
                atomsInUnitCell(1) = makeAtom(-0.5,0.0,0.5,-1)
                ! r2 = 1/4 (-1,3,2)
                atomsInUnitCell(2) = makeAtom(-0.25,3.0/4.0,0.5,-1)
                ! r3 = 1/4 (-2,3,1)
                atomsInUnitCell(3) = makeAtom(-0.5,3.0/4.0,0.25,-1)
                ! r4 = 1/4(-1,1,0)
                atomsInUnitCell(4) = makeAtom(-0.25,0.25,0.0,-1)
                ! r5 = 1/4 (-2,1,3)
                atomsInUnitCell(5) = makeAtom(-0.5,0.25,3.0/4.0,-1)
                ! r6 = 1/4 (-1,2,3)
                atomsInUnitCell(6) = makeAtom(-0.25,0.5,3.0/4.0,-1)
                ! r7 = 1/4 (-3,2,1)
                atomsInUnitCell(7) = makeAtom(-3.0/4.0,0.5,0.25,-1)
                ! r8 = 1/4 (0,2,2)
                atomsInUnitCell(8) = makeAtom(0.0,0.5,0.5,-1)
                ! r9 = 1/4 (0,1,1)
                atomsInUnitCell(9) = makeAtom(0.0,0.25,0.25,-1)
                ! r10 1/4 (-3,3,0)
                atomsInUnitCell(10) = makeAtom(-3.0/4.0,3.0/4.0,0.0,-1)
                ! r11 = (0,0,0)
                atomsInUnitCell(11) = makeAtom(0.0,0.0,0.0,-1)
                ! r12 = 1/4 (-3,0,3)
                atomsInUnitcell(12) = makeAtom(-3.0/4.0,0.0,3.0/4.0,-1)
        end subroutine hyper_kagome_unit_cell



end module unit_cells
