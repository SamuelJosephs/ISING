
module Atom 
        use omp_lib
        type Atom_t 
                real :: x,y,z, tmpx, tmpy, tmpz
                integer :: NumAtomParameters
                real , allocatable :: AtomParameters(:) ! e.g. spin
                integer :: nextAtom
                integer, allocatable :: NeighborList(:)

        end type
contains
        function makeAtom(x,y,z,atomParameters,NumAtomParameters,  NextAtom) result(res)
                real, intent(in) :: x,y,z
                integer, intent(in) :: NumAtomParameters, NextAtom
                real, intent(inout), target :: AtomParameters(NumAtomParameters)
                type(Atom_t) :: res
               res%x = x
               res%y = y
               res%z = z
               res%NumAtomParameters = NumAtomParameters
               allocate(res%atomParameters(size(atomParameters)))
               res%atomParameters = atomParameters

               allocate(res%NeighborList(0)) 
               return
        end function makeAtom

        function AtomDeepCopy(atom) result(res)
        type(Atom_t), intent(inout), pointer :: atom 
        type(Atom_t) :: res 
        res%nextAtom = atom%nextAtom 
        res%x = atom%x 
        res%y = atom%y
        res%z = atom%z
        res%numAtomParameters = atom%numAtomParameters 
        allocate(res%atomParameters(size(atom%atomPArameters)))
        res%atomParameters = atom%atomParameters 
        allocate(res%NeighborList((size(atom%NeighborList))))
        res%neighborList = atom%neighborList
        end function atomDeepCopy 
end module Atom 
