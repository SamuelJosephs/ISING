

module ChainMeshCell
        type ChainMeshCell_t 
                integer :: NumAtomsPerUnitCell
                integer :: firstAtomInMeshCell 
                real(kind=8) :: centreX, centreY, centreZ
        end type ChainMeshCell_t 
end module ChainMeshCell 

