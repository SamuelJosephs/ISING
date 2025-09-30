module cube_partition
    use constants
    use iso_fortran_env, only: dp=>real64
    implicit none
    private  ! Make everything private by default
    
    ! Define precision
    
    ! Make these types public so they can be used by calling code
    public :: bounds_type, partition_unit_cube
    
    ! Define bounds type to hold min and max points
    type :: bounds_type
        real(dp), dimension(3) :: min_point
        real(dp), dimension(3) :: max_point
    end type bounds_type
    
    ! Private types used internally
    type :: node_type
        type(bounds_type) :: bounds
        integer :: depth
        integer :: dim
        real(dp) :: split_value
        type(node_type), pointer :: left => null()
        type(node_type), pointer :: right => null()
    end type node_type
    
    type :: leaf_list_type
        type(bounds_type), allocatable :: bounds(:)
        integer :: count = 0
    end type leaf_list_type

contains
    ! Public interface
    function partition_unit_cube(num_regions) result(subdomains)
        integer, intent(in) :: num_regions
        type(bounds_type), allocatable :: subdomains(:)
        type(node_type), pointer :: root
        type(bounds_type) :: initial_bounds
        type(leaf_list_type) :: leaf_list
        logical :: debug
        debug = .false.
        ! Set up initial cube from (0,0,0) to (1,1,1)
        initial_bounds%min_point = [0.0_dp, 0.0_dp, 0.0_dp]
        initial_bounds%max_point = [1.0_dp, 1.0_dp, 1.0_dp]
        
        ! Create root node
        root => create_node(initial_bounds, 0)
        
        ! Partition space
        ! Add debugging output
        ! if (present(debug) .and. debug) then
        !     print *, "Splitting on dimension:", node%dim, " at depth:", node%depth
        !     print *, "Current bounds:"
        !     print *, "  Min:", node%bounds%min_point
        !     print *, "  Max:", node%bounds%max_point
        !     print *, "Split value:", split_value
        !     print *
        ! end if
        
        call partition_space(root, num_regions, leaf_list, debug)
        
        ! Copy results to output array
        allocate(subdomains(leaf_list%count))
        subdomains = leaf_list%bounds(1:leaf_list%count)
        
        ! Cleanup
        call cleanup_tree(root)
        if (allocated(leaf_list%bounds)) deallocate(leaf_list%bounds)
    end function partition_unit_cube
    
    ! Private implementation functions
    function create_node(bounds, depth) result(node)
        type(bounds_type), intent(in) :: bounds
        integer, intent(in) :: depth
        type(node_type), pointer :: node
        
        allocate(node)
        node%bounds = bounds
        node%depth = depth
        ! Convert depth to 0-based for cycling through dimensions, then back to 1-based for array indexing
        ! This ensures we match Python's cycling behavior while maintaining Fortran's array indexing
        node%dim = mod(depth, 3) + 1  ! Now explicitly cycles: 1->x, 2->y, 3->z
        nullify(node%left)
        nullify(node%right)
    end function create_node
    
    function calculate_volume(bounds) result(volume)
        type(bounds_type), intent(in) :: bounds
        real(dp) :: volume
        real(dp), dimension(3) :: sides
        
        sides = bounds%max_point - bounds%min_point
        volume = product(sides)
    end function calculate_volume
    
    subroutine split_bounds(original, dim, value, left, right)
        type(bounds_type), intent(in) :: original
        integer, intent(in) :: dim
        real(dp), intent(in) :: value
        type(bounds_type), intent(out) :: left, right
        
        left = original
        right = original
        left%max_point(dim) = value
        right%min_point(dim) = value
    end subroutine split_bounds
    
    recursive subroutine partition_space(node, target_regions, leaf_list, debug)
        logical, intent(in), optional :: debug
        type(node_type), pointer, intent(inout) :: node
        integer, intent(in) :: target_regions
        type(leaf_list_type), intent(inout) :: leaf_list
        type(bounds_type) :: left_bounds, right_bounds
        real(dp) :: split_value, left_vol, right_vol, total_vol
        integer :: left_regions, right_regions
        
        if (target_regions <= 1) then
            leaf_list%count = leaf_list%count + 1
            if (.not. allocated(leaf_list%bounds)) then
                allocate(leaf_list%bounds(100))
            else if (leaf_list%count > size(leaf_list%bounds)) then
                call grow_leaf_list(leaf_list)
            end if
            leaf_list%bounds(leaf_list%count) = node%bounds
            return
        end if
        
        split_value = (node%bounds%min_point(node%dim) + &
                      node%bounds%max_point(node%dim)) / 2.0_dp
        node%split_value = split_value
        
        call split_bounds(node%bounds, node%dim, split_value, left_bounds, right_bounds)
        node%left => create_node(left_bounds, node%depth + 1)
        node%right => create_node(right_bounds, node%depth + 1)
        
        left_vol = calculate_volume(left_bounds)
        right_vol = calculate_volume(right_bounds)
        total_vol = left_vol + right_vol
        
        left_regions = max(1, nint((left_vol / total_vol) * (target_regions - 1.0_dp)))
        right_regions = max(1, target_regions - left_regions)
        
        call partition_space(node%left, left_regions, leaf_list)
        call partition_space(node%right, right_regions, leaf_list)
    end subroutine partition_space
    
    subroutine grow_leaf_list(leaf_list)
        type(leaf_list_type), intent(inout) :: leaf_list
        type(bounds_type), allocatable :: temp(:)
        integer :: old_size, new_size
        
        old_size = size(leaf_list%bounds)
        new_size = old_size * 2
        
        allocate(temp(new_size))
        temp(1:old_size) = leaf_list%bounds
        call move_alloc(temp, leaf_list%bounds)
    end subroutine grow_leaf_list
    
    subroutine cleanup_tree(node)
        type(node_type), pointer, intent(inout) :: node
        
        if (associated(node)) then
            if (associated(node%left)) call cleanup_tree(node%left)
            if (associated(node%right)) call cleanup_tree(node%right)
            deallocate(node)
        end if
    end subroutine cleanup_tree
    
end module cube_partition
