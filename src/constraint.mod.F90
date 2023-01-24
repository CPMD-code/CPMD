module constraint
    USE kinds,                           ONLY: real_8

    implicit none

    !> Type used to store sparse constraint matrix
    !> @author V. Bolnykh - RWTH Aachen/Cyprus Institute
    type constraint_matrix
        ! Values stored in a sparse matrix
        real(real_8), dimension(:), allocatable :: vals
        ! Number of non-zero elements of the current row
        integer, dimension(:), allocatable :: nz_row
        ! Index of a non-zero element in the row
        integer, dimension(:), allocatable :: r_index
        ! IDs of atoms in a constraint
        integer, dimension(:, :), allocatable :: atom_id
        ! Species and atom IDs in the constraint (needed for TAU
        ! mapping)
        integer, dimension(:, :, :), allocatable :: tau_map
        ! Reference value for a constraint
        real(real_8), dimension(:), allocatable :: ref

        contains
        procedure :: mat_vec
        generic :: operator(*) => mat_vec
    end type constraint_matrix

    !> Interface of the constraint_matrix constructor
    interface constraint_matrix
        module procedure new_matrix
    end interface

contains

    !> Matrix-vector multiplication
    function mat_vec(m, v) result(r)

        !> input matrix
        class(constraint_matrix), intent(in) :: m
        !> input vector
        real(real_8), dimension(:), intent(in) :: v

        !> resulting vector
        real(real_8), dimension(:), allocatable :: r

        integer :: i, j

        allocate(r(size(v)))

        r = 0.0_real_8

        !$OMP PARALLEL DO PRIVATE(i, j)
        do i = 1, size(m%nz_row) - 1
            do j = m%nz_row(i - 1) + 1, m%nz_row(i)
                r(i) = r(i) + m%vals(j) * v(m%r_index(j))
            end do
        end do
        !$OMP END PARALLEL DO

    end function

    !> Constructor of the sparse matrix
    function new_matrix(vals, nz_row, r_index, atom_id, tau_map, ref)

        real(real_8), dimension(:), intent(in) :: vals
        integer, dimension(0:), intent(in) :: nz_row
        integer, dimension(:), intent(in) :: r_index
        integer, dimension(:,:), intent(in) :: atom_id
        integer, dimension(:,:,:), intent(in) :: tau_map
        real(real_8), dimension(:), intent(in) :: ref

        type(constraint_matrix) :: new_matrix

        integer :: i, j

        allocate(new_matrix%vals(size(vals)))
        allocate(new_matrix%nz_row(0:size(nz_row) - 1))
        allocate(new_matrix%r_index(size(r_index)))
        allocate(new_matrix%atom_id(size(atom_id, 1), size(atom_id, 2)))
        allocate(new_matrix%tau_map(size(tau_map, 1), size(tau_map, 2), &
                 size(tau_map, 3)))
        allocate(new_matrix%ref(size(ref)))

        !$OMP PARALLEL PRIVATE(i, j)
        !$OMP DO
        do i = 1, size(vals)
            new_matrix%vals(i) = vals(i)
        end do
        !$OMP END DO

        !$OMP DO
        do i = 0, size(nz_row) - 1
            new_matrix%nz_row(i) = nz_row(i)
        end do
        !$OMP END DO

        !$OMP DO
        do i = 1, size(r_index)
            new_matrix%r_index(i) = r_index(i)
        end do
        !$OMP END DO

        !$OMP DO
        do i = 1, size(atom_id, 2)
            new_matrix%atom_id(:, i) = atom_id(:, i)
        end do
        !$OMP END DO

        !$OMP DO
        do i = 1, size(tau_map, 3)
            do j = 1, size(tau_map, 2)
                new_matrix%tau_map(:, j, i) = tau_map(:, j, i)
            end do
        end do
        !$OMP END DO

        !$OMP DO
        do i = 1, size(ref)
            new_matrix%ref(i) = ref(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end function

end module constraint
