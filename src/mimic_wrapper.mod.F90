!---------------------------------------------------------------------------------------------------
! MODULE: mimic_wrapper
!> @author Viacheslav Bolnykh - RWTH Aachen/Cyprus Institute
!> @author Jogvan Magnus Haugaard Olsen - DTU (jmho@kemi.dtu.dk)
!
! DESCRIPTION:
!> Wrapper module for MiMiC main routines
!---------------------------------------------------------------------------------------------------

module mimic_wrapper

#ifdef __MIMIC
    use mimic_constants, only: mimic_bohr, mimic_covrad
    use mimic_main
    use mimic_tensors, only: initialize_tensors, &
                             terminate_tensors
#endif
    use cell, only: cell_com
    use control_pri_utils, only: control_pri
    use cnst, only: factem
    use coor, only: tau0, taup, fion, velp, lvelini
    use cotr, only: cnsval, grate, cnpar, cnsval_dest, lskcor, &
                    cotc0, ntcnst, const_nconn, const_conn
    use efld, only : textfld, extf
    use error_handling, only: stopgm
    use gvec, only: gvec_com
    use inscan_utils, only: inscan
    use ions, only: ions0, ions1
    use isos, only: isos1, isos3
    use kinds, only: real_8
    use mm_dimmod, only: clsaabox
    use mm_extrap, only: cold
    use mp_interface, only: mp_bcast, mp_sum, mp_sync
    use parac, only: parai, paral
    use readsr_utils, only: readsr, readsi, input_string_len, keyword_contains
    use rmas, only: rmass
    use system, only: maxsys, fpar, parap, cntl, spar, parm, cnti, nkpt, iatpt
    use timer, only: tihalt, tiset

    implicit none

    !> MiMiC control type, maps on the input file section
    type :: mimic_control_type
        !> vectors of the whole simulation box
        real(real_8), dimension(3,3) :: box
        !> list of covalent radii (read from external file, per species)
        real(real_8), dimension(:), allocatable :: covalent_radius
        !> list of factors to be applied to forces contributions of
        !> different code (one per code)
        real(real_8), dimension(:), allocatable :: factors
        !> paths to working folders of each of the client code
        character(len=250), dimension(:), allocatable :: paths
        !> total number of species in the simulation
        integer :: num_species
        !> maximum nuber of atoms per species
        integer :: max_atoms
        !> number of client codes
        integer :: num_client
        !> number of quantum species
        integer :: num_quantum_species
        !> maximum number of atoms per quantum species
        integer :: max_quantum_atoms
        !> total number of atoms
        integer :: num_atoms
        !> total number of quantum atoms
        integer :: num_quantum_atoms
        !> switch to turn on the computation of potential
        logical :: update_potential = .false.
        !> indicator if the dimension is switched to MM system
        logical :: dim_mm = .false.
        !> number of saved dimension status
        integer :: num_save = 0
        !> save dimension status
        logical, dimension(20) :: dim_save = .false.
        !> use total energy in SCF
        logical :: tot_scf_energy = .false.
        !> use long-range coupling scheme
        logical :: long_range_coupling = .false.
        !> distance cutoff (in bohr) defining short- and long-range regions
        real(real_8) :: cutoff_distance = huge(real_8)
        !> method used to sort fragments
        integer :: sorting_method = huge(0)
        !> how often to update the sorted fragments (default every step)
        integer :: update_sorting = huge(0)
        !> multipole expansion order
        integer :: multipole_order = huge(0)
        !> include constraints from external programs or not
        logical :: external_constraints = .true.
        !> degrees of freedom for QM subsystem
        real(real_8) :: qm_dof
        !> degrees of freedom for MM subsystem
        real(real_8) :: mm_dof
    end type mimic_control_type

    !> energy holder for QM/MM simulation
    type :: mimic_energy_type
        !> QM/MM energy
        real(real_8) :: qmmm_energy = 0.0_real_8
        !> QM energy
        real(real_8) :: qm_energy = 0.0_real_8
        !> MM energy
        real(real_8) :: mm_energy = 0.0_real_8
    end type mimic_energy_type

    ! IDs of the first and the last SR atom of the CP group
    integer, dimension(:), allocatable, private :: gr_sr_atom_start, gr_sr_atom_end
    ! IDs of the first and the last LR atom of the CP group
    integer, dimension(:), allocatable, private :: gr_lr_atom_start, gr_lr_atom_end
    ! IDs of the first and the last SR atom of the current process
    integer, dimension(:), allocatable, private :: pr_sr_atom_start, pr_sr_atom_end
    ! IDs of the first and the last LR atom of the current process
    integer, dimension(:), allocatable, private :: pr_lr_atom_start, pr_lr_atom_end

    !> internal forces array to store QM/MM and MM forces
    real(real_8), dimension(:,:,:), allocatable :: mimic_forces
    !> internal forces array to store MM forces
    real(real_8), dimension(:,:,:), allocatable :: mm_forces
    !> coordinates array storing all of the coordinates with MM box PBC applied
    real(real_8), dimension(:,:,:), allocatable :: wrapped_coordinates

    !> MiMiC options
    type(mimic_control_type), save :: mimic_control
    !> storage for the energy
    type(mimic_energy_type), save :: mimic_energy
#ifdef __MIMIC
    !> temporary storage for the allocation sizes
    type(sizes_type), save :: sizes
    !> system information
    type(system_type), save :: system_data
    !> array of subsystems associated with each client code
    type(subsystem_type), dimension(:), allocatable :: subsystems
    !> quantum subsystem
    type(quantum_fragment_type) :: quantum_fragment
#else
    character(len=*), parameter :: MIMIC_MISSING_ERROR = "CPMD is compiled without MiMiC support but MiMiC procedures are called!"
#endif

    procedure(cpmd_error_handler), pointer, save :: error_handler => cpmd_error_handler

contains

!> subroutine to handle MiMiC errors
subroutine cpmd_error_handler(err_type, message, source_file, line_num)

    !> Type of the error
    integer, intent(in) :: err_type
    !> Optional message
    character(len=*), intent(in) :: message
    !> Source file where the problem occurred
    character(len=*), intent(in) :: source_file
    !> Source line at which the error occurred
    integer, intent(in) :: line_num

    call stopgm("cpmd_error_handler", &
                "MiMiC has encountered an error, see MiMiC_error file(s) for details", &
                line_num, source_file)

end subroutine cpmd_error_handler

!> initialize overlaps within MiMiC
subroutine mimic_ifc_init_overlaps(overlaps)

    !> overlap maps
    integer, dimension(:,:), intent(in) :: overlaps

    character(*), parameter :: procedureN = 'mimic_ifc_init_overlaps'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    call mimic_init_overlaps(overlaps)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_init_overlaps

!> finalize MiMiC simulation and destroy internal objects
subroutine mimic_ifc_destroy

    character(*), parameter :: procedureN = 'mimic_ifc_destroy'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    call terminate_tensors()
    if (paral%io_parent) then
        call mimic_finalize
        call mimic_destroy
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_destroy

#ifdef __MIMIC
!> request system information: coordinates, multipoles, etc.
subroutine mimic_ifc_request_system_data(sizes, system_data)

    !> sizes of data structures
    type(sizes_type), intent(inout) :: sizes
    !> system information
    type(system_type), intent(inout) :: system_data

    character(*), parameter :: procedureN = 'mimic_ifc_request_system_data'
    integer :: isub

    call tiset(procedureN, isub)

    call mimic_request_system_data(sizes, system_data)

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_request_system_data
#endif

!> collect forces from clients
subroutine mimic_ifc_collect_forces()

    character(*), parameter :: procedureN = 'mimic_ifc_collect_forces'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    call mimic_collect_forces(subsystems)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_collect_forces

!> collect energies from clients
subroutine mimic_ifc_collect_energies()

    character(*), parameter :: procedureN = 'mimic_ifc_collect_energies'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    call mimic_collect_energies(subsystems, mimic_energy%mm_energy)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_collect_energies

!> compute QM/MM interaction energy
subroutine mimic_ifc_compute_energy

    real(real_8) :: temp_energy
    real(real_8), dimension(:,:), allocatable :: tensors, tensor_sums

    character(*), parameter :: procedureN = 'mimic_ifc_compute_energy'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (mimic_control%long_range_coupling) then
        call mimic_compute_multipoles(quantum_fragment, mimic_control%multipole_order, &
                                    parap%NRXPL(parai%mepos,1), parap%NRXPL(parai%mepos,2))
        call mp_sum(quantum_fragment%electronic_multipoles%values, &
                    size(quantum_fragment%electronic_multipoles%values), &
                    parai%allgrp)

        call mimic_compute_tensors(subsystems, &
                                 quantum_fragment, &
                                 tensors, &
                                 tensor_sums, &
                                 pr_lr_atom_start, pr_lr_atom_end)
    end if

    mimic_energy%qmmm_energy = 0.0_real_8

    call mimic_compute_energy(subsystems, &
                            quantum_fragment, &
                            mimic_energy%qmmm_energy, &
                            pr_sr_atom_start, &
                            pr_sr_atom_end, &
                            tensor_sums)
    CALL mp_sum(mimic_energy%qmmm_energy, parai%cp_grp)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_compute_energy

!> compute QM/MM forces
subroutine mimic_ifc_compute_forces()

    character(*), parameter :: procedureN = 'mimic_ifc_compute_forces'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (mimic_control%long_range_coupling) then
        call mimic_compute_multipoles(quantum_fragment, mimic_control%multipole_order, &
                                    parap%NRXPL(parai%mepos,1), parap%NRXPL(parai%mepos,2))
        call mp_sum(quantum_fragment%electronic_multipoles%values, &
                    size(quantum_fragment%electronic_multipoles%values), &
                    parai%allgrp)
    end if

    call mimic_compute_forces(subsystems, &
                            quantum_fragment, &
                            parap%NRXPL(parai%mepos,1), parap%NRXPL(parai%mepos,2), &
                            gr_sr_atom_start, gr_sr_atom_end, &
                            gr_lr_atom_start, gr_lr_atom_end, &
                            pr_sr_atom_start, pr_sr_atom_end, &
                            pr_lr_atom_start, pr_lr_atom_end, paral%parent)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_compute_forces

!> compute external potential due to MM atoms
subroutine mimic_ifc_compute_potential()

    real(real_8), dimension(:,:), allocatable :: tensors, tensor_sums

    character(*), parameter :: procedureN = 'mimic_ifc_compute_potential'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (mimic_control%long_range_coupling) then
        call mimic_compute_tensors(subsystems, &
                                 quantum_fragment, &
                                 tensors, &
                                 tensor_sums, &
                                 pr_lr_atom_start, pr_lr_atom_end)
        call mp_sum(tensor_sums, size(tensor_sums), parai%allgrp)
        tensor_sums = -tensor_sums
    end if

    call mimic_compute_potential(subsystems, &
                               quantum_fragment, &
                               parap%NRXPL(parai%mepos,1), parap%NRXPL(parai%mepos,2), &
                               gr_sr_atom_start, gr_sr_atom_end, tensor_sums)

    CALL mp_sum(extf, size(extf), parai%cp_inter_grp)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_compute_potential

!> sort fragments into short- and long-ranged shells
subroutine mimic_ifc_sort_fragments()

    integer :: i
    integer :: tot_sr_atoms
    integer :: tot_lr_atoms
    integer :: num_sr_atoms
    integer :: num_lr_atoms
    integer :: atom_start
    integer :: atom_end
    integer :: remainder

    character(*), parameter :: procedureN = 'mimic_ifc_sort_fragments'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    tot_sr_atoms = 0
    tot_lr_atoms = 0

    if (paral%io_parent) then
        write(6,'(8x,a)') 'SORTING MM ATOMS INTO SHORT- AND LONG-RANGE DOMAINS'
    end if


    do i = 1, size(subsystems)
        call mimic_sort_fragments(subsystems(i), quantum_fragment, &
                                mimic_control%cutoff_distance, mimic_control%sorting_method)

        call mimic_distribute_atoms(subsystems(i)%num_sr_atoms, parai%cp_nogrp, &
                                parai%cp_inter_me, gr_sr_atom_start(i), gr_sr_atom_end(i))
        call mimic_distribute_atoms(subsystems(i)%num_lr_atoms, parai%cp_nogrp, &
                                parai%cp_inter_me, gr_lr_atom_start(i), gr_lr_atom_end(i))
        call mimic_distribute_atoms(subsystems(i)%num_sr_atoms, parai%cp_nproc, &
                                parai%cp_me, pr_sr_atom_start(i), pr_sr_atom_end(i))
        call mimic_distribute_atoms(subsystems(i)%num_lr_atoms, parai%cp_nproc, &
                                parai%cp_me, pr_lr_atom_start(i), pr_lr_atom_end(i))
        if (paral%io_parent) then
            ! write(6,'(8x,a,i2,a)') "SUBSYSTEM ", i, ":"
            write(6,'(8x,a,t58,i8)') "SHORT-RANGE ATOMS PER PROCESS GROUP:", &
                                  gr_sr_atom_end(i) - gr_sr_atom_start(i) + 1
            write(6,'(8x,a,t58,i8)') "SHORT-RANGE ATOMS PER PROCESS:", &
                                  pr_sr_atom_end(i) - pr_sr_atom_start(i) + 1
            write(6,'(8x,a,t58,i8)') "LONG-RANGE ATOMS PER PROCESS GROUP:", &
                                  gr_lr_atom_end(i) - gr_lr_atom_start(i) + 1
            write(6,'(8x,a,t58,i8)') "LONG-RANGE ATOMS PER PROCESS:", &
                                  pr_lr_atom_end(i) - pr_lr_atom_start(i) + 1
        end if
    end do
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_sort_fragments

!> Distribute atoms across the processes of the given communicator
subroutine mimic_distribute_atoms(tot_num_atoms, comm_size, proc_id, atom_start, atom_end)

    !> total number of atoms to distribute
    integer, intent(in) :: tot_num_atoms
    !> size of the communicator, i.e., number of processes in the communicator
    integer, intent(in) :: comm_size
    !> ID of the current process (starts counting from zero)
    integer, intent(in) :: proc_id
    !> index marking the start of the chunk of the atom array to be treated by this process
    integer, intent(out) :: atom_start
    !> index marking the end of the chunk of the atom array to be treated by this process
    integer, intent(out) :: atom_end

    integer :: remainder
    integer :: n_atoms
    integer :: temp_final

    character(*), parameter :: procedureN = 'mimic_distribute_atoms'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (proc_id > comm_size) then
        CALL stopgm(procedureN, 'process id does not belong to the communicator', __LINE__, __FILE__)
    endif

    n_atoms = tot_num_atoms / comm_size
    remainder = modulo(tot_num_atoms, comm_size)
    atom_start = 0
    if (remainder > 0) then
        if (proc_id < remainder) then
            n_atoms = n_atoms + 1
        else
            atom_start = remainder
        endif
    end if
    atom_start = atom_start + proc_id * n_atoms + 1
    temp_final = atom_start + n_atoms - 1
    if (temp_final > tot_num_atoms) then
        temp_final = tot_num_atoms
    endif
    atom_end = temp_final
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_distribute_atoms

!> add internal forces array to global forces
subroutine mimic_sum_forces(fion)

    !> global forces
    real(real_8), dimension(:,:,:), intent(inout) :: fion

    real(real_8), dimension(:,:,:), allocatable :: fion_temp

    character(*), parameter :: procedureN = 'mimic_sum_forces'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    allocate(fion_temp, mold=mimic_forces)

    CALL mp_sum(mimic_forces, fion_temp, size(mimic_forces), parai%cp_grp)
    fion = fion + fion_temp

    mimic_forces = 0.0_real_8
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_sum_forces

!> save dimensions of atomic data structures
subroutine mimic_save_dim()

    character(*), parameter :: procedureN = 'mimic_save_dim'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    mimic_control%num_save = mimic_control%num_save + 1
    mimic_control%dim_save(mimic_control%num_save) = mimic_control%dim_mm
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_save_dim

!> revert any changes to dimensions of data structures
subroutine mimic_revert_dim()

    character(*), parameter :: procedureN = 'mimic_revert_dim'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    call mimic_switch_dim(go_qm= .not. mimic_control%dim_save(mimic_control%num_save))
    mimic_control%num_save = mimic_control%num_save - 1
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_revert_dim

!> switch dimensions of atomic data structures
subroutine mimic_switch_dim(go_qm)

    !> flag indicating switch MM->QM or vice versa
    logical, intent(in) :: go_qm

    character(*), parameter :: procedureN = 'mimic_switch_dim'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (.not. cntl%mimic) return

    if (go_qm) then
        maxsys%nax = mimic_control%max_quantum_atoms
        maxsys%nsx = mimic_control%num_quantum_species
        ions1%nsp = mimic_control%num_quantum_species
        ions1%nat = mimic_control%num_quantum_atoms
        mimic_control%dim_mm = .false.
    else
        maxsys%nax = mimic_control%max_atoms
        maxsys%nsx = mimic_control%num_species
        ions1%nsp = mimic_control%num_species
        ions1%nat = mimic_control%num_atoms
        mimic_control%dim_mm = .true.
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_switch_dim

!> read &MIMIC section of the input file
subroutine mimic_read_input()

    integer, parameter :: iunit = 5
    integer, parameter :: output_unit = 6
    integer, parameter :: max_unknown_lines = 250
    character(len=input_string_len) :: line, previous_line, error_message, &
                                       unknown(max_unknown_lines)
    integer :: first, last, ierr, num_unknown_lines
    logical :: error, something_went_wrong, go_on_reading
    integer :: i, j
    integer :: num_overlaps, num_factors
    integer, dimension(:, :), allocatable :: overlaps

    character(len=*), parameter :: procedureN = 'mimic_read_input'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (paral%io_parent) then
        ! write(output_unit,'(1x,a)') 'READING &MIMIC SECTION'
        num_unknown_lines = 0
        line = ' '
        previous_line = ' '
        error_message = ' '
        ierr = inscan(iunit, '&MIMIC')
        if (ierr == 0) then
            textfld = .true.
            mimic_control%tot_scf_energy = .false.
            mimic_control%long_range_coupling = .false.
            mimic_control%cutoff_distance = huge(real_8)
            mimic_control%sorting_method = CENTROID
            mimic_control%update_sorting = 1
            mimic_control%multipole_order = 2
            go_on_reading        = .true.
            something_went_wrong = .false.
            do while (go_on_reading)
                previous_line = line
                read(iunit, '(a80)', iostat=ierr) line
                if (ierr /= 0) then
                    something_went_wrong = .true.
                    go_on_reading = .false.
                else if (keyword_contains(line, '&END')) then
                    go_on_reading = .false.
                else if (keyword_contains(line, 'PATHS')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = 1
                    call readsi(line, first, last, mimic_control%num_client, error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    allocate(mimic_control%paths(mimic_control%num_client), stat=ierr)
                    if (ierr /= 0) call stopgm(procedureN, 'Error while allocating array, cf. output file', __LINE__, __FILE__)
                    allocate(mimic_control%factors(mimic_control%num_client), stat=ierr)
                    if (ierr /= 0) call stopgm(procedureN, 'Error while allocating array, cf. output file', __LINE__, __FILE__)
                    mimic_control%factors(:) = 1
                    do i = 1, mimic_control%num_client
                        read(iunit, fmt='(a)', iostat=ierr) mimic_control%paths(i)
                        if (ierr /= 0) then
                            something_went_wrong = .true.
                            go_on_reading = .false.
                        end if
                    end do
                else if (keyword_contains(line, 'DISABLE', and='CONSTRAINTS')) then
                    mimic_control%external_constraints = .false.
                else if (keyword_contains(line, 'OVERLAPS')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = 1
                    call readsi(line, first, last, num_overlaps, error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    allocate(overlaps(4, num_overlaps), stat=ierr)
                    if (ierr /= 0) call stopgm(procedureN, 'Error while allocating array, cf. output file', __LINE__, __FILE__)
                    do i = 1, num_overlaps
                        read(iunit, '(a80)', iostat=ierr) line
                        if (ierr /= 0) then
                            something_went_wrong = .true.
                            go_on_reading = .false.
                        end if
                        first = 1
                        do j = 1, 4
                            call readsi(line, first, last, overlaps(j, i), error)
                            if (error) then
                                error_message = 'ERROR WHILE READING VALUE'
                                something_went_wrong = .true.
                                go_on_reading = .false.
                            end if
                            first = last + 1
                        end do
                    end do
                else if (keyword_contains(line, 'BOX')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = 1
                    call readsr(line, first, last, mimic_control%box(1,1), error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = last + 1
                    call readsr(line, first, last, mimic_control%box(2,2), error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = last + 1
                    call readsr(line, first, last, mimic_control%box(3,3), error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                else if (keyword_contains(line, 'SUBSYSTEM', and='FACTORS')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    call readsi(line, first, last, num_factors, error)
                    if (num_factors /= mimic_control%num_client) then
                        error_message = 'NUMBER OF FACTORS SHOULD BE EQUAL TO THE NUMBER OF CLIENTS'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    endif
                    do i = 1, num_factors
                        read(iunit, '(a80)', iostat=ierr) line
                        if (ierr /= 0) then
                            something_went_wrong = .true.
                            go_on_reading = .false.
                        end if
                        first = 1
                        call readsr(line, first, last, mimic_control%factors(i), error)
                    end do
                else if (keyword_contains(line, 'TOTAL', and='SCF')) then
                    mimic_control%tot_scf_energy = .true.
                else if (keyword_contains(line, 'FRAGMENT', and='SORTING')) then
                    if (keyword_contains(line, 'CENTROID')) then
                        mimic_control%sorting_method = CENTROID
                    else if (keyword_contains(line, 'CENTER-OF-MASS')) then
                        mimic_control%sorting_method = CENTER_OF_MASS
                    else if (keyword_contains(line, 'CENTER-OF-CHARGE')) then
                        mimic_control%sorting_method = CENTER_OF_CHARGE
                    else if (keyword_contains(line, 'MINIMUM', and='DISTANCE')) then
                        mimic_control%sorting_method = MINIMUM_DISTANCE
                    else if (keyword_contains(line, 'ATOM-WISE')) then
                        mimic_control%sorting_method = ATOM_WISE
                    else
                        mimic_control%sorting_method = CENTROID
                    end if
                    if (keyword_contains(line, 'UPDATE')) then
                        read(iunit, '(a)', iostat=ierr) line
                        if (ierr /= 0) then
                            something_went_wrong = .true.
                            go_on_reading = .false.
                        end if
                        first = 1
                        call readsi(line, first, last, mimic_control%update_sorting, error)
                        if (error) then
                            error_message = 'ERROR WHILE READING VALUE'
                            something_went_wrong = .true.
                            go_on_reading = .false.
                        end if
                    end if
                else if (keyword_contains(line, 'LONG-RANGE', and='COUPLING')) then
                    mimic_control%long_range_coupling = .true.
                else if (keyword_contains(line, 'CUTOFF', and='DISTANCE')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = 1
                    call readsr(line, first, last, mimic_control%cutoff_distance, error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                else if (keyword_contains(line, 'MULTIPOLE', and='ORDER')) then
                    read(iunit, '(a)', iostat=ierr) line
                    if (ierr /= 0) then
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                    first = 1
                    call readsi(line, first, last, mimic_control%multipole_order, error)
                    if (error) then
                        error_message = 'ERROR WHILE READING VALUE'
                        something_went_wrong = .true.
                        go_on_reading = .false.
                    end if
                else
                    ! Unknown keyword
                    if (' ' /= line) then
                       if (num_unknown_lines < max_unknown_lines) then
                          num_unknown_lines = num_unknown_lines + 1
                          unknown(num_unknown_lines) = line
                       else
                          do i = 1, max_unknown_lines - 1
                             unknown(i) = unknown(i+1)
                          enddo
                          unknown(num_unknown_lines) = line
                       endif
                    endif
                end if
            end do
        else
            something_went_wrong = .true.
            error_message = 'MISSING &MIMIC SECTION - SECTION MANDATORY FOR MIMIC'
        end if
        if (something_went_wrong) THEN
            write(output_unit,'(/,1x,64("!"))')
            write(output_unit,'(1x, a, 1x, a)') 'ERROR:', 'PROBLEM WHILE READING &MIMIC SECTION:'
            write(output_unit,'(8x,a)') trim(adjustl(error_message))
            if (line /= ' ' .or. previous_line /= ' ') THEN
               write(output_unit,'(8x,a)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
               write(output_unit,'(/,1x,a)') trim(adjustl(previous_line))
               write(output_unit,'(1x,a)') trim(Adjustl(line))
            endif
            write(output_unit,'(1x,64("!"))')
            call stopgm(procedureN, 'Error while reading &MIMIC section, cf. output file', __LINE__, __FILE__)
        endif
        if (num_unknown_lines > 0) then
            write(6,'(/,1x,64("="))')
            write(6,'(1x,a,14x,a,12x,a)') '= ', 'UNKNOWN KEYWORDS IN &MIMIC SECTION', ' ='
            do i = 1, num_unknown_lines
                write(6,'(1x,a2,a60,a2)') '= ', unknown(i), ' ='
            end do
            write(6,'(1x,64("="),/)')
        end if
        ! write(output_unit,'(1x,a)') 'DONE READING &MIMIC SECTION'
    end if

    if (.not. mimic_control%long_range_coupling) then
        mimic_control%cutoff_distance = huge(real_8)
    end if

    if (cntl%mimic) then
        call mimic_init_error_handler(parai%me, error_handler)
        call mp_bcast(textfld, parai%io_source, parai%cp_grp)
        call mp_bcast(rmass%pma0, size(rmass%pma0), parai%io_source, parai%cp_grp)
        call mp_bcast(ions0%iatyp, size(ions0%iatyp), parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%num_client, parai%io_source, parai%cp_grp)
        if (.not. allocated(mimic_control%factors)) allocate(mimic_control%factors(mimic_control%num_client))
        call mp_bcast(mimic_control%factors, size(mimic_control%factors), parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%box, size(mimic_control%box), parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%long_range_coupling, parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%sorting_method, parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%update_sorting, parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%cutoff_distance, parai%io_source, parai%cp_grp)
        call mp_bcast(mimic_control%multipole_order, parai%io_source, parai%cp_grp)

        mimic_control%max_quantum_atoms = maxsys%nax
        mimic_control%num_quantum_species = maxsys%nsx

        call mp_bcast(num_overlaps, parai%io_source, parai%cp_grp)
        if (.not. allocated(overlaps)) allocate(overlaps(4, num_overlaps))
        call mp_bcast(overlaps, size(overlaps), parai%io_source, parai%cp_grp)

        call mimic_set_num_clients(mimic_control%num_client)
        allocate(gr_sr_atom_start(mimic_control%num_client))
        allocate(gr_sr_atom_end(mimic_control%num_client))
        allocate(gr_lr_atom_start(mimic_control%num_client))
        allocate(gr_lr_atom_end(mimic_control%num_client))
        allocate(pr_sr_atom_start(mimic_control%num_client))
        allocate(pr_sr_atom_end(mimic_control%num_client))
        allocate(pr_lr_atom_start(mimic_control%num_client))
        allocate(pr_lr_atom_end(mimic_control%num_client))
        allocate(subsystems(mimic_control%num_client))

        call mimic_ifc_init_overlaps(OVERLAPS)
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_read_input

!> perform handshake operation with clients
subroutine mimic_ifc_handshake

    character(len=*), parameter :: procedureN = 'mimic_ifc_handshake'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    if (paral%io_parent) then
        call mimic_handshake(mimic_control%paths)
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_handshake

!> request allocation sizes for internal data structures
subroutine mimic_ifc_request_sizes

    integer :: max_atom_pfrag
    integer, dimension(:), allocatable :: atoms_pcode
    integer, dimension(2) :: sp_map_shape
    integer :: i, j
    integer :: cp_inter_grp
    integer :: ierr

    character(len=*), parameter :: procedureN = 'mimic_ifc_request_sizes'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    allocate(atoms_pcode(mimic_control%num_client), stat=ierr)
    if (ierr /= 0) call stopgm(procedureN, 'Error while allocating array, cf. output file', __LINE__, __FILE__)

    if (paral%io_parent) then
        call mimic_request_sizes(sizes, system_data)
        max_atom_pfrag = size(sizes%atoms_pfragment, dim=1)
    else
        allocate(sizes%atoms_pcode(mimic_control%num_client))
        allocate(sizes%multipoles_order(mimic_control%num_client))
        allocate(sizes%multipoles_patom(mimic_control%num_client))
        allocate(sizes%frag_num(mimic_control%num_client))
        allocate(sizes%nbonds_pcode(mimic_control%num_client))
        allocate(sizes%nangles_pcode(mimic_control%num_client))
        allocate(sizes%types_length_pcode(mimic_control%num_client))

        sizes%nbonds_pcode(:) = 0
        sizes%nangles_pcode(:) = 0
    end if

    call mp_sync(parai%cp_grp)

    call mp_bcast(sizes%atoms_pcode, &
                  mimic_control%num_client, &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(sizes%multipoles_order, &
                  mimic_control%num_client, &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(sizes%multipoles_patom, &
                  mimic_control%num_client, &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(sizes%frag_num, &
                  mimic_control%num_client, &
                  parai%io_source, &
                  parai%cp_grp)

    call mp_bcast(sizes%types_length_pcode, &
                  mimic_control%num_client, &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(max_atom_pfrag, &
                  parai%io_source, &
                  parai%cp_grp)

    if (.not. allocated(sizes%atoms_pfragment)) &
            allocate(sizes%atoms_pfragment(max_atom_pfrag, mimic_control%num_client))
    call mp_bcast(sizes%atoms_pfragment, size(sizes%atoms_pfragment), parai%io_source, parai%cp_grp)

    call mp_bcast(sizes%num_atoms, parai%io_source, parai%cp_grp)
    call mp_bcast(sizes%num_species, parai%io_source, parai%cp_grp)
    if (.not. allocated(sizes%atoms_pspecies)) allocate(sizes%atoms_pspecies(sizes%num_species))
    call mp_bcast(sizes%atoms_pspecies, sizes%num_species, parai%io_source, parai%cp_grp)
    call mp_bcast(ions0%na, size(ions0%na), parai%io_source, parai%cp_grp)
    call mp_bcast(ions1%nat, parai%io_source, parai%cp_grp)
    call mp_bcast(ions1%nsp, parai%io_source, parai%cp_grp)
    ions1%nsp = maxsys%nsx

    mimic_control%num_species = maxsys%nsx + sizes%num_species
    mimic_control%max_atoms = max(maxsys%nax, sizes%num_atoms)
    mimic_control%num_quantum_atoms = sum(ions0%na)
    mimic_control%num_atoms = mimic_control%num_quantum_atoms + sum(sizes%atoms_pspecies)

    if (allocated(tau0)) deallocate(tau0, stat=ierr)
    if (allocated(velp)) deallocate(velp, stat=ierr)
    if (allocated(lvelini)) deallocate(lvelini, stat=ierr)
    if (allocated(taup)) deallocate(taup, stat=ierr)
    if (allocated(fion)) deallocate(fion, stat=ierr)

    allocate(tau0(3, mimic_control%max_atoms, mimic_control%num_species), stat=ierr)
    allocate(velp(3, mimic_control%max_atoms, mimic_control%num_species), stat=ierr)
    allocate(lvelini(0:mimic_control%max_atoms + 1, mimic_control%num_species), stat=ierr)
    allocate(taup(3, mimic_control%max_atoms, mimic_control%num_species), stat=ierr)
    allocate(fion(3, mimic_control%max_atoms, mimic_control%num_species), stat=ierr)
    tau0 = 0.0_real_8
    velp = 0.0_real_8
    lvelini = .false.
    fion = 0.0_real_8
    taup = 0.0_real_8

    if ( paral%io_parent ) then
        call mimic_ifc_request_system_data(sizes, system_data)
        sp_map_shape = shape(system_data%species_map)
    else
        allocate(system_data%species(maxval(sizes%atoms_pcode), &
                                     mimic_control%num_client))
        allocate(system_data%masses(sizes%num_species))
        allocate(system_data%elements(sizes%num_species))
        allocate(system_data%multipole_values(maxval(sizes%multipoles_patom), &
                    maxval(sizes%atoms_pcode), &
                    mimic_control%num_client))
        allocate(system_data%atom_fragment_ids(maxval(sizes%atoms_pfragment), &
                         maxval(sizes%frag_num), &
                         mimic_control%num_client))
    end if
    call mp_bcast(system_data%species, &
                  size(system_data%species), &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(system_data%masses, &
                  size(system_data%masses), &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(system_data%elements, &
                  size(system_data%elements), &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(system_data%multipole_values, &
                  size(system_data%multipole_values), &
                  parai%io_source, &
                  parai%cp_grp)
    call mp_bcast(sp_map_shape, &
                  size(sp_map_shape), &
                  parai%io_source, &
                  parai%cp_grp)

    if (.not. allocated(system_data%species_map)) &
        allocate(system_data%species_map(sp_map_shape(1), sp_map_shape(2)))

    call mp_bcast(system_data%species_map, &
                  size(system_data%species_map), &
                  parai%io_source, &
                  parai%cp_grp)

    do i = 1, mimic_control%num_client
        call mp_bcast(system_data%atom_fragment_ids(:,:,i), &
                      size(system_data%atom_fragment_ids(:,:,i)), &
                      parai%io_source, parai%cp_grp)
    end do

    allocate(mimic_control%covalent_radius(maxsys%nsx + sizes%num_species))
    do i = 1, sizes%num_species
        mimic_control%covalent_radius(maxsys%nsx + i) = mimic_covrad(system_data%elements(i)) * mimic_bohr
    end do
    rmass%pma0(maxsys%nsx + 1:maxsys%nsx + sizes%num_species) = system_data%masses
    ions0%iatyp(maxsys%nsx + 1:maxsys%nsx + sizes%num_species) = system_data%elements
    ions0%na(maxsys%nsx + 1 : mimic_control%num_species) = sizes%atoms_pspecies
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_request_sizes

!> initialize constraints information
subroutine mimic_init_constraints

    integer, dimension(:, :), allocatable :: temp_cnst
    real(real_8), dimension(:), allocatable :: temp_cnst_val
    real(real_8), dimension(:), allocatable :: temp_grate
    real(real_8), dimension(:), allocatable :: temp_cnsval_dest
    real(real_8), dimension(:, :), allocatable :: temp_cnpar
    integer :: tot_const_num
    integer :: n_code, n_bond, n_angle
    integer :: offset

    character(len=*), parameter :: procedureN = 'mimic_ifc_request_sizes'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    tot_const_num = cotc0%mcnstr

    if (mimic_control%external_constraints) then
         do n_code = 1, mimic_control%num_client
             tot_const_num = tot_const_num + &
                             sizes%nbonds_pcode(n_code) + &
                             sizes%nangles_pcode(n_code)
         end do
    endif

    if (tot_const_num == 0) then
        return
    end if

    if (allocated(ntcnst)) then
        call move_alloc(ntcnst, temp_cnst)
        call move_alloc(cnsval, temp_cnst_val)
        call move_alloc(cnpar, temp_cnpar)
        call move_alloc(cnsval_dest, temp_cnsval_dest)
        call move_alloc(grate, temp_grate)
    end if

    allocate(ntcnst(6, tot_const_num))
    allocate(cnsval(tot_const_num))
    allocate(cnpar(2, tot_const_num))
    allocate(cnsval_dest(tot_const_num))
    allocate(grate(tot_const_num))

    ntcnst(:,:) = 0.0_real_8
    grate(:) = 0.0_real_8
    cnpar(:, :) = 0.0_real_8
    cnsval_dest(:) = 0.0_real_8
    cnsval(:) = 0.0_real_8

    if (allocated(temp_cnst))then
        ntcnst(1:6, 1:cotc0%mcnstr) = temp_cnst(1:6, 1:cotc0%mcnstr)
        cnsval(1:cotc0%mcnstr) = temp_cnst_val
        grate(1:cotc0%mcnstr) = temp_grate
        cnpar(1:2, 1:cotc0%mcnstr) = temp_cnpar
        cnsval_dest(1:cotc0%mcnstr) = temp_cnsval_dest
    end if

    offset = cotc0%mcnstr + 1

    if (mimic_control%external_constraints) then
        do n_code = 1, mimic_control%num_client
            do n_bond = 1, sizes%nbonds_pcode(n_code)
                ntcnst(1, offset) = 1
                ntcnst(2, offset) = system_data%bonds(n_bond, n_code)%atom_i
                ntcnst(3, offset) = system_data%bonds(n_bond, n_code)%atom_j
                cnsval(offset) = system_data%bonds(n_bond, n_code)%length
                offset = offset + 1
            end do

            do n_angle = 1, sizes%nangles_pcode(n_code)
                ntcnst(1, offset) = 2
                ntcnst(2, offset) = system_data%angles(n_angle, n_code)%atom_i
                ntcnst(3, offset) = system_data%angles(n_angle, n_code)%atom_j
                ntcnst(4, offset) = system_data%angles(n_angle, n_code)%atom_k
                cnsval(offset) = system_data%angles(n_angle, n_code)%angle
                offset = offset + 1
            end do
        end do
    endif

    cotc0%mcnstr = tot_const_num

    call mimic_count_trcnst(const_nconn, const_conn, system_data%bonds, sizes%nbonds_pcode)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_init_constraints

!> intialize MiMiC's data structures
subroutine mimic_ifc_init(rhoe, extf, tau)

    !> electronic density array
    real(real_8), dimension(fpar%kr1, fpar%kr2s, fpar%kr3s), intent(in) :: rhoe
    !> external potential acting on the electronic grid
    real(real_8), dimension(fpar%kr1, fpar%kr2s, fpar%kr3s), intent(in) :: extf
    !> CPMD coordinates array
    real(real_8), dimension(:,:,:), intent(inout) :: tau

    real(real_8), dimension(3) :: rdiff
    real(real_8), dimension(3) :: origin, xm
    real(real_8), dimension(3,3) :: box
    integer, dimension(3) :: n_points, n_points_r
    integer :: remainder
    integer :: i
    integer :: num_atoms

    character(len=*), parameter :: procedureN = 'mimic_ifc_init'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    allocate(wrapped_coordinates, source=tau)
    allocate(mimic_forces, mold=tau)

    mimic_forces = 0.0_real_8

    origin(:) = 0.0_real_8
    box(:,:) = 0.0_real_8
    box(1,1) = cell_com%celldm(1)
    box(2,2) = cell_com%celldm(1) * cell_com%celldm(2)
    box(3,3) = cell_com%celldm(1) * cell_com%celldm(3)

    n_points(1) = spar%nr1s
    n_points(2) = spar%nr2s
    n_points(3) = spar%nr3s

    n_points_r(1) = fpar%kr1
    n_points_r(2) = fpar%kr2s
    n_points_r(3) = fpar%kr3s

    call mimic_allocate_mm_struct(subsystems, &
                                  mimic_control%factors, &
                                  wrapped_coordinates, &
                                  mimic_forces, &
                                  ions0%na(1:mimic_control%num_quantum_species), &
                                  mimic_control%num_species, &
                                  mimic_control%covalent_radius, &
                                  sizes, &
                                  system_data)

    do i = 1, size(subsystems)
        CALL mp_bcast(subsystems(i)%num_atoms, parai%source, parai%cp_grp)
    end do

    if (paral%io_parent) then
        call mimic_init_constraints
    end if

    call mp_bcast(cotc0%mcnstr, &
                  parai%io_source, &
                  parai%cp_grp)

    if (cotc0%mcnstr > 0) then
        if (.not. paral%io_parent) then
            if (allocated(ntcnst)) then
                deallocate(ntcnst)
                deallocate(cnsval)
                deallocate(cnpar)
                deallocate(cnsval_dest)
                deallocate(grate)
            end if
            allocate(ntcnst(6, cotc0%mcnstr))
            allocate(cnsval(cotc0%mcnstr))
            allocate(cnpar(2, cotc0%mcnstr))
            allocate(cnsval_dest(cotc0%mcnstr))
            allocate(grate(cotc0%mcnstr))
        end if

        call mp_bcast(ntcnst, &
                      size(ntcnst), &
                      parai%io_source, &
                      parai%cp_grp)
        call mp_bcast(cnsval, &
                      size(cnsval), &
                      parai%io_source, &
                      parai%cp_grp)
    end if

    if (paral%io_parent) then
        call mimic_collect_coordinates(subsystems)
    end if
    CALL mp_bcast(wrapped_coordinates, size(wrapped_coordinates), parai%io_source, parai%cp_grp)
    call mimic_allocate_qm_struct(quantum_fragment, ions0%na(1:mimic_control%num_quantum_species), &
                                ions0%zv(1:mimic_control%num_quantum_species), &
                                origin, box, n_points, n_points_r, rhoe, extf, &
                                wrapped_coordinates, mimic_forces)

    call mimic_center_qm(rdiff, subsystems, quantum_fragment)

    xm(1) = box(1,1) * 0.5_real_8
    xm(2) = box(2,2) * 0.5_real_8
    xm(3) = box(3,3) * 0.5_real_8

    call mimic_min_image(mimic_control%box, xm, subsystems)
    tau = wrapped_coordinates

    call initialize_tensors()

    if (mimic_control%long_range_coupling) then
        call mimic_compute_multipoles(quantum_fragment, mimic_control%multipole_order, &
                                    parap%NRXPL(parai%mepos,1), parap%NRXPL(parai%mepos,2))
        call mp_sum(quantum_fragment%electronic_multipoles%values, &
                    size(quantum_fragment%electronic_multipoles%values), &
                    parai%allgrp)
    else
        call mimic_ifc_sort_fragments()
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_init

!> re-center the QM subsystem and wrap the box with PBC
subroutine mimic_update_coords(tau, c0, cm, nstates, npws, inyh)

    !> CPMD coordinate array
    real(real_8), intent(inout) :: tau(:,:,:)
    !> WF positions and velocities
    complex(kind=real_8), dimension(nkpt%ngwk, nstates, 2), intent(inout) :: c0, cm
    !> number of states and PWs
    integer, intent(in) :: nstates, npws
    !> Miller indices of PWs
    integer, dimension(:,:), intent(in) :: inyh

    real(real_8), dimension(3) :: grid_center
    real(real_8), dimension(3) :: rdiff
    real(real_8), dimension(3) :: rtrans, gdiff, xm
    real(real_8) :: gdot
    complex(real_8), dimension(:), allocatable :: eig_trans
    integer :: npw, nstate
    integer :: nk, nh
    real(real_8), dimension(3,3) :: box

    character(*), parameter :: procedureN = 'mimic_update_coords'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    grid_center(1) = spar%nr1s / 2 + 1
    grid_center(2) = spar%nr2s / 2 + 1
    grid_center(3) = spar%nr3s / 2 + 1


    box(1,1) = cell_com%celldm(1)
    box(2,2) = cell_com%celldm(1) * cell_com%celldm(2)
    box(3,3) = cell_com%celldm(1) * cell_com%celldm(3)

    wrapped_coordinates = tau

    allocate(eig_trans(npws))

    call mimic_center_qm(rdiff, subsystems, quantum_fragment)
    tau = wrapped_coordinates

    !$OMP PARALLEL PRIVATE(npw, rtrans, gdiff, gdot, nstate)
    !$OMP DO
    do npw = 1, npws
        rtrans = real(inyh(1:3,npw), real_8) - grid_center

        gdiff(1) = rtrans(1) * gvec_com%b1(1) &
                 + rtrans(2) * gvec_com%b2(1) &
                 + rtrans(3) * gvec_com%b3(1)
        gdiff(2) = rtrans(1) * gvec_com%b1(2) &
                 + rtrans(2) * gvec_com%b2(2) &
                 + rtrans(3) * gvec_com%b3(2)
        gdiff(3) = rtrans(1) * gvec_com%b1(3) &
                 + rtrans(2) * gvec_com%b2(3) &
                 + rtrans(3) * gvec_com%b3(3)

        gdot = parm%tpiba * dot_product(gdiff, rdiff)
        eig_trans(npw) = cmplx(cos(gdot), -sin(gdot), real_8)
    end do
    !$OMP END DO

    !$OMP DO
    do nstate = 1, nstates
        do npw = 1, npws
            c0(npw,nstate,1) = c0(npw,nstate,1) * eig_trans(npw)
            cm(npw,nstate,1) = cm(npw,nstate,1) * eig_trans(npw)
        end do !npw
    end do !nstate
    !$OMP END DO

    if (cntl%tmdbo .and. cntl%textrap) then
        do nh = 1, cnti%mextra
            do nk = 1, nkpt%nkpnt
                !$OMP DO
                do nstate = 1, nstates
                   do npw = 1, npws
                       cold(npw,nstate,nk,nh) = cold(npw,nstate,nk,nh) * eig_trans(npw)
                   enddo
                enddo
                !$OMP END DO
            enddo
        enddo
    endif
    !$OMP END PARALLEL

    xm(1) = box(1,1) * 0.5_real_8
    xm(2) = box(2,2) * 0.5_real_8
    xm(3) = box(3,3) * 0.5_real_8

    call mimic_min_image(mimic_control%box, xm, subsystems)

    clsaabox%mm_c_trans = clsaabox%mm_c_trans + rdiff
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_update_coords

!> send atomic coordinates to clients
subroutine mimic_ifc_send_coordinates

    real(real_8), dimension(3) :: shift

    character(*), parameter :: procedureN = 'mimic_send_coordinates'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    shift(1) = mimic_control%box(1,1) - cell_com%celldm(1)
    shift(2) = mimic_control%box(2,2) - cell_com%celldm(1) * cell_com%celldm(2)
    shift(3) = mimic_control%box(3,3) - cell_com%celldm(1) * cell_com%celldm(3)
    shift = -shift * 0.5_real_8 + clsaabox%mm_c_trans

    call mimic_translate(quantum_fragment, subsystems, -shift)

    call mimic_send_coords(subsystems)

    call mimic_translate(quantum_fragment, subsystems, shift)
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_ifc_send_coordinates

!> determine QM and MM degrees of freedom
subroutine mimic_subsystem_dof()

    integer :: i, j, ia, iat, is
    integer :: num_dof
    real(real_8) :: qm_constraints
    real(real_8) :: mm_constraints
    real(real_8) :: const

    character(*), parameter :: procedureN = 'mimic_subsystem_dof'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    ! compute constraints contributions to DoFs.
    qm_constraints = 0.0_real_8
    mm_constraints = 0.0_real_8
    if (cotc0%mcnstr > 0) then
        do i = 1, cotc0%mcnstr
            ! set prefactor depending on constraint type.
            if (ntcnst(1,i) == 1) then
                const = 0.5_real_8
            else if (ntcnst(1,i) == 2) then
                const = 1.0_real_8/3.0_real_8
            else if (ntcnst(1,i) == 3) then
                const = 0.25_real_8
            else if (ntcnst(1,i) == 4) then
                const = 0.5_real_8
            else if (ntcnst(1,i) == 5) then
                const = 0.25_real_8
            else if (ntcnst(1,i) == 6) then
                const = 1.0_real_8
            else if (ntcnst(1,i) == 7) then
                const = 1.0_real_8/3.0_real_8
            else if (ntcnst(1,i) == 8) then
                const = 1.0_real_8
            end if
            ! Loop over list of atoms in constraint
            ! Unused entries are supposed to be 0
            do j = 2, 6
                iat = ntcnst(j,i)
                if (iat <= 0 .or. iat > mimic_control%num_atoms) then
                    const = 0.0_real_8
                else if (iat <= mimic_control%num_quantum_atoms) then
                    qm_constraints = qm_constraints + const
                else
                    mm_constraints = mm_constraints + const
                endif
            end do
        end do
    end if

    num_dof = 0
    do iat = 1, mimic_control%num_quantum_atoms
        ia = iatpt(1,i)
        is = iatpt(2,i)
        num_dof = num_dof + lskcor(1,iat) + lskcor(2,iat) + lskcor(3,iat)
    end do
    mimic_control%qm_dof = real(num_dof, kind=real_8) - qm_constraints

    num_dof = 0
    do iat = mimic_control%num_quantum_atoms + 1, mimic_control%num_atoms
        ia = iatpt(1,iat)
        is = iatpt(2,iat)
        num_dof = num_dof + lskcor(1,iat) + lskcor(2,iat) + lskcor(3,iat)
    end do
    mimic_control%mm_dof = real(num_dof, kind=real_8) - mm_constraints

    if (paral%io_parent) then
        write(6,'(1x,a,t58,i8)') 'DEGREES OF FREEDOM FOR QM:', nint(mimic_control%qm_dof)
        write(6,'(1x,a,t58,i8)') 'DEGREES OF FREEDOM FOR MM:', nint(mimic_control%mm_dof)
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_subsystem_dof

!> temperature of QM and MM subsystems separately (adapted from mm_localt)
subroutine mimic_subsystem_temperatures(velocities)

    !> velocities of QM and MM subsystems
    real(real_8), dimension(:,:,:) :: velocities

    integer :: iat, ia, is
    real(real_8) :: qm_kinetic_energy
    real(real_8) :: mm_kinetic_energy
    real(real_8) :: qm_temperature
    real(real_8) :: mm_temperature

    character(*), parameter :: procedureN = 'mimic_subsystem_temperatures'
    integer :: isub

    call tiset(procedureN, isub)

#ifdef __MIMIC
    qm_kinetic_energy = 0.0_real_8
    do iat = 1, mimic_control%num_quantum_atoms
        ia = iatpt(1,iat)
        is = iatpt(2,iat)
        qm_kinetic_energy = qm_kinetic_energy + 0.5_real_8 * rmass%pma(is) * &
                            (velocities(1,ia,is) * velocities(1,ia,is) + &
                             velocities(2,ia,is) * velocities(2,ia,is) + &
                             velocities(3,ia,is) * velocities(3,ia,is))
    end do
    if (mimic_control%qm_dof > 0.1_real_8) then
        qm_temperature = factem * 2.0_real_8 * qm_kinetic_energy / mimic_control%qm_dof
    else
        qm_temperature = 0.0_real_8
    end if

    mm_kinetic_energy = 0.0_real_8
    do iat = mimic_control%num_quantum_atoms + 1, mimic_control%num_atoms
        ia = iatpt(1,iat)
        is = iatpt(2,iat)
        mm_kinetic_energy = mm_kinetic_energy + 0.5_real_8 * rmass%pma(is) * &
                            (velocities(1,ia,is) * velocities(1,ia,is) + &
                             velocities(2,ia,is) * velocities(2,ia,is) + &
                             velocities(3,ia,is) * velocities(3,ia,is))
    end do
    if (mimic_control%mm_dof > 0.1_real_8) then
        mm_temperature = factem * 2.0_real_8 * mm_kinetic_energy / mimic_control%mm_dof
    else
        mm_temperature = 0.0_real_8
    endif

    if (paral%io_parent) then
        write(6,'(8x,a)') 'SUBSYSTEM TEMPERATURES:'
        write(6,'(8x,a,t58,f8.2)') 'T(QM) =', qm_temperature
        write(6,'(8x,a,t58,f8.2)') 'T(MM) =', mm_temperature
    end if
#else
    call stopgm(procedureN, MIMIC_MISSING_ERROR, __LINE__, __FILE__)
#endif

    call tihalt(procedureN, isub)

end subroutine mimic_subsystem_temperatures

end module mimic_wrapper
