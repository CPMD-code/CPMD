module constraint_utils
    USE system,                          ONLY: cntl, cnti
    USE parac,                           ONLY: paral
    USE machine,                         ONLY: m_walltime
    USE constraint
    USE constr_utils,                    ONLY: funcr, diffr
    USE dum2_utils,                      ONLY: dum2
    USE error_handling,                  ONLY: stopgm
    USE fillc_utils,                     ONLY: fillc
    USE kinds,                           ONLY: real_8
    USE timer,                           ONLY: tihalt,&
                                               tiset

    implicit none

    ! stretch constraint type
    integer, parameter :: TYPE_STRETCH = 1

    contains

    !> Initialize constraints matrix
    function build_constraints(ntcnst, cval, na, nsp) result(matrix)

        ! CPMD constraint map
        integer, dimension(:,:), intent(in) :: ntcnst
        ! Reference values for constraints
        real(real_8), dimension(:), intent(in) :: cval
        ! Number of atoms per species
        integer, dimension(:), intent(in) :: na
        ! Number of atomic species
        integer, intent(in) :: nsp
        ! Resulting sparse matrix
        type(constraint_matrix) :: matrix

        integer :: i, j, iat, isp
        integer :: nconstr, n_coupl
        integer :: ctype, ctype2, ia, ia2, ib, ib2, id
        real(real_8), dimension(:), allocatable :: vals
        integer, dimension(:), allocatable :: nz_row
        integer, dimension(:), allocatable :: nz_column, temp_nz
        integer, dimension(:,:), allocatable :: ids
        integer, dimension(:,:,:), allocatable :: tau_map
        integer :: accumulator

        ! Count constraints
        nconstr = 0
        do i = 1, size(ntcnst, 2)
            ctype = ntcnst(1, i)
            if (ctype == TYPE_STRETCH) then
                nconstr = nconstr + 1
            else
                CALL stopgm('BUILD_CONSTRAINTS','UNSUPPORTED CONSTRAINT TYPE',&
                   __LINE__,__FILE__)
            end if
        end do ! i

        allocate(nz_row(0:nconstr))
        !FIXME update for other types of constraints?
        allocate(ids(2, nconstr))
        allocate(tau_map(2, 2, nconstr))

        nz_row = 0
        ids = 0
        tau_map = 0

        ! Check coupling between constraints
        n_coupl = 0
        accumulator = 0
        do i = 1, size(ntcnst, 2)
            ctype = ntcnst(1, i)
            ia = ntcnst(2, i)
            ib = ntcnst(3, i)
            ids(1, i) = ia
            ids(2, i) = ib
            do j = 1, size(ntcnst, 2)
                ctype2 = ntcnst(1, j)
                if (ctype2 == TYPE_STRETCH) then
                    ia2 = ntcnst(2, j)
                    ib2 = ntcnst(3, j)
                    if (ia == ia2 .or. ib == ib2 .or. &
                        ia == ib2 .or. ib == ia2) then
                        accumulator = accumulator + 1
                        n_coupl = n_coupl + 1
                        allocate(temp_nz(n_coupl))
                        temp_nz = 0
                        if (allocated(nz_column)) then
                            temp_nz(1:n_coupl - 1) = nz_column
                            deallocate(nz_column)
                        end if
                        call move_alloc(temp_nz, nz_column)
                        nz_column(n_coupl) = j
                    end if ! overlap
                end if ! ctype
            end do ! j
            id = 0
            do isp = 1, nsp
                do iat = 1, na(isp)
                    id = id  + 1
                    if (id == ia) then
                        tau_map(1, 1, i) = isp
                        tau_map(2, 1, i) = iat
                    end if
                    if (id == ib) then
                        tau_map(1, 2, i) = isp
                        tau_map(2, 2, i) = iat
                    end if
                end do
            end do
            nz_row(i) = accumulator
        end do ! i

        allocate(vals(size(nz_column)))
        vals = 0.0_real_8

        matrix = new_matrix(vals, nz_row, nz_column, ids, tau_map, cval)

    end function build_constraints

    !> Update constraint values and their gradients
    subroutine update_constraints(matrix, grad, diff, tau)

        !> Constraint matrix - needed to determine connectivity
        class(constraint_matrix), intent(inout) :: matrix
        !> Constraint gradients, output
        real(real_8), dimension(:, :), intent(out) :: grad
        !> Constraint violation, output
        real(real_8), dimension(:), intent(out) :: diff
        !> Coordinate matrix
        real(real_8), dimension(:,:,:), intent(in) :: tau

        character(*), PARAMETER :: procedureN = 'update_constraints'
        integer :: isub
        integer :: i, j
        integer :: ia, ib
        integer :: cid
        real(real_8) :: dist, dist_sq
        real(real_8), dimension(3) :: r_a, r_b
        real(real_8), dimension(6) :: d_sigma

        call tiset(procedureN,isub)

        grad = 0.0_real_8

        !$OMP PARALLEL DO PRIVATE(i, ia, ib, r_a, r_b, dist_sq, d_sigma)
        do i = 1, size(matrix%nz_row) - 1
            ia = matrix%atom_id(1, i)
            ib = matrix%atom_id(2, i)
            call fillc(ia, tau, r_a)
            call fillc(ib, tau, r_b)
            call funcr(diff(i), dist_sq, matrix%ref(i), r_a, r_b)
            call diffr(d_sigma, r_a, r_b)
            grad(:, i) = grad(:, i) + d_sigma
        end do ! i
        !$OMP END PARALLEL DO

        call tihalt(procedureN,isub)

    end subroutine update_constraints

    !> Fill symmetric constraint matrix
    subroutine fill_matrix(matrix, grad, dt, dtm)

        !> Matrix to fill (SHOULD BE INITIALIZED)
        class(constraint_matrix), intent(inout) :: matrix
        !> Constraint gradients
        real(real_8), dimension(:,:), intent(in) :: grad
        !> Ionic timestep
        real(real_8), intent(in) :: dt
        !> dt / (2 * mi)
        real(real_8), dimension(:), intent(in) :: dtm

        character(*), PARAMETER :: procedureN = 'fill_matrix'
        integer :: isub
        integer :: i, j, k, is, ia, pos1, pos2, at1, at2, cid
        real(real_8), dimension(:), allocatable :: vals

        call tiset(procedureN,isub)

        allocate(vals(size(matrix%vals)))

        vals = 0.0_real_8

        !$OMP PARALLEL DO PRIVATE(i, j, cid, at1, pos1, at2, pos2) REDUCTION(+:vals)
        do i = 1, size(matrix%nz_row) - 1
            do j = matrix%nz_row(i - 1) + 1, matrix%nz_row(i)
                cid = matrix%r_index(j)
                do at1 = 1, size(matrix%atom_id, 1)
                    pos1 = (at1 - 1) * 3 + 1
                    do at2 = 1, size(matrix%atom_id, 1)
                        pos2 = (at2 - 1) * 3 + 1
                        if (matrix%atom_id(at1, i) == matrix%atom_id(at2, cid)) then
                            vals(j) = vals(j) + &
                                      dt * dtm(matrix%tau_map(1, at1, i)) * &
                                      dot_product(grad(pos1:pos1+2, i), grad(pos2:pos2+2, cid))
                        end if
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        deallocate(matrix%vals)
        call move_alloc(vals, matrix%vals)

        call tihalt(procedureN,isub)

    end subroutine fill_matrix

    !> Fill asymmetric constraint matrix
    subroutine fill_matrix_new(matrix, gradl, gradr, dt, dtm)

        !> Matrix to fill (SHOULD BE INITIALIZED)
        class(constraint_matrix), intent(inout) :: matrix
        !> Constraint gradients
        real(real_8), dimension(:,:), intent(in) :: gradl, gradr
        !> Ionic timestep
        real(real_8), intent(in) :: dt
        !> dt / (2 * mi)
        real(real_8), dimension(:), intent(in) :: dtm

        character(*), PARAMETER :: procedureN = 'fill_matrix_new'
        integer :: isub
        integer :: i, j, k, is, ia, pos1, pos2, at1, at2, cid
        real(real_8), dimension(:), allocatable :: vals

        call tiset(procedureN,isub)

        allocate(vals(size(matrix%vals)))

        vals = 0.0_real_8

        !$OMP PARALLEL DO PRIVATE(i, j, cid, at1, pos1, at2, pos2) REDUCTION(+:vals)
        do i = 1, size(matrix%nz_row) - 1
            do j = matrix%nz_row(i - 1) + 1, matrix%nz_row(i)
                cid = matrix%r_index(j)
                do at1 = 1, size(matrix%atom_id, 1)
                    pos1 = (at1 - 1) * 3 + 1
                    do at2 = 1, size(matrix%atom_id, 1)
                        pos2 = (at2 - 1) * 3 + 1
                        if (matrix%atom_id(at1, i) == matrix%atom_id(at2, cid)) then
                            vals(j) = vals(j) + &
                                      dt * dtm(matrix%tau_map(1, at1, i)) * &
                                      dot_product(gradl(pos1:pos1+2, i), gradr(pos2:pos2+2, cid))
                        end if
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        deallocate(matrix%vals)
        call move_alloc(vals, matrix%vals)

        call tihalt(procedureN,isub)

    end subroutine fill_matrix_new


    !> Perform SHAKE constraint solver
    subroutine new_shake(matrix, xlagr, tau0, taup, velp, dt, dtm)

        !> Constraint matrix
        class(constraint_matrix), intent(inout) :: matrix
        !> Lagrange multiplier array
        real(real_8), dimension(:) :: xlagr
        !> CPMD coordinate matrix
        real(real_8), dimension(:,:,:) :: tau0
        real(real_8), dimension(:,:,:) :: taup
        !> CPMD velocity matrix
        real(real_8), dimension(:,:,:) :: velp
        !> Ionic timestep
        real(real_8), intent(in) :: dt
        !> dt / (2 * mi)
        real(real_8), dimension(:), intent(in) :: dtm

        integer :: i, at, atid, spid, cat, pos
        real(real_8), dimension(:), allocatable :: sigma
        real(real_8), dimension(:, :), allocatable :: grad_sigma
        real(real_8), dimension(:, :), allocatable :: grad_current
        real(real_8), dimension(:,:,:), allocatable :: tscr
        real(real_8), parameter :: threshold = 1.0e-7
        real(real_8) :: time1, time2
        real(real_8), external :: dnrm2

        character(*), PARAMETER :: procedureN = 'new_shake'
        integer :: isub

        call tiset(procedureN,isub)

        time1 = m_walltime()

        ! Allocate temporary storage
        allocate(sigma(size(matrix%nz_row)-1))
        allocate(grad_sigma(6, size(matrix%nz_row)-1))
        allocate(grad_current(6, size(matrix%nz_row)-1))
        allocate(tscr(size(tau0, 1), size(tau0, 2), size(tau0, 3)))

        sigma = 0.0_real_8
        grad_sigma = 0.0_real_8
        grad_current = 0.0_real_8
        tscr = 0.0_real_8

        call dum2(tau0, tscr)
        call update_constraints(matrix, grad_sigma, sigma, tscr)

        do i = 1, cnti%shake_maxstep
            call dum2(taup,tscr)
            call update_constraints(matrix, grad_current, sigma, tscr)
            call fill_matrix(matrix, grad_current, dt, dtm)
            if (dnrm2(size(sigma), sigma, 1) < threshold) then
                time2 = m_walltime()
                if (paral%io_parent) then
                    write(6,'(8x,a,i3,a,t58,f8.2)') "SHAKE CONVERGED AFTER ", i, &
                                                    " ITERATIONS. TCPU = ", &
                                                    (time2 - time1) * 0.001_real_8
                end if
                exit
            end if
            if (cntl%pbicgstab) then
                call pbicgstab_solve(matrix, xlagr, sigma, threshold)
            else
                call pcg_solve(matrix, xlagr, sigma, threshold)
            endif

            do at = 1, size(matrix%nz_row) - 1
                do cat = 1, size(matrix%atom_id, 1)
                    spid = matrix%tau_map(1, cat, at)
                    atid = matrix%tau_map(2, cat, at)
                    pos = (cat - 1) * 3 + 1
                    taup(:, atid, spid) = taup(:, atid, spid) - dtm(spid) * dt * xlagr(at) * grad_sigma(pos : pos + 2, at)
                    velp(:, atid, spid) = velp(:, atid, spid) - dtm(spid) * xlagr(at) * grad_sigma(pos : pos + 2, at)
                end do
            end do

            if (i == cnti%shake_maxstep) then
                if (paral%io_parent) then
                    write(6,'(2A,I4,A,F10.8,A,F10.8,A)') ' NEW SHAKE DID NOT ', &
                                                         'CONVERGE AFTER ', i, ' ITERATIONS. ERRX=', &
                                                         dnrm2(size(sigma), sigma, 1), &
                                                         ' TOLX=', threshold, '. STOP!'
                    if (.not. cntl%pbicgstab) then
                        write(6,'(A)') " IF PCG IS NOT CONVERGING:"
                        write(6,'(A)') " 1) TRY INCREASING SHAKE_CG_ITER IN THE &CPMD SECTION. EXAMPLE:"
                        write(6,'(A)') " SHAKE_CG_ITER"
                        write(6,'(A)') " 200"
                        write(6,'(A)') " 2) TRY USING THE ALTERNATIVE PBICGSTAB SOLVER. &
                                       IN THE &CPMD SECTION, SELECT:"
                        write(6,'(A)') " NEW CONSTRAINTS PBICGSTAB "
                    else
                        write(6,'(A)') " IF PBICGSTAB IS NOT CONVERGING:"
                        write(6,'(A)') " TRY INCREASING SHAKE_CG_ITER IN THE &CPMD SECTION. EXAMPLE:"
                        write(6,'(A)') " SHAKE_CG_ITER"
                        write(6,'(A)') " 200"
                    end if
                    call stopgm('NEW_SHAKE',' ',__LINE__,__FILE__)
                end if
                exit
            end if
        end do

        tau0(:,:,:) = taup(:,:,:)

        call tihalt(procedureN,isub)

    end subroutine new_shake

    !> Perform velocity update (within RATTLE algorithm)
    subroutine new_rattle(matrix, xlagr, tau0, velp, dt, dtm)

        !> Constraint matrix
        class(constraint_matrix), intent(inout) :: matrix
        !> Lagrange multiplier array
        real(real_8), dimension(:) :: xlagr
        !> CPMD coordinate matrix
        real(real_8), dimension(:,:,:), intent(in) :: tau0
        !> CPMD velocity matrix
        real(real_8), dimension(:,:,:), intent(inout) :: velp
        !> Ionic timestep
        real(real_8), intent(in) :: dt
        real(real_8), dimension(:), intent(in) :: dtm

        real(real_8), dimension(:,:,:), allocatable :: tscr
        real(real_8), dimension(:), allocatable :: sigma
        real(real_8), dimension(:, :), allocatable :: grad_sigma
        real(real_8), dimension(:), allocatable :: dv
        integer :: i, j, at, atid, spid, cat, pos

        character(*), PARAMETER :: procedureN = 'new_rattle'
        integer :: isub

        call tiset(procedureN,isub)

        allocate(sigma(size(matrix%nz_row)-1))
        allocate(grad_sigma(6, size(matrix%nz_row)-1))
        allocate(dv(size(matrix%nz_row)-1))
        allocate(tscr(size(velp, 1), size(velp, 2), size(velp, 3)))

        sigma = 0.0_real_8
        grad_sigma = 0.0_real_8
        dv = 0.0_real_8
        tscr = 0.0_real_8

        call dum2(tau0, tscr)
        call update_constraints(matrix, grad_sigma, sigma, tscr)
        call fill_matrix(matrix, grad_sigma, dt, dtm)

        matrix%vals = matrix%vals / dt

        !$OMP PARALLEL DO PRIVATE(i, at, cat, spid, atid, pos)
        do i = 1, size(matrix%nz_row) - 1
            do cat = 1, size(matrix%atom_id, 1)
                at = matrix%atom_id(cat, i)
                spid = matrix%tau_map(1, cat, i)
                atid = matrix%tau_map(2, cat, i)
                pos = (cat - 1) * 3 + 1
                dv(i) = dv(i) + dot_product(grad_sigma(pos:pos+2, i), velp(:, atid, spid))
            end do
        end do
        !$OMP END PARALLEL DO

        call pcg_solve(matrix, xlagr, dv, 0.0000001_real_8)

        dv = 0.0_real_8

        !$OMP PARALLEL DO PRIVATE(i, at, cat, spid, atid, pos)
        do i = 1, size(matrix%nz_row) - 1
            do cat = 1, size(matrix%atom_id, 1)
                at = matrix%atom_id(cat, i)
                spid = matrix%tau_map(1, cat, i)
                atid = matrix%tau_map(2, cat, i)
                pos = (cat - 1) * 3 + 1
                velp(:, atid, spid) = velp(:, atid, spid) - dtm(spid) * xlagr(i) * grad_sigma(pos : pos + 2, i)
            end do
        end do
        !$OMP END PARALLEL DO

        call tihalt(procedureN,isub)

    end subroutine new_rattle

    !> Use the preconditioned conjugate gradient solver
    subroutine pcg_solve(matrix, x, b, tol)
        !> Matrix to solve
        class(constraint_matrix), intent(in) :: matrix
        !> Initial guss
        real(real_8), dimension(:), intent(inout) :: x
        !> right-hand side
        real(real_8), dimension(:), intent(in) :: b
        !> desired tolerance
        real(real_8) :: tol

        integer :: i, j
        real(real_8) :: denominator
        real(real_8), dimension(:), allocatable :: r
        real(real_8), dimension(:), allocatable :: precond
        real(real_8), dimension(:), allocatable :: d, s
        real(real_8), dimension(:), allocatable :: temp
        real(real_8) :: delta_start, delta_new, delta_old
        real(real_8) :: alpha, beta
        real(real_8), external :: ddot, dnrm2

        character(*), PARAMETER :: procedureN = 'pcg_solve'
        integer :: isub

        call tiset(procedureN,isub)

        allocate(r(size(matrix%nz_row)-1))
        allocate(precond(size(matrix%nz_row)-1))
        allocate(d(size(matrix%nz_row)-1))
        allocate(s(size(matrix%nz_row)-1))
        allocate(temp(size(matrix%nz_row)-1))

        r = 0.0_real_8
        precond = 0.0_real_8
        d = 0.0_real_8
        s = 0.0_real_8
        temp = 0.0_real_8

        ! Initialize preconditioner which is inverse diagonal
        !$OMP PARALLEL DO PRIVATE(i, j)
        do i = 1, size(matrix%nz_row) - 1
            do j = matrix%nz_row(i - 1) + 1, matrix%nz_row(i)
                if (matrix%r_index(j) == i) then
                    precond(i) = 1.0_real_8 / matrix%vals(j)
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        ! Fill residual and intialize CG variables
        r = b - matrix * x
        d = precond * r
        delta_new = ddot(size(r), r, 1, d, 1)
        delta_start = delta_new

        do i = 1, cnti%shake_cg_iter
            ! Check if we have solved constraints
            if (dnrm2(size(r), r, 1) < tol) then
                exit
            end if
            temp = matrix * d
            denominator = ddot(size(d), d, 1, temp, 1)
            alpha = delta_new / denominator
            call daxpy(size(d), alpha, d(1), 1, x(1), 1)
            call daxpy(size(r), -alpha, temp(1), 1, r(1), 1)
            s = precond * r
            delta_old = delta_new
            delta_new = ddot(size(r), r, 1, s, 1)
            beta = delta_new / delta_old
            call dscal(size(d), beta, d(1), 1)
            call daxpy(size(d), 1.0_real_8, s(1), 1, d(1), 1)
            if (i == cnti%shake_cg_iter) then
                if (paral%io_parent) then
                    write(6,'(A52,I6,A16)') "WARNING! PCG DID NOT CONVERGE AFTER SHAKE_CG_ITER=", cnti%shake_cg_iter, "STEPS!"
                end if
                exit
            end if
        end do

        deallocate(r)
        deallocate(precond)
        deallocate(d)
        deallocate(s)
        deallocate(temp)
        call tihalt(procedureN,isub)

    end subroutine pcg_solve

    !> Use the preconditioned biconjugate gradient stabilized solver
    subroutine pbicgstab_solve(matrix, x, b, tol)
        !> Matrix to solve
        class(constraint_matrix), intent(in) :: matrix
        !> Initial guss
        real(real_8), dimension(:), intent(inout) :: x
        !> right-hand side
        real(real_8), dimension(:), intent(in) :: b
        !> desired tolerance
        real(real_8) :: tol

        integer :: i, j
        real(real_8), dimension(:), allocatable :: precond
        real(real_8), dimension(:), allocatable :: r, r0, p, y, v, h, s, z, t, temp
        real(real_8) :: rho_n, rho, alpha, omega, beta
        real(real_8), external :: ddot, dnrm2

        character(*), PARAMETER :: procedureN = 'PBICGSTAB_solve'
        integer :: isub

        call tiset(procedureN,isub)

        ! Allocate and initialize
        allocate(precond(size(matrix%nz_row)-1))
        allocate(r(size(matrix%nz_row)-1))
        allocate(r0(size(matrix%nz_row)-1))
        allocate(p(size(matrix%nz_row)-1))
        allocate(y(size(matrix%nz_row)-1))
        allocate(v(size(matrix%nz_row)-1))
        allocate(h(size(matrix%nz_row)-1))
        allocate(s(size(matrix%nz_row)-1))
        allocate(z(size(matrix%nz_row)-1))
        allocate(t(size(matrix%nz_row)-1))
        allocate(temp(size(matrix%nz_row)-1))

        precond = 0.0_real_8
        temp = 0.0_real_8
        r = 0.0_real_8
        r0 = 0.0_real_8
        p = 0.0_real_8
        y = 0.0_real_8
        v = 0.0_real_8
        h = 0.0_real_8
        s = 0.0_real_8
        z = 0.0_real_8
        t = 0.0_real_8
        beta = 0.0_real_8
        rho = 1.0_real_8
        rho_n = 1.0_real_8
        alpha = 1.0_real_8
        omega = 1.0_real_8

        ! Initialize preconditioner which is inverse diagonal
        !$OMP PARALLEL DO PRIVATE(i, j)
        do i = 1, size(matrix%nz_row) - 1
            do j = matrix%nz_row(i - 1) + 1, matrix%nz_row(i)
                if (matrix%r_index(j) == i) then
                    precond(i) = 1.0_real_8 / matrix%vals(j)
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        ! Fill residual and intialize CG variables
        r  = b - matrix * x
        r0 = r
        ! Check if we have already solved constraints
        if (dnrm2(size(r), r, 1) < tol) then
            deallocate(precond)
            deallocate(temp)
            deallocate(r)
            deallocate(r0)
            deallocate(v)
            deallocate(y)
            deallocate(h)
            deallocate(t)
            deallocate(s)
            deallocate(z)
            deallocate(p)
            call tihalt(procedureN,isub)
            return
        end if

        ! Main loop
        do i = 1, cnti%shake_cg_iter
            rho = rho_n
            rho_n = ddot(size(r0), r0, 1, r, 1)
            beta = ( rho_n / rho ) * ( alpha / omega )
            p = r + beta * ( p - omega * v )
            y = precond * p
            v = matrix * y
            alpha = rho_n / ddot(size(r0), r0, 1, v, 1)
            h = x + alpha * y
            ! Check convergence
            temp = b - matrix * h
            if (dnrm2(size(temp), temp, 1) < tol) then
                x = h
                exit
            end if
            s = r - alpha * v
            z = precond * s
            t = matrix * z
            omega = ddot(size(t), precond * t, 1, s, 1) / ddot(size(t), precond * t, 1, t, 1)
            x = h + omega * z
            ! Check convergence
            temp = b - matrix * x
            if (dnrm2(size(temp), temp, 1) < tol) then
                exit
            end if
            r = s - omega * t
            if (i == cnti%shake_cg_iter) then
                if (paral%io_parent) then
                    write(6,'(A52,I6,A16)') "WARNING! PBICGSTAB DID NOT CONVERGE AFTER SHAKE_CG_ITER=", cnti%shake_cg_iter, " STEPS!"
                end if
                exit
            end if
        end do

        deallocate(precond)
        deallocate(r)
        deallocate(r0)
        deallocate(v)
        deallocate(y)
        deallocate(h)
        deallocate(t)
        deallocate(s)
        deallocate(z)
        deallocate(p)
        deallocate(temp)

        call tihalt(procedureN,isub)

    end subroutine pbicgstab_solve

end module constraint_utils
