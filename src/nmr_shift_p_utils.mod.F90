MODULE nmr_shift_p_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fpi,&
                                             uimag
  USE cppt,                            ONLY: gk,&
                                             hg
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nmr_position_p_utils,            ONLY: set_origin
  USE nmr_util_p_utils,                ONLY: make_123
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: chi,&
                                             chi_si_iso,&
                                             chi_si_to_shift_ppm,&
                                             shift_factor,&
                                             shift_matrix
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
!!use util_p_utils, only : eigr_dot
!!use nmr_util_p_utils, only : consolidate
!!use nmr_util_p_utils, only : apply_op_p
!!use nmr_util_p_utils, only : ffttor
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !public :: dia_shifts
  !public :: apply_op_gbgg
  !public :: compute_shift_dotproducts
  PUBLIC :: calc_p
  PUBLIC :: calc_l
  PUBLIC :: calc_l_one
  !!public :: nmr_interaction3
  !public :: nmr_interaction
  PUBLIC :: print_shift_matrix
  !public :: calc_shift_via_current

CONTAINS

  ! ==================================================================
  SUBROUTINE dia_shifts(c0,iB,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! Computes the diamagnetic ground state contribution to the shifts
    ! in the IGLO method: The current for field B along direction iB is
    ! j(r0) = e/2  (r0-R_i) x B  * |psi0_i(r0)|^2
    ! With the help of j(r), we calculate the shift contributions. The
    ! output is directly stored in the shift matrix.
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iB
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dia_shifts'

    COMPLEX(real_8), ALLOCATABLE             :: current(:,:), GvecJ(:)
    INTEGER                                  :: iBB, iBBB, ierr, ir, is
    REAL(real_8)                             :: overall_factor
    REAL(real_8), ALLOCATABLE                :: rhoi(:)

! ==--------------------------------------------------------------==

    ALLOCATE(rhoi(2*maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO check 2* - ?
    ALLOCATE(current(maxfft, 2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(GvecJ(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(current)!,1),maxfft)
    CALL zeroing(rhoi)!,        2*maxfft)
    CALL zeroing(GvecJ)!,       maxfft)

    CALL make_123(iB,iBB,iBBB)
    ! j(*,1) holds j_iBB, j(*,2) holds j_iBBB

    DO is=1,nstate
       CALL ffttor(c0(1,is),rhoi,psi,ncpw%ngw,.FALSE.)
       !$omp parallel do private(ir)
       DO ir=1,fpar%nnr1
          rhoi(ir) = rhoi(ir)**2
       ENDDO

       ! The current component j1 = j_iB+1 holds + (r-R)_iB+2 * rho_i,
       ! and                   j2 = j_iB+2 holds - (r-R)_iB+1 * rho_i.
       ! Compute it:
       CALL set_origin(is)
       CALL apply_op_rx(rhoi,+0.5_real_8,current(1,1),1._real_8,iBBB)
       CALL apply_op_rx(rhoi,-0.5_real_8,current(1,2),1._real_8,iBB)
    ENDDO                     ! istates

    CALL fft2tog(current(1,1),current(1,1),&
         current(1,2),current(1,2), psi,ncpw%nhg,.TRUE.)
    ! At this point, we have applied a field B along iB. This field has
    ! induced a current J = current(1) * e_iB+1 + current(2) * e_iB+2.
    ! This current now induces a field at location R with the strength
    ! exp iGR  G x J  / |G|^2. Compute its three (!) components:


    overall_factor = -2._real_8 * shift_factor * parm%alat
    CALL dscal(2*ncpw%nhg,overall_factor, current(1,1),1)
    CALL dscal(2*ncpw%nhg,overall_factor, current(1,2),1)


    ! First, treat the contribution...
    ! from j_iB+1 = current(1)
    ! to Bind_iB+2
    ! using G_iB
    CALL apply_op_gbgg(current(1,1),iB,GvecJ)
    CALL compute_shift_dotproducts(iB,iBBB,GvecJ,+1._real_8)

    ! from j_iB+1 = current(1)
    ! to Bind_iB
    ! using G_iBBB
    CALL apply_op_gbgg(current(1,1),iBBB,GvecJ)
    CALL compute_shift_dotproducts(iB,ib,GvecJ,-1._real_8)

    ! from j_iB+2 = current(2)
    ! to Bind_iB+1
    ! using G_iB
    CALL apply_op_gbgg(current(1,2),iB,GvecJ)
    CALL compute_shift_dotproducts(iB,iBB,GvecJ,-1._real_8)

    ! from j_iB+2 = current(2)
    ! to Bind_iB
    ! using G_iB+1
    CALL apply_op_gbgg(current(1,2),iBB,GvecJ)
    CALL compute_shift_dotproducts(iB,ib,GvecJ,+1._real_8)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(rhoi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(current,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(GvecJ,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dia_shifts
  ! ==================================================================
  SUBROUTINE apply_op_gbgg(j,ik,GvecJ)
    ! ==--------------------------------------------------------------==
    ! computes
    ! GvecJ = -i G_ik / G^2  j(G).
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: j(ncpw%nhg)
    INTEGER                                  :: ik
    COMPLEX(real_8)                          :: GvecJ(ncpw%nhg)

    INTEGER                                  :: ig, ig1

! ==--------------------------------------------------------------==

    IF (geq0) THEN
       ig1=2
       GvecJ(1) = CMPLX(0._real_8,0._real_8,kind=real_8)
    ELSE
       ig1=1
    ENDIF
    !$omp parallel do private(ig)
    DO ig=ig1,ncpw%nhg
       GvecJ(ig) = CMPLX(0._real_8,-gk(ik,ig)/hg(ig),kind=real_8) * j(ig)
    ENDDO
    RETURN
  END SUBROUTINE apply_op_gbgg
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_p(c0,c1,nstate,nu,psi)
    ! ==================================================================
    ! PARALLEL
    ! ***
    ! OPTIMIZED MAXIMALLY FOR CPU TIME,
    ! NEEDS scratch of size 12*nnr1.
    ! ***
    ! ==--------------------------------------------------------------==
    ! EXPECTS input:     c0 gnd state wfns,
    ! \                  c1 response wavefunctions onto the perturbation
    ! \                     h1 = p_j
    ! \                  j  the index of the component of p
    ! \                     needed to put the result in the right place
    ! CALCULATES OUTPUT: shift_alpha,beta  =
    ! 4 pi/G^2 exp(-iGR) * epsi(alpha,kappa,gamma) epsi(beta,mu,nu)
    ! * G_gamma
    ! * (psi0 r_mu grad_kappa psi1p_nu - psi1p_nu r_mu grad_kappa psi0)(r)
    ! 
    ! WHICH EVALUATES TO:
    ! 
    ! 4 pi/G^2 exp (-iGR) * epsi(alpha,kappa,gamma) epsi(beta,mu,nu)
    ! * (2 r_mu (grad_gamma psi1p_nu) (grad_kappa psi0)         [TERM 1]
    ! >   +  delta(gamma,mu) psi1p_nu (grad_kappa psi0)         [TERM 2]
    ! >   +  delta(kappa,mu) psi0     (grad_gamma psi1p_nu))(r) [TERM 3]
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c1(ncpw%ngw,*)
    INTEGER                                  :: nstate, nu
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_p'

    INTEGER                                  :: alpha, beta, gamma, iat, &
                                                ibeta, ierr, igamma, ikappa, &
                                                is, isa, isp, isub, kappa, &
                                                nu1, nu2
    REAL(real_8)                             :: contrib
    REAL(real_8), ALLOCATABLE                :: ppsi0(:,:), ppsi1(:,:), &
                                                produkt(:,:,:), psi0(:), &
                                                psi1(:)
    REAL(real_8), EXTERNAL                   :: eigr_dot

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset('  shift(p)',isub)
    CALL make_123(nu,nu1,nu2)
    ALLOCATE(ppsi0(fpar%nnr1, 2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ppsi1(fpar%nnr1, 2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(produkt(fpar%nnr1, 3, 2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi0(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi1(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(ppsi0)!, 2*nnr1)
    CALL zeroing(ppsi1)!, 2*nnr1)
    CALL zeroing(produkt)!, 6*nnr1)
    CALL zeroing(psi0)!, nnr1)
    CALL zeroing(psi1)!, nnr1)
    ! ==--------------------------------------------------------------==
    DO is=1,nstate
       CALL set_origin(is)
       CALL fft2tor(c0(1,is),psi0,c1(1,is),psi1,psi,ncpw%ngw,.FALSE.)
       ! ==--  ARRAYS:  grad_1,2 psi0  AND  grad_2,3 psi1:  -------------==
       CALL apply_op_p(c0(1,is),ppsi0(1,1),1,ncpw%ngw)
       CALL apply_op_p(c0(1,is),ppsi0(1,2),2,ncpw%ngw)
       CALL apply_op_p(c1(1,is),ppsi1(1,1),2,ncpw%ngw)
       CALL apply_op_p(c1(1,is),ppsi1(1,2),3,ncpw%ngw)
       CALL fft2tor(ppsi0,ppsi0,ppsi0(1,2),ppsi0(1,2),psi,ncpw%ngw,.FALSE.)
       CALL fft2tor(ppsi1,ppsi1,ppsi1(1,2),ppsi1(1,2),psi,ncpw%ngw,.FALSE.)
       ! ==--------------------------------------------------------------==
       ! TERM 1:
       CALL consolidate(ppsi0(1,1),ppsi1(1,1), psi)! kappa=1 gamma=2 -> a=3, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,3,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,3,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -

       CALL consolidate(ppsi0(1,1),ppsi1(1,2), psi)! kappa=1 gamma=3 -> a=2, -
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,2,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,2,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -

       CALL consolidate(ppsi0(1,2),ppsi1(1,2), psi)! kappa=2 gamma=3 -> a=1, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,1,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,1,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -
       ! ==--------------------------------------------------------------==
       ! TERM 2:
       DO kappa=1,2
          ikappa = kappa
          CALL consolidate(psi1,ppsi0(1,ikappa),psi)
          IF (kappa.NE.nu) THEN
             alpha  = nu
             beta   = kappa
             IF (beta.EQ.nu1) ibeta = 1
             IF (beta.EQ.nu2) ibeta = 2
             CALL daxpy(fpar%nnr1,crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
          ELSE
             DO alpha=1,3
                beta = alpha
                IF (beta.NE.nu) THEN
                   IF (beta.EQ.nu1) ibeta = 1
                   IF (beta.EQ.nu2) ibeta = 2
                   CALL daxpy(fpar%nnr1,-crge%f(is,1),psi,1,&
                        produkt(1,alpha,ibeta),1)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       ! TERM 3:
       DO gamma=2,3
          igamma = gamma-1
          CALL consolidate(psi0,ppsi1(1,igamma),psi)
          IF (gamma.EQ.nu) THEN
             DO alpha=1,3
                beta = alpha
                IF (beta.NE.nu) THEN
                   IF (beta.EQ.nu1) ibeta = 1
                   IF (beta.EQ.nu2) ibeta = 2
                   CALL daxpy(fpar%nnr1,crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
                ENDIF
             ENDDO
          ELSE
             alpha = nu
             beta  = gamma
             IF (beta.EQ.nu1) ibeta = 1
             IF (beta.EQ.nu2) ibeta = 2
             CALL daxpy(fpar%nnr1,-crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
          ENDIF
       ENDDO
       ! ==--  NEW ARRAYS: grad_3 psi0  AND  grad_1 psi1  ---------------==
       CALL apply_op_p(c0(1,is),ppsi0(1,1),3,ncpw%ngw)
       CALL apply_op_p(c1(1,is),ppsi1(1,2),1,ncpw%ngw)
       CALL fft2tor(ppsi0,ppsi0,ppsi1(1,2),ppsi1(1,2),psi,ncpw%ngw,.FALSE.)
       ! ==--  TERM 1:  -------------------------------------------------==
       CALL consolidate(ppsi0(1,2),ppsi1(1,2), psi)! kappa=2 gamma=1 -> a=3, -
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,3,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,3,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -

       CALL consolidate(ppsi0(1,1),ppsi1(1,2), psi)! kappa=3 gamma=1 -> a=2, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,2,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,2,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -

       CALL consolidate(ppsi0(1,1),ppsi1(1,1), psi)! kappa=3 gamma=2 -> a=1, -
       CALL apply_op_rx(psi,-2._real_8*crge%f(is,1),produkt(1,1,1),1._real_8,nu2)! beta=nu+1, mu=nu+2, +
       CALL apply_op_rx(psi,+2._real_8*crge%f(is,1),produkt(1,1,2),1._real_8,nu1)! beta=nu+2, mu=nu+1, -
       ! ==--------------------------------------------------------------==
       ! TERM 2:
       kappa=3
       ikappa = 1
       CALL consolidate(psi1,ppsi0(1,ikappa),psi)
       IF (kappa.NE.nu) THEN
          alpha  = nu
          beta   = kappa
          IF (beta.EQ.nu1) ibeta = 1
          IF (beta.EQ.nu2) ibeta = 2
          CALL daxpy(fpar%nnr1,crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
       ELSE
          DO alpha=1,3
             beta = alpha
             IF (beta.NE.nu) THEN
                IF (beta.EQ.nu1) ibeta = 1
                IF (beta.EQ.nu2) ibeta = 2
                CALL daxpy(fpar%nnr1,-crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
             ENDIF
          ENDDO
       ENDIF
       ! TERM 3:
       gamma  = 1
       igamma = 2
       CALL consolidate(psi0,ppsi1(1,igamma),psi)
       IF (gamma.EQ.nu) THEN
          DO alpha=1,3
             beta = alpha
             IF (beta.NE.nu) THEN
                IF (beta.EQ.nu1) ibeta = 1
                IF (beta.EQ.nu2) ibeta = 2
                CALL daxpy(fpar%nnr1,crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
             ENDIF
          ENDDO
       ELSE
          alpha = nu
          beta  = gamma
          IF (beta.EQ.nu1) ibeta = 1
          IF (beta.EQ.nu2) ibeta = 2
          CALL daxpy(fpar%nnr1,-crge%f(is,1),psi,1,produkt(1,alpha,ibeta),1)
       ENDIF
       ! ==--------------------------------------------------------------==
    ENDDO                     ! is

    CALL fft2tog(produkt(1,1,1),produkt(1,1,1),&
         produkt(1,1,2),produkt(1,1,2), psi,ncpw%nhg,.TRUE.)
    CALL fft2tog(produkt(1,2,1),produkt(1,2,1),&
         produkt(1,2,2),produkt(1,2,2), psi,ncpw%nhg,.TRUE.)
    CALL fft2tog(produkt(1,3,1),produkt(1,3,1),&
         produkt(1,3,2),produkt(1,3,2), psi,ncpw%nhg,.TRUE.)
    CALL nmr_interaction3&
         (produkt(1,1,1),produkt(1,2,1),produkt(1,3,1))
    CALL nmr_interaction3&
         (produkt(1,1,2),produkt(1,2,2),produkt(1,3,2))

    isa = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          isa = isa + 1
          DO alpha=1,3
             contrib=0.5_real_8*shift_factor&
                  *eigr_dot(isa,produkt(1,alpha,1))
             shift_matrix(nu1,alpha,iat,isp) &
                  =shift_matrix(nu1,alpha,iat,isp) + contrib

             contrib=0.5_real_8*shift_factor&
                  *eigr_dot(isa,produkt(1,alpha,2))
             shift_matrix(nu2,alpha,iat,isp) &
                  =shift_matrix(nu2,alpha,iat,isp) + contrib
             ! NB: There is a nonvanishing (but <0.0001) probability that
             ! nu1/2 and alpha indices are to be exchanged. But I think this way
             ! it is right.
          ENDDO           ! alpha
       ENDDO               ! atoms
    ENDDO                   ! atoms

    DEALLOCATE(ppsi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ppsi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(produkt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  shift(p)',isub)
    RETURN
  END SUBROUTINE calc_p
  ! ==================================================================
  SUBROUTINE calc_l(c0,c1,nstate,a,psi,ppsi)
    ! ==================================================================
    ! PARALLEL
    ! ***
    ! OPTIMIZED MAXIMALLY FOR CPU TIME,
    ! NEEDS scratch of size 9*nnr1.
    ! ***
    ! NB: ALL arrays *_contrib stay separately on the processors,
    ! only at the very end, they are GLOSUMmed. (explicitely in the
    ! main program)
    ! ==--------------------------------------------------------------==
    ! EXPECTS input:     c0 gnd state wfns,
    ! \                  c1 response wavefunctions onto the perturbation
    ! \                     h1 = L_a
    ! \                  a  the direction of the field (= the comp of L)
    ! \                     needed to put the result in the right place
    ! CALCULATES output: L_contrib(nat,nsp,l,k,a) [l,k=1..3] onto which 
    ! \                  the following contribution is added:
    ! \                  Sum_is  <phi0|  p_l  1/|r-R|  p_k  |phi1[L_a]>
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: a
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: ppsi(fpar%nnr1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_l'

    INTEGER                                  :: b, bb, bbb, iat, ierr, ir, &
                                                is, isa, isp, isub
    REAL(real_8), ALLOCATABLE                :: ppsi0(:,:), ppsi1(:,:), &
                                                produkt(:,:)
    REAL(real_8), EXTERNAL                   :: eigr_dot

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset('  shift(L)',isub)
    ALLOCATE(ppsi0(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ppsi1(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(produkt(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(ppsi0)
    CALL zeroing(ppsi1)
    CALL zeroing(produkt)!,3*nnr1)
    DO is=1,nstate
       DO b=1,3
          CALL apply_op_p(c0(1,is), ppsi0(1,b),  b, ncpw%ngw)
          CALL apply_op_p(c1(1,is), ppsi1(1,b),  b, ncpw%ngw)
          CALL fft2tor(ppsi0(1,b),ppsi0(1,b),&
               ppsi1(1,b),ppsi1(1,b), psi,ncpw%ngw,.FALSE.)
       ENDDO

       DO b=1,3
          CALL make_123(b,bb,bbb)
          DO ir=1,fpar%nnr1
             produkt(ir,b) = produkt(ir,b)&
                  +crge%f(is,1)*( ppsi1(ir,bb)*ppsi0(ir,bbb)&
                  -ppsi1(ir,bbb)*ppsi0(ir,bb) )
          ENDDO             ! ir
       ENDDO                 ! b = 1,3
    ENDDO                     ! states

    CALL fft2tog(produkt(1,1),produkt(1,1),&
         produkt(1,2),produkt(1,2), psi,ncpw%nhg,.TRUE.)
    CALL  ffttog(produkt(1,3),produkt(1,3), psi,ncpw%nhg,.TRUE.)
    CALL nmr_interaction3(produkt(1,1),produkt(1,2),produkt(1,3))

    isa = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          isa = isa + 1
          shift_matrix(a,1,iat,isp) = shift_matrix(a,1,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,1))
          shift_matrix(a,2,iat,isp) = shift_matrix(a,2,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,2))
          shift_matrix(a,3,iat,isp) = shift_matrix(a,3,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,3))
       ENDDO                 ! atoms
    ENDDO                     ! species
    DEALLOCATE(ppsi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ppsi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(produkt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  shift(L)',isub)
    RETURN
  END SUBROUTINE calc_l

  ! ==================================================================
  SUBROUTINE calc_l_one(c0,c1,nstate,a,psi,ppsi,istate)
    ! ==================================================================
    ! IDENTICAL to calc_L, except that the stuff is only done for istate
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: a
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: ppsi(fpar%nnr1)
    INTEGER                                  :: istate

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_l_one'

    INTEGER                                  :: b, bb, bbb, iat, ierr, ir, &
                                                isa, isp, isub
    REAL(real_8), ALLOCATABLE                :: ppsi0(:,:), ppsi1(:,:), &
                                                produkt(:,:)
    REAL(real_8), EXTERNAL                   :: eigr_dot

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset('  shift(L)',isub)
    ALLOCATE(ppsi0(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ppsi1(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(produkt(fpar%nnr1, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(ppsi0)
    CALL zeroing(ppsi1)
    CALL zeroing(produkt)!,3*nnr1)
    DO b=1,3
       CALL apply_op_p(c0(1,istate), ppsi0(1,b),  b, ncpw%ngw)
       CALL apply_op_p(c1(1,istate), ppsi1(1,b),  b, ncpw%ngw)
       CALL fft2tor(ppsi0(1,b),ppsi0(1,b),&
            ppsi1(1,b),ppsi1(1,b), psi,ncpw%ngw,.FALSE.)
    ENDDO

    DO b=1,3
       CALL make_123(b,bb,bbb)
       DO ir=1,fpar%nnr1
          produkt(ir,b) = produkt(ir,b)&
               +crge%f(istate,1)*( ppsi1(ir,bb)*ppsi0(ir,bbb)&
               -ppsi1(ir,bbb)*ppsi0(ir,bb) )
       ENDDO                 ! ir
    ENDDO                     ! b = 1,3

    CALL fft2tog(produkt(1,1),produkt(1,1),&
         produkt(1,2),produkt(1,2), psi,ncpw%nhg,.TRUE.)
    CALL  ffttog(produkt(1,3),produkt(1,3), psi,ncpw%nhg,.TRUE.)
    CALL nmr_interaction3(produkt(1,1),produkt(1,2),produkt(1,3))

    isa = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          isa = isa + 1
          shift_matrix(a,1,iat,isp) = shift_matrix(a,1,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,1))
          shift_matrix(a,2,iat,isp) = shift_matrix(a,2,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,2))
          shift_matrix(a,3,iat,isp) = shift_matrix(a,3,iat,isp)&
               +   shift_factor*eigr_dot(isa,produkt(1,3))
       ENDDO                 ! atoms
    ENDDO                     ! species
    DEALLOCATE(ppsi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ppsi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(produkt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  shift(L)',isub)
    RETURN
  END SUBROUTINE calc_l_one
  ! ==================================================================
  SUBROUTINE nmr_interaction(func1,func2)
    ! ==--------------------------------------------------------------==
    ! INPUT:   Two functions f(G), packed storage, dim=NHG
    ! OUTPUT:  Each function is overwritten by  4 pi / G^2 * f(G)
    ! \        The G=0 component is NOT CHANGED.
    ! PARALLEL
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: func1(2,*), func2(2,*)

    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: x, xx

    xx = fpi / parm%tpiba2
    ig1=1
    IF (geq0) ig1=2
    !$omp parallel do private(ig,x)
    DO ig=ig1,ncpw%nhg
       x = xx / hg(ig)
       func1(1,ig) = func1(1,ig) * x
       func1(2,ig) = func1(2,ig) * x
       func2(1,ig) = func2(1,ig) * x
       func2(2,ig) = func2(2,ig) * x
    ENDDO
    RETURN
  END SUBROUTINE nmr_interaction
  ! ==================================================================
  SUBROUTINE print_shift_matrix(nstate,tsusceptibility)
    ! ==--------------------------------------------------------------==
    ! PRINTS, but does NOT MODIFY the shift_matrix
    INTEGER                                  :: nstate
    LOGICAL                                  :: tsusceptibility

    CHARACTER(*), PARAMETER :: procedureN = 'print_shift_matrix'
    CHARACTER(1), DIMENSION(3), PARAMETER    :: label = (/"X","Y","Z"/)

    INTEGER                                  :: a, b, i, iat, iatom, ierr, &
                                                info, isp, j
    REAL(real_8)                             :: chksum, delta, isoshift, &
                                                scr(3,3), scr2(3,3), work(3*3)
    REAL(real_8), ALLOCATABLE                :: eigenvalues(:,:,:)

! ==--------------------------------------------------------------==

    ALLOCATE(eigenvalues(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    chksum=0.0_real_8
    ! each atom SEPARATELY:
    iatom = 0
    DO isp = 1,maxsys%nsx
       DO iat = 1, ions0%na(isp)
          iatom = iatom + 1
          IF (paral%io_parent)&
               WRITE (6, '(9("=")," ATOM ",A2,I3,'//&
               '",  NUC SHIELD MATRIX *ALL* STATES ",'//&
               '6("="),"(PV)",2("="))') elem%el(ions0%iatyp(isp)),iat
          CALL zeroing(scr)!,3*3)
          CALL zeroing(scr2)!,3*3)
          IF (tsusceptibility) THEN
             CALL daxpy(3*3,-chi_SI_to_shift_ppm,chi,1,scr,1)
             CALL daxpy(3*3,-chi_SI_to_shift_ppm,chi,1,scr2,1)
          ENDIF
          DO i=1,3
             DO j=1,3
                scr2(i,j) =  scr2(i,j) + shift_matrix(i,j,iat,isp)
                scr(i,j)  =  scr(i,j)+0.5_real_8*(shift_matrix(i,j,iat,isp)&
                     + shift_matrix(j,i,iat,isp))
             ENDDO         ! i
          ENDDO             ! j
          CALL mp_sum(scr,3*3,parai%allgrp)
          chksum=chksum+SUM(ABS(scr))
          CALL mp_sum(scr2,3*3,parai%allgrp)
          CALL dsyev('N','U', 3, scr, 3, eigenvalues(1,iat,isp),&
               work, 9, info)
          IF ((info .NE. 0).AND.paral%io_parent)&
               WRITE(6, *)'ERROR IN DSYEV, code ',info
          DO a=1,3
             IF (paral%io_parent)&
                  WRITE(6, '(A2,I3," [all] S_",A1,"* = ",3F9.2,F18.2)')&
                  elem%el(ions0%iatyp(isp)),iat,label(a),&
                  (scr2(a,b),b=1,3),&
                  eigenvalues(a,iat,isp)
          ENDDO             ! a
       ENDDO                 ! iat
    ENDDO                     ! isp
    IF (paral%io_parent)&
         WRITE(6, '(9("=")," RAW SHIELD MATRICES END ",31("="))')
    IF (paral%io_parent) WRITE(6,'(A,E12.6)') ' ChkSum(NMR) = ',chksum
    ! ==--------------------------------------------------------------==
    ! FINAL LIST, this time INCLUDING THE SUSCEPTIBILITY CORRECTION:
    IF (tsusceptibility) THEN
       delta = -chi_SI_iso*chi_SI_to_shift_ppm
       IF (paral%io_parent)&
            WRITE (6, '("===== NUC SHIELD'//&
            ' ====(PV 1)===(PV 2)===(PV 3)====(iso)==(aniso)==")')
       iatom = 0
       DO isp = 1,maxsys%nsx
          DO iat = 1, ions0%na(isp)
             iatom = iatom + 1
             isoshift = (eigenvalues(1,iat,isp)&
                  +eigenvalues(2,iat,isp)&
                  +eigenvalues(3,iat,isp)) / 3._real_8
             IF (paral%io_parent)&
                  WRITE (6, '(A2,I3,A,I3,A,5F9.2)')&
                  elem%el(ions0%iatyp(isp)),iat,'/',iatom,&
                  '   SHIELD :',&
                  (eigenvalues(i,iat,isp),i=1,3),&
                  isoshift,&
                  eigenvalues(3,iat,isp)-&
                  (eigenvalues(1,iat,isp)&
                  +eigenvalues(2,iat,isp))/2._real_8
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE (6, '(5("=")," NET SHIELD ",48("="))')
       IF (paral%io_parent)&
            WRITE (6, '(5("=")," NET SHIELD, '//&
            'susceptibility correction = ",'//&
            'F6.3,13("="))') delta
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigenvalues,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE print_shift_matrix
  ! ==================================================================
  SUBROUTINE calc_shift_via_current(j1,j2,isa,iB,shift,psi)
    ! ==--------------------------------------------------------------==
    ! PARALLEL
    ! EXPECTS: The currents j1, j2 in G-space.
    ! OUTPUT:  The shift value for atom ISA.
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: j1(ncpw%nhg), j2(ncpw%nhg)
    INTEGER                                  :: isa, iB
    REAL(real_8)                             :: shift
    COMPLEX(real_8)                          :: psi(:)

    INTEGER                                  :: iBB, iBBB, ig, isub
    REAL(real_8), EXTERNAL                   :: eigr_dot

! ==--------------------------------------------------------------==

    CALL tiset('  j/SHIFTS',isub)

    CALL make_123(iB,iBB,iBBB)

    DO ig=1,ncpw%nhg
       psi(ig) = -uimag * (j1(ig)*gk(iBBB,ig) - j2(ig)*gk(iBB,ig))
       IF (hg(ig) .GT. 1.e-4_real_8) psi(ig) = psi(ig) / hg(ig)
    ENDDO

    shift = shift - 2._real_8*parm%alat*shift_factor * eigr_dot(isa,psi)

    CALL tihalt('  j/SHIFTS',isub)
    RETURN
  END SUBROUTINE calc_shift_via_current
  ! ==================================================================

END MODULE nmr_shift_p_utils

SUBROUTINE compute_shift_dotproducts(i1,i2,GvecJ,prefactor)
  ! takes as input the array
  ! GvecJ  =  -iG_x / G^2  j_y(G)
  ! and adds the result of the expression
  ! Sum_G  exp - i G R_ion  GvecJ
  ! to the shift matrix (i1,i2)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac, ONLY : paral,parai
  USE response_pmod , ONLY:shift_matrix
  USE ions , ONLY:ions0,ions1
  IMPLICIT NONE
  INTEGER                                    :: i1, i2
  REAL(real_8)                               :: GvecJ(*), prefactor

  INTEGER                                    :: iat, isa, isp
  REAL(real_8), EXTERNAL                     :: eigr_dot

! ==--------------------------------------------------------------==

  isa = 0
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        isa = isa + 1
        shift_matrix(i1,i2,iat,isp) &
             =shift_matrix(i1,i2,iat,isp)&
             + prefactor*eigr_dot(isa,GvecJ)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE compute_shift_dotproducts

SUBROUTINE nmr_interaction3(func1,func2,func3)
  ! ==--------------------------------------------------------------==
  ! PARALLEL
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw,parm
  USE parac, ONLY : paral,parai
  USE cnst , ONLY:fpi
  USE geq0mod , ONLY:geq0
  USE cppt , ONLY:hg
  IMPLICIT NONE
  REAL(real_8)                               :: func1(2,*), func2(2,*), &
                                                func3(2,*)

  INTEGER                                    :: ig, ig1
  REAL(real_8)                               :: x, xx

! ==--------------------------------------------------------------==

  xx = fpi / parm%tpiba2
  ig1=1
  IF (geq0) ig1=2
  !$omp parallel do private(ig,x)
  DO ig=ig1,ncpw%nhg
     x = xx / hg(ig)
     func1(1,ig) = func1(1,ig) * x
     func1(2,ig) = func1(2,ig) * x
     func2(1,ig) = func2(1,ig) * x
     func2(2,ig) = func2(2,ig) * x
     func3(1,ig) = func3(1,ig) * x
     func3(2,ig) = func3(2,ig) * x
  ENDDO
  RETURN
END SUBROUTINE nmr_interaction3
! ==================================================================
