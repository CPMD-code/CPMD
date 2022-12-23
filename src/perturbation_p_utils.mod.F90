MODULE perturbation_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE forces_driver,                   ONLY: forces
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: dmbi,&
                                             nmr_options,&
                                             response1
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ortho_p
  PUBLIC :: restrain_ngw_zero
  PUBLIC :: lag_mult
  PUBLIC :: energy

CONTAINS

  ! ==================================================================
  SUBROUTINE ortho_p(nstate,c0,c1)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate), &
                                                c1(nkpt%ngwk,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ortho_p'

    INTEGER                                  :: ierr, isub
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8), ALLOCATABLE, SAVE          :: overlap(:)

    CALL tiset('ortho_p    ',isub)

    IF (ifirst .EQ. 0) THEN
       ifirst = 1
       ALLOCATE(overlap(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(overlap)!,nstate*nstate)
    ENDIF

    CALL ovlap(nstate,overlap,c0,c1)

    CALL mp_sum(overlap,nstate*nstate,parai%allgrp)
    CALL dgemm('N','N',2*ncpw%ngw,nstate,nstate,-1._real_8,c0,2*ncpw%ngw,overlap,&
         nstate,+1._real_8,c1,2*ncpw%ngw)


    CALL tihalt('ortho_p    ',isub)
    RETURN
  END SUBROUTINE ortho_p
  ! ==================================================================
  SUBROUTINE restrain_ngw_zero(cx,ngw_zero,ngw,nstate)
    INTEGER                                  :: ngw_zero, ngw, nstate
    COMPLEX(real_8)                          :: cx(ngw,nstate)

    INTEGER                                  :: ig, is

    IF (ngw_zero .LT. ngw) THEN
       !$omp parallel do private(is,ig)
       DO is=1,nstate
          DO ig=ngw_zero+1,ngw
             cx(ig,is) = CMPLX(0._real_8,0._real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE restrain_ngw_zero
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE lag_mult(c0,c1,psi,rhoe,z11,nstate)
    ! ==--------------------------------------------------------------==
    ! Calculates <psi_i| H^0 |psi_j> and diagonalizes it, such that
    ! the |psi_i> array contains the eigenstates of H^0.
    ! The diagonalization is NOT done if TNMR and TLOCALIZE are set.
    ! C1 is overwritten.
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate,1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lag_mult'

    INTEGER                                  :: ierr, info, is, lwork
    REAL(real_8), ALLOCATABLE                :: eigv(:), work(:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    ALLOCATE(eigv(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! calculate the potential and apply it to c0
    CALL forces(c0,c1,tau0,fion,rhoe,psi,nstate,1,.FALSE.,&
         .FALSE.)

    ! calculate the hamiltonian matrix <c0|H|c0> = <c0|c1>:
    CALL ovlap(nstate,z11,c0(:,:,1),c1)
    CALL mp_sum(z11,nstate*nstate,parai%allgrp)


    ! diagonalize the hamiltonian matrix <c0|H|c0> = z11.
    ! IF TNMR, then we need localized (instead of canonical) orbitals.
    IF (((response1%tnmr.OR.response1%tepr).AND.nmr_options%tlocalize).OR. response1%tinteraction) THEN
       CALL dscal(nstate*nstate,-1._real_8,z11,1)
    ELSE
       IF (paral%parent) THEN
          lwork = MAX(1, 3*nstate - 1)
          ALLOCATE(work(lwork),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)! TODO optimal lwork
          IF (cntl%tlsd) THEN
             CALL dsyev('V','U',spin_mod%nsup,z11(1,1),nstate,eigv(1),work,lwork,&
                  info)
             CALL dsyev('V','U',spin_mod%nsdown,z11(spin_mod%nsup+1,spin_mod%nsup+1),nstate,&
                  eigv(spin_mod%nsup+1),work,lwork,info)
          ELSE
             CALL dsyev('V','U',nstate,z11,nstate,eigv,work,lwork,info)
          ENDIF
          DEALLOCATE(work,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          IF (info.NE.0) CALL stopgm('lag_mult','dsyev',& 
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(z11,SIZE(z11),parai%source,parai%allgrp)
       CALL mp_bcast(eigv,SIZE(eigv),parai%source,parai%allgrp)
       ! rotate c0 in canonical form
       CALL rotate(1._real_8,c0(:,:,1),0._real_8,c1,z11,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       CALL dcopy(2*nstate*ncpw%ngw,c1,1,c0,1)
       CALL zeroing(z11)!,nstate*nstate)
       !$omp parallel do private(is)
       DO is=1,nstate
          z11(is,is)=-eigv(is)
       ENDDO
    ENDIF

    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lag_mult
  ! ==================================================================
  FUNCTION energy(a,b,nstate) RESULT(my_res)
    ! INPUT:  A,B: Wavefunction arrays
    ! OUTPUT:  = <A|B>
    ! NB: GLOSUM must be done elsewhere!!
    COMPLEX(real_8)                          :: a(ncpw%ngw,*), b(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: my_res

    INTEGER                                  :: is, isub, ngw1

    CALL tiset('ENERGY       ',isub)
    ngw1=ncpw%ngw
    IF (dmbi%cutoff_restr) ngw1=dmbi%ngw_zero
    my_res = 0._real_8
    DO is=1,nstate
       my_res = my_res + dotp(ngw1,a(:,is),b(:,is))
    ENDDO
    CALL tihalt('ENERGY       ',isub)
  END FUNCTION energy
  ! ==================================================================

END MODULE perturbation_p_utils
