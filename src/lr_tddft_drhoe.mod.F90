MODULE lr_tddft_drhoe
  USE cppt,                            ONLY: nzh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: nolr,&
                                             td01,&
                                             td03
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: ovlap,&
                                             ovlap_add
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: ipcurr
  USE poin,                            ONLY: rhoo
  USE readsr_utils,                    ONLY: xstring
  USE rho1ofr_utils,                   ONLY: rho1ofr,&
                                             rhoabofr
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: relax_rho
  PUBLIC :: print_ex_rho

CONTAINS

  SUBROUTINE relax_rho(c0,c1,c2,sc0,rhoo,rhoe,psi,nstate,nroot)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: sc0(*)
    REAL(real_8)                             :: rhoo(*), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: nroot
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nolr,nroot,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'relax_rho'

    INTEGER                                  :: i, ierr
    REAL(real_8)                             :: rsum0, rsum1
    REAL(real_8), ALLOCATABLE                :: r1mat(:)

    rsum0=0._real_8
    !$omp parallel do private(I) shared(fpar,RHOO) &
    !$omp  reduction(+:RSUM0)
    DO i=1,fpar%nnr1
       rsum0=rsum0+rhoo(i)
    ENDDO
    rsum0=rsum0*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ! 
    ALLOCATE(r1mat(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Z*C0
    CALL rho1ofr(c0,c2,crge%f(:,1),rhoe,psi,nstate)
    ! C0*C0
    CALL daxpy(fpar%nnr1*clsd%nlsd,1._real_8,rhoo,1,rhoe,1)
    ! R1mat=x*x
    CALL ovlap(nstate,r1mat,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1))
    IF (.NOT.td03%tda) THEN
       CALL ovlap_add(nstate,r1mat,c1(1,1,td01%fstate,2),c1(1,1,td01%fstate,2))
    ENDIF
    CALL mp_sum(r1mat,nstate*nstate,parai%allgrp)
    ! c*R1mat
    CALL rotate(-1._real_8,c0,0._real_8,c2,r1mat,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    DEALLOCATE(r1mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) THEN
       CALL daxpy(fpar%nnr1,-1._real_8,rhoe(1,2),1,rhoe(1,1),1)
       ! X*X alpha spin
       CALL rhoabofr(spin_mod%nsup,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1),rhoe(:,1),psi)
       IF (.NOT.td03%tda) THEN
          CALL rhoabofr(spin_mod%nsup,c1(:,:,td01%fstate,2),&
               c1(:,:,td01%fstate,2),rhoe(:,1),psi)
       ENDIF
       ! C-*C alpha spin
       CALL rhoabofr(spin_mod%nsup,c2(:,:),c0(:,:),rhoe(:,1),psi)
       ! X*X beta spin
       CALL rhoabofr(spin_mod%nsdown,c1(:,spin_mod%nsup+1:,td01%fstate,1),&
            c1(:,spin_mod%nsup+1:,td01%fstate,1),rhoe(:,2),psi)
       IF (.NOT.td03%tda) THEN
          CALL rhoabofr(spin_mod%nsdown,c1(:,spin_mod%nsup+1:,td01%fstate,2),&
               c1(:,spin_mod%nsup+1:,td01%fstate,2),rhoe(:,2),psi)
       ENDIF
       ! C-*C beta spin
       CALL rhoabofr(spin_mod%nsdown,c2(:,spin_mod%nsup+1:),c0(:,spin_mod%nsup+1:),rhoe(:,2),psi)
       CALL daxpy(fpar%nnr1,1._real_8,rhoe(1,2),1,rhoe(1,1),1)
    ELSE
       ! X*X
       CALL rhoabofr(nstate,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1),rhoe(:,1),psi)
       IF (.NOT.td03%tda) THEN
          CALL rhoabofr(nstate,c1(:,:,td01%fstate,2),&
               c1(:,:,td01%fstate,2),rhoe(:,1),psi)
       ENDIF
       ! C-*C
       CALL rhoabofr(nstate,c2(:,:),c0(:,:),rhoe(:,1),psi)
    ENDIF
    ! Check number of electrons
    rsum1=0._real_8
    !$omp parallel do private(I) shared(fpar,RHOE) &
    !$omp  reduction(+:RSUM1)
    DO i=1,fpar%nnr1
       rsum1=rsum1+rhoe(i,1)
    ENDDO
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    CALL mp_sum(rsum0,parai%allgrp)
    CALL mp_sum(rsum1,parai%allgrp)
    IF (paral%parent) THEN
       IF (ABS(rsum0-rsum1).GT.1.e-10_real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               " WARNING : INTEGRAL OVER RELAXED DENSITY IS INCORRECT "
          IF (paral%io_parent)&
               WRITE(6,'(A,T46,F20.10)')&
               " WARNING : GROUND STATE DENSITY",rsum0
          IF (paral%io_parent)&
               WRITE(6,'(A,T46,F20.10)')&
               " WARNING : RELAXED DENSITY",rsum1
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE relax_rho
  ! ==================================================================
  SUBROUTINE print_ex_rho(rhoe,psi,tau0)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'print_ex_rho'

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, i1, i2, ierr, n1, n2

    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    !$omp parallel do private(I)
    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO i=1,ncpw%nhg
       vtemp(i) = psi(nzh(i))
    ENDDO
    IF (cntl%tpath) THEN
       IF (ipcurr.EQ.0) THEN
          filen='DENSITY_ES'
       ELSE
          cflbod='DENSITY_ES_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
       ENDIF
    ELSE
       filen='DENSITY_ES'
    ENDIF
    ! 
    CALL densto(vtemp,tau0,filen)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE print_ex_rho
  ! ==========================================================
END MODULE lr_tddft_drhoe
