MODULE lr_force_utils
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE hfx_drivers,                     ONLY: hfxpsi,&
                                             hfxrpa
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE poin,                            ONLY: potr
  USE rho1ofr_utils,                   ONLY: rho1ofr
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: clsd,&
                                             spin_mod,&
                                             tdsp1
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho1_utils,                  ONLY: give_scr_v1ofrho1,&
                                             give_scr_vtdofrho1,&
                                             v1ofrho1,&
                                             vtdofrho1
  USE vpsi_utils,                      ONLY: vpsimt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_force
  PUBLIC :: give_scr_lr_force
  PUBLIC :: lr_force_sub
  PUBLIC :: lru_force
  PUBLIC :: give_scr_lru_force

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_force(c0,c1,c2,sc0,e2,ddxc,psi,drhoe,eigv,&
       nstate,TYPE,orbital)
    ! ==--------------------------------------------------------------==
    ! computes E(2) the second order functional and the corresponding
    ! electronic forces 
    COMPLEX(real_8)                          :: c0(:,:), c1(:,:), c2(:,:), &
                                                sc0(*)
    REAL(real_8)                             :: e2(:), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: drhoe(:,:), eigv(:)
    INTEGER                                  :: nstate
    CHARACTER(len=*)                         :: TYPE, orbital

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_force'

    INTEGER                                  :: i, ierr, isub, k
    LOGICAL                                  :: tcis
    REAL(real_8), ALLOCATABLE                :: focc(:)

    CALL tiset('  LR_FORCE',isub)
    ! ==--------------------------------------------------------------==
    ! == C2 = f*(H-e)*C1  for canonical orbitals                      ==
    ! == C2 = f*(H-F)*C1  else                                        ==
    ! ==--------------------------------------------------------------==
    CALL hpsi(c1(:,1:nstate),c2,sc0,potr,psi(:,1),nstate,1,clsd%nlsd)
    IF (cntl%tlsd) THEN
       CALL hfxpsi(c0(:,1:tdsp1%nupel),c1(:,1:tdsp1%nupel),c2(:,1:tdsp1%nupel),&
            crge%f(1:tdsp1%nupel,1),-1._real_8,psi(:,1),tdsp1%nupel,&
            tdsp1%nupel)
       k=tdsp1%nupel+1
       CALL hfxpsi(c0(:,k:k+tdsp1%ndoel-1),c1(:,k:k+tdsp1%ndoel-1),c2(:,k:k+tdsp1%ndoel-1), &
            & crge%f(k:k+tdsp1%ndoel-1,1),-1._real_8,psi(:,1),tdsp1%ndoel,tdsp1%ndoel)
    ELSE
       CALL hfxpsi(c0(:,1:nstate),c1,c2,crge%f(:,1),-1._real_8,psi(:,1),nstate,nstate)
    ENDIF
    IF (INDEX(orbital,"CANON").NE.0) THEN
       DO i=1,nstate
          CALL daxpy(2*ncpw%ngw,eigv(i),c1(1,i),1,c2(1,i),1)
       ENDDO
    ELSE
       CALL rotate(-1._real_8,c1,1._real_8,c2,eigv,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ENDIF
    ! ==--  Phonon Code  ---------------------------------------------==
    ! ==--  Orbital Hardness -----------------------------------------==
    IF (INDEX(TYPE,"PHONON").NE.0 .OR.INDEX(TYPE,"ORBHARD").NE.0) THEN
       DO i=1,nstate
          CALL dscal(2*ncpw%ngw,crge%f(i,1),c2(1,i),1)
          e2(3)=e2(3)-dotp(ncpw%ngw,c2(:,i),c1(:,i))
       ENDDO
       CALL rho1ofr(c0,c1,crge%f(:,1),drhoe,psi(:,1),nstate)
       CALL v1ofrho1(e2,drhoe,ddxc,psi)
       CALL vpsimt(c0,c2,crge%f(:,1),drhoe,psi(:,1),nstate,clsd%nlsd,.TRUE.)
       ! ==--  Time Dependent DFT ---------------------------------------==
    ELSEIF (INDEX(TYPE,"TDAMB").NE.0) THEN
       ! Nothing more to do
    ELSEIF((INDEX(TYPE,"TDA").NE.0)&
         .OR.(INDEX(TYPE,"TDAPB").NE.0) ) THEN

       ALLOCATE(focc(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(I)
       DO i=1,nstate
          focc(i)=1._real_8
       ENDDO
       CALL rho1ofr(c0,c1,focc,drhoe,psi(:,1),nstate)
       CALL vtdofrho1(e2,drhoe,ddxc,psi,.TRUE.)
       CALL vpsimt(c0,c2,focc,drhoe,psi(:,1),nstate,clsd%nlsd,.TRUE.)
       tcis=INDEX(TYPE,"TDAPB").EQ.0
       CALL switch_functionals
       CALL hfxrpa(c0,c1,c2,psi(:,1),nstate,tcis)
       CALL switch_functionals
       DEALLOCATE(focc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSEIF (INDEX(TYPE,"ZFORCE").NE.0) THEN
       DO i=1,nstate
          CALL dscal(2*ncpw%ngw,crge%f(i,1),c2(1,i),1)
          e2(3)=e2(3)-dotp(ncpw%ngw,c2(:,i),c1(:,i))
       ENDDO
       CALL rho1ofr(c0,c1,crge%f(:,1),drhoe,psi(:,1),nstate)
       CALL v1ofrho1(e2,drhoe,ddxc,psi)
       CALL vpsimt(c0,c2,crge%f(:,1),drhoe,psi(:,1),nstate,clsd%nlsd,.TRUE.)
       ! ==--------------------------------------------------------------==
    ELSE
       CALL stopgm('LR_FORCE','NO OPTION DEFINED.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('  LR_FORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_force
  ! ==================================================================
  SUBROUTINE give_scr_lr_force(llr_force,TYPE,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_force
    CHARACTER(len=*)                         :: TYPE
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lfff, lhpsi, lv1ofrho1

    CALL give_scr_hpsi(lhpsi,tag,crge%n)
    IF (INDEX(TYPE,"PHONON").NE.0 .OR.INDEX(TYPE,"ORBHARD").NE.0) THEN
       CALL give_scr_v1ofrho1(lv1ofrho1,tag)
       llr_force=crge%n*crge%n+100
       lfff=0
    ELSEIF (INDEX(TYPE,"cntl%tddft").NE.0) THEN
       CALL give_scr_vtdofrho1(lv1ofrho1,tag)
       llr_force=2*ncpw%ngw*crge%n+lhpsi
       lfff=crge%n+100
    ELSEIF(INDEX(TYPE,"TDA").NE.0 .OR. INDEX(TYPE,"TDAPB").NE.0 .OR.&
         INDEX(TYPE,"TDAMB").NE.0) THEN
       CALL give_scr_vtdofrho1(lv1ofrho1,tag)
       llr_force=0
       lfff=crge%n+100
    ELSEIF (INDEX(TYPE,"ZFORCE").NE.0) THEN
       CALL give_scr_v1ofrho1(lv1ofrho1,tag)
       llr_force=0
       lfff=crge%n+100
    ENDIF
    llr_force=MAX(llr_force,lhpsi,lv1ofrho1+lfff)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lr_force
  ! ==================================================================
  SUBROUTINE lr_force_sub(c0,c1,c2,sc0,e2,ddxc,psi,drhoe,eigv,&
       nstate,msub,TYPE,orbital)
    ! ==--------------------------------------------------------------==
    ! C0 only span the relevant subspace!
    ! EIGV is also only within the subspace
    COMPLEX(real_8)                          :: sc0(*)
    REAL(real_8)                             :: e2(3), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: drhoe(:,:), eigv(*)
    INTEGER                                  :: nstate, msub
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,msub), &
                                                c1(nkpt%ngwk,msub), &
                                                c0(nkpt%ngwk,msub)
    CHARACTER(len=*)                         :: TYPE, orbital

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_force_sub'

    INTEGER                                  :: i, ierr, isub
    REAL(real_8), ALLOCATABLE                :: focc(:)

! VARIABLES
! ==--------------------------------------------------------------==

    CALL tiset('  LR_FORCE',isub)
    IF (cntl%tlsd) CALL stopgm('LR_FORCE_SUB','LSD NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == C2 = f*(H-e)*C1  for canonical orbitals                      ==
    ! == C2 = f*(H-F)*C1  else                                        ==
    ! ==--------------------------------------------------------------==
    CALL hpsi(c1(:,1:msub),c2,sc0,potr,psi(:,1),msub,1,clsd%nlsd)
    IF (INDEX(orbital,"CANON").NE.0) THEN
       DO i=1,msub
          CALL daxpy(2*ncpw%ngw,eigv(i),c1(1,i),1,c2(1,i),1)
       ENDDO
    ELSE
       CALL rotate(-1._real_8,c1,1._real_8,c2,eigv,msub,2*ncpw%ngw,cntl%tlsd,msub,msub)
    ENDIF
    IF (INDEX(TYPE,"TDA").NE.0) THEN
       ALLOCATE(focc(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(I)
       DO i=1,nstate
          focc(i)=1._real_8
       ENDDO
       ! WARNING
       ! 
       IF (cntl%tlsd) CALL stopgm("LR_FORCE_SUB","CORRECT THIS BUG",& 
            __LINE__,__FILE__)
       ! 
       ! RHO1OFR MAKES AN ASSUMPTION ON THE STATES IN LSD THAT
       ! IS NOT CORRECT - I THINK WE SHOULD USE RHOABOFR HERE
       ! 
       ! WARNING
       CALL rho1ofr(c0,c1,focc,drhoe,psi(:,1),msub)
       CALL vtdofrho1(e2,drhoe,ddxc,psi,.TRUE.)
       CALL vpsimt(c0,c2,focc,drhoe,psi(:,1),msub,clsd%nlsd,.TRUE.)
       DEALLOCATE(focc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       CALL stopgm('LR_FORCE','NO OPTION DEFINED.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('  LR_FORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_force_sub
  ! ==================================================================
  SUBROUTINE lru_force(c0,c1,c2,cgrad,psi,drhoe,ufmat,&
       nstate,msub)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), psi(:)
    REAL(real_8)                             :: drhoe(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: cgrad(nstate,*)
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate), &
                                                c1(nkpt%ngwk,nstate)
    INTEGER                                  :: msub
    REAL(real_8)                             :: ufmat(msub,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lru_force'

    INTEGER                                  :: i, ierr, isub
    REAL(real_8), ALLOCATABLE                :: focc(:), xmat(:)

! VARIABLES
! ==--------------------------------------------------------------==

    CALL tiset(' LRU_FORCE',isub)
    IF (cntl%tlsd) CALL stopgm('LRU_FORCE','LSD NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ALLOCATE(focc(msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(xmat(msub * msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(I)
    DO i=1,msub
       focc(i)=1._real_8
    ENDDO
    CALL vpsimt(c1,c2,focc,drhoe,psi,msub,clsd%nlsd,.TRUE.)
    CALL ovlap2(ncpw%ngw,nstate,msub,cgrad,c0,c2,.TRUE.)
    CALL ovlap(msub,xmat,c1,c1)
    CALL dgemm('T','N',nstate,msub,msub,-1._real_8,ufmat,msub,xmat,msub,&
         1.0_real_8,cgrad,nstate)
    DEALLOCATE(focc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(xmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL mp_sum(cgrad,nstate*msub,parai%allgrp)
    CALL dscal(nstate*msub,2._real_8,cgrad,1)
    ! ==--------------------------------------------------------------==
    CALL tihalt(' LRU_FORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lru_force
  ! ==================================================================
  SUBROUTINE give_scr_lru_force(llr_force,msub,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_force, msub
    CHARACTER(len=30)                        :: tag

    llr_force=msub+msub*msub+100
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lru_force
  ! ==================================================================

END MODULE lr_force_utils
