MODULE canon_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: canon
  PUBLIC :: give_scr_canon

CONTAINS

  ! ==================================================================
  SUBROUTINE canon(c0,c2,f,nstate,eigv)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: eigv(nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'canon'

    INTEGER                                  :: ierr, info, is, naux
    LOGICAL                                  :: debug
    REAL(real_8), ALLOCATABLE                :: aux(:), hmat(:,:)

    debug=.FALSE.
    IF (lspin2%tlse) CALL stopgm('CANON','TLSE NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    naux=MAX(nstate*nstate,25*nstate)
    ALLOCATE(hmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(naux),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL ovlap(nstate,hmat,c0,c2)
    CALL summat(hmat,nstate)
    CALL dscal(nstate*nstate,-1._real_8,hmat,1)
    IF (paral%parent) THEN
       IF (cntl%tlsd) THEN
          info=0
          CALL dsyev('V','U',spin_mod%nsup,hmat(1,1),nstate,eigv(1),aux,naux,&
               info)
          IF (info.NE.0) CALL stopgm('CANON','DSYEV:UP:INFO',& 
               __LINE__,__FILE__)
          info=0
          CALL dsyev('V','U',spin_mod%nsdown,hmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate,eigv(&
               spin_mod%nsup+1),aux,naux,info)
          IF (info.NE.0) CALL stopgm('CANON','DSYEV:DOWN:INFO',& 
               __LINE__,__FILE__)
       ELSE
          info=0
          CALL dsyev('V','U',nstate,hmat,nstate,eigv,aux,naux,info)
          IF (info.NE.0) CALL stopgm('CANON','DSYEV:INFO',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    CALL mp_bcast(hmat,SIZE(hmat),parai%source,parai%allgrp)
    CALL mp_bcast(eigv,SIZE(eigv),parai%source,parai%allgrp)
    DO is=1,nstate
       IF (f(is).GT.1.e-6_real_8) eigv(is)=eigv(is)/f(is)
    ENDDO
    CALL rotate(1._real_8,c0,0._real_8,c2,hmat,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,c2,1,c0,1)
    DEALLOCATE(hmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE canon
  ! ==================================================================
  SUBROUTINE give_scr_canon(lcanon,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcanon
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lcanon=0
    lcanon=nstate*nstate+MAX(nstate*nstate,25*nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_canon
  ! ==================================================================


END MODULE canon_utils
