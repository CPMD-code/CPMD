MODULE lr_ortho_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr03
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_ortho
  PUBLIC :: lr_gauge
  PUBLIC :: lr_orthos
  PUBLIC :: give_scr_lr_ortho

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_ortho(nstate,c0,c1)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), &
                                                c1(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_ortho'

    INTEGER                                  :: ierr, isub, llr_ortho
    REAL(real_8), ALLOCATABLE                :: aux(:)

! variables
! ==--------------------------------------------------------------==

    CALL tiset('   LR_ORTHO',isub)
    CALL give_scr_lr_ortho(llr_ortho,nstate)
    ALLOCATE(aux(llr_ortho),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL ovlap(nstate,aux,c0,c1)
    CALL mp_sum(aux,nstate*nstate,parai%allgrp)
    CALL lr_gauge(nstate,aux)
    CALL rotate(-1._real_8,c0,1._real_8,c1,aux,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    IF (geq0) CALL zclean(c1,nstate,ncpw%ngw)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('   LR_ORTHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_ortho
  ! ==================================================================
  SUBROUTINE lr_gauge(nstate,ovlp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ovlp(nstate,nstate)

    REAL(real_8), PARAMETER                  :: w_ij = 0.5_real_8, &
                                                w_ji = 0.5_real_8 

    INTEGER                                  :: i, j
    REAL(real_8)                             :: f_ij, f_ji

! ==--------------------------------------------------------------==

    IF (lr03%tpara_gauge) RETURN
    IF (lr03%tgauge_all.AND.cntl%tlsd) THEN
       DO i=1,spin_mod%nsup-1
          DO j=i+1,spin_mod%nsup
             f_ij=ovlp(i,j)
             f_ji=ovlp(j,i)
             ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
             ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
          ENDDO
       ENDDO
       DO i=spin_mod%nsup+1,nstate-1
          DO j=i+1,nstate
             f_ij=ovlp(i,j)
             f_ji=ovlp(j,i)
             ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
             ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
          ENDDO
       ENDDO
    ELSE IF (lr03%tgauge_all) THEN
       DO i=1,nstate-1
          DO j=i+1,nstate
             f_ij=ovlp(i,j)
             f_ji=ovlp(j,i)
             ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
             ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
          ENDDO
       ENDDO
    ELSE IF (lspin2%tlse) THEN
       DO i=1,nstate-2
          DO j=nstate-1,nstate
             f_ij=ovlp(i,j)
             f_ji=ovlp(j,i)
             ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
             ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
          ENDDO
       ENDDO
       i=nstate-1
       j=nstate
       f_ij=ovlp(i,j)
       f_ji=ovlp(j,i)
       ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
       ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
    ELSE IF (cntl%tlsd) THEN
       DO i=1,spin_mod%nsup-1
          DO j=i+1,spin_mod%nsup
             IF (crge%f(i,1).NE.crge%f(j,1)) THEN
                f_ij=ovlp(i,j)
                f_ji=ovlp(j,i)
                ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
                ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
             ENDIF
          ENDDO
       ENDDO
       DO i=spin_mod%nsup+1,nstate-1
          DO j=i+1,nstate
             IF (crge%f(i,1).NE.crge%f(j,1)) THEN
                f_ij=ovlp(i,j)
                f_ji=ovlp(j,i)
                ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
                ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO i=1,nstate-1
          DO j=i+1,nstate
             IF (crge%f(i,1).NE.crge%f(j,1)) THEN
                f_ij=ovlp(i,j)
                f_ji=ovlp(j,i)
                ovlp(i,j)=w_ij*f_ij + (1.0_real_8-w_ji)*f_ji
                ovlp(j,i)=w_ji*f_ji + (1.0_real_8-w_ij)*f_ij
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_gauge
  ! ==================================================================
  SUBROUTINE lr_orthos(nstate,n1,c0,c1)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, n1
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), &
                                                c1(ncpw%ngw,n1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_orthos'

    INTEGER                                  :: ierr, isub, llr_ortho, ns
    REAL(real_8), ALLOCATABLE                :: aux(:)

! variables
! ==--------------------------------------------------------------==

    CALL tiset('  LR_ORTHOS',isub)
    CALL give_scr_lr_ortho(llr_ortho,n1)
    ALLOCATE(aux(llr_ortho),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tlsd) THEN
       ns=n1/2
       CALL ovlap2(ncpw%ngw,spin_mod%nsup,ns,aux,c0,c1,.TRUE.)
       CALL mp_sum(aux,spin_mod%nsup*ns,parai%allgrp)
       CALL dgemm('N','N',2*ncpw%ngw,ns,spin_mod%nsup,-1._real_8,c0,2*ncpw%ngw,&
            aux,spin_mod%nsup,1._real_8,c1,2*ncpw%ngw)
       CALL ovlap2(ncpw%ngw,spin_mod%nsdown,ns,aux,c0(1,spin_mod%nsup+1),c1(1,ns+1),.TRUE.)
       CALL mp_sum(aux,spin_mod%nsdown*ns,parai%allgrp)
       CALL dgemm('N','N',2*ncpw%ngw,ns,spin_mod%nsdown,-1._real_8,c0(1,spin_mod%nsup+1),2*ncpw%ngw,&
            aux,spin_mod%nsdown,1._real_8,c1(1,ns+1),2*ncpw%ngw)
    ELSE
       CALL ovlap2(ncpw%ngw,nstate,n1,aux,c0,c1,.TRUE.)
       CALL mp_sum(aux,nstate*n1,parai%allgrp)
       CALL dgemm('N','N',2*ncpw%ngw,n1,nstate,-1._real_8,c0,2*ncpw%ngw,&
            aux,nstate,1._real_8,c1,2*ncpw%ngw)
    ENDIF
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (geq0) CALL zclean(c1,n1,ncpw%ngw)
    CALL tihalt('  LR_ORTHOS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_orthos
  ! ==================================================================
  SUBROUTINE give_scr_lr_ortho(llr_ortho,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_ortho, nstate

! ==--------------------------------------------------------------==

    llr_ortho=nstate*nstate
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lr_ortho
  ! ==================================================================


END MODULE lr_ortho_utils
