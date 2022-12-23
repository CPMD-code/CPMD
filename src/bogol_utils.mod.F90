MODULE bogol_utils
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE tauf,                            ONLY: itaur

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bogol
  PUBLIC :: give_scr_bogol

CONTAINS

  ! ================================================================
  SUBROUTINE bogol(c0,c2,sc0,we,focc,&
       rhoe,psi,nstate,bogo)
    ! ==------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,*), sc0(*)
    REAL(real_8)                             :: rhoe(fpar%nnr1,*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: focc(nstate,*), we(nstate,*)
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,*)
    REAL(real_8)                             :: bogo

    CHARACTER(*), PARAMETER                  :: procedureN = 'bogol'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, ib, ierr, ikind, ikk, &
                                                ikpt, kbeg, kend, kinc, &
                                                lbogol, nkpoint
    REAL(real_8), ALLOCATABLE                :: scr(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL give_scr_bogol(lbogol,tag,nstate)
    ALLOCATE(scr(lbogol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    bogo = 0._real_8
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock)CALL rkpt_swap(c0,nstate,ikpt,&
            'HGKP HGKM MASKGW EIGKR TWNL C0')
       DO ikind=1,nkpoint
          ikk=kpbeg(ikpt)+ikind
          IF (cntl%tlsd) THEN
             ! .. alpha spin
             CALL hpsi(c0(:,1:spin_mod%nsup,ikind),c2,sc0,&
                  rhoe(1,1),psi,spin_mod%nsup,ikind,1)
             DO i=1,spin_mod%nsup
                IF (tkpts%tkpnt) THEN
                   scr(i)=-ddot(2*nkpt%ngwk,c0(1,i,ikind),1,c2(1,i),1)
                ELSE
                   scr(i)=-dotp(ncpw%ngw,c0(:,i,ikind),c2(:,i))
                ENDIF
             ENDDO
             CALL mp_sum(scr,spin_mod%nsup,parai%allgrp)
             DO i=1,spin_mod%nsup
                bogo=bogo+wk(ikk)*focc(i,ikk)*(scr(i)-we(i,ikk))
             ENDDO
             ! .. beta spin
             ! use correct tau-potential
             itaur=2
             ib=spin_mod%nsup
             CALL hpsi(c0(:,ib+1:ib+spin_mod%nsdown,ikind),c2(1,ib+1),sc0,&
                  rhoe(1,2),psi,spin_mod%nsdown,ikind,1)
             itaur=1
             DO i=1,spin_mod%nsdown
                IF (tkpts%tkpnt) THEN
                   scr(i)=-ddot(2*nkpt%ngwk,c0(1,ib+i,ikind),1,c2(1,ib+i),1)
                ELSE
                   scr(i)=-dotp(ncpw%ngw,c0(:,ib+i,ikind),c2(:,ib+i))
                ENDIF
             ENDDO
             CALL mp_sum(scr,spin_mod%nsdown,parai%allgrp)
             DO i=1,spin_mod%nsdown
                bogo=bogo+wk(ikk)*focc(ib+i,ikk)*(scr(i)-we(ib+i,ikk))
             ENDDO
          ELSE
             CALL hpsi(c0(:,:,ikind),c2,sc0,&
                  rhoe,psi,nstate,ikind,clsd%nlsd)
             DO i=1,nstate
                IF (tkpts%tkpnt) THEN
                   scr(i)=-ddot(2*nkpt%ngwk,c0(1,i,ikind),1,c2(1,i),1)
                ELSE
                   scr(i)=-dotp(ncpw%ngw,c0(:,i,ikind),c2(:,i))
                ENDIF
             ENDDO
             CALL mp_sum(scr,nstate,parai%allgrp)
             DO i=1,nstate
                bogo=bogo+wk(ikk)*focc(i,ikk)*(scr(i)-we(i,ikk))
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! ==------------------------------------------------------------==
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    RETURN
  END SUBROUTINE bogol
  ! ================================================================
  SUBROUTINE give_scr_bogol(lbogol,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lbogol
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lhpsi

    CALL give_scr_hpsi(lhpsi,tag,nstate)
    lbogol=MAX(nstate,lhpsi)
    tag='MAX(NSTATE,'//tag(1:19)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_bogol
  ! ==================================================================

END MODULE bogol_utils
