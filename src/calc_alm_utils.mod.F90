MODULE calc_alm_utils
  USE cppt,                            ONLY: twnl
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: afnl,&
                                             alm,&
                                             biln,&
                                             cnl,&
                                             nlptr
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean_k
!!use ovlap_utils, only : ovlap_h

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_alm
  PUBLIC :: calc_k_alm
  PUBLIC :: give_scr_calc_alm

CONTAINS

  ! ==================================================================
  SUBROUTINE calc_alm
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATES OVERLAP MATRIX BETWEEN NON-LOCAL PSEUDOPOTENTIAL ==
    ! ==  PROJECTORS                                                  ==
    ! ==  K POINTS VERSION (ALM BECOMES COMPLEX -- IMAGP=2)           ==
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_alm'

    COMPLEX(real_8)                          :: cp1
    INTEGER :: i1, i2, ia, ia1, ierr, ig, ik, ikpt, is, is1, isa, isa1, isub, &
      iv, kbeg, kend, kinc, ncnl, nkpoint, nli, nlm1, nlmx2, nni, nnl
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: tlsdtrue
    REAL(real_8)                             :: dummy

! dummy parameter for calling _swap functions
! ==--------------------------------------------------------------==

    IF (ifirst.EQ.0) THEN
       ! Allocate Array ALM and CNL
       ifirst=1
       ! NLM is initialized in SETSYS. 
       nlmx2=MAX(1,imagp*nlm*nlm*kpts_com%nkptall)
       ncnl=MAX(1,2*nkpt%ngwk*nlm)
       nnl=MAX(1,3*nlm)
       nni=MAX(1,crge%n*nlm*imagp)
       ALLOCATE(alm(imagp,nlm,nlm,kpts_com%nkptall),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cnl(nkpt%ngwk,nlm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(nlptr(3,nlm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(afnl(imagp,nlm,nni/(imagp*nlm)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(biln(imagp,nlm,nni/(imagp*nlm)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%parent.AND.nlm.GT.0) CALL prmem('  CALC_ALM')
       nli=0
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             DO iv=1,nlps_com%ngh(is)
                nli=nli+1
                nlptr(1,nli)=is
                nlptr(2,nli)=isa
                nlptr(3,nli)=iv
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  CALC_ALM',isub)
    ! ==--------------------------------------------------------------==
    ! Calculate overlap
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       ! TODO ask if swap is still needed
       IF (tkpts%tkblock) THEN
          CALL rkpt_swap(dummy,1,ikpt,'EIGKR TWNL')
       ENDIF
       DO ik=1,nkpoint
          isa1=0
          nlm1=0
          DO is1=1,ions1%nsp
             DO ia1=1,ions0%na(is1)
                isa1=isa1+1
                DO iv=1,nlps_com%ngh(is1)
                   cp1=(0.0_real_8,-1.0_real_8)**nghtol(iv,is1)
                   nlm1=nlm1+1
                   ! CNL is temporarily used as scr space
                   DO ig=1,nkpt%ngwk
                      cnl(ig,nlm1)=cp1*twnl(ig,iv,is1,ik)*eigkr(ig,isa1,ik)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          ! Check (always true except when debugging)
          IF (nlm1.NE.nlm) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'NLM1=',nlm1,'NLM=',nlm
             CALL stopgm('CALC_ALM','WRONG NON-LOCAL NUMBER',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ..aa
          IF (cntl%tlsd) THEN
             tlsdtrue=.TRUE.
             cntl%tlsd=.FALSE.
          ELSE
             tlsdtrue=.FALSE.
          ENDIF
          IF (tkpts%tkpnt) THEN
             IF (geq0) CALL zclean_k(cnl,nlm,ncpw%ngw)
             CALL ovlap_h(nlm,alm(1,1,1,ik),cnl,cnl)
          ELSE
             CALL ovlap  (nlm,alm(1,1,1,ik),cnl,cnl)
          ENDIF
          ! ..aa
          IF (tlsdtrue) cntl%tlsd=.TRUE.
          ! ALM is symmetric or hermitian.
          DO i1=1,nlm
             DO i2=i1+1,nlm
                alm(1,i1,i2,ik)=alm(1,i2,i1,ik)
             ENDDO
          ENDDO
          IF (tkpts%tkpnt) THEN
             DO i1=1,nlm
                DO i2=i1+1,nlm
                   alm(2,i1,i2,ik)=-alm(2,i2,i1,ik)
                ENDDO
             ENDDO
          ELSE
          ENDIF
          IF (tkpts%tkpnt) THEN
             CALL sumhmat(alm(1,1,1,ik),nlm)
          ELSE
             CALL summat (alm(:,:,:,ik),nlm)
          ENDIF
       ENDDO
       IF (tkpts%tkblock) CALL wkpt_swap(dummy,1,ikpt,'ALM')
    ENDDO
    CALL tihalt('  CALC_ALM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_alm
  ! ==================================================================
  SUBROUTINE calc_k_alm(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE ALM FOR ONE SET OF KPOINTS                         ==
    ! == USED WITH TKBLOCK                                            ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ikpt

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_k_alm'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: cp1
    INTEGER                                  :: i1, i2, ia1, ierr, ig, ik, &
                                                is1, isa1, isub, iv, lscralm, &
                                                nlm1
    LOGICAL                                  :: tlsdtrue
    REAL(real_8), ALLOCATABLE                :: scralm(:)

    CALL tiset('CALC_K_ALM',isub)
    CALL give_scr_calc_alm(lscralm,tag)
    ALLOCATE(scralm(lscralm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DO ik=1,nkpbl(ikpt)
       isa1=0
       nlm1=0
       DO is1=1,ions1%nsp
          DO ia1=1,ions0%na(is1)
             isa1=isa1+1
             DO iv=1,nlps_com%ngh(is1)
                cp1=(0.0_real_8,-1.0_real_8)**nghtol(iv,is1)
                nlm1=nlm1+1
                ! CNL is temporarily used as scr space
                DO ig=1,nkpt%ngwk
                   cnl(ig,nlm1)=cp1*twnl(ig,iv,is1,ik)*eigkr(ig,isa1,ik)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ! Check (always true except when debugging)
       IF (nlm1.NE.nlm) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'NLM1=',nlm1,'NLM=',nlm
          CALL stopgm('CALC_ALM','WRONG NON-LOCAL NUMBER',& 
               __LINE__,__FILE__)
       ENDIF
       ! ..     aa
       IF (cntl%tlsd) THEN
          tlsdtrue=.TRUE.
          cntl%tlsd=.FALSE.
       ELSE
          tlsdtrue=.FALSE.
       ENDIF
       IF (tkpts%tkpnt) THEN
          IF (geq0) CALL zclean_k(cnl,nlm,ncpw%ngw)
          CALL ovlap_h(nlm,alm(1,1,1,ik),cnl,cnl)
       ELSE
          CALL ovlap  (nlm,alm(1,1,1,ik),cnl,cnl)
       ENDIF
       ! ..     aa
       IF (tlsdtrue) cntl%tlsd=.TRUE.
       ! ALM is symmetric or hermitian.
       DO i1=1,nlm
          DO i2=i1+1,nlm
             alm(1,i1,i2,ik)=alm(1,i2,i1,ik)
          ENDDO
       ENDDO
       IF (tkpts%tkpnt) THEN
          DO i1=1,nlm
             DO i2=i1+1,nlm
                alm(2,i1,i2,ik)=-alm(2,i2,i1,ik)
             ENDDO
          ENDDO
       ELSE
       ENDIF
       IF (tkpts%tkpnt) THEN
          CALL sumhmat(alm(1,1,1,ik),nlm)
       ELSE
          CALL summat (alm(:,:,:,ik),nlm)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scralm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('CALC_K_ALM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_k_alm
  ! ==================================================================
  SUBROUTINE give_scr_calc_alm(lcalc_alm,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lcalc_alm
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lcalc_alm=imagp*nlm*(nlm+1)/2
    tag='IMAGP*NLM*(NLM+1)/2'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_calc_alm
  ! ==================================================================

END MODULE calc_alm_utils
