MODULE lr_xcpot_utils
  USE density_functionals_utils,       ONLY: pade_lda
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: func1,&
                                             func2,&
                                             func3,&
                                             mfxcc_is_pade,&
                                             mfxcc_is_skipped,&
                                             mfxcx_is_skipped
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01,&
                                             lr03,&
                                             lrf1,&
                                             lrf2,&
                                             lrf3
  USE nlcc,                            ONLY: corel,&
                                             roct
  USE parac,                           ONLY: paral
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: cntl,&
                                             fpar
  USE tbxc,                            ONLY: toldcode
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_xcpot

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_xcpot(ddxc,rhoe,swfun)
    ! ==--------------------------------------------------------------==
    ! calculate the 2nd functional derivative of the exchange-
    ! correlation potential used in perturbation theory
    ! ==--------------------------------------------------------------==
    ! ..arguments
    CHARACTER(len=*), PARAMETER              :: procedureN = 'lr_xcpot'
    REAL(real_8)                             :: ddxc(fpar%nnr1,*), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    LOGICAL                                  :: swfun

    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: dens

! ==--------------------------------------------------------------==

    IF (cntl%use_xc_driver .and. lr03%txc_analytic) CALL stopgm(procedureN,'dmu/dn unavailable with the new xc driver; &
                                                   &please use XC_DD_ANALYTIC or numerical derivatives',&
                                                    __LINE__,__FILE__)
    IF (lspin2%tlse.AND.lr03%txc_analytic) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'LR_XCPOT: use XC_DD_ANALYTIC or numeric with LSE'
       CALL stopgm(procedureN,'ROKS not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    IF (swfun) CALL switch_functionals
    IF (lr03%txc_analytic) THEN
       IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
          IF (cntl%tlsd) THEN
             CALL zeroing(ddxc(:,1:3))!,3*nnr1)
          ELSE
             CALL zeroing(ddxc(:,1))!,nnr1)
          ENDIF
       ELSEIF (func1%mfxcc == mfxcc_is_pade.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
          IF (cntl%tlsd) THEN
             CALL stopgm(procedureN,&
                  'dmu/dn for this functional not available',& 
                  __LINE__,__FILE__)
          ELSE
             DO ir=1,fpar%nnr1
                dens=rhoe(ir,1)
                IF (corel%tinlc) dens=dens+roct(ir)
                IF (dens.LT.small)THEN
                   ddxc(ir,1)=0._real_8
                ELSE
                   ddxc(ir,1)=pade_lda(dens,2)
                ENDIF
             ENDDO
          ENDIF
       ELSE
          CALL stopgm(procedureN,&
               'dmu/dn for this functional not available',& 
               __LINE__,__FILE__)
       ENDIF
    ELSEIF (lr01%lopti.EQ.0 .OR. lr01%lopti.EQ.2) THEN
       ! ..this is needed for quadratic line search in cntl%pcg
       IF (cntl%tlsd) THEN
          CALL zeroing(ddxc(:,1:3))!,3*nnr1)
       ELSE
          DO ir=1,fpar%nnr1
             dens=rhoe(ir,1)
             IF (corel%tinlc) dens=dens+roct(ir)
             IF (dens.LT.small)THEN
                ddxc(ir,1)=0._real_8
             ELSE
                ddxc(ir,1)=pade_lda(dens,2)
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    IF (swfun) CALL switch_functionals
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_xcpot
  ! ==================================================================
END MODULE lr_xcpot_utils
