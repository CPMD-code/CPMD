MODULE rnlsm2_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE mm_dimmod,                       ONLY: mmdim
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: dfnl,&
                                             eigr
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm2
  PUBLIC :: give_scr_rnlsm2

CONTAINS

#if defined(__SR8000) || defined(_vpp_) || defined (__ES) 

  ! ==================================================================
  SUBROUTINE RNLSM2(C0,NSTATE,IKPT,IKIND)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DFNL WHICH IS USED IN THE SUBROUTINE RNLFOR       ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==          NOT IMPLEMENTED FOR TSHEL(IS)=TRUE                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: NSTATE
    COMPLEX(real_8)                          :: C0(nkpt%ngwk,NSTATE)
    INTEGER                                  :: IKPT, IKIND

    CHARACTER(*), PARAMETER                  :: procedureN = 'RNLSM2'
    COMPLEX(real_8), PARAMETER               :: ZONE = (1._real_8,0._real_8) ,&
                                                ZZERO = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: CI, ZFAC
    COMPLEX(real_8), ALLOCATABLE             :: EISCR(:,:,:)
    INTEGER                                  :: I, IA, IG, II, IS, ISA, ISA0, &
                                                ISUB, IV
    LOGICAL                                  :: MASK
    REAL(real_8)                             :: ARG1, ARG2, ARG3, CII, CIR, &
                                                EI, ER, TFAC
    REAL(real_8), ALLOCATABLE                :: DAI(:,:,:,:)

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF(NLM.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL TISET('    RNLSM2',ISUB)
    ALLOCATE(eiscr(nkpt%ngwk, 3, mmdim%naxq),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dai(imagp, 3, mmdim%naxq, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(DAI)!,IMAGP*3*mmdim%naxq*NSTATE)
    IF(tkpts%tkpnt) THEN
       ZFAC=parm%tpiba*ZONE
    ELSE
       TFAC=2._real_8*parm%tpiba
    ENDIF
    MASK=IMAGP.EQ.2
    ISA0=0
    DO IS=1,ions1%nsp
       DO IV=1,nlps_com%ngh(IS)
          CI=(0.0_real_8,-1.0_real_8)**(NGHTOL(IV,IS)+1)
          CIR=REAL(CI)
          CII=AIMAG(CI)
          IF(tkpts%tkpnt) THEN
#ifdef __SR8000
             !poption parallel, tlocal(ISA,ARG1,ARG2,ARG3,ER,EI)
             !voption indep(EISCR)
#endif
             !$omp parallel do private(IA,IG,ISA,ARG1,ARG2,ARG3,ER,EI)
             DO IA=1,ions0%na(IS)
                ISA=ISA0+IA
                ! Make use of the special structure of CI
                IF(ABS(CIR).GT.0.5_real_8) THEN
                   ! CI is real
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   DO IG=1,NGW
                      ARG1=(GK(1,IG)+&
                           &                 RK(1,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR
                      ARG2=(GK(2,IG)+&
                           &                 RK(2,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR
                      ARG3=(GK(3,IG)+&
                           &                 RK(3,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR

                      ER=REAL(EIGKR(IG,ISA,IKIND))
                      EI=AIMAG(EIGKR(IG,ISA,IKIND))
                      EISCR(IG,1,IA) = CMPLX(ARG1*ER,ARG1*EI,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(ARG2*ER,ARG2*EI,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(ARG3*ER,ARG3*EI,kind=real_8)
                   ENDDO
#ifdef __SR8000
                   !poption parallel
#endif
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   DO IG=NGW+1,nkpt%ngwk
                      ARG1=(-GK(1,IG-NGW)+&
                           &                 RK(1,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR         
                      ARG2=(-GK(2,IG-NGW)+&
                           &                 RK(2,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR        
                      ARG3=(-GK(3,IG-NGW)+&
                           &                 RK(3,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR

                      ER=REAL(EIGKR(IG,ISA,IKIND))
                      EI=AIMAG(EIGKR(IG,ISA,IKIND))
                      EISCR(IG,1,IA) = CMPLX(ARG1*ER,ARG1*EI,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(ARG2*ER,ARG2*EI,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(ARG3*ER,ARG3*EI,kind=real_8)
                   ENDDO
                ELSE
                   ! CI is imaginary
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   DO IG=1,NGW
                      ARG1=(GK(1,IG)+&
                           &                 RK(1,IKIND))*TWNL(IG,IV,IS,IKIND)*CII           
                      ARG2=(GK(2,IG)+&
                           &                 RK(2,IKIND))*TWNL(IG,IV,IS,IKIND)*CII          
                      ARG3=(GK(3,IG)+&
                           &                 RK(3,IKIND))*TWNL(IG,IV,IS,IKIND)*CII

                      ER=REAL(EIGKR(IG,ISA,IKIND))
                      EI=AIMAG(EIGKR(IG,ISA,IKIND))
                      EISCR(IG,1,IA) = CMPLX(-ARG1*EI,ARG1*ER,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(-ARG2*EI,ARG2*ER,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(-ARG3*EI,ARG3*ER,kind=real_8)
                   ENDDO
#ifdef __SR8000
                   !poption parallel
#endif
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   DO IG=NGW+1,nkpt%ngwk
                      ARG1=(-GK(1,IG-NGW)+&
                           &                 RK(1,IKIND))*TWNL(IG,IV,IS,IKIND)*CII            
                      ARG2=(-GK(2,IG-NGW)+&
                           &                 RK(2,IKIND))*TWNL(IG,IV,IS,IKIND)*CII           
                      ARG3=(-GK(3,IG-NGW)+&
                           &                 RK(3,IKIND))*TWNL(IG,IV,IS,IKIND)*CII

                      ER=REAL(EIGKR(IG,ISA,IKIND))
                      EI=AIMAG(EIGKR(IG,ISA,IKIND))
                      EISCR(IG,1,IA) = CMPLX(-ARG1*EI,ARG1*ER,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(-ARG2*EI,ARG2*ER,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(-ARG3*EI,ARG3*ER,kind=real_8)
                   ENDDO
                ENDIF
                IF(GEQ0) EISCR(NGW+1,1,IA)=CMPLX(0._real_8,0._real_8,kind=real_8)
                IF(GEQ0) EISCR(NGW+1,2,IA)=CMPLX(0._real_8,0._real_8,kind=real_8)
                IF(GEQ0) EISCR(NGW+1,3,IA)=CMPLX(0._real_8,0._real_8,kind=real_8)
             ENDDO
             CALL ZGEMM('C','N',3*ions0%na(IS),NSTATE,nkpt%ngwk,ZFAC,&
                  &             EISCR(1,1,1),nkpt%ngwk,C0(1,1),nkpt%ngwk,ZZERO,&
                  &             DAI(1,1,1,1),3*mmdim%naxq)
          ELSE
#ifdef __SR8000
             !poption parallel, tlocal(ISA,ARG1,ARG2,ARG3,ER,EI) 
             !voption indep(EISCR)
#endif 
             !$omp parallel do private(IA,IG,ISA,ARG1,ARG2,ARG3,ER,EI)
             DO IA=1,ions0%na(IS)
                ISA=ISA0+IA
                ! Make use of the special structure of CI
                IF(ABS(CIR).GT.0.5_real_8) THEN
                   ! CI is real
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   DO IG=1,NGW
                      ! !                ARG1=GK(1,IG)*TWNL(IG,IV,IS,1)*CIR
                      ! !                ARG2=GK(2,IG)*TWNL(IG,IV,IS,1)*CIR
                      ! !                ARG3=GK(3,IG)*TWNL(IG,IV,IS,1)*CIR
                      ! !                ER=real(EIGR(IG,ISA))
                      ! !                EI=aimag(EIGR(IG,ISA))
                      ARG1=GK(1,IG)
                      ARG2=GK(2,IG)
                      ARG3=GK(3,IG)
                      ER=CIR*TWNL(IG,IV,IS,1)*REAL(EIGR(IG,ISA))
                      EI=CIR*TWNL(IG,IV,IS,1)*AIMAG(EIGR(IG,ISA))
                      EISCR(IG,1,IA) = CMPLX(ARG1*ER,ARG1*EI,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(ARG2*ER,ARG2*EI,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(ARG3*ER,ARG3*EI,kind=real_8)
                   ENDDO
                ELSE
                   ! CI is imaginary
#ifdef _vpp_
                   !OCL NOALIAS
#endif 
                   !dir$ prefervector
                   DO IG=1,NGW
                      ! !                ARG1=GK(1,IG)*TWNL(IG,IV,IS,1)*CII
                      ! !                ARG2=GK(2,IG)*TWNL(IG,IV,IS,1)*CII
                      ! !                ARG3=GK(3,IG)*TWNL(IG,IV,IS,1)*CII
                      ! !                ER=real(EIGR(IG,ISA))
                      ! !                EI=aimag(EIGR(IG,ISA))
                      ARG1=GK(1,IG)
                      ARG2=GK(2,IG)
                      ARG3=GK(3,IG)
                      ER=CII*TWNL(IG,IV,IS,1)*REAL(EIGR(IG,ISA))
                      EI=CII*TWNL(IG,IV,IS,1)*AIMAG(EIGR(IG,ISA))
                      EISCR(IG,1,IA) = CMPLX(-ARG1*EI,ARG1*ER,kind=real_8)
                      EISCR(IG,2,IA) = CMPLX(-ARG2*EI,ARG2*ER,kind=real_8)
                      EISCR(IG,3,IA) = CMPLX(-ARG3*EI,ARG3*ER,kind=real_8)
                   ENDDO
                ENDIF
                IF(GEQ0) THEN 
                   EISCR(1,1,IA)=0.5_real_8*EISCR(1,1,IA)
                   EISCR(1,2,IA)=0.5_real_8*EISCR(1,2,IA)
                   EISCR(1,3,IA)=0.5_real_8*EISCR(1,3,IA)
                ENDIF
             ENDDO
             IF(nkpt%ngwk.GT.0) THEN
                IF (ions0%na(IS).GT.1) THEN
                   CALL DGEMM('T','N',3*ions0%na(IS),NSTATE,2*nkpt%ngwk,TFAC,&
                        &             EISCR(1,1,1),2*nkpt%ngwk,C0(1,1),2*nkpt%ngwk,0.0_real_8,&
                        &             DAI(1,1,1,1),3*mmdim%naxq)
                ELSE
                   CALL DGEMM('T','N',3,NSTATE,2*nkpt%ngwk,TFAC,&
                        &             EISCR(1,1,1),2*nkpt%ngwk,C0(1,1),2*nkpt%ngwk,0.0_real_8,&
                        &             DAI(1,1,1,1),3*mmdim%naxq)
                ENDIF
             ENDIF
          ENDIF
          CALL mp_sum(DAI,IMAGP*3*mmdim%naxq*NSTATE,parai%allgrp)
#ifdef __SR8000
          !poption parallel, tlocal(II)  
#endif 
          !$omp parallel do private(I,IA,IG,II)
          DO I=parap%NST12(parai%mepos,1),parap%NST12(parai%mepos,2)
             II=I-parap%NST12(parai%mepos,1)+1
#ifdef _vpp_
             !OCL NOALIAS
             DO IA=1,ions0%na(IS)
                DFNL(1,ISA0+IA,IV,1,II,IKIND)=DAI(1,1,IA,I)
                DFNL(1,ISA0+IA,IV,2,II,IKIND)=DAI(1,2,IA,I)
                DFNL(1,ISA0+IA,IV,3,II,IKIND)=DAI(1,3,IA,I)
                IF(MASK) THEN 
                   DFNL(2,ISA0+IA,IV,1,II,IKIND)=DAI(2,1,IA,I)
                   DFNL(2,ISA0+IA,IV,2,II,IKIND)=DAI(2,2,IA,I)
                   DFNL(2,ISA0+IA,IV,3,II,IKIND)=DAI(2,3,IA,I)
                ENDIF
             ENDDO
#else 
             !dir$ prefervector
             DO IG=1,IMAGP
                DO IA=1,ions0%na(IS)
                   DFNL(IG,ISA0+IA,IV,1,II,IKIND)=DAI(IG,1,IA,I)
                   DFNL(IG,ISA0+IA,IV,2,II,IKIND)=DAI(IG,2,IA,I)
                   DFNL(IG,ISA0+IA,IV,3,II,IKIND)=DAI(IG,3,IA,I)
                ENDDO
             ENDDO
#endif 
          ENDDO
       ENDDO
       ISA0=ISA0+ions0%na(IS)
    ENDDO

    DEALLOCATE(eiscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dai,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL TIHALT('    RNLSM2',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE RNLSM2
  ! ==================================================================
  SUBROUTINE GIVE_SCR_RNLSM2(LRNLSM2,TAG,NSTATE)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: LRNLSM2
    CHARACTER(LEN=30)                        :: TAG
    INTEGER                                  :: NSTATE

! ==--------------------------------------------------------------==

    IF(NLM.EQ.0) THEN
       LRNLSM2=0
    ELSE
       LRNLSM2=2*nkpt%ngwk*mmdim%naxq*3+IMAGP*mmdim%naxq*NSTATE*3
       TAG=   '2*NGWK*NAXq*3+IMAGP*NAXq*...'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE GIVE_SCR_RNLSM2
  ! ==================================================================
#else 
  ! ==================================================================
  SUBROUTINE RNLSM2(C0,NSTATE,IKPT,IKIND)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DFNL WHICH IS USED IN THE SUBROUTINE RNLFOR       ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==          NOT IMPLEMENTED FOR TSHEL(IS)=TRUE                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: NSTATE
    COMPLEX(real_8)                          :: C0(nkpt%ngwk,NSTATE)
    INTEGER                                  :: IKPT, IKIND

    CHARACTER(*), PARAMETER                  :: procedureN = 'RNLSM2'
    COMPLEX(real_8), PARAMETER               :: ZONE = (1._real_8,0._real_8) ,&
                                                ZZERO = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: CI, ZFAC
    COMPLEX(real_8), ALLOCATABLE             :: EISCR(:,:)
    INTEGER                                  :: I, IA, ierr, IG, II, IS, ISA, &
                                                ISA0, ISUB, IV, K
    REAL(real_8)                             :: ARG, CII, CIR, EI, ER, TFAC
    REAL(real_8), ALLOCATABLE                :: DAI(:,:,:)

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF(NLM.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL TISET('    RNLSM2',ISUB)
    ALLOCATE(eiscr(nkpt%ngwk, mmdim%naxq),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dai(imagp, mmdim%naxq, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(DAI)!,IMAGP*mmdim%naxq*NSTATE)
    IF(tkpts%tkpnt) THEN
       ZFAC=parm%tpiba*ZONE
    ELSE
       TFAC=2._real_8*parm%tpiba
    ENDIF
    DO K=1,3
       ISA0=0
       DO IS=1,ions1%nsp
          DO IV=1,nlps_com%ngh(IS)
             CI=(0.0_real_8,-1.0_real_8)**(NGHTOL(IV,IS)+1)
             CIR=REAL(CI)
             CII=AIMAG(CI)
             IF(tkpts%tkpnt) THEN
                !$omp parallel do private(ia,isa,ig,arg,er,ei)
#ifdef __SR8000 
                !poption parallel, tlocal(ia,isa,ig,arg,er,ei)
#endif 
                DO IA=1,ions0%na(IS)
                   ISA=ISA0+IA
                   ! Make use of the special structure of CI
                   IF(ABS(CIR).GT.0.5_real_8) THEN
                      ! CI is real
                      DO IG=1,ncpw%ngw
                         ARG=(GK(K,IG)+&
                                ! TODO refactor this 
                              &                   RK(K,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR

                         ER=REAL(EIGKR(IG,ISA,IKIND))
                         EI=AIMAG(EIGKR(IG,ISA,IKIND))
                         EISCR(IG,IA) = CMPLX(ARG*ER,ARG*EI,kind=real_8)
                      ENDDO
                      DO IG=ncpw%ngw+1,nkpt%ngwk
                         ARG=(-GK(K,IG-ncpw%ngw)+&
                                ! TODO refactor this
                              &                   RK(K,IKIND))*TWNL(IG,IV,IS,IKIND)*CIR

                         ER=REAL(EIGKR(IG,ISA,IKIND))
                         EI=AIMAG(EIGKR(IG,ISA,IKIND))
                         EISCR(IG,IA) = CMPLX(ARG*ER,ARG*EI,kind=real_8)
                      ENDDO
                   ELSE
                      ! CI is imaginary
                      DO IG=1,ncpw%ngw
                         ARG=(GK(K,IG)+&
                                ! TODO refactor this
                              &                   RK(K,IKIND))*TWNL(IG,IV,IS, IKIND) *CII

                         ER=REAL(EIGKR(IG,ISA,IKIND))
                         EI=AIMAG(EIGKR(IG,ISA,IKIND))
                         EISCR(IG,IA) = CMPLX(-ARG*EI,ARG*ER,kind=real_8)
                      ENDDO
                      DO IG=ncpw%ngw+1,nkpt%ngwk
                         ARG=(-GK(K,IG-ncpw%ngw)+&
                                ! TODO refactor this
                              &                   RK(K,IKIND))*TWNL(IG,IV,IS,IKIND)*CII

                         ER=REAL(EIGKR(IG,ISA,IKIND))
                         EI=AIMAG(EIGKR(IG,ISA,IKIND))
                         EISCR(IG,IA) = CMPLX(-ARG*EI,ARG*ER,kind=real_8)
                      ENDDO
                   ENDIF
                   IF(GEQ0) EISCR(ncpw%ngw+1,IA)=CMPLX(0._real_8,0._real_8,kind=real_8)
                ENDDO
                CALL ZGEMM('C','N',ions0%na(IS),NSTATE,nkpt%ngwk,ZFAC,&
                     &               EISCR(1,1),nkpt%ngwk,C0(1,1),nkpt%ngwk,ZZERO,&
                     &               DAI(1,1,1),mmdim%naxq)
             ELSE
                !$omp parallel do private(ia,isa,ig,arg,er,ei)
#ifdef __SR8000
                !poption parallel, tlocal(ia,isa,ig,arg,er,ei)
#endif 
                DO IA=1,ions0%na(IS)
                   ISA=ISA0+IA
                   ! Make use of the special structure of CI
                   IF(ABS(CIR).GT.0.5_real_8) THEN
                      ! CI is real
                      DO IG=1,ncpw%ngw
                         ARG=GK(K,IG)*TWNL(IG,IV,IS,1)*CIR
                         ER=REAL(EIGR(IG,ISA,1))
                         EI=AIMAG(EIGR(IG,ISA,1))
                         EISCR(IG,IA) = CMPLX(ARG*ER,ARG*EI,kind=real_8)
                      ENDDO
                   ELSE
                      ! CI is imaginary
                      DO IG=1,ncpw%ngw
                         ARG=GK(K,IG)*TWNL(IG,IV,IS,1)*CII
                         ER=REAL(EIGR(IG,ISA,1))
                         EI=AIMAG(EIGR(IG,ISA,1))
                         EISCR(IG,IA) = CMPLX(-ARG*EI,ARG*ER,kind=real_8)
                      ENDDO
                   ENDIF
                   IF(GEQ0) EISCR(1,IA)=0.5_real_8*EISCR(1,IA)
                ENDDO
                IF(nkpt%ngwk.GT.0) THEN
                   IF (ions0%na(IS).GT.1) THEN
                      CALL DGEMM('T','N',ions0%na(IS),NSTATE,2*nkpt%ngwk,TFAC,&
                           &               EISCR(1,1),2*nkpt%ngwk,C0(1,1),2*nkpt%ngwk,0.0_real_8,&
                           &               DAI(1,1,1),mmdim%naxq)
                   ELSE
                      CALL DGEMV('T',2*nkpt%ngwk,NSTATE,TFAC,C0(1,1),2*nkpt%ngwk,&
                           &               EISCR(1,1),1,0.0_real_8,DAI(1,1,1),mmdim%naxq)
                   ENDIF
                ENDIF
             ENDIF
             CALL mp_sum(DAI,IMAGP*mmdim%naxq*NSTATE,parai%allgrp)
             DO I=parap%NST12(parai%mepos,1),parap%NST12(parai%mepos,2)
                II=I-parap%NST12(parai%mepos,1)+1
                CALL DCOPY(IMAGP*ions0%na(IS),DAI(1,1,I),1,&
                     &                   DFNL(1,ISA0+1,IV,K,II,IKIND),1)
             ENDDO
          ENDDO
          ISA0=ISA0+ions0%na(IS)
       ENDDO
    ENDDO

    DEALLOCATE(eiscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dai,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL TIHALT('    RNLSM2',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE RNLSM2
  ! ==================================================================
  SUBROUTINE GIVE_SCR_RNLSM2(LRNLSM2,TAG,NSTATE)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: LRNLSM2
    CHARACTER(LEN=30)                        :: TAG
    INTEGER                                  :: NSTATE

    IF(NLM.EQ.0) THEN
       LRNLSM2=0
    ELSE
       LRNLSM2=2*nkpt%ngwk*mmdim%naxq+IMAGP*mmdim%naxq*NSTATE
       TAG=   '2*NGWK*NAXq+IMAGP*NAXq*NSTATE'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE GIVE_SCR_RNLSM2
  ! ==================================================================
#endif 

END MODULE rnlsm2_utils
