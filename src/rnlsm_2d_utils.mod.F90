MODULE rnlsm_2d_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlm,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ddfnl,&
                                             eigr
  USE system,                          ONLY: ncpw,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm_2d
  PUBLIC :: give_scr_rnlsm_2d

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm_2d(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DDFNL WHICH IS USED IN THE SUBROUTINE SD_NL2      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm_2d'

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: eiscr(:,:)
    INTEGER                                  :: i, ia, ierr, ig, ii, is, isa, &
                                                isa0, isub, iv
    REAL(real_8)                             :: arg(6), cii, cir, ei, er, tfac
    REAL(real_8), ALLOCATABLE                :: dai(:,:)

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  RNLSM_2D',isub)
    ALLOCATE(eiscr(ncpw%ngw,6),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dai(6,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tfac=2._real_8*parm%tpiba2
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          ci=(0.0_real_8,-1.0_real_8)**(nghtol(iv,is)+2)
          cir=REAL(ci)
          cii=AIMAG(ci)
          DO ia=1,ions0%na(is)
             isa=isa0+ia
             ! ..Make use of the special structure of CI
             IF (ABS(cir).GT.0.5_real_8) THEN
                ! ..CI is real
                DO ig=1,ncpw%ngw
                   arg(1)=gk(1,ig)*gk(1,ig)*twnl(ig,iv,is,1)*cir
                   arg(2)=gk(1,ig)*gk(2,ig)*twnl(ig,iv,is,1)*cir
                   arg(3)=gk(1,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cir
                   arg(4)=gk(2,ig)*gk(2,ig)*twnl(ig,iv,is,1)*cir
                   arg(5)=gk(2,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cir
                   arg(6)=gk(3,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cir
                   er=REAL(eigr(ig,isa,1))
                   ei=AIMAG(eigr(ig,isa,1))
                   eiscr(ig,1) = CMPLX(arg(1)*er,arg(1)*ei,kind=real_8)
                   eiscr(ig,2) = CMPLX(arg(2)*er,arg(2)*ei,kind=real_8)
                   eiscr(ig,3) = CMPLX(arg(3)*er,arg(3)*ei,kind=real_8)
                   eiscr(ig,4) = CMPLX(arg(4)*er,arg(4)*ei,kind=real_8)
                   eiscr(ig,5) = CMPLX(arg(5)*er,arg(5)*ei,kind=real_8)
                   eiscr(ig,6) = CMPLX(arg(6)*er,arg(6)*ei,kind=real_8)
                ENDDO
             ELSE
                ! ..CI is imaginary
                DO ig=1,ncpw%ngw
                   arg(1)=gk(1,ig)*gk(1,ig)*twnl(ig,iv,is,1)*cii
                   arg(2)=gk(1,ig)*gk(2,ig)*twnl(ig,iv,is,1)*cii
                   arg(3)=gk(1,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cii
                   arg(4)=gk(2,ig)*gk(2,ig)*twnl(ig,iv,is,1)*cii
                   arg(5)=gk(2,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cii
                   arg(6)=gk(3,ig)*gk(3,ig)*twnl(ig,iv,is,1)*cii
                   er=REAL(eigr(ig,isa,1))
                   ei=AIMAG(eigr(ig,isa,1))
                   eiscr(ig,1) = CMPLX(-arg(1)*ei,arg(1)*er,kind=real_8)
                   eiscr(ig,2) = CMPLX(-arg(2)*ei,arg(2)*er,kind=real_8)
                   eiscr(ig,3) = CMPLX(-arg(3)*ei,arg(3)*er,kind=real_8)
                   eiscr(ig,4) = CMPLX(-arg(4)*ei,arg(4)*er,kind=real_8)
                   eiscr(ig,5) = CMPLX(-arg(5)*ei,arg(5)*er,kind=real_8)
                   eiscr(ig,6) = CMPLX(-arg(6)*ei,arg(6)*er,kind=real_8)
                ENDDO
             ENDIF
             IF (geq0) THEN
                DO i=1,6
                   eiscr(1,i)=0.5_real_8*eiscr(1,i)
                ENDDO
             ENDIF
             CALL dgemm('T','N',6,nstate,2*ncpw%ngw,tfac,eiscr(1,1),2*ncpw%ngw,c0(&
                  1,1),2*ncpw%ngw,0.0_real_8,dai(1,1),6)
             CALL mp_sum(dai,6*nstate,parai%allgrp)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ii=i-parap%nst12(parai%mepos,1)+1
                ddfnl(isa,iv,1,ii)=dai(1,i)
                ddfnl(isa,iv,2,ii)=dai(2,i)
                ddfnl(isa,iv,3,ii)=dai(3,i)
                ddfnl(isa,iv,4,ii)=dai(4,i)
                ddfnl(isa,iv,5,ii)=dai(5,i)
                ddfnl(isa,iv,6,ii)=dai(6,i)
             ENDDO
          ENDDO
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO
    DEALLOCATE(eiscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dai,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  RNLSM_2D',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlsm_2d
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm_2d(lrnlsm,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    IF (nlm.EQ.0) THEN
       lrnlsm=0
    ELSE
       lrnlsm=2*ncpw%ngw*6+6*nstate+100
       tag=   '2*NGW*6+6*NSTATE'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm_2d
  ! ==================================================================

END MODULE rnlsm_2d_utils
