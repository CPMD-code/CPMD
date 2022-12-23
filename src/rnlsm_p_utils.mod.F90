MODULE rnlsm_p_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             ndfnl,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: eigr
  USE system,                          ONLY: maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm3
  PUBLIC :: give_scr_rnlsm3

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm3(c0,nstate,dd_fnl)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DDFNL WHICH IS USED IN PERTURBATION THEORY FOR    ==
    ! ==   PHONON
    ! ==  K-POINT VERSION IS NOT IMPLEMENTED                              ==
    ! ==          NOT IMPLEMENTED FOR TSHEL(IS)=TRUE                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    REAL(real_8) :: dd_fnl(imagp,ions1%nat,maxsys%nhxs,3,3,ndfnl,nkpt%nkpnt)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm3'

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: eiscr(:,:)
    INTEGER                                  :: i, ia, ierr, ig, ii, &
                                                il_eiscr, is, isa, isa0, &
                                                isub, iv, k, k2
    REAL(real_8)                             :: arg, cii, cir, ei, er, tfac
    REAL(real_8), ALLOCATABLE                :: dai(:,:,:)

    CALL tiset('    RNLSM3',isub)
    ! ==--------------------------------------------------------------==
    ! SCR partition
    il_eiscr = 2*nkpt%ngwk*maxsys%nax
    ALLOCATE(eiscr(nkpt%ngwk, maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dai(imagp, maxsys%nax, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(dai)!,imagp*maxsys%nax*nstate)
    IF (tkpts%tkpnt) THEN
       CALL stopgm(' RNLSM3', 'KPOINT NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSE
       tfac=2._real_8*parm%tpiba*parm%tpiba
    ENDIF
    DO k=1,3
       DO k2=1,3
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                ci=(0.0_real_8,-1.0_real_8)**(nghtol(iv,is)+2)
                cir=REAL(ci)
                cii=AIMAG(ci)
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   ! Make use of the special structure of CI
                   IF (ABS(cir).GT.0.5_real_8) THEN
                      ! CI is real
                      DO ig=1,ncpw%ngw
                         arg=gk(k,ig)*gk(k2,ig)*twnl(ig,iv,is,1)*cir
                         er=REAL(eigr(ig,isa,1))
                         ei=AIMAG(eigr(ig,isa,1))
                         eiscr(ig,ia) = CMPLX(arg*er,arg*ei,kind=real_8)
                      ENDDO
                   ELSE
                      ! CI is imaginary
                      DO ig=1,ncpw%ngw
                         arg=gk(k,ig)*gk(k2,ig)*twnl(ig,iv,is,1)*cii
                         er=REAL(eigr(ig,isa,1))
                         ei=AIMAG(eigr(ig,isa,1))
                         eiscr(ig,ia) = CMPLX(-arg*ei,arg*er,kind=real_8)
                      ENDDO
                   ENDIF
                   IF (geq0) eiscr(1,ia)=0.5_real_8*eiscr(1,ia)
                ENDDO
                CALL dgemm('T','N',ions0%na(is),nstate,2*nkpt%ngwk,tfac,eiscr(1,1),2*&
                     nkpt%ngwk,c0(1,1),2*nkpt%ngwk,0.0_real_8,dai(1,1,1),maxsys%nax)
                CALL mp_sum(dai,imagp*maxsys%nax*nstate,parai%allgrp)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   CALL dcopy(imagp*ions0%na(is),dai(1,1,i),1,dd_fnl(1,isa0+1,iv,&
                        k,k2,ii,1),1)
                ENDDO
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eiscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dai,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    RNLSM3',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlsm3
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm3(lrnlsm3,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm3
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lrnlsm3=2*nkpt%ngwk*maxsys%nax+imagp*maxsys%nax*nstate
    tag=   '2*NGWK*maxsys%nax+IMAGP*maxsys%nax*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm3
  ! ==================================================================

END MODULE rnlsm_p_utils
