MODULE rnl_dk_p_utils
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  USE response_pmod,                   ONLY: ddfnl_ddk,&
                                             ddtwnl_ddk,&
                                             dfnl_dk,&
                                             dtwnl_dk
  USE sumfnl_utils,                    ONLY: sumfnl
  USE system,                          ONLY: maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnl_dk_p
  PUBLIC :: give_scr_rnl_dk_p

CONTAINS

  ! ==================================================================
  SUBROUTINE rnl_dk_p(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                    THE ARRAY DFNL_DK                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnl_dk_p'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8), &
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: c00(:,:), scr(:,:)
    INTEGER                                  :: ia, idir, idir2, ierr, ig, &
                                                is, isa, isa0, istate, isub, &
                                                iv
    REAL(real_8)                             :: cii, cir, ei, er, t

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  RNL_DK_P',isub)
    ! SCR test.
    ALLOCATE(scr(2*ncpw%ngw,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ALLOCATE(c00(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c00)!,SIZE(c00))

    CALL zeroing(scr)!,SIZE(scr))
    !$omp parallel do private(ISTATE,IG)
    DO istate = 1, nstate
       DO ig = 1, ncpw%ngw
          c00(ig,istate) = c0(ig,istate)
          c00(ig+ncpw%ngw,istate) = CONJG(c0(ig,istate))
       ENDDO
    ENDDO

    imagp = 2

    IF (ncpw%ngw.EQ.0) GOTO 1000
    DO idir = 1,3
       isa0=0
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             cir=REAL(ci)
             cii=AIMAG(ci)
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                ! Make use of the special structure of CI
                IF (ABS(cir).GT.0.5_real_8) THEN
                   ! CI is real
                   !$omp parallel do private(IG,ER,EI,T)
                   DO ig=1,ncpw%ngw
                      er=REAL(eigkr(ig,isa,1))
                      ei=AIMAG(eigkr(ig,isa,1))
                      t=dtwnl_dk(ig,iv,is,idir)*cir
                      scr(ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                   ENDDO
                   ! =============== and the part from NGW+1 to 2*NGW
                   !$omp parallel do private(IG,ER,EI,T)
                   DO ig=1,ncpw%ngw
                      er=REAL(eigkr(ig,isa,1))
                      ei=-AIMAG(eigkr(ig,isa,1))
                      t=dtwnl_dk(ncpw%ngw+ig,iv,is,idir)*cir
                      scr(ncpw%ngw+ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                   ENDDO
                ELSE
                   ! CI is imaginary
                   !$omp parallel do private(IG,ER,EI,T)
                   DO ig=1,ncpw%ngw
                      er=REAL(eigkr(ig,isa,1))
                      ei=AIMAG(eigkr(ig,isa,1))
                      t=dtwnl_dk(ig,iv,is,idir)*cii
                      scr(ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                   ENDDO

                   ! =============== and the part from NGW+1 to 2*NGW
                   !$omp parallel do private(IG,ER,EI,T)
                   DO ig=1,ncpw%ngw
                      er=REAL(eigkr(ig,isa,1))
                      ei=-AIMAG(eigkr(ig,isa,1))
                      t=dtwnl_dk(ncpw%ngw+ig,iv,is,idir)*cii
                      scr(ncpw%ngw+ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                   ENDDO

                ENDIF
                IF (geq0) THEN
                   ! SCR(1,IA)=0.5_real_8*SCR(1,IA)
                   scr(ncpw%ngw+1,ia)=CMPLX(0._real_8,0._real_8,kind=real_8)
                ENDIF
             ENDDO

             CALL zgemm('C','N',ions0%na(is),nstate,2*ncpw%ngw,zone,scr(1,1),2* ncpw%ngw,&
                  c00(1,1),2*ncpw%ngw,zzero,dfnl_dk(1,isa0+1,iv,1,idir),ions1%nat*&
                  maxsys%nhxs)

          ENDDO
          isa0=isa0+ions0%na(is)
       ENDDO

       CALL sumfnl(dfnl_dk(1,1,1,1,idir),nstate)

    ENDDO       ! IDIR


    ! ****************************
    ! SECOND DERIVATIVES
    ! ****************************


    ! Form Factor relative to the second derivatives
    DO idir = 1,3
       DO idir2 = idir,3
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                cir=REAL(ci)
                cii=AIMAG(ci)
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   ! Make use of the special structure of CI
                   IF (ABS(cir).GT.0.5_real_8) THEN
                      ! CI is real
                      !$omp parallel do private(IG,ER,EI,T)
                      DO ig=1,ncpw%ngw
                         er=REAL(eigkr(ig,isa,1))
                         ei=AIMAG(eigkr(ig,isa,1))
                         t=ddtwnl_ddk(ig,iv,is,idir,idir2)*cir
                         scr(ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                      ENDDO
                      ! =============== and the part from NGW+1 to 2*NGW
                      !$omp parallel do private(IG,ER,EI,T)
                      DO ig=1,ncpw%ngw
                         er=REAL(eigkr(ig,isa,1))
                         ei=-AIMAG(eigkr(ig,isa,1))
                         t=ddtwnl_ddk(ncpw%ngw+ig,iv,is,idir,idir2)*cir
                         scr(ncpw%ngw+ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                      ENDDO
                   ELSE
                      ! CI is imaginary
                      !$omp parallel do private(IG,ER,EI,T)
                      DO ig=1,ncpw%ngw
                         er=REAL(eigkr(ig,isa,1))
                         ei=AIMAG(eigkr(ig,isa,1))
                         t=ddtwnl_ddk(ig,iv,is,idir,idir2)*cii
                         scr(ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                      ENDDO
                      ! =============== and the part from NGW+1 to 2*NGW
                      !$omp parallel do private(IG,ER,EI,T)
                      DO ig=1,ncpw%ngw
                         er=REAL(eigkr(ig,isa,1))
                         ei=-AIMAG(eigkr(ig,isa,1))
                         t=ddtwnl_ddk(ncpw%ngw+ig,iv,is,idir,idir2)*cii
                         scr(ncpw%ngw+ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                      ENDDO

                   ENDIF
                   IF (geq0) THEN
                      ! SCR(1,IA)=0.5_real_8*SCR(1,IA)
                      scr(ncpw%ngw+1,ia)=CMPLX(0._real_8,0._real_8,kind=real_8)
                   ENDIF
                ENDDO

                CALL zgemm('C','N',ions0%na(is),nstate,2*ncpw%ngw,zone,scr(1,1), 2*&
                     ncpw%ngw,c00(1,1),2*ncpw%ngw,zzero,ddfnl_ddk(1,isa0+1,iv,1,&
                     idir, idir2),ions1%nat*maxsys%nhxs)


             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO

          CALL sumfnl(ddfnl_ddk(1,1,1,1,idir,idir2),nstate)
          CALL dcopy(2*ions1%nat*maxsys%nhxs*nstate,ddfnl_ddk(1,1,1,1,idir,idir2),&
               1,ddfnl_ddk(1,1,1,1,idir2,idir),1)
       ENDDO
    ENDDO
1000 CONTINUE

    ! ==-------------------------------------------------------------==      
    imagp =1

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(c00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('  RNL_DK_P',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnl_dk_p
  ! ==================================================================
  SUBROUTINE give_scr_rnl_dk_p(lrnl_dk_p,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnl_dk_p
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lsumfnl

    IF (nlm.EQ.0) THEN
       lrnl_dk_p=0
    ELSE
       lsumfnl = 1
       lrnl_dk_p=MAX(4*nkpt%ngwk*maxsys%nax,lsumfnl)
       tag   ='MAX(2*NGWK*maxsys%nax,LSUMFNL)'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnl_dk_p
  ! ==================================================================


END MODULE rnl_dk_p_utils
