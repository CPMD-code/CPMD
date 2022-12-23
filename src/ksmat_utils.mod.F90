MODULE ksmat_utils
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             xxmat,&
                                             zxmat
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE spin,                            ONLY: lspin2
  USE system,                          ONLY: cntl,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ksmat
  PUBLIC :: give_scr_ksmat

CONTAINS

  ! ==================================================================
  SUBROUTINE ksmat(c2,vpot,psi,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE KOHN-SHAM MATRIX IN THE ATOMIC WAVEFUNCTION BASIS       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(:,:)
    REAL(real_8)                             :: vpot(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, ikind

    COMPLEX(real_8)                          :: pab(1)
    INTEGER                                  :: ia, il_auxc, il_ddia, is, &
                                                ist, isub, natst
    LOGICAL                                  :: tfdist2
    REAL(real_8), DIMENSION(20)              :: foc = 1.0_real_8

! ==--------------------------------------------------------------==

    IF (lspin2%tlse) CALL stopgm('KSMAT','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    ! AK      CALL zeroing(PAB)!,NNR1)
    pab(1)=0.0_real_8
    CALL tiset('     KSMAT',isub)
    ! ==--------------------------------------------------------------==
    tfdist2=cntl%tfdist
    cntl%tfdist=.FALSE.
    CALL fnl_set('SAVE')
    CALL fnlalloc(atwp%numaormax,.FALSE.,.FALSE.)
    CALL give_scr_fnonloc(il_auxc,il_ddia,atwp%numaormax)
    ist=1
    DO is=1,ions1%nsp
       natst=atwf_mod%numaor(is)
       DO ia=1,ions0%na(is)
          CALL rnlsm(catom(:,ist:ist+natst-1),natst,1,ikind,.FALSE.)
          CALL zeroing(c2(:,1:natst))!,nkpt%ngwk*natst)
          CALL vpsi(catom(:,ist:ist+natst-1),c2,foc,vpot,psi,natst,ikind,1,.TRUE.)
          CALL fnonloc(c2,foc,natst,ikind,1,.TRUE.)
          IF (tkpts%tkpnt) THEN
             CALL ovlap2_c(nkpt%ngwk,atwp%nattot,natst,zxmat(1,ist),catom,c2)
          ELSE
             CALL ovlap2(nkpt%ngwk,atwp%nattot,natst,xxmat(1,ist),catom,c2,.TRUE.)
          ENDIF
          ist=ist+natst
       ENDDO
    ENDDO
    CALL fnldealloc(.FALSE.,.FALSE.)
    CALL fnl_set('RECV')
    cntl%tfdist=tfdist2
    IF (tkpts%tkpnt) THEN
       CALL dscal(2*atwp%nattot*atwp%nattot,-1._real_8,zxmat(1,1),1)
    ELSE
       CALL dscal(  atwp%nattot*atwp%nattot,-1._real_8,xxmat(1,1),1)
    ENDIF
    CALL tihalt('     KSMAT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ksmat
  ! ==================================================================
  SUBROUTINE give_scr_ksmat(lksmat,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lksmat
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_auxc, il_ddia, lfnonloc, &
                                                lrnlsm

    CALL give_scr_rnlsm(lrnlsm,tag,atwp%nattot,.FALSE.)
    CALL give_scr_fnonloc(il_auxc,il_ddia,atwp%numaormax)
    lfnonloc = il_auxc + il_ddia
    lksmat=MAX(lfnonloc,lrnlsm)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ksmat
  ! ==================================================================

END MODULE ksmat_utils
