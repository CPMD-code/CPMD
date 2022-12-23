MODULE k_forces_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE fnonloc_utils,                   ONLY: give_scr_fnonloc
  USE nlforce_utils,                   ONLY: give_scr_nlforce
  USE nlps,                            ONLY: imagp
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE spin,                            ONLY: lspin2
  USE summat_utils,                    ONLY: give_scr_summat
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: cnti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_kforces

CONTAINS 
  ! ==================================================================
  SUBROUTINE give_scr_kforces(lforces,tag,nstate,lproj,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: lproj, tfor

    INTEGER :: il_auxc, il_ddia, il_fsc, il_gam, il_ksener, il_psiab, il_scr, &
      lortho, lrnlsm, lrscpot, lsummat, lsymvec

    CALL give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    IF (tfor) THEN
       CALL give_scr_symvec(lsymvec,tag)
    ELSE
       lsymvec=0
    ENDIF
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    IF (pslo_com%tivan .AND. lproj .AND. cnti%iproj.NE.0) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       il_gam = imagp*nstate*nstate
    ENDIF
    CALL give_scr_summat(lsummat,tag,nstate)
    il_auxc=MAX(il_auxc,2*lsummat)
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    il_scr=il_gam+il_auxc+il_ddia+il_fsc+il_psiab+100
    il_ksener = 2*nstate*nstate + 3*2*nstate + nstate/2 +&
         MAX(1,2*2*nstate-1)
    CALL give_scr_ortho(lortho,tag,nstate)
    lforces = MAX(lrnlsm,lrscpot,lsymvec,il_scr,il_ksener,lortho)
    tag='LRNLSM,LRSCPOT,LSYMVEC,..'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_kforces
  ! ==================================================================
END MODULE k_forces_utils
