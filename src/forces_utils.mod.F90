MODULE forces_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE fnonloc_utils,                   ONLY: give_scr_fnonloc
  USE jrotation_utils,                 ONLY: set_orbdist
  USE nlforce_utils,                   ONLY: give_scr_nlforce
  USE nlps,                            ONLY: imagp
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE spin,                            ONLY: lspin2
  USE summat_utils,                    ONLY: give_scr_summat
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: cnti,&
                                             cntl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_forces

CONTAINS
  ! ==================================================================
  SUBROUTINE give_scr_forces(lforces,tag,nstate,lproj,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: lproj, tfor

    INTEGER :: il_amat, il_auxc, il_ddia, il_fsc, il_gam, il_psiab, il_scr, &
      il_scrdip, lrnlsm, lrscpot, lsummat, lsymvec, nstx

    lsummat=0
    lsymvec=0
    lrnlsm=0
    lrscpot=0
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    IF (tfor) THEN
       CALL give_scr_symvec(lsymvec,tag)
    ENDIF
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    il_auxc=0
    IF (pslo_com%tivan .AND. lproj .AND. cnti%iproj.NE.0) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
       CALL give_scr_summat(lsummat,tag,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       IF (cntl%tdmal) THEN
          il_gam = 1
          lsummat=0
       ELSE
          il_gam = imagp*nstate*nstate
          CALL give_scr_summat(lsummat,tag,nstate)
       ENDIF
    ENDIF
    il_scrdip=0
    IF (cntl%tfield) CALL give_scr_opeigr(il_scrdip,tag,nstate)
    IF (cntl%tdmal) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       il_amat=3*nstate*nstx
    ELSE
       il_amat=0
       il_auxc=MAX(il_auxc,2*lsummat)
       il_auxc=MAX(il_auxc,nstate**2)! AUXC space for OVLAP (Laio A.)
    ENDIF
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    il_scr=il_gam+il_auxc+il_ddia+il_fsc+il_psiab+il_scrdip+100
    il_scr=il_scr+il_amat
    lforces = MAX(lrnlsm,lrscpot,lsymvec,il_scr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_forces
  ! ==================================================================

END MODULE forces_utils
