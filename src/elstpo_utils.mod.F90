MODULE elstpo_utils
  USE cppt,                            ONLY: hipz,&
                                             indz,&
                                             nzh,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE hip_utils,                       ONLY: hip
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE system,                          ONLY: ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: elstpo

CONTAINS

  ! ==================================================================
  SUBROUTINE elstpo(rhog,eirop,eivps,vscr,qphi)
    ! ==--------------------------------------------------------------==
    ! ==            COMPUTES ELECTROSTATIC POTENTIAL                  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: rhog(:), eirop(:), eivps(:), &
                                                vscr(:), qphi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'elstpo'

    COMPLEX(real_8)                          :: rtot
    INTEGER                                  :: ig, ig1, isub

    CALL tiset(procedureN,isub)

    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          CALL setfftn(0)
          ! ..Hockney method
          CALL zeroing(vscr)!,maxfft)
          ig1=1
          IF (geq0) ig1=2
          !CDIR NODEP
          !$omp parallel do private(IG,RTOT)
          DO ig=ig1,ncpw%nhg
             rtot=rhog(ig)+eirop(ig)
             vscr(nzh(ig))=rtot
             vscr(indz(ig))=CONJG(rtot)
          ENDDO
          IF (geq0) vscr(nzh(1))=rhog(1)+eirop(1)
          CALL  invfftn(vscr,.FALSE.,parai%allgrp)
          CALL hip(vscr,qphi)
          CALL  fwfftn(vscr,.FALSE.,parai%allgrp)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%nhg
             eivps(ig)=eivps(ig)+vscr(nzh(ig))! bugfix
          ENDDO
          !$omp parallel do private(IG,RTOT)
          DO ig=1,ncpw%nhg
             rtot=rhog(ig)+eirop(ig)
             vscr(ig)=hipz(ig)*rtot+eivps(ig)! Be consistent with VOFRHOH
          ENDDO
       ELSEIF (isos3%ps_type.EQ.2.OR.isos3%ps_type.EQ.3) THEN
          ! ..Martyna & Tuckerman method
          IF (geq0) THEN
             ! !!         VSCR(1)=EIVPS(1)
             rtot=rhog(1)+eirop(1)
             vscr(1)=rtot*scg(1)! Be consistent with VOFRHOT
             ig1=2
          ELSE
             ig1=1
          ENDIF
          !$omp parallel do private(IG,RTOT)
          DO ig=ig1,ncpw%nhg
             rtot=rhog(ig)+eirop(ig)
             vscr(ig)=rtot*scg(ig)+eivps(ig)! bugfix
          ENDDO
       ELSE
          CALL stopgm('ELSTPO','PS_TYPE not implemented',&
               __LINE__,__FILE__)
       ENDIF
    ELSE
       IF (geq0) THEN
          ! !!       VSCR(1)=EIVPS(1)
          rtot=rhog(1)+eirop(1)
          vscr(1)=rtot*scg(1)! Be consistent with VOFRHOA
          ig1=2
       ELSE
          ig1=1
       ENDIF
       !$omp parallel do private(IG,RTOT)
       DO ig=ig1,ncpw%nhg
          rtot=rhog(ig)+eirop(ig)
          vscr(ig)=rtot*scg(ig)+eivps(ig)! bugfix
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE elstpo
  ! ==================================================================

END MODULE elstpo_utils
