MODULE rho1pri_utils
  !!use densto_utils, only : densto
  USE cppt,                            ONLY: nzh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE rho1ofr_utils,                   ONLY: rho1ofr
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rho1pri

CONTAINS

  ! ==================================================================
  SUBROUTINE rho1pri(c0,c1,tau0,nstate,nres,label)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE THE LINEAR RESPONSE DENSITY AND DUMP              ==
    ! ==  ALL DATA TO FILE DENSITY                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: nres
    CHARACTER(len=1)                         :: label

    CHARACTER(*), PARAMETER                  :: procedureN = 'rho1pri'

    CHARACTER(len=1)                         :: cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:), rhog(:)
    INTEGER                                  :: ierr, ig, il_psi_1d, &
                                                il_psi_2d, il_rhoe_1d, &
                                                il_rhoe_2d, ir, is
    REAL(real_8), ALLOCATABLE                :: rhoe(:,:)

    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhog(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    DO is=1,nres
       CALL rho1ofr(c0,c1(:,:,is),crge%f(:,1),rhoe,psi(:,1),nstate)
       ! 
       IF (cntl%tlsd) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%nhg
             rhog(ig) = psi(nzh(ig),1)
          ENDDO
          IF (paral%io_parent)&
               WRITE(cipnum,'(I1)') is
          filen=label//'LRDENS'//cipnum//'_ALPHA'
          CALL densto(rhog,tau0,filen)
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             psi(ir,1)=CMPLX(rhoe(ir,2),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%nhg
             rhog(ig) = psi(nzh(ig),1)
          ENDDO
          IF (paral%io_parent)&
               WRITE(cipnum,'(I1)') is
          filen=label//'LRDENS'//cipnum//'_BETA'
          CALL densto(rhog,tau0,filen)
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%nhg
             rhog(ig) = psi(nzh(ig),1)
          ENDDO
          IF (paral%io_parent)&
               WRITE(cipnum,'(I1)') is
          filen=label//'LRDENS'//cipnuM
          CALL densto(rhog,tau0,filen)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rho1pri
  ! ==================================================================

END MODULE rho1pri_utils
