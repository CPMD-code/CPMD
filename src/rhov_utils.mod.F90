MODULE rhov_utils
  USE cppt,                            ONLY: indz,&
                                             inyh,&
                                             nzh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhov
  PUBLIC :: give_scr_rhov

CONTAINS

  ! ==================================================================
  SUBROUTINE rhov(nstate,is1,is2,rsumv,psi)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE VANDERBILT DENSITY       ==
    ! ==                                                              ==
    ! == N_V(G) = SUM_I,IJ RHO_I,IJ Q_I,JI(G) E^-IG.R_I               ==
    ! == RHO_I,IJ = SUM_N < BETA_I,I | PSI_N >< PSI_N | BETA_I,J >    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, is1, is2
    REAL(real_8)                             :: rsumv
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhov'

    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:), deltar(:), qg(:)
    INTEGER                                  :: i, ia, ierr, ig, il_ctmp, &
                                                il_ddia, il_deltar, il_qg, &
                                                is, isa, isa0, isub, iv, jv
    REAL(real_8)                             :: fac, sum
    REAL(real_8), ALLOCATABLE                :: ddia(:)

    CALL tiset('      RHOV',isub)
    CALL setfftn(0)
    IF (cntl%tfdist) CALL stopgm('RHOV','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('RHOV','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! IF(NOGRP.GT.1) CALL STOPGM('RHOV','NO GROUPS ALLOWED')
    ! SCR partition
    il_deltar = 2*ncpw%nhg
    il_qg     = 2*ncpw%nhg
    il_ddia   = maxsys%nax
    il_ctmp   = 2*ncpw%nhg
    ALLOCATE(deltar(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ddia(maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ctmp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(deltar)!,nhg)
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                CALL qvan2(iv,jv,is,qg)
                IF (cntl%bigmem) THEN
                   CALL zeroing(ddia(1:ions0%na(is)))
                   DO i=is1,is2
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         ddia(ia)=ddia(ia)+&
                              crge%f(i,1)*fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                      ENDDO
                   ENDDO
                   fac=1.0_real_8
                   IF (iv.NE.jv) fac=2.0_real_8
                   CALL dgemv('N',2*ncpw%nhg,ions0%na(is),fac,eigrb(1,isa0+1),2*ncpw%nhg,&
                        ddia,1,0.0_real_8,ctmp,1)
                   !$omp parallel do private(IG)
                   !CDIR NODEP
                   DO ig=1,ncpw%nhg
                      deltar(ig)=deltar(ig)+qg(ig)*ctmp(ig)
                   ENDDO
                ELSE
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      !$omp parallel do private(IG)
#ifdef __SR8000
                      !poption parallel
#endif
                      DO ig=1,ncpw%nhg
                         ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                      ENDDO
                      sum=0.0_real_8
                      !$omp parallel do private(I) reduction(+:SUM)
                      DO i=is1,is2
                         sum = sum +&
                              crge%f(i,1)*fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                      ENDDO
                      IF (iv.NE.jv) sum=2._real_8*sum
                      !$omp parallel do private(IG)
                      !CDIR NODEP
                      DO ig=1,ncpw%nhg
                         deltar(ig)=deltar(ig)+qg(ig)*sum*ctmp(ig)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    IF (geq0) THEN
       rsumv=REAL(deltar(1))
    ELSE
       rsumv=0.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(psi)!,maxfft)
    !$omp parallel do private(IG) shared(PSI)
    !CDIR NODEP
#ifdef __SR8000
    !poption parallel
#endif
    DO ig=1,ncpw%nhg
       psi(nzh(ig))=deltar(ig)
       psi(indz(ig))=CONJG(deltar(ig))
    ENDDO
    CALL invfftn(psi, .FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deltar,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ctmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('      RHOV',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhov
  ! ==================================================================
  SUBROUTINE give_scr_rhov(lrhov,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhov
    CHARACTER(len=30)                        :: tag

    lrhov=6*ncpw%nhg+maxsys%nax
    tag='6*NHG+maxsys%nax'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhov
  ! ==================================================================

END MODULE rhov_utils
