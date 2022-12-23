MODULE rhov1_utils
  USE cppt,                            ONLY: inyh
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhov1
  PUBLIC :: give_scr_rhov1

CONTAINS

  ! ==================================================================
  SUBROUTINE rhov1(rg,nstate,ist)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE VANDERBILT DENSITY       ==
    ! ==                FOR A SINGLE STATE                            ==
    ! == N_V(G) = SUM_I,IJ RHO_I,IJ Q_I,JI(G) E^-IG.R_I               ==
    ! == RHO_I,IJ = SUM_N < BETA_I,I | PSI_N >< PSI_N | BETA_I,J >    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: rg(ncpw%nhg)
    INTEGER                                  :: nstate, ist

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhov1'

    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:), deltar(:), qg(:)
    INTEGER                                  :: ia, ierr, ig, il_ctmp, &
                                                il_deltar, il_qg, is, isa, &
                                                isa0, isub, iv, jv, length
    REAL(real_8)                             :: sum

    CALL tiset('     RHOV1',isub)
    IF (cntl%tfdist) CALL stopgm('RHOV1','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('RHOV1','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! SCR partition
    il_deltar = 2*ncpw%nhg
    il_qg     = 2*ncpw%nhg
    il_ctmp   = 2*ncpw%nhg
    length    = il_deltar+il_qg+il_ctmp
    ALLOCATE(deltar(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ctmp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(deltar)!,nhg)
    ! FIXME: check out, whether this can be threaded with OpenMP. AK 2005/05/14
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                CALL qvan2(iv,jv,is,qg)
                IF (cntl%bigmem) THEN
                   CALL zeroing(ctmp)!,nhg)
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      sum = fnl(1,isa,iv,ist,1)*fnl(1,isa,jv,ist,1)
                      IF (iv.NE.jv) sum=2._real_8*sum
                      CALL daxpy(2*ncpw%nhg,sum,eigrb(1,isa),1,ctmp(1),1)
                   ENDDO
                   !$omp parallel do private(IG)
                   !CDIR NODEP
                   DO ig=1,ncpw%nhg
                      deltar(ig)=deltar(ig)+qg(ig)*ctmp(ig)
                   ENDDO
                ELSE
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      !$omp parallel do private(IG)
                      DO ig=1,ncpw%nhg
                         ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                      ENDDO
                      sum = fnl(1,isa,iv,ist,1)*fnl(1,isa,jv,ist,1)
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
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(IG)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       rg(ig)=rg(ig)+deltar(ig)
    ENDDO
    CALL tihalt('     RHOV1',isub)
    ! ==--------------------------------------------------------------==
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
    RETURN
  END SUBROUTINE rhov1
  ! ==================================================================
  SUBROUTINE give_scr_rhov1(lrhov1,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhov1
    CHARACTER(len=30)                        :: tag

    lrhov1=6*ncpw%nhg
    tag='6*NHG'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhov1
  ! ==================================================================

END MODULE rhov1_utils
