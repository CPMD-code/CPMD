MODULE rnlfl_utils
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE system,                          ONLY: cntl,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlfl

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlfl(fion,gamma,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE FORCE ON THE IONIC DEGREES OF FREEDOM DUE TO THE        ==
    ! ==  ORTHOGONALITY CONSTRAINED IN VANDERBILT PSEUDO-POTENTIALS   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: GAMMA(nstate,*)
    INTEGER                                  :: nkpoint

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                isub, iv, j, jv, k
    REAL(real_8)                             :: temp, tk

    CALL tiset('     RNLFL',isub)
    IF (cntl%tfdist) CALL stopgm('FNL_RSPACE','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('RNLFL','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO ik=1,nkpoint
       DO k=1,3
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      temp=qq(iv,jv,is)
                      IF (ABS(temp).GT.1.e-5_real_8) THEN
                         DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            DO j=1,nstate
                               tk=temp*GAMMA(i,j)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  fion(k,ia,is)=fion(k,ia,is)-tk*&
                                       dfnl(1,isa,jv,k,ii,ik)*fnl(1,isa,iv,j,ik)
                               ENDDO
                            ENDDO
                         ENDDO
                         DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            DO j=1,nstate
                               tk=temp*GAMMA(j,i)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  fion(k,ia,is)=fion(k,ia,is)-tk*&
                                       dfnl(1,isa,iv,k,ii,ik)*fnl(1,isa,jv,j,ik)
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt('     RNLFL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlfl
  ! ==================================================================

END MODULE rnlfl_utils
