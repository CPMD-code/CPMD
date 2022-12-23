MODULE nlsl_utils
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl
  USE str2,                            ONLY: becs
  USE strs,                            ONLY: dlam
  USE system,                          ONLY: cntl,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlsl

CONTAINS

  ! ==================================================================
  SUBROUTINE nlsl(gamma,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==        THE STRESS TENSOR CONTRIBUTION DUE TO THE             ==
    ! ==  ORTHOGONALITY CONSTRAINED IN VANDERBILT PSEUDO-POTENTIALS   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: GAMMA(nstate,*)

    INTEGER                                  :: i, ia, ii, is, isa, isa0, &
                                                isub, iv, j, jv, kk
    REAL(real_8)                             :: temp, tk, tt

    CALL tiset('      NLSL',isub)
    IF (cntl%tfdist) CALL stopgm('NLSL','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('NLSL','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO kk=1,6
       tt=0.0_real_8
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   temp=qq(iv,jv,is)
                   IF (ABS(temp).GT.1.e-5_real_8) THEN
                      IF (iv.NE.jv) temp=2.0_real_8*temp
                      ! This loops are split to get better cash usage
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         DO j=i,nstate
                            tk=temp*GAMMA(i,j)
                            IF (i.EQ.j) tk=0.5_real_8*tk
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               tt=tt +&
                                    tk*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,j,1)
                            ENDDO
                         ENDDO
                      ENDDO
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         DO j=i,nstate
                            tk=temp*GAMMA(j,i)
                            IF (i.EQ.j) tk=0.5_real_8*tk
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               tt=tt +&
                                    tk*becs(1,isa,jv,kk,ii,1)*fnl(1,isa,iv,j,1)
                            ENDDO
                         ENDDO
                      ENDDO
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         DO j=i,nstate
                            tk=temp*GAMMA(i,j)
                            IF (i.EQ.j) tk=0.5_real_8*tk
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               tt=tt +&
                                    tk*becs(1,isa,jv,kk,ii,1)*fnl(1,isa,iv,j,1)
                            ENDDO
                         ENDDO
                      ENDDO
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         DO j=i,nstate
                            tk=temp*GAMMA(j,i)
                            IF (i.EQ.j) tk=0.5_real_8*tk
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               tt=tt +&
                                    tk*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,j,1)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
       dlam(kk) = -tt
    ENDDO
    CALL tihalt('      NLSL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlsl
  ! ==================================================================

END MODULE nlsl_utils
