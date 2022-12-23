MODULE augchg_utils
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: cntl,&
                                             maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: augchg

CONTAINS

  ! ==================================================================
  ! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- AUGCHG_C WILL BE WRITTEN  ==
  ! ==================================================================
  SUBROUTINE augchg(fnl,f,qa,numorb)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: qa(ions1%nat)
    INTEGER                                  :: numorb
    REAL(real_8) :: f(numorb), fnl(ions1%nat,maxsys%nhxs,numorb)

    INTEGER                                  :: i, ia, is, isa, isa0, iv, jv
    REAL(real_8)                             :: sum

    IF (cntl%tfdist) CALL stopgm('AUGCHG','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   sum=0.0_real_8
                   DO i=1,numorb
                      sum=sum+f(i)*fnl(isa,iv,i)*fnl(isa,jv,i)
                   ENDDO
                   IF (iv.NE.jv) sum=2.0_real_8*sum
                   qa(isa)=qa(isa)+sum*qq(iv,jv,is)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE augchg
  ! ==================================================================


END MODULE augchg_utils
