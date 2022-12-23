MODULE teststore_utils
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cnti,&
                                             restf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: teststore

CONTAINS

  ! ==================================================================
  FUNCTION teststore(infi)
    ! ==--------------------------------------------------------------==
    ! ==  DO WE HAVE TO WRITE A RESTART FILE                          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: infi
    LOGICAL                                  :: teststore

    INTEGER                                  :: i

! ==--------------------------------------------------------------==

    restf%nrcopy=0
    IF (restf%nstepwr(1).LT.0) THEN
       teststore=(MOD(infi,store1%istore).EQ.0) .OR. (infi.EQ.cnti%nomore)
    ELSE
       teststore=(infi.EQ.restf%nstepwr(restf%nfnow))
       IF (teststore) THEN
          restf%nrcopy=1
          DO i=restf%nfnow+1,cnti%nrestf
             IF (restf%nstepwr(restf%nfnow).NE.restf%nstepwr(i)) GOTO 10
             restf%nrcopy=restf%nrcopy+1
          ENDDO
10        CONTINUE
          restf%nfnow=restf%nfnow+restf%nrcopy
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION teststore
  ! ==================================================================

END MODULE teststore_utils
