MODULE fillc_utils
  USE cotr,                            ONLY: duat
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fillc

CONTAINS

  SUBROUTINE fillc(iat,tau,x)
    ! ==--------------------------------------------------------------==
    ! == Extract the coordinates X(1:3) of IAT index in TAU (TSCR)    ==
    ! == Works also if dummy atoms                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iat
    REAL(real_8)                             :: tau(3,maxsys%nax*maxsys%nsx), &
                                                x(3)

    INTEGER                                  :: id, ityp, naa

! ==--------------------------------------------------------------==

    IF (iat.LE.0) THEN
       ! Do nothing
       RETURN
    ELSEIF (iat.LE.ions1%nat) THEN
       ! Real atoms
       x(1)=tau(1,iat)
       x(2)=tau(2,iat)
       x(3)=tau(3,iat)
    ELSEIF (iat.LE.ions1%nat+duat%ndat) THEN
       ! Dummy atoms (type 1-4).
       naa=iat-ions1%nat
       ityp=duat%listda(naa,1)
       id=duat%listda(naa,2)
       IF (ityp.EQ.4) THEN
          x(1)=duat%dummy4(1,id)
          x(2)=duat%dummy4(2,id)
          x(3)=duat%dummy4(3,id)
       ELSEIF (ityp.EQ.3) THEN
          x(1)=duat%dummy3(1,id)
          x(2)=duat%dummy3(2,id)
          x(3)=duat%dummy3(3,id)
       ELSEIF (ityp.EQ.2) THEN
          x(1)=duat%dummy2(1,id)
          x(2)=duat%dummy2(2,id)
          x(3)=duat%dummy2(3,id)
       ELSEIF (ityp.EQ.1) THEN
          x(1)=duat%dummy1(1,id)
          x(2)=duat%dummy1(2,id)
          x(3)=duat%dummy1(3,id)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fillc

END MODULE fillc_utils
