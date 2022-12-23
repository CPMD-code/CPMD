MODULE exdipo_utils
  USE dipo_utils,                      ONLY: dipo
  USE dipomod,                         ONLY: moment
  USE kinds,                           ONLY: real_8
  USE lodipo_utils,                    ONLY: lodipo
  USE lodp,                            ONLY: dmomlo,&
                                             exd,&
                                             numbld,&
                                             trmom
  USE prop,                            ONLY: prop1
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: exdipo
  PUBLIC :: transm

CONTAINS

  ! ==================================================================
  SUBROUTINE exdipo(ist,tau0,v,eirop)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES EXCITED STATE DIPOLE MOMENT                       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ist
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: v(fpar%nnr1), eirop(ncpw%nhg)

    INTEGER                                  :: i, isub

    CALL tiset('    EXDIPO',isub)
    IF (prop1%locd) THEN
       CALL lodipo(eirop,v)
       DO i=1,numbld
          exd(1,i,ist)=dmomlo(1,i)
          exd(2,i,ist)=dmomlo(2,i)
          exd(3,i,ist)=dmomlo(3,i)
       ENDDO
    ELSE
       CALL dipo(tau0,eirop,v)
       exd(1,1,ist)=moment%dmom(1)
       exd(2,1,ist)=moment%dmom(2)
       exd(3,1,ist)=moment%dmom(3)
    ENDIF
    CALL tihalt('    EXDIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE exdipo
  ! ==================================================================
  SUBROUTINE transm(ist,tau0,v,eirop)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES TRANSITION MOMENT                                 ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ist
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: v(fpar%nnr1), eirop(ncpw%nhg)

    INTEGER                                  :: i, isub

    CALL tiset('    TRANSM',isub)
    IF (prop1%locd) THEN
       CALL lodipo(eirop,v)
       DO i=1,numbld
          trmom(1,i,ist)=dmomlo(1,i)
          trmom(2,i,ist)=dmomlo(2,i)
          trmom(3,i,ist)=dmomlo(3,i)
       ENDDO
    ELSE
       CALL dipo(tau0,eirop,v)
       trmom(1,1,ist)=moment%dmom(1)
       trmom(2,1,ist)=moment%dmom(2)
       trmom(3,1,ist)=moment%dmom(3)
    ENDIF
    CALL tihalt('    TRANSM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE transm
  ! ==================================================================


END MODULE exdipo_utils
