MODULE prtgyr_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fbohr
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: avcor,&
                                             avgyr,&
                                             avsus,&
                                             fpcor,&
                                             fpgyr,&
                                             fpsus

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prtgyr

CONTAINS

  ! ==================================================================
  SUBROUTINE prtgyr
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: is

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(1X,63("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(" *",19X,A,20X,"*")') ' CHARACTERISTIC RADII '
    IF (paral%io_parent)&
         WRITE(6,'(1X,63("*"))')
    ! ==--------------------------------------------------------------==
    ! ==  INFO ON RADII OF GYRATION PER SPECIES                       ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(A,A,A)') '   SPECIES  R OF GYRATION',&
         '  FREE PARTICLE',&
         '  (IN ANGSTROM)'
    DO is=1,ions1%nsp
       IF (paral%io_parent)&
            WRITE(6,'(I4,2X,A2,3X,F8.5,7X,F8.5)')&
            is,elem%el(ions0%iatyp(is)),&
            avgyr(is)/fbohr,fpgyr(is)/fbohr
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  INFO ON SUSCEPTIBILITY RADII PER SPECIES                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(A)') '            R OF SUSCEPTIBILITY '
    DO is=1,ions1%nsp
       IF (paral%io_parent)&
            WRITE(6,'(I4,2X,A2,3X,F8.5,7X,F8.5)')&
            is,elem%el(ions0%iatyp(is)),&
            avsus(is)/fbohr,fpsus(is)/fbohr
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  INFO ON CORRELATION RADII PER SPECIES                       ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(A)') '            R OF CORRELATION '
    DO is=1,ions1%nsp
       IF (paral%io_parent)&
            WRITE(6,'(I4,2X,A2,3X,F8.5,7X,F8.5)')&
            is,elem%el(ions0%iatyp(is)),&
            avcor(is)/fbohr,fpcor(is)/fbohr
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,63("*"))')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prtgyr
  ! ==================================================================

END MODULE prtgyr_utils
