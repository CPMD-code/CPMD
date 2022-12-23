MODULE wrccfl_utils
  USE bsym,                            ONLY: autocm,&
                                             rtsasb
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: infi,&
                                             iteropt
  USE system,                          ONLY: cnti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wrccfl

CONTAINS

  ! ==================================================================
  SUBROUTINE wrccfl(etotls,etotbs,etoths)
    ! CB: Writes BS dynamics results into a file      
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: etotls, etotbs, etoths

    INTEGER, PARAMETER                       :: iunit = 33 

    CHARACTER(len=100)                       :: fname
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: couplj, jwn

! ==--------------------------------------------------------------==

    fname='BS_ENERG'
    IF (paral%io_parent) CALL fileopen(iunit,fname,fo_app,ferror)
    ! ==--------------------------------------------------------------==
    couplj = (etoths - etotbs) * rtsasb
    jwn=couplj*autocm
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       IF ((infi.LE.1).OR.(infi.EQ.cnti%nomore)) THEN
          WRITE(6,'(/,A47,F17.8,A7)')&
               'The projected energy of the low spin state is:',&
               etotls, 'a.u.'
          WRITE(6,'(A47,F17.8,A7)')&
               'The coupling constant J is                   :',&
               couplj, 'a.u.'
          WRITE(6,'(A47,F17.3,A7)') '~', jwn, 'cm^-1.'
       ELSE
          WRITE(6,'(A30,F17.8,A7,/)')&
               'The coupling constant J is ~ :',&
               jwn, 'cm^-1.'
       ENDIF
       WRITE(iunit,'(I10,F18.10,F18.10,F18.10,F18.10,F9.2)')&
            iteropt%nfi,etotls,etotbs,etoths,couplj,jwn
    END IF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CALL fileclose(iunit)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE wrccfl

END MODULE wrccfl_utils
