MODULE td_os_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: td01,&
                                             td03
  USE opeigr_utils,                    ONLY: give_scr_opeigr,&
                                             opeigr
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE symm,                            ONLY: symmi,&
                                             symmt
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE td_os_berry_utils,               ONLY: td_os_berry

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_os

CONTAINS

  ! ======================================================================
  SUBROUTINE td_os(nstate, c0, nlr, eigv, c1)
    ! ======================================================================
    ! These from raman_p:
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw, nstate)
    INTEGER                                  :: nlr
    REAL(real_8)                             :: eigv(nlr)
    COMPLEX(real_8)                          :: c1(ncpw%ngw, nstate, nlr)

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_os'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ierr, ilrv, length
    REAL(real_8), ALLOCATABLE                :: scros(:)

! Input:
! n. of occupied states
! occupied KS orbitals
! n. of linear-response vectors
! transition energies
! linear-response orbitals
! -------------------------------------------------------------
! Sometimes the scratch space SCR defined in specpt is not big
! enough for opeigr. A new scratch array SCROS is defined here
! if needed.
! -------------------------------------------------------------
! Local:
! LB   ----------------------------------------------------------------------

    IF (.NOT.td03%tda) RETURN
    IF (cntl%tlsd) RETURN
    IF (td01%ns_tri.NE.0) RETURN
    IF (td01%msubta.GT.0) RETURN

    ! ==--------------------------------------------------------------==
    IF (symmt%tpgauto .OR. symmt%tmsym .OR. symmi%indpg.NE.0) THEN
       CALL stopgm('td_os','POINT GROUP symmetry not implemented.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Allocate SCROS
    length=0
    CALL give_scr_opeigr(length,tag,nstate)

    ! disabling scratch usage
    ALLOCATE(scros(length),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) CALL prmem('     TD_OS')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"LB:",61("-"))')
    DO ilrv = 1, nlr
       CALL td_os_berry(c0,c1(1,1,ilrv),eigv(ilrv),nstate)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,"LB:",61("-"))')

    DEALLOCATE(scros,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE td_os
  ! ======================================================================

END MODULE td_os_utils
