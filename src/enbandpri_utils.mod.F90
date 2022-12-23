MODULE enbandpri_utils
  USE cnst,                            ONLY: ry
  USE elct,                            ONLY: crge
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE parac,                           ONLY: paral
  USE system,                          ONLY: parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: enbandpri

CONTAINS

  ! ==================================================================
  SUBROUTINE enbandpri(eigen,focc,efermi,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==  PRINT THE ENERGY BANDS IN THE FILE ENERGYBANDS              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: efermi
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: eigen(nstate,nkpoint), &
                                                focc(nstate,nkpoint)

    INTEGER                                  :: ii, ikind, istate, isub, iunit
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: hartree

! ==--------------------------------------------------------------==

    CALL tiset (' ENBANDPRI',isub)
    hartree=2._real_8*ry
    iunit=30
    IF (paral%io_parent)&
         CALL fileopen(iunit,'ENERGYBANDS',fo_def,ferror)
    IF (paral%io_parent)&
         WRITE(iunit,'(1X,A8,I12,F20.15)') 'SYMMETRY', parm%ibrav,parm%alat
    IF (paral%io_parent)&
         WRITE(iunit,'(1X,2I12,F20.15,I12)')&
         nkpoint,nstate,efermi*hartree,NINT(crge%nel)
    DO ikind=1,nkpoint
       IF (paral%io_parent)&
            WRITE(iunit,'(I5,4(F20.15))')&
            ikind,&
            (rk(ii,ikind),ii=1,3),&
            wk(ikind)
       DO istate=1,nstate
          IF (paral%io_parent)&
               WRITE(iunit,'(5X,I5,2F20.15)')&
               istate,&
               eigen(istate,ikind)*hartree,&
               focc(istate,ikind)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(iunit)
    CALL tihalt(' ENBANDPRI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE enbandpri
  ! ==================================================================

END MODULE enbandpri_utils
