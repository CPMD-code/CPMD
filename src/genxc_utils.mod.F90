MODULE genxc_utils
  USE error_handling,                  ONLY: stopgm
  USE functionals_utils,               ONLY: xc
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE tbxc,                            ONLY: eexc,&
                                             tabx,&
                                             toldcode,&
                                             vvxc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: genxc

CONTAINS

  ! ==================================================================
  SUBROUTINE genxc
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'genxc'

    INTEGER                                  :: i, ierr
    REAL(real_8)                             :: ec, ex, rr, vc, vx

! ==--------------------------------------------------------------==

    IF (toldcode) THEN
       IF (tabx%narray.GT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,10) tabx%rmaxxc,tabx%narray
10        FORMAT(' TABLE FOR EXCHANGE-CORRELATION ENERGY',/,&
               '    RMAX:  ',t56,1pe10.4,/,&
               '    NARRAY:',t61,i5,/)
          ! ..Allocate memory for XC-tables

          ! these arrays must be allocated as 0:
          ! TODO implement via memory90
          ALLOCATE(eexc(0:tabx%narray),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(vvxc(0:tabx%narray),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          tabx%ddro=tabx%rmaxxc/REAL(tabx%narray-1,kind=real_8)
          DO i=0,tabx%narray
             rr=REAL(i,kind=real_8)*tabx%ddro
             CALL xc(rr,ex,ec,vx,vc)
             eexc(i)=ex+ec
             vvxc(i)=vx+vc
          ENDDO
       ENDIF
       IF (paral%parent) CALL prmem('     GENXC')
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE genxc
  ! ==================================================================

END MODULE genxc_utils
