MODULE gs_disortho_utils
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gs_disortho

CONTAINS

  ! ==================================================================
  SUBROUTINE gs_disortho(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
    ! ==     USING A DISTRIBUTED OVERLAP MATRIX ALGORITHM             ==
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gs_disortho'

    INTEGER                                  :: ierr, index0, index1, isub
    REAL(real_8)                             :: local_nrm, nrm
    REAL(real_8), ALLOCATABLE                :: local_opw(:), opw(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('GSDISORTHO',isub)
    ! 
    ! ==================================================================
    ALLOCATE(opw(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(local_opw(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! NORMALIZE FIRST VECTOR
    local_nrm = 2 * ddot(2*nkpt%ngwk, c0(1,1), 1, c0(1,1), 1)
    IF (geq0) THEN
       local_nrm = local_nrm - REAL(c0(1,1))**2
    ENDIF
    CALL mp_sum(local_nrm,nrm,parai%allgrp)
    nrm = SQRT(nrm)
    CALL dscal(2*nkpt%ngwk, 1/nrm, c0(1,1), 1 )
    DO index0 = 2, nstate
       DO index1=1,index0-1
          local_opw(index1) = 2* ddot(2*nkpt%ngwk, c0(1,index1),1, c0(1,&
               index0), 1)
          IF (geq0) THEN
             local_opw(index1) = local_opw(index1) -REAL(c0(1,index1))&
                  * REAL(c0(1,index0))
          ENDIF
       ENDDO
       CALL mp_sum(local_opw,opw,index0-1,parai%allgrp)
       DO index1=1, index0-1
          CALL daxpy(2*nkpt%ngwk, -opw(index1),c0(1,index1), 1, c0(1,&
               index0), 1)
       ENDDO
       ! NORMALIZE TO ONE
       local_nrm = 2 * ddot(2*nkpt%ngwk,c0(1,index0),1,c0(1,index0),1)
       IF (geq0) THEN
          local_nrm = local_nrm -REAL(c0(1,index0))* REAL(c0(1,index0)&
               )
       ENDIF
       CALL mp_sum(local_nrm,nrm,parai%allgrp)
       nrm = SQRT(nrm)
       CALL dscal(2*nkpt%ngwk, 1/nrm, c0(1,index0), 1 )
    ENDDO
    DEALLOCATE(opw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(local_opw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('GSDISORTHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gs_disortho
  ! ==================================================================

END MODULE gs_disortho_utils
