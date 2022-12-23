MODULE response_p_utils
  USE coor,                            ONLY: tau0,&
                                             taup
  USE do_perturbation_p_utils,         ONLY: do_perturbation
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai
  USE response_pmod,                   ONLY: nmr_para
  USE ropt,                            ONLY: iteropt
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: maxsys,&
                                             ncpw
  USE testex_utils,                    ONLY: testex
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: response_p

CONTAINS

  ! ==================================================================
  SUBROUTINE response_p
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'response_p'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), cs(:,:)
    INTEGER                                  :: filenum, ierr, irec(100), &
                                                msglen, nstate
    REAL(real_8), ALLOCATABLE                :: eigv(:)

! ==--------------------------------------------------------------==

    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (soft_com%exsoft) RETURN
    nstate=crge%n
    IF (tkpts%tkpnt) CALL stopgm('RESPONSE_P','k-points not implemented',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOC MEMORY                                                ==
    ! ==--------------------------------------------------------------==
    CALL fnlalloc(nstate,.FALSE.,.FALSE.)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   taup)!,3*maxsys%nax*maxsys%nsx)
    ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cs)!,SIZE(cs))
    ALLOCATE(c0(ncpw%ngw,nstate,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0)!,SIZE(c0))
    ALLOCATE(eigv(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   eigv)!,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  LOAD DATA                                                   ==
    ! ==--------------------------------------------------------------==
    msglen = 20
    CALL mp_bcast(filn,nmr_para%nmr_supersource,&
         nmr_para%nmr_supergroup)
    CALL mp_bcast(filbod,nmr_para%nmr_supersource,&
         nmr_para%nmr_supergroup)

    filenum = 42+nmr_para%nmr_mygroup
    IF (restart1%restart) THEN
       CALL read_irec(irec)
       CALL zhrwf(filenum,irec,c0,cs,nstate,eigv,tau0,taup,taup,iteropt%nfi)
       CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL do_perturbation(c0,cs,nstate)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.FALSE.,.FALSE.)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE response_p
  ! ==================================================================

END MODULE response_p_utils
