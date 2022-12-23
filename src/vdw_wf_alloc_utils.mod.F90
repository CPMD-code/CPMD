MODULE vdw_wf_alloc_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1
  USE pimd,                            ONLY: np_local
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE vdwcmod,                         ONLY: &
       icontfragw, icontfragwi, ifragw, iwfcref, natwfcx, nfrags, nfragx, &
       nwfcx, rwann, rwfc, spr, swann, taufrag, tauref, trwanncx, twannupx, &
       vdwl, vdwwfl, wwfcref

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vdw_wf_alloc

CONTAINS

  ! ==================================================================
  SUBROUTINE vdw_wf_alloc
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'vdw_wf_alloc'

    INTEGER                                  :: ierr, numx

    IF (.NOT.vdwl%vdwd) RETURN

    nwfcx=crge%n
    nfragx=ions1%nat
    IF (cntl%tpath) THEN
       numx=np_local
    ELSE
       numx=1
    ENDIF
    ALLOCATE(rwfc(3,nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(spr(nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taufrag(3,nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tauref(3,maxsys%nax,maxsys%nsx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rwann(3,nwfcx*nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(swann(1*nwfcx*nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(icontfragw(nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(icontfragwi(nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iwfcref(natwfcx,nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ifragw(nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(wwfcref(natwfcx,nwfcx,nfragx,numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nfrags(numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(twannupx(numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    twannupx(:)=vdwwfl%twannup
    ALLOCATE(trwanncx(numx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    trwanncx(:)=vdwwfl%trwannc
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vdw_wf_alloc
  ! ==================================================================

END MODULE vdw_wf_alloc_utils
