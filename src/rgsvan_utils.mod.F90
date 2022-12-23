MODULE rgsvan_utils
  USE csmat_utils,                     ONLY: csmat
  USE kinds,                           ONLY: real_8
  USE rgs_utils,                       ONLY: uinv
  USE sfac,                            ONLY: ldf1
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgsvan

CONTAINS

  ! ==================================================================
  SUBROUTINE rgsvan(c0,fnla,nstate,smat)
    ! ==--------------------------------------------------------------==
    ! ==  Gram-Schmidt orthogonalization for Vanderbilt pp            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8)                             :: fnla(:,:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: smat(:)

    INTEGER                                  :: isub
    LOGICAL                                  :: tlsd_back

!(ncpw%ngw,*)
!(ldf1,maxsys%nhxs,*)
! ==--------------------------------------------------------------==

    IF (nstate.LE.0) RETURN
    CALL tiset('    RGSVAN',isub)
    tlsd_back=cntl%tlsd
    cntl%tlsd=.FALSE.
    CALL csmat(smat,c0,fnla,nstate,1)
    CALL uinv('U',smat,nstate,nstate)
    IF (ncpw%ngw.GT.0)&
         CALL dtrmm('R','U','N','N',2*ncpw%ngw,nstate,1._real_8,smat,nstate,c0,&
         2*ncpw%ngw)
    IF (ldf1.GT.0)&
         CALL dtrmm('R','U','N','N',ldf1*maxsys%nhxs,nstate,1._real_8,smat,nstate,&
         fnla,ldf1*maxsys%nhxs)
    cntl%tlsd=tlsd_back
    CALL tihalt('    RGSVAN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgsvan
  ! ==================================================================

END MODULE rgsvan_utils
