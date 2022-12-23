MODULE rrane_utils
  USE cppt,                            ONLY: hg
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpts,                            ONLY: tkpts
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: cntr,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE utils,                           ONLY: zclean,&
                                             zclean_k

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rrane

CONTAINS

  ! ==================================================================
  SUBROUTINE rrane(c0,cm,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cm(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,*)

    INTEGER                                  :: i, ikind, ikpt, j, kbeg, &
                                                kend, kinc, nkpoint
    REAL(real_8)                             :: ampex, gcutwon

! Variables
!  real(real_8), pointer :: rcm(:)
! ==--------------------------------------------------------------==
! ==  Randomize electronic wavefunctions around loaded ones       ==
! ==--------------------------------------------------------------==

    IF (pslo_com%tivan) CALL stopgm('RRANE','RANDOMIZE WAVEFUNCTION NOT ' //&
         'SUPPORTED FOR VANDERBILT PP',& 
         __LINE__,__FILE__)

    !  call reshape_inplace(cm, (/2*nkpt%ngwk/), rcm)

    gcutwon=4._real_8/gvec_com%gcutw
    IF (tkpts%tkpnt) THEN
       CALL inq_swap(kbeg,kend,kinc)
       DO ikpt=kbeg,kend,kinc
          nkpoint=nkpbl(ikpt)
          IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,'MASKGW C0')
          DO ikind = 1, nkpoint
             DO i=1,nstate
                CALL repprngu_vec_cmplx(nkpt%ngwk,cm)
                DO j=1,ncpw%ngw
                   ampex=cntr%ampre*EXP(-(hg(j))*gcutwon)
                   c0(j,i,ikind)=c0(j,i,ikind)+ampex*cm(j)
                   c0(ncpw%ngw+j,i,ikind)=c0(ncpw%ngw+j,i,ikind)+ampex*cm(ncpw%ngw+j)
                ENDDO
             ENDDO
             IF (geq0) CALL zclean_k(c0(1,1,ikind),nstate,ncpw%ngw)
             CALL c_clean(c0(1,1,ikind),nstate,ikind)
          ENDDO
          IF (tkpts%tkblock) CALL wkpt_swap(c0,nstate,ikpt,'C0')
       ENDDO
    ELSE
       DO i=1,nstate
          CALL repprngu_vec_cmplx(ncpw%ngw,cm)
          DO j=1,ncpw%ngw
             ampex=cntr%ampre*EXP(-(hg(j))*gcutwon)
             c0(j,i,1)=c0(j,i,1)+ampex*cm(j)
          ENDDO
       ENDDO
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rrane
  ! ==================================================================

END MODULE rrane_utils
