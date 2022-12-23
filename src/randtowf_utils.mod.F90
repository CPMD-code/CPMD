MODULE randtowf_utils
  USE cppt,                            ONLY: hg
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpnt,                            ONLY: rk
  USE kpts,                            ONLY: tkpts
  USE prng_utils,                      ONLY: repprngu_mat_cmplx
  USE system,                          ONLY: cntr,&
                                             ncpw,&
                                             nkpt,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: randtowf

CONTAINS

  ! ==================================================================
  SUBROUTINE randtowf(c0,nstate,ikind,ikk)
    ! ==--------------------------------------------------------------==
    ! == Build a wavefunction:                                        ==
    ! == Symmetric of hermitian                                       ==
    ! == C0(NGWK,NSTATE) input:  Randomized values                    ==
    ! ==          output: Wavefunctions                               ==
    ! == NSTATE   Number of states                                    ==
    ! == IKIND    Index of kpoints (input) in the set of 1:NKPOINT    ==
    ! == IKK      Index of kpoints (input) in the set of 1:NKPTS      ==
    ! ==          if 0 do not use                                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate, ikind, ikk

    CHARACTER(*), PARAMETER                  :: procedureN = 'randtowf'

    INTEGER                                  :: i, ig, ig1, isub
    LOGICAL                                  :: trk0
    REAL(real_8)                             :: ampexp

!(nkpt%ngwk,nstate)
! Variables
!  real(real_8), pointer :: rc0(:)
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    !  call reshape_inplace(c0, (/nkpt%ngwk * nstate/), rc0)

    CALL repprngu_mat_cmplx(nkpt%ngwk,nstate,c0)
    IF (geq0) THEN
       ! The first component is for G=0.
       ig1=2
    ELSE
       ig1=1
    ENDIF
    IF (.NOT.tkpts%tkpnt.OR.ikk.EQ.0.OR.ikind.EQ.0) THEN
       trk0=.TRUE.
       ! This is evaluated if IF(1) is false.
       ! Do not group (2) with 91) to avoid evaluation of (2) if (1)=.TRUE. 
    ELSEIF((rk(1,ikk).EQ.0._real_8).AND.(rk(2,ikk).EQ.0._real_8).AND.&
         (rk(3,ikk).EQ.0._real_8)) THEN
       trk0=.TRUE.
    ELSE
       trk0=.FALSE.
    ENDIF
    IF (trk0) THEN
       ! Gamma point.
       DO i=1,nstate
          IF (geq0) c0(1,i) = CMPLX(REAL(c0(1,i)),0.0_real_8,kind=real_8)
          DO ig=ig1,ncpw%ngw
             ampexp=cntr%ampre*EXP((-4._real_8*SQRT(hg(ig)))/REAL(spar%ngwks,kind=real_8))
             c0(ig,i) = ampexp*(c0(ig,i)-CMPLX(0.5_real_8,0.5_real_8,kind=real_8))
          ENDDO
       ENDDO
       IF (tkpts%tkpnt) THEN
          DO i=1,nstate
             IF (geq0) c0(ncpw%ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
             DO ig=ig1,ncpw%ngw
                c0(ig+ncpw%ngw,i) = CONJG(c0(ig,i))
             ENDDO
          ENDDO
          IF (ikind.NE.0) CALL c_clean(c0,nstate,ikind)
       ENDIF
    ELSE
       DO i=1,nstate
          IF (geq0) c0(ncpw%ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
          DO ig=ig1,ncpw%ngw
             ampexp=cntr%ampre*EXP((-4._real_8*SQRT(hg(ig)))/REAL(spar%ngwks,kind=real_8))
             c0(ig,    i) = ampexp*(c0(ig,    i)-CMPLX(0.5_real_8,0.5_real_8,kind=real_8))
             c0(ig+ncpw%ngw,i) = ampexp*(c0(ig+ncpw%ngw,i)-CMPLX(0.5_real_8,0.5_real_8,kind=real_8))
          ENDDO
       ENDDO
       ! Clean coefficient G/|G+RK|^2 > Ecutoff if spherical cutoff.
       CALL c_clean(c0,nstate,ikind)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE randtowf

END MODULE randtowf_utils
