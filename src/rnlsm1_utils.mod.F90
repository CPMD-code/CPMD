#include "cpmd_global.h"

MODULE rnlsm1_utils
  USE cppt,                            ONLY: twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: tkpts
  USE mm_dimmod,                       ONLY: mmdim
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE sfac,                            ONLY: fnl
  USE sumfnl_utils,                    ONLY: give_scr_sumfnl,&
                                             sumfnl
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm1
  PUBLIC :: give_scr_rnlsm1

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm1(c0,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY FNL (ALSO CALLED BEC IN SOME VANDERBILT ROUTINES) ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm1'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: scr(:,:)
    INTEGER                                  :: first_state, ia, ierr, ig, &
                                                is, isa0, isub, iv, &
                                                last_state, n_state
    REAL(real_8)                             :: cii, cir, ei, er, t

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )
    ! SCR test.

    ALLOCATE(scr(nkpt%ngwk, mmdim%naxq), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate SCR',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(fnl(:,:,:,:,ikind))!,imagp*ions1%nat*maxsys%nhxs*nstate)
    IF (nkpt%ngwk.EQ.0) GOTO 1000
    ! ==--------------------------------------------------------------==
    CALL part_1d_get_blk_bounds(nstate,parai%cp_inter_me,parai%cp_nogrp,&
         first_state,last_state)
    n_state = last_state - first_state + 1
    ! ==--------------------------------------------------------------==
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
          cir=REAL(ci)
          cii=AIMAG(ci)
          ! Make use of the special structure of CI
          IF (ABS(cir).GT.0.5_real_8) THEN
             ! CI is real
#ifdef __SR8000 
             !poption parallel
             ! soption unroll(2)
#endif 
             !$OMP parallel do private(IA,IG,ER,EI,T) __COLLAPSE2
             DO ia=1,ions0%na(is)
                !mb              isa=isa0+ia
#ifdef __SR8000
                ! soption unroll(8)
#endif 
                DO ig=1,nkpt%ngwk
                   !mb              er=dreal(eigkr(ig,isa,ikind))
                   !mb              ei=dimag(eigkr(ig,isa,ikind))
                   er=REAL(eigkr(ig,isa0+ia,ikind),kind=real_8)
                   ei=AIMAG(eigkr(ig,isa0+ia,ikind))
                   t=twnl(ig,iv,is,ikind)*cir
                   scr(ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                ENDDO
             ENDDO
          ELSE
             ! CI is imaginary
#ifdef __SR8000 
             !poption parallel
             ! soption unroll(2)
#endif 
             !$OMP parallel do private(IA,IG,ER,EI,T) __COLLAPSE2
             DO ia=1,ions0%na(is)
                !mb            isa=isa0+ia
#ifdef __SR8000
                ! soption unroll(8)
#endif 
                DO ig=1,nkpt%ngwk
                   !mb              er=dreal(eigkr(ig,isa,ikind))
                   !mb              ei=dimag(eigkr(ig,isa,ikind))
                   er=REAL(eigkr(ig,isa0+ia,ikind),kind=real_8)
                   ei=AIMAG(eigkr(ig,isa0+ia,ikind))
                   t=twnl(ig,iv,is,ikind)*cii
                   scr(ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                ENDDO
             ENDDO
          ENDIF
          IF (geq0) THEN
             IF (tkpts%tkpnt) THEN
#ifdef __SR8000
                !poption parallel
#endif
                DO ia=1,ions0%na(is)
                   scr(ncpw%ngw+1,ia)=CMPLX(0._real_8,0._real_8,kind=real_8)
                ENDDO
             ELSE
#ifdef __SR8000
                !poption parallel
#endif
                DO ia=1,ions0%na(is)
                   scr(1,ia)=0.5_real_8*scr(1,ia)
                ENDDO
             ENDIF
          ENDIF
          ! 
          ! This is the part that cost the most
          ! 
          IF (tkpts%tkpnt) THEN
             CALL zgemm('C','N',ions0%na(is),n_state,nkpt%ngwk,zone,scr(1,1),&
                  nkpt%ngwk,c0(1,first_state),nkpt%ngwk,zzero,&
                  fnl(1,isa0+1,iv,first_state,ikind),ions1%nat*maxsys%nhxs)
          ELSE
             IF (ions0%na(is).GT.1) THEN
                CALL dgemm('T','N',ions0%na(is),n_state,2*nkpt%ngwk,2._real_8,scr(1,1),&
                     2*nkpt%ngwk,c0(1,first_state),2*nkpt%ngwk,0.0_real_8,&
                     fnl(1,isa0+1,iv,first_state,ikind),ions1%nat*maxsys%nhxs)
             ELSE
                CALL dgemv('T',2*nkpt%ngwk,n_state,2._real_8,c0(1,first_state),&
                     2*nkpt%ngwk,scr(1,1),1,0.0_real_8,&
                     fnl(1,isa0+1,iv,first_state,ikind),ions1%nat*maxsys%nhxs)
             ENDIF
          ENDIF

       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO

    ! ==--------------------------------------------------------------==
    IF (parai%cp_nogrp.GT.1) THEN

       CALL mp_sum(fnl(:,:,:,:,ikind),imagp*ions1%nat*maxsys%nhxs*nstate,&
            parai%cp_inter_grp)
    ENDIF
    ! ==--------------------------------------------------------------==      

1000 CONTINUE

    CALL sumfnl(fnl(:,:,:,:,ikind),nstate)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

  END SUBROUTINE rnlsm1
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm1(lrnlsm1,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm1
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lsumfnl

    IF (nlm.EQ.0) THEN
       lrnlsm1=0
    ELSE
       IF (cntl%tfdist) THEN
          lrnlsm1=2*nkpt%ngwk*mmdim%naxq+2
          tag   ='2*NGWK*NAXq+2'
       ELSE
          CALL give_scr_sumfnl(lsumfnl,tag,nstate)
          lrnlsm1=MAX(2*nkpt%ngwk*mmdim%naxq,lsumfnl)
          tag   ='MAX(2*NGWK*NAXq,LSUMFNL)'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm1
  ! ==================================================================

END MODULE rnlsm1_utils
