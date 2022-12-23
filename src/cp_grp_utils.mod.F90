MODULE cp_grp_utils
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cp_grp_get_sizes
  PUBLIC :: cp_grp_redist

  INTERFACE cp_grp_redist
     MODULE PROCEDURE cp_grp_redist_d
     MODULE PROCEDURE cp_grp_redist_z
  END INTERFACE cp_grp_redist

  INTEGER,DIMENSION(:,:),ALLOCATABLE,PUBLIC :: cp_grp_get_cp_rank

CONTAINS
  ! ==================================================================
  SUBROUTINE cp_grp_get_sizes(ngwk_l,ngw_l,&
       geq0_l,igeq0_l,firstk_g,lastk_g,first_g,last_g,i_g0_l)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(out), OPTIONAL           :: ngwk_l, ngw_l
    LOGICAL, INTENT(out), OPTIONAL           :: geq0_l
    INTEGER, INTENT(out), OPTIONAL           :: igeq0_l, firstk_g, lastk_g, &
                                                first_g, last_g, i_g0_l

    INTEGER                                  :: first, firstk, i_g0, last, &
                                                lastk, n, nk

    nk = CEILING(REAL(nkpt%ngwk,kind=real_8)/REAL(parai%cp_nogrp,kind=real_8))
    firstk = parai%cp_inter_me*nk + 1
    lastk = MIN((parai%cp_inter_me+1)*nk,nkpt%ngwk)
    ! protect for OOB
    IF (firstk>nkpt%ngwk) THEN
       firstk = 1
       lastk = 0
    ENDIF
    IF (PRESENT(firstk_g)) firstk_g = firstk
    IF (PRESENT(lastk_g)) lastk_g = lastk
    IF (PRESENT(ngwk_l)) ngwk_l = lastk - firstk + 1
    ! 
    ! index of the G=0 element
    i_g0 = 0
    IF (tkpts%tkpnt) THEN
       IF (GEq0) i_g0 = 1 + ncpw%ngw
    ELSE
       IF (GEq0) i_g0 = 1
    ENDIF
    ! 
    ! 
    n = CEILING(REAL(ncpw%ngw,kind=real_8)/REAL(parai%cp_nogrp,kind=real_8))
    first = parai%cp_inter_me*n + 1
    last = MIN((parai%cp_inter_me+1)*n,ncpw%ngw)
    ! protect for OOB
    IF (first>ncpw%ngw) THEN
       first = 1
       last = 0
    ENDIF
    IF (PRESENT(first_g)) first_g = first
    IF (PRESENT(last_g)) last_g = last
    IF (PRESENT(ngw_l)) ngw_l = last - first + 1
    ! 
    ! 
    IF (PRESENT(geq0_l)) THEN
       geq0_l = i_g0.GE.firstk.AND.i_g0.LE.lastk
    ENDIF

    IF (PRESENT(igeq0_l)) THEN
       igeq0_l = 0
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) igeq0_l = parai%cp_me
       CALL mp_sum(igeq0_l,parai%cp_grp)
    ENDIF

    IF (PRESENT(i_g0_l)) THEN
       i_g0_l = -HUGE(0)
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) i_g0_l = i_g0 - firstk + &
            1
    ENDIF

       ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_get_sizes
  ! ==================================================================
  SUBROUTINE cp_grp_redist_z(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_d(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_d

  ! ==================================================================
END MODULE cp_grp_utils
! ==================================================================

SUBROUTINE cp_grp_zero_g(c0,n,m,first_g,last_g)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: first_g, last_g

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: i, j

  IF (first_g.EQ.1.AND.last_g.EQ.n) RETURN
  !$omp parallel do &
  !$omp            private(i) &
  !$omp            shared(first_g,last_g,N)
  DO j = 1,m
     DO i = 1,first_g-1
        c0(i,j) = zzero
     ENDDO
     DO i = last_g+1,n
        c0(i,j) = zzero
     ENDDO
  ENDDO
  !$omp end parallel do
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_g
! ==================================================================

SUBROUTINE cp_grp_zero_states(c0,n,m,ibeg,iend)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac , ONLY:parai
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: ibeg(0:parai%cp_nproc-1), &
                                                iend(0:parai%cp_nproc-1)

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: nn

  IF (ibeg(parai%cp_me).GT.1) THEN
     nn = n * ( ibeg(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,1),1)
  ENDIF
  IF (iend(parai%cp_me).LT.m) THEN
     nn = n * ( m - iend(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,iend(parai%cp_me)+1),1)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_states
! ==================================================================
SUBROUTINE cp_grp_copy_wfn_to_local(c0,ld,C0_l,LD_l,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy the right C0 block into the local array
     DO i=1,n
        CALL dcopy(2*M_l,c0(ibeg,i),1,C0_l(1,i),1)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_wfn_to_local
! ==================================================================
SUBROUTINE cp_grp_copy_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     DO i=1,n
        CALL dcopy(2*M_l,C0_l(1,i),1,c0(ibeg,i),1)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_local_to_wfn
! ==================================================================
SUBROUTINE cp_grp_add_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i, j

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     !$omp parallel do private(i,j) shared(N,M_l,ibeg)
     DO i=1,n
        DO j=1,M_l
           c0(ibeg+j-1,i) = c0(ibeg+j-1,i) + C0_l(j,i)
        ENDDO
     ENDDO
     !$omp end parallel do
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_add_local_to_wfn
