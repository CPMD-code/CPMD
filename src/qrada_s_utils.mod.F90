MODULE qrada_s_utils
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: gl
  USE dylmr_utils,                     ONLY: dbess
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggnh,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo,&
                                             qspl1
  USE radin_utils,                     ONLY: radin,&
                                             radlg
  USE str2,                            ONLY: qrada
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdbp,                            ONLY: ncpr1,&
                                             qfunc,&
                                             qrl,&
                                             r,&
                                             rab
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: qrada_s
  !public :: forqrad

CONTAINS

  ! ==================================================================
  SUBROUTINE qrada_s(jl,fint)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE INITALIZES ARRAYS   BETAA QRADA                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: jl(*), fint(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'qrada_s'

    INTEGER                                  :: ierr, il, ir, is, istep, &
                                                isub, iv, jv, l, lqx, lval, &
                                                mmax, nt
    REAL(real_8)                             :: fqrad, ggl, xg
    REAL(real_8), ALLOCATABLE                :: qsp1(:), qsp2(:), qsp3(:)

    CALL tiset('   QRADA_S',isub)
    ALLOCATE(qsp1(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qsp2(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qsp3(4*nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          lval=ncpr1%nvales(is)
          nt=ncpr1%nvales(is)
          mmax=ncpr1%kkbeta(is)
          istep=nt*nt
          ! ==--------------------------------------------------------------==
          ! ==  CALCULATION OF ARRAY  QRADA(IGL,IV,JV,IS)
          ! ==--------------------------------------------------------------==
          lqx=2*(lval-1)+1
          IF (qspl1%qspline) THEN
             ! ..QRADA is calculated on the fly from the spline function for qrad
          ELSE
             ! ..This is not fully correct! Qrada is calculated from a spline in 
             ! G space, qrada should be calculated as the derivative of the same 
             ! function. But this is only for static calculations (Parrinello-Rahman
             ! runs are always with QSPLINE=.TRUE.) and the introduced error is
             ! very small
             CALL dcopy(nsplpo,ggnh,1,qsp1,1)
             DO l=1,lqx
                DO iv=1,ncpr1%nbeta(is)
                   DO jv=iv,ncpr1%nbeta(is)
                      CALL zeroing(qsp2)!,nsplpo)
                      DO il=nsplpa,nsplpe
                         ggl=ggnh(il)
                         xg=SQRT(ggl)*parm%tpiba
                         CALL dbess(xg,l,mmax,r(1,is),jl)
                         IF (pslo_com%tpseu(is)) THEN
                            DO ir=1,mmax
                               fint(ir)=qrl(ir,iv,jv,l,is)*jl(ir)
                            ENDDO
                         ELSE
                            DO ir=1,mmax
                               fint(ir)=qfunc(ir,iv,jv,is)*jl(ir)
                            ENDDO
                         ENDIF
                         IF (pslo_com%tlog(is)) THEN
                            CALL radlg(mmax,fint,rab(1,is),fqrad)
                         ELSE
                            CALL radin(mmax,ncpr1%cmesh(is),fint,fqrad)
                         ENDIF
                         qsp2(il)=fpi*fqrad
                      ENDDO
                      CALL mp_sum(qsp2,nsplpo,parai%allgrp)
                      CALL curv1(nsplpo,qsp1(1),qsp2(1),0.0_real_8,0.0_real_8,3,&
                           qsp3(1),fint,0.0_real_8,ierr)
                      DO il=1,ncpw%nhgl
                         ggl=gl(il)
                         qrada(il,iv,jv,l,is)=curv2(ggl,nsplpo,qsp1(1),&
                              qsp2(1),qsp3(1),0.0_real_8)
                         qrada(il,jv,iv,l,is)=qrada(il,iv,jv,l,is)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    DEALLOCATE(qsp1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qsp2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qsp3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('   QRADA_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE qrada_s
  ! ==================================================================

END MODULE qrada_s_utils
