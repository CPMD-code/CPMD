MODULE nlccset_utils
  USE bessm_utils,                     ONLY: bessov
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE cppt,                            ONLY: hg
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlcc,                            ONLY: corecg,&
                                             corei,&
                                             corel,&
                                             corer,&
                                             rcgrid,&
                                             rhoc,&
                                             rhocspl
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggnh,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo
  USE radin_utils,                     ONLY: radin,&
                                             radlg
  USE system,                          ONLY: maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdbp,                            ONLY: ncpr1,&
                                             r,&
                                             rab,&
                                             rscore
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlccset

CONTAINS

  ! ==================================================================
  SUBROUTINE nlccset
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==           THE FORM FACTORS OF THE CORE CHARGES               ==
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'nlccset'

    INTEGER                                  :: ierr, ig, il, ir, is, isub, m
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: arg, cpi, tmp, vol, xg
    REAL(real_8), ALLOCATABLE                :: work(:)

! ==--------------------------------------------------------------==

    CALL tiset('   NLCCSET',isub)
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       ALLOCATE(rhocspl(nsplpo,2,ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(rhocspl)!,nsplpo*2*ions1%nsp)
       ifirst=1
       ALLOCATE(work(maxsys%mmaxx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
       DO is=1,ions1%nsp
          IF (corel%tnlcc(is)) THEN
             IF (pslo_com%tvan(is)) THEN
                ! ..VDBT
                work(1) = 0._real_8
                DO ir=2,ncpr1%meshva(is)
                   work(ir)=rscore(ir,is)
                ENDDO
                IF (pslo_com%tlog(is)) THEN
                   CALL radlg(ncpr1%meshva(is),work,rab(1,is),tmp)
                ELSE
                   CALL radin(ncpr1%meshva(is),ncpr1%cmesh(is),work,tmp)
                ENDIF
                DO il=nsplpa,nsplpe
                   IF (il.NE.1) THEN
                      xg=SQRT(ggnh(il))*parm%tpiba
                      work(1) = 0._real_8
                      DO ir=2,ncpr1%meshva(is)
                         work(ir)=rscore(ir,is)*SIN(xg*r(ir,is))/r(ir,is)
                      ENDDO
                      IF (pslo_com%tlog(is)) THEN
                         CALL radlg(ncpr1%meshva(is),work,rab(1,is),tmp)
                      ELSE
                         CALL radin(ncpr1%meshva(is),ncpr1%cmesh(is),work,tmp)
                      ENDIF
                      rhocspl(il,1,is) = tmp/xg
                   ELSE
                      rhocspl(1,1,is) = tmp
                   ENDIF
                ENDDO
                CALL mp_sum(rhocspl(:,:,is),nsplpo,parai%allgrp)
                CALL curv1(nsplpo,ggnh(1),rhocspl(1,1,is),0._real_8,0._real_8,3,&
                     rhocspl(1,2,is),work,0._real_8,ierr)
             ELSEIF (corei%nlcct(is).EQ.2) THEN
                ! NUMERICAL FORM
                tmp = 0._real_8
                DO il = corei%meshcc(is),1,-1
                   tmp=tmp+ABS(corecg(il,is))
                   IF (tmp.GT.1.e-8_real_8) THEN
                      m = il
                      GOTO 100
                   ENDIF
                ENDDO
                m=corei%meshcc(is)
100             CONTINUE
                DO il=nsplpa,nsplpe
                   xg=SQRT(ggnh(il))*parm%tpiba
                   CALL bessov(rcgrid(1,is),corer%clogcc(is),m,corecg(1,is),0,xg,&
                        rcgrid(m,is),tmp)
                   rhocspl(il,1,is)=tmp*fpi
                ENDDO
                CALL mp_sum(rhocspl(:,:,is),nsplpo,parai%allgrp)
                CALL curv1(nsplpo,ggnh(1),rhocspl(1,1,is),0._real_8,0._real_8,3,&
                     rhocspl(1,2,is),work,0._real_8,ierr)
             ELSEIF (corei%nlcct(is).EQ.1) THEN
                ! ANALYTICAL FORM
                cpi=(pi/corer%enlcc(is))**1.5_real_8
                !$omp parallel do private(IG,ARG)
#ifdef __SR8000
                !poption parallel, tlocal(IG,ARG)
#endif 
                DO ig=1,nsplpo
                   arg=parm%tpiba2*ggnh(ig)/corer%enlcc(is)/4.0_real_8
                   rhocspl(ig,1,is)=(corer%anlcc(is)+corer%bnlcc(is)/corer%enlcc(is)*(1.5_real_8-&
                        arg))*EXP(-arg)*cpi
                ENDDO
             ELSE
                WRITE(6,*) ' UNKNOWN TYPE OF NLCC NLCCT(IS)=',corei%nlcct(is)
                CALL stopgm('NLCCSET',' ',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDDO
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    vol=1._real_8/parm%omega
    DO is=1,ions1%nsp
       IF (corel%tnlcc(is)) THEN
          !$omp parallel do private(IG) shared(NSPLPO)
#ifdef __SR8000
          !poption parallel, tlocal(IG)
#endif 
          DO ig=1,ncpw%nhg
             rhoc(ig,is)=vol*curv2(hg(ig),nsplpo,ggnh(1),rhocspl(1,1,is),&
                  rhocspl(1,2,is),0._real_8)
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL tihalt('   NLCCSET',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlccset
  ! ==================================================================

END MODULE nlccset_utils
