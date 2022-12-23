MODULE formf_utils
  USE atom,                            ONLY: atom_common,&
                                             patom1,&
                                             patom2,&
                                             rv,&
                                             rw,&
                                             vr
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggnh,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo,&
                                             voo
  USE radin_utils,                     ONLY: radin,&
                                             radlg
  USE ragg,                            ONLY: raggio
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE special_functions,               ONLY: cp_erf
  USE system,                          ONLY: maxsys,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: simpsn
  USE vdbp,                            ONLY: ncpr1,&
                                             r,&
                                             rab,&
                                             rucore
  USE zeroing_utils,                   ONLY: zeroing

  !$ USE omp_lib, ONLY : omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: formfn
  PUBLIC :: formsg
  PUBLIC :: formfv
  PUBLIC :: formup
  !public :: sinc

CONTAINS

  ! ==================================================================
  SUBROUTINE formfn(is,dfint,fint)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE FORM FACTORS OF PSEUDOPOTENTIAL (VPS)          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: dfint(maxsys%mmaxx,*), &
                                                fint(maxsys%mmaxx)

    INTEGER                                  :: ierr, il, ir, isub, mepe
    LOGICAL                                  :: zer
    REAL(real_8)                             :: check, flx, gg, r2, rcpf, xg, &
                                                xrg

! ==--------------------------------------------------------------==

    CALL tiset('    FORMFN',isub)
    CALL zeroing(fint)!,maxsys%mmaxx)
    DO ir=1,atom_common%meshvp(is)
       fint(ir)=(rv(ir,is)*rv(ir,is)*rv(ir,is)*vr(ir,is,dpot_mod%lloc(is))+&
            ions0%zv(is)*cp_erf(rv(ir,is)/raggio(is))*rv(ir,is)*rv(ir,is))
       zer = ABS(rv(ir,is)).LT.1.e-8_real_8
       IF (.NOT.zer) THEN
          check = vr(ir,is,dpot_mod%lloc(is))+ions0%zv(is)*&
               cp_erf(rv(ir,is)/raggio(is))/rv(ir,is)
       ELSE
          check = vr(ir,is,dpot_mod%lloc(is))+2._real_8*&
               ions0%zv(is)/(SQRT(pi)*raggio(is))
       ENDIF
       IF (ABS(check).LT.1.e-8_real_8) fint(ir)=0._real_8
    ENDDO
    CALL zeroing(voo(:,1,is))!,nsplpo)
    mepe=1
    !$omp parallel do private(IL,XG,MEPE,IR,XRG,FLX)
    DO il=nsplpa,nsplpe
       !$  MEPE=OMP_GET_THREAD_NUM()+1
       xg=SQRT(ggnh(il))*parm%tpiba
       IF (xg.GT.1.e-6_real_8) THEN
          DO ir=1,atom_common%meshvp(is)
             xrg = rv(ir,is)*xg
             dfint (ir,mepe) = fint(ir)* SIN(xrg)/xrg
          ENDDO
       ELSE
          CALL dcopy(atom_common%meshvp(is),fint(1),1,dfint(1,mepe),1)
       ENDIF
       CALL simpsn (atom_common%meshvp(is), dfint(1,mepe), flx)
       voo(il,1,is) = atom_common%clogvp(is)*fpi* flx
    ENDDO
    IF (patom1%pconf) THEN
       r2=patom2%prc(is)**2
       rcpf=SQRT(pi**3)*patom2%prc(is)**3*patom2%palpha(is)
       !$omp parallel do private(IL,GG)
       DO il=nsplpa,nsplpe
          gg=ggnh(il)*parm%tpiba2
          IF (gg.GT.1.e-8_real_8) THEN
             voo(il,1,is)=voo(il,1,is)-rcpf*SQRT(gg)*EXP(-0.25_real_8*gg*r2)
          ENDIF
       ENDDO
    ENDIF
    CALL mp_sum(voo(:,:,is),nsplpo,parai%allgrp)
    CALL curv1(nsplpo,ggnh,voo(1,1,is),0.0_real_8,0.0_real_8,3,&
         voo(1,2,is),fint,0.0_real_8,ierr)
    CALL tihalt('    FORMFN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE formfn
  ! ==================================================================
  SUBROUTINE formsg(is,rs1,rs2)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE FORM FACTORS OF PSEUDOPOTENTIAL (VPS)          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: rs1(*), rs2(*)

    INTEGER                                  :: ic, ierr, il, isub
    REAL(real_8)                             :: gg, ggl, r2, rcpf, rcrc

    CALL tiset('    FORMSG',isub)
    rcrc=sgpp2%rcsg(is)*sgpp2%rcsg(is)
    r2=raggio(is)**2
    CALL zeroing(voo(:,1,is))!,nsplpo)
    !$omp parallel do private(IL,GG)
    DO il=nsplpa,nsplpe
       gg=ggnh(il)*parm%tpiba2
       IF (gg.LT.1.e-8_real_8) THEN
          voo(il,1,is) = pi*ions0%zv(is)*(2._real_8*rcrc-r2)
       ELSE
          voo(il,1,is) = -fpi*ions0%zv(is)/gg*&
               (EXP(-0.5_real_8*gg*rcrc)-EXP(-0.25_real_8*gg*r2))
       ENDIF
    ENDDO
    rcpf=SQRT((2._real_8*pi)**3)*sgpp2%rcsg(is)**3
    !$omp parallel do private(IL,GG)
    DO il=nsplpa,nsplpe
       gg=ggnh(il)*parm%tpiba2
       rs1(il) = rcpf*EXP(-0.5_real_8*gg*rcrc)
    ENDDO
    DO ic=1,sgpp1%nclsg(is)
       IF (ic.EQ.1) THEN
          !$omp parallel do private(IL,GGL)
          DO il=nsplpa,nsplpe
             ggl=ggnh(il)
             voo(il,1,is) = voo(il,1,is)+sgpp2%clsg(1,is)*rs1(il)
          ENDDO
       ELSEIF (ic.EQ.2) THEN
          !$omp parallel do private(IL,GG,GGL)
          DO il=nsplpa,nsplpe
             ggl=ggnh(il)
             gg=ggl*parm%tpiba2
             voo(il,1,is) = voo(il,1,is)+&
                  sgpp2%clsg(2,is)*(3._real_8-rcrc*gg)*rs1(il)
          ENDDO
       ELSEIF (ic.EQ.3) THEN
          !$omp parallel do private(IL,GG,GGL)
          DO il=nsplpa,nsplpe
             ggl=ggnh(il)
             gg=ggl*parm%tpiba2
             voo(il,1,is) = voo(il,1,is)+&
                  sgpp2%clsg(3,is)*(15._real_8-10._real_8*rcrc*gg+rcrc*rcrc*gg*gg)*rs1(il)
          ENDDO
       ELSEIF (ic.EQ.4) THEN
          !$omp parallel do private(IL,GG,GGL)
          DO il=nsplpa,nsplpe
             ggl=ggnh(il)
             gg=ggl*parm%tpiba2
             voo(il,1,is) = voo(il,1,is)+&
                  sgpp2%clsg(4,is)*(105._real_8-105._real_8*rcrc*gg+21._real_8*rcrc*rcrc*gg*gg-&
                  rcrc*rcrc*rcrc*gg*gg*gg)*rs1(il)
          ENDDO
       ELSE
          CALL stopgm('FORMSG','NCLSG TOO BIG',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    IF (patom1%pconf) THEN
       r2=patom2%prc(is)**2
       rcpf=SQRT(pi**3)*patom2%prc(is)**3*patom2%palpha(is)
       !$omp parallel do private(IL,GG)
       DO il=nsplpa,nsplpe
          gg=ggnh(il)*parm%tpiba2
          IF (gg.GT.1.e-8_real_8) THEN
             voo(il,1,is)=voo(il,1,is)-rcpf*SQRT(gg)*EXP(-0.25_real_8*gg*r2)
          ENDIF
       ENDDO
    ENDIF
    CALL mp_sum(voo(:,:,is),nsplpo,parai%allgrp)
    CALL curv1(nsplpo,ggnh,voo(1,1,is),0.0_real_8,0.0_real_8,3,&
         voo(1,2,is),rs1,0.0_real_8,ierr)
    CALL tihalt('    FORMSG',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE formsg
  ! ==================================================================
  SUBROUTINE formfv(is,fx,vscr)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE FORM FACTORS OF PSEUDOPOTENTIAL (VPS)          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: fx(maxsys%mmaxx,*), vscr(*)

    REAL(real_8), PARAMETER                  :: eps = 1.e-5_real_8 

    INTEGER                                  :: ierr, il, intx, ir, isub, iz, &
                                                mepe
    REAL(real_8)                             :: figlx, fint, gg, r2, rcpf, &
                                                test, xg

    CALL tiset('    FORMFV',isub)
    test=0._real_8
    iz=1
    fx(1,1)=0._real_8
    vscr(1)=0._real_8
    DO ir=2,ncpr1%meshva(is)
       test=0.5_real_8*rucore(ir,1,is)+ions0%zv(is)*cp_erf(r(ir,is)/raggio(is))
       IF (ABS(test).LT.eps) THEN
          vscr(ir)=0._real_8
       ELSE
          vscr(ir)=test
          iz=ir
       ENDIF
       ! VSCR =   SCREENED PSEUDOPOTENTIAL
       fx(ir,1)=vscr(ir)*r(ir,is)
    ENDDO
    intx=MIN(40*(iz/40+1)+1,ncpr1%meshva(is))
    IF ( pslo_com%tlog(is) ) THEN
       CALL radlg(intx,fx(1,1),rab(1,is),fint)
    ELSE
       CALL radin(intx,ncpr1%cmesh(is),fx(1,1),fint)
    ENDIF
    CALL zeroing(voo(:,1,is))!,nsplpo)
    mepe=1
    !$omp parallel do private(IL,XG,IR,MEPE,FIGLX) shared(INTX)
    DO il=nsplpa,nsplpe
       !$      MEPE=OMP_GET_THREAD_NUM()+1
       xg=SQRT(ggnh(il))*parm%tpiba
       IF (xg.GT.1.e-6_real_8) THEN
          DO ir=1,intx
             fx(ir,mepe)=vscr(ir)*SIN(r(ir,is)*xg)
          ENDDO
          IF ( pslo_com%tlog(is) ) THEN
             CALL radlg(intx,fx(1,mepe),rab(1,is),figlx)
          ELSE
             CALL radin(intx,ncpr1%cmesh(is),fx(1,mepe),figlx)
          ENDIF
          voo(il,1,is)=fpi*figlx/xg
       ELSE
          voo(il,1,is)=fpi*fint
       ENDIF
    ENDDO
    IF (patom1%pconf) THEN
       r2=patom2%prc(is)**2
       rcpf=SQRT(pi**3)*patom2%prc(is)**3*patom2%palpha(is)
       !$omp parallel do private(IL,GG)
       DO il=nsplpa,nsplpe
          gg=ggnh(il)*parm%tpiba2
          IF (gg.GT.1.e-8_real_8) THEN
             voo(il,1,is)=voo(il,1,is)-rcpf*SQRT(gg)*EXP(-0.25_real_8*gg*r2)
          ENDIF
       ENDDO
    ENDIF
    CALL mp_sum(voo(:,:,is),nsplpo,parai%allgrp)
    CALL curv1(nsplpo,ggnh,voo(1,1,is),0.0_real_8,0.0_real_8,3,&
         voo(1,2,is),fx,0.0_real_8,ierr)
    CALL tihalt('    FORMFV',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE formfv
  ! ==================================================================
  REAL(real_8) FUNCTION sinc(x)
    REAL(real_8) :: x
    REAL(real_8), PARAMETER :: eps=EPSILON(1.0_real_8)
    sinc=1.0_real_8
    IF (ABS(x).GT.eps) sinc=SIN(x)/x
  END FUNCTION sinc
  ! ==================================================================
  SUBROUTINE formup(is,dfint,fint)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE FORM FACTORS OF PSEUDOPOTENTIAL (VPS)          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: dfint(maxsys%mmaxx,*), &
                                                fint(maxsys%mmaxx)

    INTEGER                                  :: ierr, il, ir, isub, mepe
    REAL(real_8)                             :: flx, gg, r2, rcpf, xg, xrg

! ==--------------------------------------------------------------==

    CALL tiset('    FORMUP',isub)
    CALL zeroing(fint)!,maxsys%mmaxx)
    DO ir=1,atom_common%meshvp(is)
       fint(ir)=rv(ir,is)*rv(ir,is)*vr(ir,is,1)+&
            ions0%zv(is)*cp_erf(rv(ir,is)/raggio(is))*rv(ir,is)
    ENDDO
    CALL zeroing(voo(:,1,is))!,nsplpo)
    mepe=1
    !$omp parallel do private(IL,XG,MEPE,IR,XRG,FLX)
    DO il=nsplpa,nsplpe
       !$      MEPE=OMP_GET_THREAD_NUM()+1
       xg=SQRT(ggnh(il))*parm%tpiba
       IF (xg.GT.1.e-6_real_8) THEN
          DO ir=1,atom_common%meshvp(is)
             xrg = rv(ir,is)*xg
             dfint (ir,mepe) = fint(ir) * sinc(xrg) * rw(ir,is)
          ENDDO
       ELSE
          DO ir=1,atom_common%meshvp(is)
             dfint (ir,mepe) = fint(ir) * rw(ir,is)
          ENDDO
       ENDIF
       CALL simpsn (atom_common%meshvp(is), dfint(1,mepe), flx)
       voo(il,1,is) = fpi* flx
    ENDDO
    IF (patom1%pconf) THEN
       r2=patom2%prc(is)**2
       rcpf=SQRT(pi**3)*patom2%prc(is)**3*patom2%palpha(is)
       !$omp parallel do private(IL,GG)
       DO il=nsplpa,nsplpe
          gg=ggnh(il)*parm%tpiba2
          IF (gg.GT.1.e-8_real_8) THEN
             voo(il,1,is)=voo(il,1,is)-rcpf*SQRT(gg)*EXP(-0.25_real_8*gg*r2)
          ENDIF
       ENDDO
    ENDIF
    CALL mp_sum(voo(:,:,is),nsplpo,parai%allgrp)
    CALL curv1(nsplpo,ggnh,voo(1,1,is),0.0_real_8,0.0_real_8,3,&
         voo(1,2,is),fint,0.0_real_8,ierr)
    CALL tihalt('    FORMUP',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE formup
  ! ==================================================================

END MODULE formf_utils
