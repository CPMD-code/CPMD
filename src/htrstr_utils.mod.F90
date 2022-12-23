MODULE htrstr_utils
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: hg,&
                                             inyh,&
                                             nzh,&
                                             rhops
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE str2,                            ONLY: drhovg,&
                                             gagk
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             dehc,&
                                             deht,&
                                             delta
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ncpw,&
                                             parm,&
                                             cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: htrstr

CONTAINS

  ! ==================================================================
  SUBROUTINE htrstr(eht,v,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == HARTREE ENERGY CONTRIBUTION TO THE STRESS TENSOR             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eht
    COMPLEX(real_8)                          :: v(*), eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: cfpibg, chgm1, ei123, rhet, &
                                                rhog, rhogs, rhopr
    INTEGER                                  :: ia, ig, ig1, is, isa, isub, kk
    REAL(real_8)                             :: ei, er, er2, fpibg, r2

! Variables

#ifdef __VECTOR 
    COMPLEX(real_8) :: t1,t2,t3,t4,t5,t6,fac1,fac2
#endif 
    ! ==--------------------------------------------------------------==
    CALL tiset('    HTRSTR',isub)
    IF (geq0) THEN
       ig1=2
    ELSE
       ig1=1
    ENDIF
#if defined(__VECTOR) && (! (__SR8000))
    t1=(0.0_real_8,0.0_real_8)
    t2=(0.0_real_8,0.0_real_8)
    t3=(0.0_real_8,0.0_real_8)
    t4=(0.0_real_8,0.0_real_8)
    t5=(0.0_real_8,0.0_real_8)
    t6=(0.0_real_8,0.0_real_8)
    IF (cntl%bigmem) THEN
       IF (pslo_com%tivan) THEN
          !$omp parallel do private(IG,ISA,IA,IS,R2,RHOPR, &
          !$omp  RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC1,FAC2) &
          !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
          DO ig=ig1,ncpw%nhg
             rhopr=(0.0_real_8,0.0_real_8)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                r2=raggio(is)*raggio(is)
                rhopr=rhopr+rhops(is,ig)*r2*0.5_real_8*eigrb(ig,isa)
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac1=cfpibg*rhogs
             fac2=fac1*(rhog*chgm1+rhopr)
             t1=t1+(fac2*gagk(ig,1)+fac1*drhovg(ig,1))
             t2=t2+(fac2*gagk(ig,2)+fac1*drhovg(ig,2))
             t3=t3+(fac2*gagk(ig,3)+fac1*drhovg(ig,3))
             t4=t4+(fac2*gagk(ig,4)+fac1*drhovg(ig,4))
             t5=t5+(fac2*gagk(ig,5)+fac1*drhovg(ig,5))
             t6=t6+(fac2*gagk(ig,6)+fac1*drhovg(ig,6))
             IF (cntl%tlsd) THEN
                t1=t1+(fac2*gagk(ig,1)+fac1*drhovg(ig,7))
                t2=t2+(fac2*gagk(ig,2)+fac1*drhovg(ig,8))
                t3=t3+(fac2*gagk(ig,3)+fac1*drhovg(ig,9))
                t4=t4+(fac2*gagk(ig,4)+fac1*drhovg(ig,10))
                t5=t5+(fac2*gagk(ig,5)+fac1*drhovg(ig,11))
                t6=t6+(fac2*gagk(ig,6)+fac1*drhovg(ig,12))
             ENDIF
          ENDDO
          dehc(1)=dehc(1)+t1
          dehc(2)=dehc(2)+t2
          dehc(3)=dehc(3)+t3
          dehc(4)=dehc(4)+t4
          dehc(5)=dehc(5)+t5
          dehc(6)=dehc(6)+t6
       ELSE
          !$omp parallel do private(IG,ISA,IA,IS,R2,RHOPR, &
          !$omp  RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC2) &
          !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
          DO ig=ig1,ncpw%nhg
             rhopr=(0.0_real_8,0.0_real_8)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                r2=raggio(is)*raggio(is)
                rhopr=rhopr+rhops(is,ig)*r2*0.5_real_8*eigrb(ig,isa)
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac2=cfpibg*rhogs*(rhog*chgm1+rhopr)
             t1=t1+fac2*gagk(ig,1)
             t2=t2+fac2*gagk(ig,2)
             t3=t3+fac2*gagk(ig,3)
             t4=t4+fac2*gagk(ig,4)
             t5=t5+fac2*gagk(ig,5)
             t6=t6+fac2*gagk(ig,6)
          ENDDO
          dehc(1)=dehc(1)+t1
          dehc(2)=dehc(2)+t2
          dehc(3)=dehc(3)+t3
          dehc(4)=dehc(4)+t4
          dehc(5)=dehc(5)+t5
          dehc(6)=dehc(6)+t6
       ENDIF
    ELSE
       IF (pslo_com%tivan) THEN
          !$omp parallel do private(IG,ISA,IA,IS,R2,EI123,ER,EI,ER2,RHOPR, &
          !$omp  RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC1,FAC2) &
          !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
          DO ig=ig1,ncpw%nhg
             rhopr=(0.0_real_8,0.0_real_8)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                r2=raggio(is)*raggio(is)
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                er2=0.5_real_8*r2*rhops(is,ig)
                rhopr=rhopr+CMPLX(er*er2,ei*er2,kind=real_8)
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac1=cfpibg*rhogs
             fac2=fac1*(rhog*chgm1+rhopr)
             t1=t1+(fac2*gagk(ig,1)+fac1*drhovg(ig,1))
             t2=t2+(fac2*gagk(ig,2)+fac1*drhovg(ig,2))
             t3=t3+(fac2*gagk(ig,3)+fac1*drhovg(ig,3))
             t4=t4+(fac2*gagk(ig,4)+fac1*drhovg(ig,4))
             t5=t5+(fac2*gagk(ig,5)+fac1*drhovg(ig,5))
             t6=t6+(fac2*gagk(ig,6)+fac1*drhovg(ig,6))
             IF (cntl%tlsd) THEN
                t1=t1+(fac2*gagk(ig,1)+fac1*drhovg(ig,7))
                t2=t2+(fac2*gagk(ig,2)+fac1*drhovg(ig,8))
                t3=t3+(fac2*gagk(ig,3)+fac1*drhovg(ig,9))
                t4=t4+(fac2*gagk(ig,4)+fac1*drhovg(ig,10))
                t5=t5+(fac2*gagk(ig,5)+fac1*drhovg(ig,11))
                t6=t6+(fac2*gagk(ig,6)+fac1*drhovg(ig,12))
             ENDIF
          ENDDO
          dehc(1)=dehc(1)+t1
          dehc(2)=dehc(2)+t2
          dehc(3)=dehc(3)+t3
          dehc(4)=dehc(4)+t4
          dehc(5)=dehc(5)+t5
          dehc(6)=dehc(6)+t6
       ELSE
          !$omp parallel do private(IG,ISA,IA,IS,R2,EI123,ER,EI,ER2,RHOPR, &
          !$omp  RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC2) &
          !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
          DO ig=ig1,ncpw%nhg
             rhopr=(0.0_real_8,0.0_real_8)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                r2=raggio(is)*raggio(is)
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                er2=0.5_real_8*r2*rhops(is,ig)
                rhopr=rhopr+CMPLX(er*er2,ei*er2,kind=real_8)
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac2=cfpibg*rhogs*(rhog*chgm1+rhopr)
             t1=t1+fac2*gagk(ig,1)
             t2=t2+fac2*gagk(ig,2)
             t3=t3+fac2*gagk(ig,3)
             t4=t4+fac2*gagk(ig,4)
             t5=t5+fac2*gagk(ig,5)
             t6=t6+fac2*gagk(ig,6)
          ENDDO
          dehc(1)=dehc(1)+t1
          dehc(2)=dehc(2)+t2
          dehc(3)=dehc(3)+t3
          dehc(4)=dehc(4)+t4
          dehc(5)=dehc(5)+t5
          dehc(6)=dehc(6)+t6
       ENDIF
    ENDIF
#elif defined(__SR8000) 
    IF (cntl%bigmem) THEN
       IF (pslo_com%tivan) THEN
          !poption parallel
          !poption tlocal(IG,ISA,RHOPR,IA,IS,R2)
          !poption tlocal(RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC1,FAC2)
          !poption psum(DEHC)
          DO ig=ig1,nhg
             rhopr=(0.0_real_8,0.0_real_8)
             isa=0
             DO is=1,ions1%nsp
                r2=raggio(is)*raggio(is)
                DO ia=1,ions0%na(is)
                   isa=isa+1
                   rhopr=rhopr+rhops(is,ig)*r2*0.5_real_8*eigrb(ig,isa)
                ENDDO
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac1=cfpibg*rhogs
             fac2=fac1*(rhog*chgm1+rhopr)
             dehc(1)=dehc(1)+(fac2*gagk(ig,1)+fac1*drhovg(ig,1))
             dehc(2)=dehc(2)+(fac2*gagk(ig,2)+fac1*drhovg(ig,2))
             dehc(3)=dehc(3)+(fac2*gagk(ig,3)+fac1*drhovg(ig,3))
             dehc(4)=dehc(4)+(fac2*gagk(ig,4)+fac1*drhovg(ig,4))
             dehc(5)=dehc(5)+(fac2*gagk(ig,5)+fac1*drhovg(ig,5))
             dehc(6)=dehc(6)+(fac2*gagk(ig,6)+fac1*drhovg(ig,6))
             IF (cntl%tlsd) THEN
                dehc(1)=dehc(1)+(fac2*gagk(ig,1)+fac1*drhovg(ig,7))
                dehc(2)=dehc(2)+(fac2*gagk(ig,2)+fac1*drhovg(ig,8))
                dehc(3)=dehc(3)+(fac2*gagk(ig,3)+fac1*drhovg(ig,9))
                dehc(4)=dehc(4)+(fac2*gagk(ig,4)+fac1*drhovg(ig,10))
                dehc(5)=dehc(5)+(fac2*gagk(ig,5)+fac1*drhovg(ig,11))
                dehc(6)=dehc(6)+(fac2*gagk(ig,6)+fac1*drhovg(ig,12))
             ENDIF
          ENDDO
       ELSE
          !poption parallel
          !poption tlocal(IG,ISA,RHOPR,IA,IS,R2)
          !poption tlocal(RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC2)
          !poption psum(DEHC)
          DO ig=ig1,nhg
             rhopr=(0.0_real_8,0.0_real_8)
             isa=0
             DO is=1,ions1%nsp
                r2=raggio(is)*raggio(is)
                DO ia=1,ions0%na(is)
                   isa=isa+1
                   rhopr=rhopr+rhops(is,ig)*r2*0.5_real_8*eigrb(ig,isa)
                ENDDO
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac2=cfpibg*rhogs*(rhog*chgm1+rhopr)
             dehc(1)=dehc(1)+fac2*gagk(ig,1)
             dehc(2)=dehc(2)+fac2*gagk(ig,2)
             dehc(3)=dehc(3)+fac2*gagk(ig,3)
             dehc(4)=dehc(4)+fac2*gagk(ig,4)
             dehc(5)=dehc(5)+fac2*gagk(ig,5)
             dehc(6)=dehc(6)+fac2*gagk(ig,6)
          ENDDO
       ENDIF
    ELSE
       IF (pslo_com%tivan) THEN
          !poption parallel
          !poption tlocal(IG,ISA,RHOPR,IA,IS,R2,EI123,ER,EI,ER2)
          !poption tlocal(RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC1,FAC2)
          !poption psum(DEHC)
          DO ig=ig1,nhg
             rhopr=(0.0_real_8,0.0_real_8)
             isa=0
             DO is=1,ions1%nsp
                r2=raggio(is)*raggio(is)
                DO ia=1,ions0%na(is)
                   isa=isa+1
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                   er=REAL(ei123)
                   ei=AIMAG(ei123)
                   er2=0.5_real_8*r2*rhops(is,ig)
                   rhopr=rhopr+CMPLX(er*er2,ei*er2,kind=real_8)
                ENDDO
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac1=cfpibg*rhogs
             fac2=fac1*(rhog*chgm1+rhopr)
             dehc(1)=dehc(1)+(fac2*gagk(ig,1)+fac1*drhovg(ig,1))
             dehc(2)=dehc(2)+(fac2*gagk(ig,2)+fac1*drhovg(ig,2))
             dehc(3)=dehc(3)+(fac2*gagk(ig,3)+fac1*drhovg(ig,3))
             dehc(4)=dehc(4)+(fac2*gagk(ig,4)+fac1*drhovg(ig,4))
             dehc(5)=dehc(5)+(fac2*gagk(ig,5)+fac1*drhovg(ig,5))
             dehc(6)=dehc(6)+(fac2*gagk(ig,6)+fac1*drhovg(ig,6))
             IF (cntl%tlsd) THEN
                dehc(1)=dehc(1)+(fac2*gagk(ig,1)+fac1*drhovg(ig,7))
                dehc(2)=dehc(2)+(fac2*gagk(ig,2)+fac1*drhovg(ig,8))
                dehc(3)=dehc(3)+(fac2*gagk(ig,3)+fac1*drhovg(ig,9))
                dehc(4)=dehc(4)+(fac2*gagk(ig,4)+fac1*drhovg(ig,10))
                dehc(5)=dehc(5)+(fac2*gagk(ig,5)+fac1*drhovg(ig,11))
                dehc(6)=dehc(6)+(fac2*gagk(ig,6)+fac1*drhovg(ig,12))
             ENDIF
          ENDDO
       ELSE
          !poption parallel
          !poption tlocal(IG,ISA,RHOPR,IA,IS,R2,EI123,ER,EI,ER2)
          !poption tlocal(RHET,RHOG,RHOGS,FPIBG,CFPIBG,CHGM1,FAC2)
          !poption psum(DEHC)
          DO ig=ig1,nhg
             rhopr=(0.0_real_8,0.0_real_8)
             isa=0
             DO is=1,ions1%nsp
                r2=raggio(is)*raggio(is)
                DO ia=1,ions0%na(is)
                   isa=isa+1
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                   er=REAL(ei123)
                   ei=AIMAG(ei123)
                   er2=0.5_real_8*r2*rhops(is,ig)
                   rhopr=rhopr+CMPLX(er*er2,ei*er2,kind=real_8)
                ENDDO
             ENDDO
             rhet=v(nzh(ig))
             rhog=rhet+eirop(ig)
             rhogs=CONJG(rhog)
             fpibg=fpi/parm%tpiba2/hg(ig)
             cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
             chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
             fac2=cfpibg*rhogs*(rhog*chgm1+rhopr)
             dehc(1)=dehc(1)+fac2*gagk(ig,1)
             dehc(2)=dehc(2)+fac2*gagk(ig,2)
             dehc(3)=dehc(3)+fac2*gagk(ig,3)
             dehc(4)=dehc(4)+fac2*gagk(ig,4)
             dehc(5)=dehc(5)+fac2*gagk(ig,5)
             dehc(6)=dehc(6)+fac2*gagk(ig,6)
          ENDDO
       ENDIF
    ENDIF
#else 
    IF (cntl%bigmem) THEN
       !$omp parallel do private(ig,isa,is,ia,kk,rhopr,r2, &
       !$omp rhet,rhog,rhogs,fpibg,cfpibg,chgm1) &
       !$omp reduction(+:dehc)
       DO ig=ig1,ncpw%nhg
          rhopr = (0.0_real_8,0.0_real_8)
          isa=0
          DO is=1,ions1%nsp
             r2 = raggio(is)*raggio(is)
             DO ia=1,ions0%na(is)
                isa=isa+1
                rhopr = rhopr + rhops(is,ig)*r2*0.5_real_8*eigrb(ig,isa)
             ENDDO
          ENDDO
          rhet=v(nzh(ig))
          rhog=rhet+eirop(ig)
          rhogs=CONJG(rhog)
          fpibg=fpi/parm%tpiba2/hg(ig)
          cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
          chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
          DO kk=1,6
             dehc(kk) = dehc(kk) +&
                  cfpibg*rhogs*(rhog*chgm1+rhopr)*gagk(ig,kk)
          ENDDO
          IF (pslo_com%tivan) THEN
             DO kk=1,6
                dehc(kk) = dehc(kk) + cfpibg*rhogs*drhovg(ig,kk)
                IF (cntl%tlsd) dehc(kk) = dehc(kk) + cfpibg*rhogs*drhovg(ig,6+kk)
             ENDDO
          ENDIF
       ENDDO
    ELSE
       !$omp parallel do private(ig,isa,is,ia,kk,rhopr,r2,ei123,er,ei,er2, &
       !$omp rhet,rhog,rhogs,fpibg,cfpibg,chgm1) &
       !$omp reduction(+:dehc)
       DO ig=ig1,ncpw%nhg
          rhopr = (0.0_real_8,0.0_real_8)
          isa=0
          DO is=1,ions1%nsp
             r2 = raggio(is)*raggio(is)
             DO ia=1,ions0%na(is)
                isa=isa+1
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                er2=0.5_real_8*r2*rhops(is,ig)
                rhopr = rhopr + CMPLX(er*er2,ei*er2,kind=real_8)
             ENDDO
          ENDDO
          rhet=v(nzh(ig))
          rhog=rhet+eirop(ig)
          rhogs=CONJG(rhog)
          fpibg=fpi/parm%tpiba2/hg(ig)
          cfpibg=CMPLX(fpibg,0._real_8,kind=real_8)
          chgm1=CMPLX(1._real_8/hg(ig)/parm%tpiba2,0.0_real_8,kind=real_8)
          DO kk=1,6
             dehc(kk) = dehc(kk) +&
                  cfpibg*rhogs*(rhog*chgm1+rhopr)*gagk(ig,kk)
          ENDDO
          IF (pslo_com%tivan) THEN
             DO kk=1,6
                dehc(kk) = dehc(kk) + cfpibg*rhogs*drhovg(ig,kk)
                IF (cntl%tlsd) dehc(kk) = dehc(kk) + cfpibg*rhogs*drhovg(ig,6+kk)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
#endif 
    DO kk=1,6
       deht(kk) = -eht*delta(alpha(kk),beta(kk)) +&
            2.0_real_8*parm%omega*REAL(dehc(kk))
    ENDDO
    CALL tihalt('    HTRSTR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE htrstr
  ! ==================================================================

END MODULE htrstr_utils
