MODULE vofrhoh_utils
  USE cppt,                            ONLY: gk,&
                                             hipz,&
                                             indz,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             vps
  USE eextern_utils,                   ONLY: eextern
  USE efld,                            ONLY: textfld
  USE eicalc_utils,                    ONLY: eicalc
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE hip_utils,                       ONLY: give_qphi,&
                                             hip
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE td_input,                        ONLY: td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhoh
  PUBLIC :: give_scr_vofrhoh

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhoh(tau0,fion,rhoe,v,vtemp,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==--------------------------------------------------------------==
    ! ==  RHOE: in  electronic density in real space                  ==
    ! ==  V:    out electronic density in G space                     ==
    ! ==--------------------------------------------------------------==
    ! EHR[
    ! EHR]
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    REAL(real_8), INTENT(in)                 :: rhoe(:)
    COMPLEX(real_8)                          :: v(:), vtemp(:)
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhoh'

    COMPLEX(real_8)                          :: ei123, eps, gx, gy, gz, rhet, &
                                                rhets, rhog, rhot, rp, txx, &
                                                tyy, tzz, vp, vzcp
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), qphi(:)
    INTEGER                                  :: ia, ierr, ig, ig1, il_qphi, &
                                                ir, is, isa, isub, nnrs
    REAL(real_8)                             :: ehip1, ehip2, eze, omtp

    CALL tiset(procedureN,isub)
    CALL setfftn(0)
    ! SCR partition
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (isos1%tclust) THEN
       CALL give_qphi(il_qphi)
       ALLOCATE(qphi(il_qphi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Arrays inside SCR.
    nnrs  = spar%nr1s*spar%nr2s*spar%nr3s
    CALL eicalc(eivps,eirop)
    ! ..External Field or forces of external potential
    ener_com%eext=0._real_8
    ! EHR[
    IF ((textfld.OR.cntl%texadd).AND.(.NOT.td_prop%td_extpot)) THEN
       ! EHR]
       CALL eextern(rhoe,v,eirop,tau0,fion,ener_com%eext,tfor)
       CALL mp_sum(ener_com%eext,parai%allgrp)
    ENDIF
    ! TRANSFORM THE DENSITY TO G SPACE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(rhoe(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==                       ISOLATED SYSTEMS                       ==
    ! ==          CALCULATE THE HARTREE POTENTIAL AND ENERGY          ==
    ! ==--------------------------------------------------------------==
    ! ZERO CORRECTION AND PSEUDO-POTENTIAL
    IF (geq0) THEN
       vp=eivps(1)
       rp=eirop(1)
       vploc=REAL(vp)
       rhog=v(nzh(1))
       rhot=rhog+rp
       vtemp(1)=hipz(1)*rhot+vp
       eze=0.5_real_8*hipz(1)*REAL(rhot)*REAL(rhot)
       eps=0.5_real_8*REAL(vp*rhog)
       ig1=2
    ELSE
       eze=0.0_real_8
       eps=(0.0_real_8,0.0_real_8)
       ig1=1
    ENDIF
    !$omp parallel do private(IG,VP,RP,RHOG,RHOT,VZCP) &
    !$omp  reduction(+:EZE,EPS)
    DO ig=ig1,ncpw%nhg
       vp=eivps(ig)
       rp=eirop(ig)
       rhog=v(nzh(ig))
       rhot=rhog+rp
       vzcp=hipz(ig)*rhot
       vtemp(ig)=vp+vzcp
       eze=eze+(REAL(rhot)*REAL(vzcp)+AIMAG(rhot)*AIMAG(vzcp))
       eps=eps+(REAL(vp)*REAL(rhog)+AIMAG(vp)*AIMAG(rhog))
    ENDDO
    ener_com%epseu=2._real_8*REAL(eps)*parm%omega
    ! ==------------------------------------------------------------==
    ! FORCES ON ATOMS DUE TO ZERO CORRECTION AND PSEUDO-POTENTIAL
    IF (tfor) THEN
       omtp=2._real_8*parm%omega*parm%tpiba
       IF (geq0) THEN
          ig1=2
       ELSE
          ig1=1
       ENDIF
#if defined (__VECTOR)
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             IF (cntl%bigmem) THEN
                DO ig=ig1,ncpw%nhg
                   ei123=eigrb(ig,isa)
                   rp=eirop(ig)
                   rhet=v(nzh(ig))
                   rhog=rhet+rp
                   rhets=CONJG(rhet)
                   gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                   gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                   gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                   vzcp=CONJG(hipz(ig)*rhog)
                   txx=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gx
                   tyy=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gy
                   tzz=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gz
                   fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                   fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                   fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
                ENDDO
             ELSE
                DO ig=ig1,ncpw%nhg
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                   rp=eirop(ig)
                   rhet=v(nzh(ig))
                   rhog=rhet+rp
                   rhets=CONJG(rhet)
                   gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                   gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                   gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                   vzcp=CONJG(hipz(ig)*rhog)
                   txx=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gx
                   tyy=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gy
                   tzz=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gz
                   fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                   fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                   fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
                ENDDO
             ENDIF
          ENDDO
       ENDDO
#else
       DO ig=ig1,ncpw%nhg
          rp=eirop(ig)
          rhet=v(nzh(ig))
          rhog=rhet+rp
          rhets=CONJG(rhet)
          gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
          gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
          gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
          vzcp=CONJG(hipz(ig)*rhog)
          isa=0
          DO is=1,ions1%nsp
             txx=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gx
             tyy=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gy
             tzz=(rhops(is,ig)*vzcp+vps(is,ig)*rhets)*gz
             DO ia=1,ions0%na(is)
                isa=isa+1
                IF (cntl%bigmem) THEN
                   ei123=eigrb(ig,isa)
                ELSE
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                ENDIF
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ENDDO
       ENDDO
#endif
    ENDIF
    ! ==------------------------------------------------------------==
    ! HARTREE POTENTIAL
    ig1=1
    IF (geq0) ig1=2
    !CDIR NODEP
    !$omp parallel do private(IG)
    DO ig=ig1,ncpw%nhg
       v(nzh(ig))=v(nzh(ig))+eirop(ig)
       v(indz(ig))=v(indz(ig))+CONJG(eirop(ig))
    ENDDO
    IF (geq0) v(nzh(1))=v(nzh(1))+eirop(1)
    CALL invfftn(v,.FALSE.,parai%allgrp)
    CALL hip(v,qphi)
    ehip1=0.0_real_8
#if defined(__VECTOR)
    !$omp parallel do private(IR) reduction(+:EHIP1)
#else
    !$omp parallel do private(IR) reduction(+:EHIP1) schedule(static)
#endif
    DO ir=1,fpar%nnr1
       ehip1=ehip1+REAL(v(ir))*rhoe(ir)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ehip2=0.0_real_8
    !$omp parallel do private(IG,VZCP) reduction(+:EHIP2)
    DO ig=1,ncpw%nhg
       vzcp=v(nzh(ig))
       ehip2=ehip2+(REAL(vzcp)*REAL(eirop(ig))+&
            AIMAG(vzcp)*AIMAG(eirop(ig)))
       vtemp(ig)=vtemp(ig)+vzcp
    ENDDO
    IF (geq0) ehip2=ehip2-0.5_real_8*REAL(v(nzh(1)))*REAL(eirop(1))
    ener_com%ehep = 0.5_real_8*ehip1*parm%omega/REAL(nnrs,kind=real_8)&
         + ehip2*parm%omega&
         + eze*parm%omega
    ener_com%ehee = 0._real_8
    ener_com%ehii = 0._real_8
    ener_com%eht = 0.5_real_8*ehip1*parm%omega/REAL(nnrs,kind=real_8)&
         + ehip2*parm%omega&
         + eze*parm%omega&
         + ener_com%esr&
         - ener_com%eself/REAL(parai%nproc,kind=real_8)
    ! FORCES ON ATOMS DUE TO HARTREE TERM
    IF (tfor) THEN
       omtp=2._real_8*parm%omega*parm%tpiba
       ig1=1
       IF (geq0) ig1=2
#if defined (__VECTOR)
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             IF (cntl%bigmem) THEN
                DO ig=ig1,ncpw%nhg
                   ei123=eigrb(ig,isa)
                   gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                   gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                   gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                   vzcp=CONJG(v(nzh(ig)))
                   txx=rhops(is,ig)*vzcp*gx
                   tyy=rhops(is,ig)*vzcp*gy
                   tzz=rhops(is,ig)*vzcp*gz
                   fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                   fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                   fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
                ENDDO
             ELSE
                DO ig=ig1,ncpw%nhg
                   gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                   gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                   gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                   gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                   gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                   gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                   vzcp=CONJG(v(nzh(ig)))
                   txx=rhops(is,ig)*vzcp*gx
                   tyy=rhops(is,ig)*vzcp*gy
                   tzz=rhops(is,ig)*vzcp*gz
                   fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                   fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                   fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
                ENDDO
             ENDIF
          ENDDO
       ENDDO
#else
       DO ig=ig1,ncpw%nhg
          gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
          gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
          gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
          vzcp=CONJG(v(nzh(ig)))
          isa=0
          DO is=1,ions1%nsp
             txx=rhops(is,ig)*vzcp*gx
             tyy=rhops(is,ig)*vzcp*gy
             tzz=rhops(is,ig)*vzcp*gz
             DO ia=1,ions0%na(is)
                isa=isa+1
                IF (cntl%bigmem) THEN
                   ei123=eigrb(ig,isa)
                ELSE
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                ENDIF
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ENDDO
       ENDDO
#endif
    ENDIF
    ! ..calculate energy for diagonalization schemes
    IF (cntl%tdiag) THEN
       ! TRANSFORM THE DENSITY TO G SPACE
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir)=CMPLX(rhoe(ir),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(v,.FALSE.,parai%allgrp)
       IF (geq0) THEN
          rp=eirop(1)
          rhog=v(nzh(1))
          eze=-0.5_real_8*hipz(1)*REAL(rhog)*REAL(rhog)
          eze=eze+0.5_real_8*hipz(1)*REAL(rp)*REAL(rp)
          ig1=2
       ELSE
          eze=0.0_real_8
          eps=(0.0_real_8,0.0_real_8)
          ig1=1
       ENDIF
       !$omp parallel do private(IG,RHOG,VZCP,RP) reduction(+:EZE)
       DO ig=ig1,ncpw%nhg
          rhog=v(nzh(ig))
          vzcp=hipz(ig)*rhog
          eze=eze+(-REAL(rhog)*REAL(vzcp)-AIMAG(rhog)*AIMAG(vzcp))
          rp=eirop(ig)
          vzcp=hipz(ig)*rp
          eze=eze+(REAL(rp)*REAL(vzcp)+AIMAG(rp)*AIMAG(vzcp))
       ENDDO
       eze=eze*parm%omega
       ehip1 = 0.5_real_8*ehip1*parm%omega/REAL(nnrs,kind=real_8)
       ehip2 = ehip2*parm%omega
       ener_com%ehep = ehip2-ehip1 + eze
       ener_com%ehee = 0._real_8
       ener_com%ehii = 0._real_8
       ener_com%eht =  ehip2-ehip1&
            + eze&
            + ener_com%esr&
            - ener_com%eself/REAL(parai%nproc,kind=real_8)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (isos1%tclust) THEN
       DEALLOCATE(qphi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vofrhoh
  ! ==================================================================
  SUBROUTINE give_scr_vofrhoh(lvofrhoh,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhoh
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_qphi

    CALL give_qphi(il_qphi)
    ! EIVPS(2*NHG) EIROP(2*NHG)
    lvofrhoh=4*ncpw%nhg+il_qphi
    tag      ='2*NHG*2+2*(NR1H+1)*NR2H*NR3PL'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhoh
  ! ==================================================================

END MODULE vofrhoh_utils
