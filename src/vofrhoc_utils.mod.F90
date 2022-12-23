MODULE vofrhoc_utils
  USE cnst,                            ONLY: pi
  USE cofor_utils,                     ONLY: cofor
  USE cppt,                            ONLY: gk,&
                                             indz,&
                                             inyh,&
                                             nzh,&
                                             scg,&
                                             vps
  USE eextern_utils,                   ONLY: eextern
  USE efld,                            ONLY: extf,&
                                             textfld
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com,&
                                             ener_d
  USE error_handling,                  ONLY: stopgm
  USE extpotmod,                       ONLY: extpot
  USE fcas,                            ONLY: fion_a,&
                                             fion_ab,&
                                             init_flag
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE gcener_utils,                    ONLY: gcener
  USE geq0mod,                         ONLY: geq0
  USE graden_utils,                    ONLY: graden
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlcc,                            ONLY: corel,&
                                             corer,&
                                             vnlt
  USE parac,                           ONLY: parai,&
                                             paral
  USE potfor_utils,                    ONLY: potfor
  USE ppener_utils,                    ONLY: ppener
  USE rnlfor_utils,                    ONLY: rnlfor
  USE ropt,                            ONLY: iteropt
  USE rpiiint_utils,                   ONLY: rpiiint
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE spin,                            ONLY: lspin4
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhoc
  PUBLIC :: give_scr_vofrhoc
  !public :: hvl_ab
  !public :: interxc
  PUBLIC :: vcomp
  PUBLIC :: ctfor
  !public :: pcasfor
  PUBLIC :: potab

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhoc(tau0,rhoe,v,vtemp,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE HARTREE AND XC ENERGY AND POTENTIALS (+FORCES) FOR      ==
    ! ==  THE GROUND STATE AND OFF-DIAGONAL ELEMENT IN THE CAS22 METH.==
    ! ==--------------------------------------------------------------==
    ! ==  RHOE(.,1) Density of all 2 occ states                       ==
    ! ==  RHOE(.,2) RHO(ab)                                           ==
    ! ==  RHOE(.,3) "Ground state" density                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhoc'

    COMPLEX(real_8)                          :: ee, eh, ehab, ei, eps, epsab
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), vtmp(:)
    INTEGER                                  :: ierr, ig, il_vtmp, ir, isub, &
                                                nnrs
    REAL(real_8)                             :: egca, esx, fide(1,1,1), &
                                                rint(4), sgcc, sgcx, sxc, &
                                                sxcab, wk(1)
    REAL(real_8), ALLOCATABLE                :: focc(:), grad(:,:), rvtmp(:)

    CALL tiset('   VOFRHOC',isub)
    CALL setfftn(0)
    IF (tstress) CALL stopgm("VOFRHOC","STRESS NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    IF (cntl%ttau) THEN
       CALL stopgm('VOFRHOC','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (tfor) THEN
       IF (init_flag.EQ.0) THEN
          ALLOCATE(fion_a(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(fion_ab(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          init_flag=1
       ENDIF
       CALL zeroing(fion_a)!,3*maxsys%nax*maxsys%nsx)
       CALL rpiiint(esx,tau0,fion_a,iteropt%iesr,tfor)
       CALL zeroing(fion_ab)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Test integral of densities
    CALL zeroing(rint)!,4)
    DO ir=1,fpar%nnr1
       rint(1)=rint(1)+rhoe(ir,1)
       rint(2)=rint(2)+rhoe(ir,2)
       rint(3)=rint(3)+rhoe(ir,3)
       rint(4)=rint(4)+rhoe(ir,4)
    ENDDO
    CALL mp_sum(rint,4,parai%allgrp)
    DO ir=1,4
       rint(ir)=rint(ir)*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ENDDO
    IF (paral%parent) THEN
       IF ((ABS(rint(1)-(crge%nel-2._real_8)).GT.1.e-5_real_8) .AND.paral%io_parent)&
            WRITE(6,*) " WARNING : DENSITY RHO(1..N) ",rint(1)
       IF ((ABS(rint(2)).GT.1.e-5_real_8) .AND.paral%io_parent)&
            WRITE(6,*) " WARNING : DENSITY RHO(A..B) ",rint(2)
       IF ((ABS(rint(3)-crge%nel).GT.1.e-5_real_8) .AND.paral%io_parent)&
            WRITE(6,*) " WARNING : DENSITY RHO(1..A) ",rint(3)
       IF ((ABS(rint(4)-crge%nel).GT.1.e-5_real_8) .AND.paral%io_parent)&
            WRITE(6,*) " WARNING : DENSITY RHO(1..B) ",rint(4)
    ENDIF
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL eicalc(eivps,eirop)
    ! External Field
    ener_c%eext_a=0._real_8
    ener_c%eext_ab=0._real_8
    IF (textfld) THEN
       CALL eextern(rhoe(:,2),v,eirop,tau0,fion_ab,ener_c%eext_ab,tfor)
       CALL eextern(rhoe(:,3),v,eirop,tau0,fion_a,ener_c%eext_a,tfor)
       CALL eextern(rhoe(:,4),v,eirop,tau0,fide,ener_c%eext_2,.FALSE.)
    ENDIF
    ! ==----------------HARTREE AND PSEUOPOTENTIAL--------------------==
    ! TRANSFORM THE DENSITY TO G SPACE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
       v(ir,2)=CMPLX(rhoe(ir,2),0._real_8,kind=real_8)
       v(ir,3)=CMPLX(rhoe(ir,3),0._real_8,kind=real_8)
       v(ir,4)=CMPLX(rhoe(ir,4),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
    CALL  fwfftn(v(:,3),.FALSE.,parai%allgrp)
    CALL  fwfftn(v(:,4),.FALSE.,parai%allgrp)
    ! COMPUTE CONTRIBUTION TO FORCES ON IONS FROM GROUND STATE
    IF (tfor) THEN
       CALL potfor(fion_a,v(:,3),eirop)
       CALL pcasfor(fion_ab,v(:,2),eirop)
       wk=1._real_8
       ALLOCATE(focc(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(IR)
       DO ir=1,crge%n-1
          focc(ir)=2._real_8
       ENDDO
       focc(ir)=0._real_8
       CALL rnlfor(fion_a,focc,wk,crge%n,1)
       DEALLOCATE(focc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! COMPUTE ELECTROSTATIC AND PSEUDOPOTENTIAL ENERGIES.
    CALL ppener(eh,ei,ee,eps,v(:,3),vtemp(:,3),eivps,eirop)
    CALL hvl_ab(ehab,epsab,v(:,1),v(:,2),vtemp(:,1),vtemp(:,2),&
         eivps,eirop)
    ener_c%eht_a = REAL(eh)*parm%omega
    ener_c%epseu_a = 2._real_8*REAL(eps)*parm%omega
    ener_c%eht_ab = REAL(ehab)*parm%omega
    ener_c%epseu_ab = 2._real_8*REAL(epsab)*parm%omega
    CALL ppener(eh,ei,ee,eps,v(:,4),vtemp(:,4),eivps,eirop)
    ener_c%eht_2 = REAL(eh)*parm%omega
    ener_c%epseu_2 = 2._real_8*REAL(eps)*parm%omega
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ADD CORE CHARGE TO RHOE
    IF (corel%tinlc) THEN
       CALL stopgm("VOFRHOC","NLCC not implemented",& 
            __LINE__,__FILE__)
       !vw commented out next line as vtmp is not yet allocated at that point
       !vw CALL corec(rhoe(:,3),vtmp,v(:,3))
       CALL dcopy(2*ncpw%nhg,vtemp(1,3),1,vnlt(1),1)
    ENDIF
    ! VTEMP (Potential in G-Space) -FFT-> V(R)
    CALL zeroing(v(:,1))
    DO ig=1,ncpw%nhg
       v(indz(ig),1) = CONJG(vtemp(ig,1))
       v(nzh(ig),1)  = vtemp(ig,1)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),1) = vtemp(1,1)
    ENDIF
    CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
    ! 
    CALL zeroing(v(:,2))
    DO ig=1,ncpw%nhg
       v(indz(ig),2) = CONJG(vtemp(ig,2))
       v(nzh(ig),2)  = vtemp(ig,2)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),2) = vtemp(1,2)
    ENDIF
    CALL  invfftn(v(:,2),.FALSE.,parai%allgrp)
    ! 
    CALL zeroing(v(:,3))
    DO ig=1,ncpw%nhg
       v(indz(ig),3) = CONJG(vtemp(ig,3))
       v(nzh(ig),3)  = vtemp(ig,3)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),3) = vtemp(1,3)
    ENDIF
    CALL  invfftn(v(:,3),.FALSE.,parai%allgrp)
    ! ADD EXTERNAL POTENTIAL TO V
    IF (textfld) THEN
       ! potential from classical interface
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extf(ir),0._real_8,kind=real_8)
          v(ir,2)=v(ir,2)+CMPLX(extf(ir),0._real_8,kind=real_8)
          v(ir,3)=v(ir,3)+CMPLX(extf(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    IF (cntl%texpot) THEN
       ! static external potential
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extpot(ir),0._real_8,kind=real_8)
          v(ir,2)=v(ir,2)+CMPLX(extpot(ir),0._real_8,kind=real_8)
          v(ir,3)=v(ir,3)+CMPLX(extpot(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    ! XC ENERGY FOR INTERACTION TERM
    CALL interxc(sxcab,rhoe(:,1),rhoe(:,2),v(:,1),v(:,2))
    ! COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)
    CALL xcener(sxc,ener_com%vxc,rhoe(:,3:),rhoe(:,3:),v(:,3:))
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    IF (cntl%tgc) THEN
       ! CALCULATE THE GRADIENT OF THE DENSITY
       ALLOCATE(grad(fpar%nnr1,4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       il_vtmp = MAX(fpar%nnr1,ncpw%nhg*2)
       ALLOCATE(vtmp(il_vtmp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(rvtmp(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL fwfftn(v(:,3),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,3),vtemp(:,3),nzh)
       CALL graden(rhoe(:,3),v(:,3),grad,vtmp)
       ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
       CALL gcener(sgcx,sgcc,rhoe(:,3:),v(:,3:),vtemp(:,3),rvtmp,grad,&
            tstress)
       DEALLOCATE(grad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(vtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(rvtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! V CONTAINS THE TOTAL POTENTIAL IN R-SPACE, MOVE IT TO RHOE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,1)=REAL(v(ir,1))
       rhoe(ir,2)=REAL(v(ir,2))
       rhoe(ir,3)=REAL(v(ir,3))
    ENDDO
    ! CORRECT ENERGY FOR CORE CHARGE
    IF (corel%tinlc) THEN
       sxc=sxc-corer%excco
       sgcx=sgcx-corer%egcxco
       sgcc=sgcc-corer%egccco
    ENDIF
    ! SINGLE PIECES OF THE ENERGY:
    nnrs = spar%nr1s*spar%nr2s*spar%nr3s
    egca=(sgcc+sgcx)*parm%omega/REAL(nnrs,kind=real_8)
    ener_c%exc_a=sxc*parm%omega/REAL(nnrs,kind=real_8) + egca
    ener_c%exc_ab=sxcab*parm%omega/REAL(nnrs,kind=real_8)
    ! COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)
    CALL xcener(sxc,ener_com%vxc,rhoe(:,4:),rhoe(:,4:),v(:,4:))
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    IF (cntl%tgc) THEN
       ! CALCULATE THE GRADIENT OF THE DENSITY
       ALLOCATE(grad(fpar%nnr1,4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       il_vtmp = MAX(fpar%nnr1,ncpw%nhg*2)
       ALLOCATE(vtmp(il_vtmp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)! TODO /2?
       ALLOCATE(rvtmp(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL fwfftn(v(:,4),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,4),vtemp(:,4),nzh)
       CALL graden(rhoe(:,4),v(:,4),grad,vtmp)
       ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
       CALL gcener(sgcx,sgcc,rhoe(:,4:),v(:,4:),vtemp(:,4),rvtmp,grad,&
            tstress)
       DEALLOCATE(grad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(vtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(rvtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! CORRECT ENERGY FOR CORE CHARGE
    IF (corel%tinlc) THEN
       sxc=sxc-corer%excco
       sgcx=sgcx-corer%egcxco
       sgcc=sgcc-corer%egccco
    ENDIF
    ! SINGLE PIECES OF THE ENERGY:
    nnrs = spar%nr1s*spar%nr2s*spar%nr3s
    egca=(sgcc+sgcx)*parm%omega/REAL(nnrs,kind=real_8)
    ener_c%exc_2=sxc*parm%omega/REAL(nnrs,kind=real_8) + egca
    IF (corel%tinlc.AND.tfor) THEN
       ! TRANSF. TOTAL POTENTIAL IN G-SPACE
       ! PUT THE TOTAL POTENTIAL IN G-SPACE INTO VTEMP
       CALL  fwfftn(v(:,3),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,3),vtemp(:,3),nzh)
    ENDIF
    ! FORCE ON IONS DUE TO NLCC
    IF (corel%tinlc.AND.tfor) CALL cofor(fion_a,vtemp(:,3))
    ! 
    CALL tihalt('   VOFRHOC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vofrhoc
  ! ==================================================================
  SUBROUTINE give_scr_vofrhoc(lvofrhoc,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhoc
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lht, lxc

    lht      = 2*ncpw%nhg*2 + crge%n
    IF (cntl%tgc) THEN
       lxc    = 4*fpar%nnr1 + MAX(fpar%nnr1,2*ncpw%nhg)
    ELSE
       lxc    = 0
    ENDIF
    lvofrhoc = MAX(lht,lxc)
    tag      ='MAX(LHT,LXC)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhoc
  ! ==================================================================
  SUBROUTINE hvl_ab(eh,eps,cr0,crab,v0,vab,eivps,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == ELECTROSTATIC AND PSEUDOPOTENTIAL ENERGIES FOR INTERACTION   ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eh, eps, cr0(*), crab(*), &
                                                v0(ncpw%nhg), vab(ncpw%nhg), &
                                                eivps(ncpw%nhg), &
                                                eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: rabs, vcg, vp
    INTEGER                                  :: ig, ig1

    IF (geq0) THEN
       ig1=2
       v0(1)=CMPLX(0._real_8,0._real_8,kind=real_8)
       vab(1)=scg(1)*(eirop(1)+cr0(nzh(1)))+eivps(1)
    ELSE
       ig1=1
    ENDIF
    eh=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    eps=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    DO ig=ig1,ncpw%nhg
       rabs=CONJG(crab(nzh(ig)))
       v0(ig)=scg(ig)*crab(nzh(ig))
       vcg=scg(ig)*(eirop(ig)+cr0(nzh(ig)))
       vp=eivps(ig)
       vab(ig)=vcg+vp
       ! Hartree term
       eh=eh+2._real_8*vcg*rabs
       ! Pseudopotential energy
       eps=eps+rabs*vp
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hvl_ab
  ! ==================================================================
  SUBROUTINE interxc(sxcab,r0,rab,v0,vab)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sxcab, r0(*), rab(*)
    COMPLEX(real_8)                          :: v0(*), vab(*)

    REAL(real_8), PARAMETER :: small = 1.e-12_real_8, &
      third = 0.33333333333333_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: cf, r3

! ==--------------------------------------------------------------==

    cf=4._real_8*(3._real_8/pi)**third
    sxcab=0._real_8
    DO ir=1,fpar%nnr1
       IF (ABS(rab(ir)).GT.small.AND.r0(ir).GT.small) THEN
          r3=r0(ir)**third
          sxcab=sxcab-cf*rab(ir)*r3
          vab(ir)=vab(ir)-cf*r3
          v0(ir)=v0(ir)-cf*third*rab(ir)/(r3*r3)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interxc
  ! ==================================================================
  SUBROUTINE vcomp(vcas,vab)
    ! ==--------------------------------------------------------------==
    ! ==  COMPILES THE TOTAL POTENTIALS FOR THE CAS22 METHOD          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vcas(fpar%nnr1,3), &
                                                vab(fpar%nnr1,3)

    INTEGER                                  :: ir
    REAL(real_8)                             :: ca, sa, x, y, z

! ==--------------------------------------------------------------==

    ca=COS(ener_d%casang)
    sa=SIN(ener_d%casang)
    x=sa*sa
    y=1._real_8-x
    z=-2._real_8*ca*sa
    !$omp parallel do private(IR) shared(X,Y,Z)
    DO ir=1,fpar%nnr1
       vcas(ir,1)=y*vcas(ir,1)+x*vab(ir,3)+z*vab(ir,1)
       vcas(ir,2)=y*vcas(ir,2)+2._real_8*x*vab(ir,3)
       vcas(ir,3)=y*vcas(ir,3)
    ENDDO
    !$omp parallel do private(IR) shared(Z)
    DO ir=1,fpar%nnr1
       vab(ir,1)=0.5_real_8*z*vab(ir,2)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vcomp
  ! ==================================================================
  SUBROUTINE ctfor(fion)
    ! ==--------------------------------------------------------------==
    ! ==  COMPILES THE TOTAL FORCES FOR THE CAS22 METHOD              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: ca, sa, x, y, z

! ==--------------------------------------------------------------==

    ca=COS(ener_d%casang)
    sa=SIN(ener_d%casang)
    x=sa*sa
    y=1._real_8-x
    z=-2._real_8*ca*sa
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          fion(1,ia,is)=y*fion(1,ia,is)+x*fion_a(1,ia,is)&
               +z*fion_ab(1,ia,is)
          fion(2,ia,is)=y*fion(2,ia,is)+x*fion_a(2,ia,is)&
               +z*fion_ab(2,ia,is)
          fion(3,ia,is)=y*fion(3,ia,is)+x*fion_a(3,ia,is)&
               +z*fion_ab(3,ia,is)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ctfor
  ! ==================================================================
  SUBROUTINE pcasfor(fion,v,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: v(fpar%nnr1), eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhets, &
                                                txx, tyy, tzz, vcgs, vtot
    INTEGER                                  :: ia, ig, ig1, is, isa
    REAL(real_8)                             :: omtp

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
                rhets=CONJG(v(nzh(ig)))
                gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                vtot=(scg(ig)*CONJG(eirop(ig))+vps(is,ig))*rhets
                txx=vtot*gx
                tyy=vtot*gy
                tzz=vtot*gz
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ELSE
             DO ig=ig1,ncpw%nhg
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                rhets=CONJG(v(nzh(ig)))
                gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
                gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
                gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
                vtot=(scg(ig)*CONJG(eirop(ig))+vps(is,ig))*rhets
                txx=vtot*gx
                tyy=vtot*gy
                tzz=vtot*gz
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ENDIF
       ENDDO
    ENDDO
#else
    DO ig=ig1,ncpw%nhg
       rhets=CONJG(v(nzh(ig)))
       vcgs=scg(ig)*CONJG(eirop(ig))
       gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
       gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
       gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
       isa=0
       DO is=1,ions1%nsp
          vtot=(vcgs+vps(is,ig))*rhets
          txx=vtot*gx
          tyy=vtot*gy
          tzz=vtot*gz
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
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pcasfor
  ! ==================================================================
  SUBROUTINE potab(epotab,vab)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==        POTENTIAL CONTRIBUTION TO THE PENALTY FUNCTION        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: epotab, vab(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'potab'

    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:)
    INTEGER                                  :: ierr, ig
    LOGICAL                                  :: debug

    debug=.FALSE.
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL eicalc(eivps,eirop)
    ! 
    epotab=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    DO ig=1,ncpw%nhg
       epotab=epotab+2._real_8*lspin4%apenal*CONJG(vab(ig))*eivps(ig)
       vab(ig)=0.5_real_8*lspin4%apenal*eivps(ig)
    ENDDO
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE potab
  ! ==================================================================

END MODULE vofrhoc_utils
