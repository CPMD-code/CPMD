MODULE td_force_utils
  USE afbdr_utils,                     ONLY: afbdr,&
                                             give_scr_afbdr
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: kr1m
  USE fftmain_utils,                   ONLY: fwfftn
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03,&
                                             lrsym,&
                                             nolr,&
                                             td01,&
                                             td03,&
                                             tshl
  USE lr_ortho_utils,                  ONLY: lr_ortho
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE opt_lr_utils,                    ONLY: give_scr_opt_lr,&
                                             opt_lr
  USE ovlap_utils,                     ONLY: ovlap,&
                                             ovlap_add
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: rhoo
  USE rho1ofr_utils,                   ONLY: rho1ofr,&
                                             rhoabofr
  USE rnlfor_utils,                    ONLY: rnlfor
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rotate
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             group,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho1_utils,                  ONLY: give_scr_vtdofrho1,&
                                             v1ofrho1,&
                                             vtdofrho1
  USE vpsi_utils,                      ONLY: vpsimt
  USE vtd2_utils,                      ONLY: give_scr_vtd2,&
                                             vtd2,&
                                             vtd2t
  USE wrgeo_utils,                     ONLY: wrgeof
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_force
  PUBLIC :: give_scr_td_force
  !public :: tdrhof

CONTAINS

  ! ==================================================================
  SUBROUTINE td_force(c0,c1,c2,sc0,eigv,rhoe,psi,eigt,nstate,nroot,&
       tau0,fion,orbital,kprint)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8)                             :: eigv(:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigt(*)
    INTEGER                                  :: nstate, nroot
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nolr,nroot,*)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_force'

    COMPLEX(real_8)                          :: edummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: cr(:,:), cz(:,:), gde(:), &
                                                pme(:)
    INTEGER                                  :: i, ierr, ig, il_rhoe, isub, &
                                                it, lddxc_1d, lddxc_2d, ngde, &
                                                npme, nvpp
    REAL(real_8)                             :: e2(5), eind, t1, t2, trace, &
                                                vp, wk_temp(1)
    REAL(real_8), ALLOCATABLE                :: ddxc(:,:), focc(:), &
                                                r1mat(:,:), tdfi(:,:,:), &
                                                vp2(:,:), vpp(:)
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset(procedureN,isub)
    ! ==================================================================
    t1=m_walltime()
    wk_temp(1)=1.0_real_8
    ! ALLOCATE MEMORY
    ALLOCATE(r1mat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cz(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tdfi(3,maxsys%nax,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%parent .AND. kprint.GE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,/)') " CALCULATION OF TDDFT FORCES"
    ENDIF
    ! ==--------------------------------------------------------------==
    lr01%lopti=lr01%lopti+10
    IF (lr03%txc_analytic) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=2*(2*clsd%nlsd-1)
       it=2!nnr1+1
       IF (cntl%tlsd) it=4!3*nnr1+1
    ELSEIF (lr01%lopti.EQ.0 .OR. lr01%lopti.EQ.2) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=(2*clsd%nlsd-1)
       it=1
    ELSE
       lddxc_1d=1
       lddxc_2d=1
       it=1
    ENDIF
    ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL lr_xcpot(ddxc(:,1:it-1) ,rhoo,.FALSE.)
    CALL lr_xcpot(ddxc(:,it:),rhoo,.TRUE. )
    ALLOCATE(focc(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO i=1,nstate
       focc(i)=1._real_8
    ENDDO
    ! R1mat=x*x
    CALL ovlap(nstate,r1mat,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1))
    IF (.NOT.td03%tda) CALL ovlap_add(nstate,r1mat,c1(1,1,td01%fstate,2),&
         c1(1,1,td01%fstate,2))
    CALL mp_sum(r1mat,nstate*nstate,parai%allgrp)
    ! c*R1mat
    CALL rotate(-1._real_8,c0,0._real_8,c2,r1mat,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ! domega/dR
    CALL zeroing(tdfi)!,3*maxsys%nax*ions1%nsp)
    CALL afbdr(c0,c2,focc,psi(:,1),nstate,tdfi)
    CALL rnlsm(c1(:,:,td01%fstate,1),nstate,1,1,.TRUE.)
    CALL rnlfor(tdfi,focc,wk_temp,nstate,1)
    IF (.NOT.td03%tda) THEN
       CALL rnlsm(c1(:,:,td01%fstate,2),nstate,1,1,.TRUE.)
       CALL rnlfor(tdfi,focc,wk_temp,nstate,1)
    ENDIF
    CALL zeroing(rhoe)!,clsd%nlsd*nnr1)
    IF (cntl%tlsd) THEN
       CALL rhoabofr(spin_mod%nsup,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1),rhoe(:,1),psi(:,1))
       IF (.NOT.td03%tda) CALL rhoabofr(spin_mod%nsup,c1(:,:,td01%fstate,2),&
            c1(:,:,td01%fstate,2),rhoe(:,1),psi(:,1))
       CALL rhoabofr(spin_mod%nsup,c2(:,:),c0(:,:),rhoe(:,1),psi(:,1))
       CALL rhoabofr(spin_mod%nsdown,c1(:,spin_mod%nsup+1:,td01%fstate,1),&
            c1(:,spin_mod%nsup+1:,td01%fstate,1),rhoe(:,2),psi(:,1))
       IF (.NOT.td03%tda) CALL rhoabofr(spin_mod%nsdown,c1(:,spin_mod%nsup+1:,td01%fstate,2),&
            c1(:,spin_mod%nsup+1:,td01%fstate,2),rhoe(:,2),psi(:,1))
       CALL rhoabofr(spin_mod%nsdown,c2(:,spin_mod%nsup+1:),c0(:,spin_mod%nsup+1:),rhoe(:,2),psi(:,1))
       CALL daxpy(fpar%nnr1,1._real_8,rhoe(1,2),1,rhoe,1)
    ELSE
       CALL rhoabofr(nstate,c1(:,:,td01%fstate,1),c1(:,:,td01%fstate,1),rhoe(:,1),psi(:,1))
       IF (.NOT.td03%tda) CALL rhoabofr(nstate,c1(:,:,td01%fstate,2),&
            c1(:,:,td01%fstate,2),rhoe(:,1),psi(:,1))
       CALL rhoabofr(nstate,c2(:,:),c0(:,:),rhoe(:,1),psi(:,1))
    ENDIF
    CALL tdrhof(tdfi,rhoe(:,1),psi(:,1))
    IF (.NOT.td03%tda) THEN
       CALL dscal(3*maxsys%nax*ions1%nsp,-0.5_real_8,tdfi,1)
    ELSE
       CALL dscal(3*maxsys%nax*ions1%nsp,-1.0_real_8,tdfi,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! RHS of Z-Matrix equation (domega/dc)
    ! ==--------------------------------------------------------------==
    CALL v1ofrho1(e2,rhoe,ddxc,psi)
    eind = -e2(1)-e2(2)
    CALL zeroing(cr)!,ngw*nstate)
    CALL vpsimt(c0,cr,crge%f(:,1),rhoe,psi(:,1),nstate,clsd%nlsd,.FALSE.)
    IF (.NOT.td03%tda) CALL dscal(fpar%nnr1,0.5_real_8,rhoe,1)
    ! 
    ! d(XWC)/dC
    il_rhoe   = MAX(fpar%krx,kr1m*group%nogrp)*fpar%kr2s*fpar%kr3s !* clsd%nlsd
    ALLOCATE(vp2(il_rhoe,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rho1ofr(c0,c1(:,:,td01%fstate,1),focc,rhoe,psi(:,1),nstate)
    IF (lrsym.EQ.3) THEN
       CALL vtd2t(vp2,rhoe,psi,.TRUE.)
    ELSE
       CALL vtd2(vp2,rhoe,psi,.TRUE.)
    ENDIF
    CALL vpsimt(c0,cr,focc,vp2,psi(:,1),nstate,clsd%nlsd,.TRUE.)
    DEALLOCATE(vp2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL vtdofrho1(e2,rhoe,ddxc(:,it:),psi,.TRUE.)
    eind = eind-e2(1)-e2(2)
    CALL vpsimt(c1(:,:,td01%fstate,1),cr,focc,rhoe,psi(:,1),nstate,clsd%nlsd,.TRUE.)
    CALL zeroing(c2)!,ngw*nstate)
    CALL vpsimt(c0,c2,focc,rhoe,psi(:,1),nstate,clsd%nlsd,.TRUE.)
    CALL ovlap(nstate,r1mat,c0,c2)
    CALL mp_sum(r1mat,nstate*nstate,parai%allgrp)
    CALL rotate(-1._real_8,c1(:,:,td01%fstate,1),1._real_8,cr,r1mat,nstate,2*ncpw%ngw,&
         cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ! 
    CALL dscal(2*ncpw%ngw*nstate,-1._real_8,cr,1)
    CALL lr_ortho(nstate,c0,cr)
    ! 
    ! Solve Z-Matrix equation
    nvpp = ncpw%ngw
    IF (lr01%lopti.EQ.0) THEN
       npme = MAX((ncpw%ngw*nstate+8)*cnti%mdiis/2,ncpw%ngw*nstate) !vw rm *2
       npme = MAX(ncpw%ngw*nstate*lr01%mldiis,npme) !vw rm *2
       ngde = MAX(((ncpw%ngw*nstate+8)*cnti%mdiis)/2,1) !vw rm *2
       ngde = MAX(ncpw%ngw*nstate*lr01%mldiis,ngde)!vw rm *2
    ELSEIF (lr01%lopti.EQ.1) THEN
       npme = 1
       ngde = 1
    ELSE IF (lr01%lopti.EQ.2) THEN
       npme = ncpw%ngw*nstate !vw rm *2
       ngde = 1
    ELSE IF (lr01%lopti.EQ.3) THEN
       npme = (ncpw%ngw*nstate+8)*cnti%mdiis/2 !vw rm *2
       ngde = ((ncpw%ngw*nstate+8)*cnti%mdiis)/2 !vw rm *2
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' WRONG OPTION FOR LINEAR RESPONSE OPTIMIZATION'
       CALL stopgm('TD_FORCE',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(pme(npme),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(ngde),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nvpp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..always use preconditioner
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    cntl%prec=.TRUE.
    CALL ksdiag(vpp)
    trace=dasum(nstate,eigv,1)/REAL(nstate,kind=real_8)
    DO ig=1,ncpw%ngw
       vp=vpp(ig)+trace
       vpp(ig)=ABS(vp/(vp**2+lr02%lr_hthrs**2))
    ENDDO
    CALL opt_lr(c0,cz,c2,sc0,eind,eigv,rhoe,cr,edummy,edummy,ddxc,vpp,&
         pme,gde,psi,nstate,"ZFORCE",orbital)
    lr01%lopti=lr01%lopti-10
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! Z dF/dR C
    CALL rho1ofr(c0,cz,crge%f(:,1),rhoe,psi(:,1),nstate)
    CALL tdrhof(tdfi,rhoe(:,1),psi(:,1))
    CALL afbdr(cz,c0,crge%f,psi(:,1),nstate,tdfi)
    CALL afbdr(c0,cz,crge%f,psi(:,1),nstate,tdfi)
    CALL mp_sum(tdfi,3*maxsys%nax*ions1%nsp,parai%allgrp)
    ! ..TSH[
    IF (tshl%tdtully.AND.(tshl%s0_sh)) THEN
       ! do nothing
    ELSE
       CALL daxpy(3*maxsys%nax*ions1%nsp,-1.0_real_8,tdfi,1,fion,1)
    ENDIF
    ! ..TSH]
    ! 
    ! MAKE Z Vectors available through C2
    CALL dcopy(2*ncpw%ngw*nstate,cz,1,c2,1)
    ! 
    DEALLOCATE(focc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tdfi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(r1mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.kprint.GE.1) THEN
       CALL wrgeof(tau0,fion)
    ENDIF
    t2=m_walltime()
    IF (paral%parent.AND.kprint.GE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,F10.3)')&
            " TIME FOR TDDFT FORCE CALCULATION [s]",(t2-t1)*0.001_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_force
  ! ==================================================================
  SUBROUTINE give_scr_td_force(ltd_force,nstate,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltd_force, nstate
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lafbdr, lopt, lvtd2, &
                                                lvtdofrho1

    CALL give_scr_afbdr(lafbdr,nstate,tag)
    CALL give_scr_opt_lr(lopt,"ZFORCE",tag)
    CALL give_scr_vtdofrho1(lvtdofrho1,tag)
    CALL  give_scr_vtd2(lvtd2,tag)
    ltd_force=0
    ltd_force=MAX(ltd_force,lafbdr,lopt,lvtdofrho1,lvtd2)+200
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_td_force
  ! ==================================================================
  SUBROUTINE tdrhof(fion,rho,psi)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), rho(:)
    COMPLEX(real_8)                          :: psi(:)

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhets, &
                                                txx, tyy, tzz, vcgs
    INTEGER                                  :: ia, ig, ir, is, isa
    REAL(real_8)                             :: omtp

    DO ir=1,fpar%nnr1
       psi(ir)=CMPLX(rho(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    ! 
    omtp=2._real_8*parm%omega*parm%tpiba
    DO ig=1,ncpw%nhg
       rhets=CONJG(psi(nzh(ig)))
       gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
       gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
       gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
       vcgs=scg(ig)*rhets
       isa=0
       DO is=1,ions1%nsp
          txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gx
          tyy=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gy
          tzz=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gz
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
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tdrhof
  ! ==================================================================

END MODULE td_force_utils
