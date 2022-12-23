MODULE k_forces_driver
  USE andp,                            ONLY: rin0,&
                                             rmix
  USE cnst,                            ONLY: ry
  USE csize_utils,                     ONLY: csize
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE forces_driver,                   ONLY: gscal,&
                                             gscal_c
  USE frsblk_c_utils,                  ONLY: reorder_c
  USE frsblk_utils,                    ONLY: reorder
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hesele_utils,                    ONLY: hesele
  USE hnlmat_utils,                    ONLY: hnlmat
  USE k_hesele_utils,                  ONLY: k_hesele
  USE k_pcgrad_utils,                  ONLY: k_pcgrad
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: kcnorm,&
                                             kgemax,&
                                             rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE nlforce_utils,                   ONLY: give_scr_nlforce,&
                                             nlforce
  USE nlps,                            ONLY: imagp
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE ortho_utils,                     ONLY: ortho,&
                                             preortho
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pcgrad_driver,                   ONLY: pcgrad
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE reshaper,                        ONLY: reshape_inplace
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rnlfl_utils,                     ONLY: rnlfl
  USE rnlfor_utils,                    ONLY: rnlfor
  USE rnlrh_utils,                     ONLY: rnlrh
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE rotate_utils,                    ONLY: rotate
  USE rpiiint_utils,                   ONLY: rpiiint
  USE rswfmod,                         ONLY: rsactive
  USE sfac,                            ONLY: fnl,&
                                             ldf1
  USE sort_utils,                      ONLY: sort2
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE stress_utils,                    ONLY: stress
  USE summat_utils,                    ONLY: give_scr_summat
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE utils,                           ONLY: dgive,&
                                             zclean,&
                                             zclean_k
  USE vdw_utils,                       ONLY: vdw,&
                                             vdw_wf
  USE vdwcmod,                         ONLY: &
       empvdwi, empvdwr, idvdw, ivdw, jvdw, vdwbe, vdwi, vdwl, vdwr, vdwrm, &
       vdwst
  USE vofrho_utils,                    ONLY: vofrho
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k_forces

CONTAINS 
  ! ==================================================================
  SUBROUTINE k_forces(c0,c2,sc0,tau0,pme,gde,vpp,eigv,fion,rhoe,&
       psi,f_save,nstate,nkpoint,lproj,tfor,&
       drhomax,drhonorm,calcrho,istep)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                     THE TOTAL ENERGY                         ==
    ! ==                  THE ELECTRONIC FORCES                       ==
    ! ==                  THE FORCES ON THE IONS                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: pme(*), gde(*)
    REAL(real_8)                             :: vpp(:), fion(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: f_save(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,*)
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,nstate), &
                                                c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate,*)
    INTEGER                                  :: nkpoint
    LOGICAL                                  :: lproj, tfor
    REAL(real_8)                             :: drhomax, drhonorm
    LOGICAL                                  :: calcrho
    INTEGER                                  :: istep

    CHARACTER(*), PARAMETER                  :: procedureN = 'k_forces'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: zee
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:), zgam(:,:)
    COMPLEX(real_8), EXTERNAL                :: zdotc
    INTEGER :: i, ib, idamax, ierr, ik, il_auxc, il_ddia, il_fsc, il_gam, &
      imax, isub, j, lsummat, nfto, nleft, nnx, nnxs, nwfc
    LOGICAL                                  :: debug
    REAL(real_8)                             :: bottom, ee, enband
    REAL(real_8), ALLOCATABLE                :: ddia(:), fsc(:), gam(:,:)
    REAL(real_8), EXTERNAL                   :: dasum, ddot

    CALL tiset('  K_FORCES',isub)

    gnmax=0.0_real_8
    gnorm=0.0_real_8
    ener_com%amu    = -2000.0_real_8
    bottom = 2000.0_real_8
    IF (istep.EQ.1) THEN
       ropt_mod%sdiis=.TRUE.
       ropt_mod%spcg=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    IF (pslo_com%tivan .AND. lproj .AND. cnti%iproj.NE.0) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       il_gam = imagp*nstate*nstate
    ENDIF
    CALL give_scr_summat(lsummat,tag,nstate)
    il_auxc=MAX(il_auxc,lsummat)
    il_auxc=MAX(il_auxc,nstate**2)  ! AUXC space for OVLAP (Laio A.)
    il_fsc=nstate
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE POTENTIAL AND THE FORCE ON THE IONS          ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! == Compute the residual part of the ion-ion interaction (ESR)   ==
    ! == due to the overlap of the smeared ionic charge densities     ==
    ! == corresponding to different atoms.                            ==
    ! ==--------------------------------------------------------------==
    nfto=3*maxsys%nax*maxsys%nsx
    CALL zeroing(fion)!,nfto)
    CALL rpiiint(ener_com%esr,tau0,fion,iteropt%iesr,tfor)
    ! ==--------------------------------------------------------------==
    ! == Calculate the Stress Tensor                                  ==
    ! ==--------------------------------------------------------------==
    IF (ropt_mod%calste) CALL stress(c0,tau0,crge%f,psi(:,1),nstate)
    ! ==--------------------------------------------------------------==
    ! == Compute the LOCAL PARTS of the energy and the forces, which  ==
    ! == depend exclusively on the electronic density.                ==
    ! ==--------------------------------------------------------------==
    nnx=fpar%nnr1*clsd%nlsd
    ! 
    IF (calcrho) THEN
       CALL dcopy(nnx,rin0,1,rhoe,1)
       CALL vofrho(tau0,fion,rhoe,psi(:,:),tfor,ropt_mod%calste)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE ELECTRONIC FORCE                             ==
    ! ==--------------------------------------------------------------==

    IF (.NOT. ropt_mod%convwf) THEN
       gemax = 0.0_real_8
       cnorm = 0.0_real_8
       ener_com%etot  = 0.0_real_8
    ENDIF

    DO ik=1,nkpoint

       CALL rnlsm(c0(:,:,ik),nstate,1,ik,tfor)

       ALLOCATE(gam(imagp*nstate, nstate), stat=ierr)
       IF (ierr /= 0) CALL stopgm('K_FORCES', 'ERROR ALLOCATING GAM',& 
            __LINE__,__FILE__)

       ALLOCATE(auxc(il_auxc), stat=ierr)
       IF (ierr /= 0) CALL stopgm('K_FORCES', 'ERROR ALLOCATING AUXC',& 
            __LINE__,__FILE__)

       ALLOCATE(ddia(il_ddia), stat=ierr)
       IF (ierr /= 0) CALL stopgm('K_FORCES', 'ERROR ALLOCATING DDIA',& 
            __LINE__,__FILE__)

       ALLOCATE(fsc(il_fsc), stat=ierr)
       IF (ierr /= 0) CALL stopgm('K_FORCES', 'ERROR ALLOCATING FSC',& 
            __LINE__,__FILE__)

       CALL zeroing(c2)!,nkpt%ngwk*nstate)

       ! ==--------------------------------------------------------------==
       ! == COMPUTE THE FORCE ON THE ELECTRONIC DEGREES OF FREEDOM DUE   ==
       ! == TO THE LOCAL POTENTIAL (STORED IN RHOE)                      ==
       ! ==--------------------------------------------------------------==

       ! IF(TIVAN) THEN
       ! CALL DCOPY(NNX,RHOE,1,PSI,1) 
       ! ==--------------------------------------------------------------==
       ! == TRANSF. TOTAL POTENTIAL IN G-SPACE                           ==
       ! == PUT THE TOTAL POTENTIAL IN G-SPACE INTO VTEMP                ==
       ! ==--------------------------------------------------------------==
       ! CALL FWFFT(PSI)
       ! CALL ZGTHR(NHG,PSI,VTEMP,NZH)
       ! IF(cntl%tlsd) THEN
       ! CALL FWFFT(PSI(1,2))
       ! CALL ZGTHR(NHG,PSI(1,2),VTEMP(1,2),NZH)
       ! ENDIF
       ! ==--------------------------------------------------------------==
       ! == CALCULATE DEEQ FOR VANDERBILT PP                             ==
       ! == FORCE ON IONS DUE TO THE "VANDERBILT CHARGE"                 ==
       ! ==--------------------------------------------------------------==
       ! IF(cntl%tlsd) THEN
       ! CALL NEWD(FNL(1,1,1,1),DEEQ(1,1,1,1),F(1,1),
       ! *              VTEMP(1,1),FION,NSUP,TFOR)
       ! CALL NEWD(FNL(1,1,1,NSUP+1),DEEQ(1,1,1,2),F(NSUP+1,1),
       ! &              VTEMP(1,2),FION,NSDOWN,TFOR)
       ! ELSE
       ! CALL NEWD(FNL,DEEQ,F,VTEMP,FION,N,TFOR)
       ! ENDIF
       ! ENDIF

       CALL vpsi(c0(:,:,ik),c2,crge%f(:,ik),rhoe,psi(:,1),nstate,ik,clsd%nlsd,.TRUE.)
       rsactive=.FALSE.

       IF (pslo_com%tivan) THEN
          IF (tkpts%tkpnt)  CALL stopgm('K_FORCES',&
               'VANDERBILT AND K-POINTS NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
          IF (ropt_mod%convwf) THEN
             CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)

             IF (geq0) THEN
                CALL zclean(c2,nstate,ncpw%ngw)
             ENDIF

             DEALLOCATE(gam,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(auxc,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(ddia,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(fsc,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)

             IF (cntl%tlsd) THEN
                CALL ev_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),fnl,&
                     ik,spin_mod%nsup,bottom,ener_com%amu,calcrho)
                ib = spin_mod%nsup+1
                CALL ev_ksener(c0(1,ib,ik),c2(1,ib),sc0,eigv(ib,ik),fnl,&
                     ik,spin_mod%nsdown,bottom,ener_com%amu,calcrho)
             ELSE
                CALL ev_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),fnl,&
                     ik,nstate,bottom,ener_com%amu,calcrho)
             ENDIF
             GOTO 100
          ENDIF

          IF (lproj.AND.cnti%iproj.NE.0) THEN
             CALL ovlap(nstate,gam,c2,c0(:,:,ik))
             CALL hnlmat(gam,crge%f,nstate)
             CALL mp_sum(gam,nstate*nstate,parai%allgrp)

             ! H TRACE
             ee  = 0.0_real_8
             !$omp parallel do private(I) reduction(+:EE)
             DO i=1,nstate
                ee=ee+(-gam(i,i))
             ENDDO
             ener_com%etot = ener_com%etot + ee

             IF (tfor .AND. calcrho)&
                  CALL rnlfl(fion,gam,nstate,nkpoint)
             IF (cnti%iproj.EQ.1 .AND. .NOT. ropt_mod%convwf) THEN
                DO i=1,nstate
                   CALL daxpy(2*nkpt%ngwk,-gam(i,i),c0(1,i,ik),1,c2(1,i),1)
                ENDDO
             ELSEIF (cnti%iproj.EQ.2 .AND. .NOT. ropt_mod%convwf) THEN
                CALL rotate(-1.0_real_8,c0(:,:,ik),1.0_real_8,c2,gam,&
                     nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
             ENDIF
             CALL nlforce(c2,crge%f,gam,auxc,ddia,nstate)
          ELSE
             CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)
          ENDIF
       ELSE
          CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)
          IF (ropt_mod%convwf) THEN
             DEALLOCATE(gam,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(auxc,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(ddia,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(fsc,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             IF (geq0) THEN
                IF (tkpts%tkpnt) THEN
                   CALL zclean_k(c2,nstate,ncpw%ngw)
                ELSE
                   CALL zclean(c2,nstate,ncpw%ngw)
                ENDIF
             ENDIF
             IF (cntl%tlsd) THEN
                IF (tkpts%tkpnt)THEN
                   CALL k_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),&
                        ik,spin_mod%nsup,bottom,ener_com%amu,calcrho)
                   ib = spin_mod%nsup+1
                   CALL k_ksener(c0(1,ib,ik),c2(1,ib),sc0,eigv(ib,ik),&
                        ik,spin_mod%nsdown,bottom,ener_com%amu,calcrho)
                ELSE
                   CALL ev_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),fnl,&
                        ik,spin_mod%nsup,bottom,ener_com%amu,calcrho)
                   ib = spin_mod%nsup+1
                   CALL ev_ksener(c0(1,ib,ik),c2(1,ib),sc0,eigv(ib,ik),fnl,&
                        ik,spin_mod%nsdown,bottom,ener_com%amu,calcrho)
                ENDIF
             ELSE
                IF (tkpts%tkpnt)THEN
                   CALL k_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),&
                        ik,nstate,bottom,ener_com%amu,calcrho)
                ELSE
                   CALL ev_ksener(c0(1,1,ik),c2,sc0,eigv(1,ik),fnl,&
                        ik,nstate,bottom,ener_com%amu,calcrho)
                ENDIF
             ENDIF
             GOTO 100
          ENDIF

          ! ETOT
          IF (.NOT. calcrho) THEN
             ee  = 0.0_real_8
             IF (tkpts%tkpnt) THEN
                DO i=1,nstate
                   ee=ee-ddot(2*nkpt%ngwk,c0(1,i,ik),1,c2(1,i),1)*wk(ik)
                ENDDO
                CALL mp_sum(ee,parai%allgrp)
                ener_com%etot = ener_com%etot + ee
             ELSE
                DO i=1,nstate
                   ee=ee-dotp(ncpw%ngw,c0(:,i,1),c2(:,i))
                ENDDO
                CALL mp_sum(ee,parai%allgrp)
                ener_com%etot = ener_com%etot + ee
             ENDIF
          ENDIF
          IF (lproj) THEN
             DO i=1,nstate
                IF (crge%f(i,ik).EQ.0._real_8) THEN
                   fsc(i)=0._real_8
                ELSE
                   fsc(i)=1._real_8
                ENDIF
             ENDDO
             IF (cnti%iproj.EQ.0) THEN
             ELSEIF (cnti%iproj.EQ.1) THEN
                IF (tkpts%tkpnt) THEN
                   DO i=1,nstate
                      zee=fsc(i)*zdotc(nkpt%ngwk,c0(1,i,ik),1,c2(1,i),1)
                      CALL mp_sum(zee,parai%allgrp)
                      CALL zaxpy(nkpt%ngwk,-zee,c0(1,i,ik),1,c2(1,i),1)
                   ENDDO
                ELSE
                   DO i=1,nstate
                      ee=fsc(i)*dotp(ncpw%ngw,c0(:,i,ik),c2(:,1))
                      CALL mp_sum(ee,parai%allgrp)
                      CALL daxpy(2*ncpw%ngw,-ee,c0(1,i,ik),1,c2(1,1),1)
                   ENDDO
                ENDIF
             ELSEIF (cnti%iproj.EQ.2) THEN
                IF (tkpts%tkpnt) THEN
                   ALLOCATE(zgam(nstate,nstate), stat=ierr)
                   IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem',&
                        __LINE__,__FILE__)
                   CALL ovlap_c(nstate,zgam,c2,c0(1,1,ik))
                   CALL mp_sum(zgam,nstate*nstate,parai%allgrp)
                   CALL gscal_c(nstate,zgam,fsc)
                   CALL rotate_c(CMPLX(-1.0_real_8,0._real_8,kind=real_8),c0(1,1,ik),&
                        CMPLX( 1.0_real_8,0._real_8,kind=real_8),c2,zgam,nstate)
                   DEALLOCATE(zgam, stat=ierr)
                   IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem',&
                        __LINE__,__FILE__)
                ELSE
                   CALL ovlap(nstate,gam,c2,c0(:,:,ik))
                   CALL mp_sum(gam,nstate*nstate,parai%allgrp)
                   CALL gscal(nstate,gam,fsc)
                   CALL rotate(-1.0_real_8,c0(:,:,ik),1.0_real_8,c2,gam,&
                        nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       IF (geq0) THEN
          IF (tkpts%tkpnt) THEN
             CALL zclean_k(c2,nstate,ncpw%ngw)
          ELSE
             CALL zclean(c2,nstate,ncpw%ngw)
          ENDIF
       ENDIF

       DEALLOCATE(gam,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(auxc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(ddia,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(fsc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       IF (tkpts%tkpnt) THEN
          CALL csize(c2,nstate,kgemax(ik),kcnorm(ik))
          gemax = MAX(gemax,kgemax(ik))
          cnorm = MAX(cnorm,kcnorm(ik))
       ELSE
          CALL csize(c2,nstate,gemax,cnorm)
       ENDIF

       IF (tkpts%tkpnt)THEN
          IF (cntl%diis) THEN
             CALL k_hesele(dt2bye,vpp,ik)
             CALL k_odiis(c0(:,:,ik),c2,vpp,ik,nstate,pme,gde,&
                  dt2bye,ropt_mod%sdiis,calcrho)
          ELSE IF (cntl%pcg) THEN
             CALL k_hesele(dt2bye,vpp,ik)
             CALL k_pcgrad(c0(1,1,ik),c2,sc0,vpp,pme,&
                  rhoe,psi(:,1),tau0,nstate,ropt_mod%spcg,ik)
          ENDIF
       ELSE
          IF (cntl%diis) THEN
             IF (ropt_mod%sdiis) CALL hesele(dt2bye,vpp)
             CALL odiis(c0,c2,vpp,nstate,pme,gde,dt2bye,ropt_mod%sdiis)
          ELSE IF (cntl%pcg) THEN
             IF (ropt_mod%spcg) CALL hesele(dt2bye,vpp)
             CALL pcgrad(c0(:,:,1),c2,sc0,vpp,pme,&
                  rhoe,psi,tau0,nstate,ropt_mod%spcg)
          ENDIF
       ENDIF

       ! ==--------------------------------------------------------------==
       ! ==  ORTHOGONALIZATION                                           ==
       ! ==--------------------------------------------------------------==
       IF (cntl%nonort)  THEN
          IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
       ELSE
          CALL preortho(c0(1,1,ik),nstate)
          IF (pslo_com%tivan) CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
          CALL ortho(nstate,c0(:,:,ik),c2)
       ENDIF
100    CONTINUE
    ENDDO    ! KPOINT LOOP
    ! do i = 1,10
    ! write(6,'(i3,10f12.6)') i,(real(c0(ig,i,1)),ig=1,10)
    ! enddo
    ! stop 'c0'
    ! ==--------------------------------------------------------------==
    ! == Calculate NON-LOCAL PARTS of the energy and forces, which    ==
    ! == depend on the wavefunctions themselves, and not on the       ==
    ! == electronic density.                                          ==
    ! ==--------------------------------------------------------------==

    CALL rnlrh(ener_com%enl,nstate,nkpoint)
    IF (tfor .AND. calcrho) THEN
       CALL rnlfor(fion,crge%f,wk,nstate,nkpoint)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions from van der Waals interaction  ==
    ! ==--------------------------------------------------------------==
    vdwr%evdw=0._real_8
    CALL zeroing(vdwr%devdw)!,6)
    IF (vdwl%vdwc) CALL vdw(tau0,empvdwi%nvdw,idvdw,ivdw,jvdw,vdwst,vdwrm,vdwbe,&
         empvdwr%vdweps,empvdwr%s6grim,vdwi%nxvdw,vdwi%nyvdw,vdwi%nzvdw,vdwr%evdw,fion,vdwr%devdw)
    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions for Wannier-based van der Waals ==
    ! ==--------------------------------------------------------------==
    IF (vdwl%vdwd) THEN
       vdwr%evdw=0._real_8
       CALL zeroing(vdwr%devdw)!,6)
       nwfc=nstate
       IF (tfor) CALL vdw_wf(tau0,fion,nwfc)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == TOTAL ENERGY                                                 ==
    ! ==--------------------------------------------------------------==
    IF (calcrho) THEN
       ener_com%eeig = 0._real_8
       ener_com%etot = ener_com%ekin + ener_com%eht + ener_com%epseu + ener_com%enl + ener_com%exc + vdwr%evdw

       ! Sum ETOT,EKIN,EPSEU,ENL,EHT,EHEP,EHEE,EHII,EXC,VXC,EGC
       CALL mp_sum(ener_com%etot,parai%allgrp)
       CALL mp_sum(ener_com%ekin,parai%allgrp)
       CALL mp_sum(ener_com%epseu,parai%allgrp)
       CALL mp_sum(ener_com%enl,parai%allgrp)
       CALL mp_sum(ener_com%eht,parai%allgrp)
       CALL mp_sum(ener_com%ehep,parai%allgrp)
       CALL mp_sum(ener_com%ehee,parai%allgrp)
       CALL mp_sum(ener_com%ehii,parai%allgrp)
       CALL mp_sum(ener_com%exc,parai%allgrp)
       CALL mp_sum(ener_com%vxc,parai%allgrp)
       CALL mp_sum(ener_com%egc,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  BAND ENERGY (only at the end)                               ==
    ! ==--------------------------------------------------------------==
    IF (ropt_mod%convwf .AND. paral%parent .AND. calcrho) THEN
       IF (paral%io_parent)&
            WRITE(6,'("BOTTOM LEVEL",F14.7)') bottom*2*ry
       IF (paral%io_parent)&
            WRITE(6,'("MINIMUM CB",F14.7)') ener_com%amu*2*ry
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') ' EIGENVALUES(EV) AND OCCUPATION:'
       enband = 0.0_real_8
       IF (cntl%tlsd) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)') '_____ SPIN UP'
          DO ik=1,nkpt%nkpts
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I8,4F12.6)')&
                  'K POINT:',ik,(rk(i,ik),i=1,3),wk(ik)
             DO i=1,spin_mod%nsup,2
                nleft=MIN(spin_mod%nsup-i,1)
                IF (paral%io_parent)&
                     WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                     (i+j,eigv(i+j,ik)*(2*ry),&
                     crge%f(i+j,ik),j=0,nleft)
                enband = enband +  eigv(i,ik)*crge%f(i,ik)*wk(ik)
                IF (i+1 .LE. nstate)&
                     enband = enband +  eigv(i+1,ik)*crge%f(i+1,ik)*wk(ik)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(/,A)') '_____ SPIN DOWN'
          ib = spin_mod%nsup
          DO ik=1,nkpt%nkpts
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I8,4F12.6)')&
                  'K POINT:',ik,(rk(i,ik),i=1,3),wk(ik)
             DO i=1,spin_mod%nsdown,2
                nleft=MIN(spin_mod%nsdown-i,1)
                IF (paral%io_parent)&
                     WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                     (i+j,eigv(ib+i+j,ik)*(2*ry),&
                     crge%f(ib+i+j,ik),j=0,nleft)
                enband = enband +  eigv(ib+i,ik)*crge%f(ib+i,ik)*wk(ik)
                IF (i+1 .LE. nstate)&
                     enband = enband +  eigv(ib+i+1,ik)*crge%f(ib+i+1,ik)*wk(ik)
             ENDDO
          ENDDO
       ELSE
          DO ik=1,nkpt%nkpts
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I8,4F12.6)')&
                  'K POINT:',ik,(rk(i,ik),i=1,3),wk(ik)
             DO i=1,nstate,2
                nleft=MIN(nstate-i,1)
                IF (paral%io_parent)&
                     WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                     (i+j,eigv(i+j,ik)*(2*ry),&
                     crge%f(i+j,ik),j=0,nleft)
                enband = enband +  eigv(i,ik)*crge%f(i,ik)*wk(ik)
                IF (i+1 .LE. nstate)&
                     enband = enband +  eigv(i+1,ik)*crge%f(i+1,ik)*wk(ik)
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(" >>>>> BAND  ENERGY",F14.7)') enband
    ENDIF

    IF (calcrho) THEN
       rsactive=cntl%krwfn
       CALL rhoofr_c(c0,rhoe,psi(:,1),nstate)
       CALL dcopy(nnx,rhoe,1,rmix,1)! TODO check this
       CALL daxpy(nnx,-1.0_real_8,rin0,1,rmix,1)! TODO check this
       imax=idamax(nnx,rmix(1,1),1)! TODO check this
       drhomax=ABS(dgive(rmix(1,1),imax))! TODO check this
       drhonorm=dasum(nnx,rmix(1,1),1)! TODO check this
       CALL mp_sum(drhonorm,parai%allgrp)
       CALL mp_max(drhomax,parai%allgrp)
       nnxs=fpar%kr1s*fpar%kr2s*fpar%kr3s
       drhomax=SQRT(parm%omega*drhomax/REAL(nnxs,kind=real_8))
       drhonorm=parm%omega*drhonorm/REAL(nnxs,kind=real_8)
       CALL dcopy(nnx,rhoe,1,rin0,1)
    ENDIF

    IF (tfor .AND. calcrho) THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
       IF (paral%parent) THEN
          CALL symvec(fion)
          CALL taucl(fion)
          CALL gsize(fion,gnmax,gnorm)
       ENDIF
    ENDIF
    CALL tihalt('  K_FORCES',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_forces
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE ev_ksener(c0,c2,sc0,eigv,fnla,ik,nstate,&
       bottom,amu,calcrho)
    ! ==--------------------------------------------------------------==
    ! ==                          CALCULATES:                         ==
    ! ==       Eigenvalues and eigenvectors of the Hamiltonian        ==
    ! ==  Rotates the optimized wfns in order to obtain the KS states ==
    ! ==--------------------------------------------------------------==

    ! Arguments:
    REAL(real_8)                             :: fnla(ldf1,maxsys%nhxs,*)
    INTEGER                                  :: ik, nstate
    REAL(real_8)                             :: eigv(nstate)
    COMPLEX(real_8), TARGET                  :: sc0(nkpt%ngwk,nstate)
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    REAL(real_8)                             :: bottom, amu
    LOGICAL                                  :: calcrho

    CHARACTER(*), PARAMETER                  :: procedureN = 'ev_ksener'

    COMPLEX(real_8), ALLOCATABLE             :: work(:)
    INTEGER                                  :: i, ierr, j, jj, lwork
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: TLSD_tmp
    REAL(real_8), ALLOCATABLE                :: aux(:), hamilt(:,:)
    REAL(real_8), POINTER                    :: fnlb(:,:,:)

    lwork   = MAX(1,2*nstate-1)

    ! TODO align for BG
    ! TODO deallocate these arrays 
    ALLOCATE(hamilt(nstate, nstate), stat=ierr)
    ALLOCATE(aux(3*nstate), stat=ierr)
    ALLOCATE(work(lwork), stat=ierr)
    ALLOCATE(INDEX(nstate), stat=ierr)

    CALL zeroing(index)!,nstate)
    ! ==--------------------------------------------------------------==
    CALL zeroing(hamilt)!,nstate*nstate)

    TLSD_tmp = cntl%tlsd
    cntl%tlsd = .FALSE.
    CALL ovlap(nstate,hamilt,c0,c2)
    cntl%tlsd = TLSD_tmp

    CALL mp_sum(hamilt,SIZE(hamilt),parai%allgrp)

    CALL dsyev('V','U',nstate,hamilt,nstate,eigv,&
         work,lwork,ierr)
    IF (lwork.LT.INT(work(1)).AND.paral%io_parent) THEN
       WRITE(6,*) ' WARNING| LWORK SUBOPTIMAL FOR ZHEEV'
       WRITE(6,*) '   LWORK:',lworK
       WRITE(6,*) ' OPTIMAL:',INT(work(1))
    ENDIF
    CALL reorder(nstate,nstate,eigv,hamilt)

    CALL dgemm('N','N',2*nkpt%ngwk,nstate,nstate,1.0_real_8,c0,2*nkpt%ngwk,hamilt,&
         nstate,0.0_real_8,sc0,2*nkpt%ngwk)
    CALL dcopy(2*nkpt%ngwk*nstate,sc0,1,c0,1)
    CALL dgemm('N','N',2*nkpt%ngwk,nstate,nstate,1.0_real_8,c2,2*nkpt%ngwk,hamilt,&
         nstate,0.0_real_8,sc0,2*nkpt%ngwk)
    CALL dcopy(2*nkpt%ngwk*nstate,sc0,1,c2,1)

    CALL sort2(eigv,nstate,index)
    IF (calcrho) THEN
       CALL dscal(nstate,-0.5_real_8,eigv,1)
    ELSE
       CALL dscal(nstate,-1.0_real_8,eigv,1)
    ENDIF
    DO i = 1 , nstate/2
       aux(1) = eigv(nstate-i+1)
       eigv(nstate-i+1) = eigv(i)
       eigv(i) = aux(1)
       jj = INDEX(nstate-i+1)
       INDEX(nstate-i+1) = INDEX(i)
       INDEX(i) = jj
    ENDDO
    DO i=1,nstate
       j=INDEX(i)
       CALL dcopy(2*nkpt%ngwk,sc0(1,j),1,c2(1,i),1)
    ENDDO
    CALL  dcopy(2*nkpt%ngwk*nstate,c0,1,sc0,1)
    DO i=1,nstate
       j=INDEX(i)
       CALL dcopy(2*nkpt%ngwk,sc0(1,j),1,c0(1,i),1)
       amu = MAX(amu,eigv(i))
       bottom = MIN(bottom,eigv(i))
    ENDDO
    IF (pslo_com%tivan) THEN
       CALL reshape_inplace(sc0, (/ldf1,maxsys%nhxs,nstate/), fnlb)
       CALL dgemm('N','N',ldf1*maxsys%nhxs,nstate,nstate,1.0_real_8,fnla,&
            ldf1*maxsys%nhxs,hamilt,nstate,0.0_real_8,fnlb,ldf1*maxsys%nhxs)
       DO i=1,nstate
          j=INDEX(i)
          CALL dcopy(ldf1*maxsys%nhxs,fnlb(1,1,j),1,fnla(1,1,i),1)
       ENDDO
    ENDIF

    DEALLOCATE(hamilt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ev_ksener
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE k_ksener(c0,c2,sc0,eigv,ik,nstate,bottom,amu,&
       calcrho)
    ! ==--------------------------------------------------------------==
    ! ==                          CALCULATES:                         ==
    ! ==       Eigenvalues and eigenvectors of the Hamiltonian        ==
    ! ==  Rotates the optimized wfns in order to obtain the KS states ==
    ! ==--------------------------------------------------------------==

    ! Arguments:
    INTEGER                                  :: ik, nstate
    REAL(real_8)                             :: eigv(nstate)
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,nstate), &
                                                c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    REAL(real_8)                             :: bottom, amu
    LOGICAL                                  :: calcrho

    CHARACTER(*), PARAMETER                  :: procedureN = 'k_ksener'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: hamilt(:,:), work(:)
    INTEGER                                  :: i, ierr, j, jj, lwork
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    REAL(real_8), ALLOCATABLE                :: aux(:)

    lwork   = MAX(1,2*nstate-1)

    ALLOCATE(hamilt(nstate, nstate), stat=ierr)
    ALLOCATE(aux(3*nstate), stat=ierr)
    ALLOCATE(work(lwork), stat=ierr)
    ALLOCATE(INDEX(nstate), stat=ierr)

    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)

    CALL zeroing(index)!,nstate)
    ! ==--------------------------------------------------------------==

    CALL ovlap_c(nstate,hamilt,c0,c2)

    CALL mp_sum(hamilt,SIZE(hamilt),parai%allgrp)

    CALL zheev('V','U',nstate,hamilt,nstate,eigv,&
         work,lwork,aux,ierr)
    IF (lwork.LT.INT(work(1)).AND.paral%io_parent) THEN
       WRITE(6,*) ' WARNING| LWORK SUBOPTIMAL FOR ZHEEV'
       WRITE(6,*) '   LWORK:',lworK
       WRITE(6,*) ' OPTIMAL:',INT(work(1))
    ENDIF
    CALL reorder_c(nstate,nstate,eigv,hamilt)
    CALL zgemm('N','N',nkpt%ngwk,nstate,nstate,zone,c0,nkpt%ngwk,hamilt,&
         nstate,zzero,sc0,nkpt%ngwk)
    CALL dcopy(2*nkpt%ngwk*nstate,sc0,1,c0,1)
    CALL zgemm('N','N',nkpt%ngwk,nstate,nstate,zone,c2,nkpt%ngwk,hamilt,&
         nstate,zzero,sc0,nkpt%ngwk)
    CALL dcopy(2*nkpt%ngwk*nstate,sc0,1,c2,1)

    CALL sort2(eigv,nstate,index)
    IF (calcrho) THEN
       CALL dscal(nstate,-0.5_real_8,eigv,1)
    ELSE
       CALL dscal(nstate,-1.0_real_8,eigv,1)
    ENDIF
    DO i = 1 , nstate/2
       aux(1) = eigv(nstate-i+1)
       eigv(nstate-i+1) = eigv(i)
       eigv(i) = aux(1)
       jj = INDEX(nstate-i+1)
       INDEX(nstate-i+1) = INDEX(i)
       INDEX(i) = jj
    ENDDO
    DO i=1,nstate
       j=INDEX(i)
       CALL dcopy(2*nkpt%ngwk,sc0(1,j),1,c2(1,i),1)
    ENDDO
    CALL  dcopy(2*nkpt%ngwk*nstate,c0,1,sc0,1)
    DO i=1,nstate
       j=INDEX(i)
       CALL dcopy(2*nkpt%ngwk,sc0(1,j),1,c0(1,i),1)
       amu = MAX(amu,eigv(i))
       bottom = MIN(bottom,eigv(i))
    ENDDO

    DEALLOCATE(hamilt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_ksener
  ! ==================================================================
END MODULE k_forces_driver
