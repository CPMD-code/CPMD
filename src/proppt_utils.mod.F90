MODULE proppt_utils
  USE adat,                            ONLY: elem
  USE atomc_utils,                     ONLY: atomc
  USE cnst,                            ONLY: au_deb
  USE condu,                           ONLY: condpa,&
                                             conduct,&
                                             conduct2,&
                                             iucond,&
                                             normcon
  USE conduct_utils,                   ONLY: conductivity
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE core_spect_utils,                ONLY: core_spectra
  USE cores,                           ONLY: coresl
  USE cppt,                            ONLY: indz,&
                                             nzh
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: ddipo,&
                                             give_scr_ddipo
  USE difrho_utils,                    ONLY: difrho
  USE dipo_utils,                      ONLY: dipo,&
                                             rsdipo
  USE dipomod,                         ONLY: moment
  USE dist_prowfn_utils,               ONLY: dist_prowfn
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE elstpo_utils,                    ONLY: elstpo
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE espchg_utils,                    ONLY: atfield,&
                                             espsolv,&
                                             printcg,&
                                             selectp
  USE exdipo_utils,                    ONLY: exdipo,&
                                             transm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_mark,&
                                             fo_new,&
                                             fo_verb
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_utils,                    ONLY: give_scr_forces
  USE g_loc,                           ONLY: glocal
  USE g_loc_dr_utils,                  ONLY: g_loc_dr_comp,&
                                             g_loc_dr_real
  USE geq0mod,                         ONLY: geq0
  USE hip_utils,                       ONLY: give_qphi
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kddipo_utils,                    ONLY: kddipo
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ldos_utils,                      ONLY: ldos
  USE ldosmod,                         ONLY: cldos
  USE localize_utils,                  ONLY: localize
  USE lodipo_utils,                    ONLY: lodipo
  USE lodp,                            ONLY: &
       dmomlo, exd, extd, focc, nsdip, numbld, rcc, trmom, xmaxld, xminld, &
       ymaxld, yminld, zmaxld, zminld
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE molorb_utils,                    ONLY: molorb
  USE mp_interface,                    ONLY: mp_bcast
  USE ortho_utils,                     ONLY: give_scr_ortho,&
                                             ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: rhoo
  USE pola,                            ONLY: alphap,&
                                             ipolarise,&
                                             iupola,&
                                             nfpolar,&
                                             tpolarb,&
                                             tzeff,&
                                             zeff
  USE polarise_utils,                  ONLY: give_scr_polarproj,&
                                             polarproj
  USE potmed_utils,                    ONLY: potmed,&
                                             printpot
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: icubeorb,&
                                             prop1,&
                                             prop2,&
                                             prop4,&
                                             prop7,&
                                             rcut,&
                                             z_11
  USE propin_utils,                    ONLY: propin
  USE prowfn_utils,                    ONLY: prowfn
  USE proylm_utils,                    ONLY: proylm
  USE pslo,                            ONLY: pslo_com
  USE readsr_utils,                    ONLY: xstring
  USE rho1ofr_utils,                   ONLY: rhoabofr
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt
  USE rscpot_utils,                    ONLY: rscpot
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: isetone,&
                                             read_irec,&
                                             write_irec
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE store_types,                     ONLY: irec_pwv,&
                                             restart1,&
                                             rout1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             group,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE vofrho_utils,                    ONLY: vofrho
  USE wann,                            ONLY: wannl
  USE wannier_print_utils,             ONLY: wannier_print
  USE wc_dos_utils,                    ONLY: wc_dos
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: proppt
  PUBLIC :: propcal
  PUBLIC :: give_scr_propcal
  PUBLIC :: prdip
  PUBLIC :: give_scr_espc

CONTAINS

  ! ==================================================================
  SUBROUTINE proppt
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'proppt'

    CHARACTER(len=128)                       :: cubefilename, filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:), cs(:,:), &
                                                eirop(:), eivps(:), psi(:,:), &
                                                qphi(:), sc0(:,:), vtemp(:)
    INTEGER :: i, ia, iat, ie, ierr, ig, il_psi_1d, il_psi_2d, il_qphi, &
      il_rhoe_1d, il_rhoe_2d, ippc, ir, irec(100), isp, j, ji, lo, lr, lrho, &
      lscr, mlen, n_cubefiles, nc, nl, nnat, nstate, nxx
    INTEGER, ALLOCATABLE                     :: isel(:)
    LOGICAL                                  :: oldstatus, statusdummy, tinfo
    REAL(real_8)                             :: cubecenter(3)
    REAL(real_8), ALLOCATABLE :: achrg(:), center(:,:), echrg(:), efield(:), &
      eigv(:), fback(:), ichrg(:), reirop(:), rhoe(:,:), rpsi(:), scr(:)

! ==--------------------------------------------------------------==

    IF ((cntl%tqmmm.AND.paral%parent).AND.paral%io_parent)&
         WRITE(6,*) 'PROPPT| WARNING QM/MM UNTESTED'
    ! ==--------------------------------------------------------------==
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (soft_com%exsoft) RETURN
    nstate=crge%n
    tinfo=.FALSE.
    condpa%tconduct=.FALSE.
    tpolarb=.FALSE.
    CALL propin(tinfo)
    ! ==--------------------------------------------------------------==
    ! If not required property to calculate -> Return
    IF (.NOT.tinfo)RETURN
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",19X,A,20X,"*")') ' PROPERTY CALCULATIONS '
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"),/)')
    ENDIF
    IF (prop1%pwfn.AND.(wannl%twann.OR.prop1%locl)) CALL stopgm('PROPPT','PROJECT '&
         // 'WAVEFUNCTION AND LOCALIZATION ARE INCOMPATIBLE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL setfftn(0)
    prop2%numorb = MAX(crge%n,prop2%numorb)
    IF (tpolarb) prop2%numorb=MAX(prop2%numorb,3)
    ALLOCATE(c0(nkpt%ngwk,prop2%numorb,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    lenbk=0
    nc=2*nkpt%ngwk*prop2%numorb
    IF (prop1%locl) THEN
       lenbk=nxxfun(prop2%numorb)
       nc=MAX(2*lenbk*parai%cp_nproc,nc)
    ELSEIF (glocal%tgwannier.OR.prop1%dberry) THEN
       lenbk=nxxfun(prop2%numorb)
       nc=MAX(2*2*lenbk*parai%cp_nproc,nc)
    ENDIF
    ALLOCATE(cs(nkpt%ngwk,prop2%numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mm_dim(mm_go_mm,oldstatus)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (cntl%tdiag.OR.coresl%tcores.OR.condpa%tconduct.OR.cldos%tldos) THEN
       ALLOCATE(rhoo(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigv(nstate*nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! Avoid 'not allocated' runtime error
       ALLOCATE(rhoo(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eigv(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL fnlalloc(prop2%numorb,.FALSE.,.FALSE.)
    ! ==--------------------------------------------------------------==
    ! ==  LOAD DATA                                                   ==
    ! ==--------------------------------------------------------------==
    IF (prop1%ceig) restart1%reigv=.TRUE.
    CALL read_irec(irec)
    CALL mm_dim(mm_go_mm,statusdummy)
    CALL zhrwf(1,irec,c0,cs,prop2%numorb,eigv,tau0,taup,taup,iteropt%nfi)
    CALL mp_bcast(tau0,SIZE(tau0),parai%io_source,parai%cp_grp)
    CALL mm_dim(mm_go_qm,statusdummy)
    ALLOCATE(fback(prop2%numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (prop2%numorb.GT.crge%n) THEN
       CALL zeroing(fback)!,prop2%numorb)
       IF (cntl%tlsd) THEN
          nl=(prop2%numorb-crge%n)/2
          CALL dcopy(spin_mod%nsup,crge%f(1,1),1,fback(1),1)
          CALL dcopy(spin_mod%nsdown,crge%f(spin_mod%nsup+1,1),1,fback(spin_mod%nsup+nl+1),1)
          CALL dcopy(2*ncpw%ngw*spin_mod%nsup,c0(1,1,1),1,cs(1,1),1)
          CALL dcopy(2*ncpw%ngw*nl,c0(1,spin_mod%nsup+spin_mod%nsdown+1,1),1,cs(1,spin_mod%nsup+1),1)
          CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,c0(1,spin_mod%nsup+1,1),1,cs(1,spin_mod%nsup+nl+1),1)
          CALL dcopy(2*ncpw%ngw*nl,c0(1,spin_mod%nsup+spin_mod%nsdown+nl+1,1),1,&
               cs(1,spin_mod%nsup+spin_mod%nsdown+nl+1),1)
          CALL dcopy(2*ncpw%ngw*prop2%numorb,cs(1,1),1,c0(1,1,1),1)
          spin_mod%nsup=spin_mod%nsup+nl
          spin_mod%nsdown=spin_mod%nsdown+nl
       ELSE
          CALL dcopy(crge%n,crge%f(1,1),1,fback(1),1)
       ENDIF
       DEALLOCATE(crge%f,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(crge%f(prop2%numorb,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(prop2%numorb,fback(1),1,crge%f(1,1),1)
    ENDIF
    DEALLOCATE(fback,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE DIPOLE MOMENTS                                    ==
    ! ==--------------------------------------------------------------==
    IF (prop1%dberry) THEN
       IF (tkpts%tkpnt) THEN
          ! TODO move C2, SC0 into KDDIPO
          ALLOCATE(c2(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(c2)
          ALLOCATE(sc0(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL kddipo(tau0,c0,c2,cs,sc0,prop2%numorb)

          DEALLOCATE(c2,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(sc0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(c2(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(c2)
          ALLOCATE(sc0(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(sc0)
          ALLOCATE(center(4,nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(center)
          CALL ddipo(tau0,c0(:,:,1),c2,cs,sc0,prop2%numorb,center)
          IF (paral%parent) CALL prdip
          DEALLOCATE(c2,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(sc0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(center,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (prop1%ldip.OR.prop1%locd.OR.prop1%lext.OR.prop1%ltrm) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ==------------------------------------------------------------==
       CALL phfac(tau0)
       CALL rnlsm(c0(:,:,1),prop2%numorb,1,1,.FALSE.)
       IF (prop1%lext) THEN
          ALLOCATE(fback(prop2%numorb),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(prop2%numorb,crge%f(1,1),1,fback(1),1)
          DO i=1,extd%nexdip
             DO j=1,prop2%numorb
                ji=(i-1)*prop2%numorb+j
                crge%f(j,1)=fback(j)-focc(ji)
             ENDDO
             CALL difrho(c0,rhoe,psi,prop2%numorb)
             CALL zeroing(eirop)!,nhg)
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IR)
             DO ir=1,fpar%nnr1
                psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
             ENDDO
             CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
             CALL exdipo(i,tau0,psi,eirop)
          ENDDO
          CALL dcopy(prop2%numorb,fback(1),1,crge%f(1,1),1)
       ENDIF
       IF (prop1%ltrm) THEN
          DO i=1,extd%ntrmom
             CALL zeroing(rhoe(:,1))!,nnr1)
             CALL zeroing(eirop)!,nhg)
             CALL rhoabofr(1,c0(:,nsdip(1,i):nsdip(1,i),1),c0(:,nsdip(2,i):nsdip(2,i),1),rhoe(:,1),psi(:,1))
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IR)
             DO ir=1,fpar%nnr1
                psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
             ENDDO
             CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
             CALL transm(i,tau0,psi,eirop)
          ENDDO
       ENDIF
       IF (prop1%ldip.OR.prop1%locd) THEN
          IF (tkpts%tkpnt) THEN
             CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
          ELSE
             CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
          ENDIF
          CALL eicalc(eivps,eirop)
          IF (prop1%lrsdip) THEN
             IF (prop1%ldip) THEN
                IF (paral%parent .AND. .NOT.isos1%tclust) THEN
                   IF (paral%io_parent)&
                        WRITE(6,*)&
                        'WARNING: DIPOLE MOMENT IN REAL SPACE WITH PBC'
                ENDIF
                CALL rsdipo(tau0,eirop,psi,rhoe)
             ENDIF
          ELSE
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IR)
             DO ir=1,fpar%nnr1
                psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
             ENDDO
             CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
             IF (prop1%ldip) CALL dipo(tau0,eirop,psi)
             IF (prop1%locd) CALL lodipo(eirop,psi)
          ENDIF
       ENDIF
       IF (paral%parent) CALL prdip
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Calculates conductivity                                      ==
    ! ==--------------------------------------------------------------==
    IF (condpa%tconduct) THEN
       CALL conductivity(c0,eigv,crge%f,crge%n,nkpt%nkpnt)
       IF (paral%io_parent) THEN    ! cmb
          WRITE(6,'(A)')&
               ' CONDUCTIVITIES(OHM^-1.CM^-1)  VERSUS FREQUENCY(EV):'
          DO i=1,condpa%nconduct
             WRITE(6,'(1X,F10.3,5X,2(1PE18.5),5X,I8)')&
                  (REAL(i-1,kind=real_8)+0.5_real_8)*condpa%condstep,conduct(i),&
                  conduct2(i)-conduct(i)**2,normcon(i)
          ENDDO
          WRITE(6,'(A,E12.6)') ' ChkSum(CONDUCTIVITY) = ',&
               SUM(ABS(conduct(1:condpa%nconduct)))
       ENDIF! cmb
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Calculates polarisability                                    ==
    ! ==--------------------------------------------------------------==
    IF (tpolarb) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL phfac(tau0)
       ! Compute electronic density in real space
       IF (tkpts%tkpnt) THEN
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
       ELSE
          CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       ENDIF
       ! Calculate Potential
       CALL vofrho(tau0,fion,rhoe,psi,.FALSE.,.FALSE.)
       ! Calculate polarisability
       CALL polarproj(c0,eigv,crge%f,ener_com%amu,rhoe,psi,crge%n,nkpt%nkpnt)

       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Calculates core optical spectra                              ==
    ! ==--------------------------------------------------------------==
    IF (coresl%tcores) THEN
       CALL phfac(tau0)
       CALL core_spectra(tau0,c0,eigv,crge%f,crge%n,nkpt%nkpnt)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  LOCALIZE WAVEFUNCTION                                       ==
    ! ==--------------------------------------------------------------==
    IF (prop1%locl) THEN
       wannl%twann=.TRUE.
       lenbk=nxxfun(crge%n)
       nxx=MAX(2*lenbk*parai%nproc,2*nkpt%ngwk*crge%n)
       ALLOCATE(sc0(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL phfac(tau0)
       CALL localize(tau0,c0,cs,sc0,crge%n)
       CALL write_irec(irec)
       irec(irec_pwv) = isetone(.FALSE.)
       DEALLOCATE(sc0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  PROJECT WAVEFUNCTIONS ON AO                                 ==
    ! ==--------------------------------------------------------------==
    IF (prop1%pwfn) THEN
       ! Orthogonalize
       CALL phfac(tau0)
       IF (pslo_com%tivan) THEN
          CALL rnlsm(c0(:,:,1),prop2%numorb,1,1,.FALSE.)
       ENDIF
       CALL ortho(prop2%numorb,c0(:,:,1),cs)
       ! Project wfn
       IF (tkpts%tkpnt) CALL stopgm('PROWFN','K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (.NOT.cntl%tdmal) THEN
          CALL prowfn(c0,tau0,prop2%numorb)
       ELSE
          CALL dist_prowfn(c0,tau0,prop2%numorb)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  PROJECT WAVEFUNCTIONS ON  SPHERICAL HARMONICS               ==
    ! ==--------------------------------------------------------------==
    IF (prop1%pylm) THEN
       lr=1
       CALL give_scr_ortho(lo,tag,prop2%numorb)
       lscr=MAX(lo,lr)
       CALL ortho(prop2%numorb,c0(:,:,1),cs)
       CALL proylm(c0,prop7%centylm,prop7%numylm,prop2%numorb,prop7%rylmax,prop7%nylmax)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE LAYER PROJECTED DENSITY OF STATES                 ==
    ! ==--------------------------------------------------------------==
    IF (cldos%tldos) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL phfac(tau0)
       CALL ldos(c0,psi,eigv,crge%n)
    ENDIF
    IF (prop1%tavgp) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(vtemp(ncpw%nhg), stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL give_qphi(il_qphi)
       ALLOCATE(qphi(il_qphi), stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Compute electronic density in real space
       CALL phfac(tau0)
       IF (pslo_com%tivan) THEN
          CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
       ENDIF
       IF (tkpts%tkpnt) THEN
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
       ELSE
          CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       ENDIF
       ALLOCATE(achrg(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ichrg(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(achrg)!,ions1%nat+1)
       CALL zeroing(ichrg)!,ions1%nat+1)
       CALL eicalc(eivps,eirop)
       DO i=1,fpar%nnr1
          psi(i,1)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
       DO ig=1,ncpw%nhg
          vtemp(ig) = psi(nzh(ig),1)
       ENDDO
       CALL elstpo(vtemp,eirop,eivps,psi(:,1),qphi)
       CALL dcopy(2*ncpw%nhg,psi,1,vtemp,1)
       IF (rout1%rhoout) THEN
          filen='ELECT_POT'
          CALL densto(vtemp,tau0,filen)
       ENDIF
       CALL zeroing(psi)!,maxfft)
       DO ig=1,ncpw%nhg
          psi(indz(ig),1) = CONJG(vtemp(ig))
          psi(nzh(ig),1)  = vtemp(ig)
       ENDDO
       IF (geq0.AND.isos1%tclust) THEN
          psi(nzh(1),1) = vtemp(1)
       ELSEIF (geq0) THEN
          psi(nzh(1),1) = CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDIF
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
       DO i=1,fpar%nnr1
          rhoe(i,1)=REAL(psi(i,1))
       ENDDO
       CALL potmed(rhoe,rcut,tau0,achrg,ichrg)
       IF (paral%parent) THEN
          CALL printpot(tau0,achrg)
       ENDIF
       IF (paral%parent) CALL prmem('    PROPPT')
       DEALLOCATE(vtemp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(qphi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       DEALLOCATE(achrg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ichrg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF                     ! if (TAVGP)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE ATOMIC CHARGES                                    ==
    ! ==--------------------------------------------------------------==
    IF (prop1%espc) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vtemp(ncpw%nhg), stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL give_qphi(il_qphi)
       ALLOCATE(qphi(il_qphi), stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(reirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Compute electronic density in real space
       CALL phfac(tau0)
       IF (pslo_com%tivan) THEN
          CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
       ENDIF
       IF (tkpts%tkpnt) THEN
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
          CALL stopgm('PROPPT','k-points not implemented',& 
               __LINE__,__FILE__)
       ELSE
          CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       ENDIF
       ALLOCATE(achrg(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(echrg(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL eicalc(eivps,eirop)

       ALLOCATE(rpsi(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL atomc(rhoe(:,1),rpsi,tau0,achrg)

       DEALLOCATE(rpsi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

#ifdef __SR8000
       !poption parallel
#endif
       DO i=1,fpar%nnr1
          psi(i,1)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          vtemp(ig) = psi(nzh(ig),1)
       ENDDO
       CALL elstpo(vtemp,eirop,eivps,psi(:,1),qphi)
       CALL dcopy(2*ncpw%nhg,psi,1,vtemp,1)
       CALL zeroing(psi)!,maxfft)
       !CDIR NODEP
       DO ig=1,ncpw%nhg
          psi(indz(ig),1) = CONJG(vtemp(ig))
          psi(nzh(ig),1)  = vtemp(ig)
       ENDDO
       IF (geq0.AND.isos1%tclust) THEN
          psi(nzh(1),1) = vtemp(1)
       ELSEIF (geq0) THEN
          psi(nzh(1),1) = CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDIF
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(I)
       DO i=1,fpar%nnr1
          rhoe(i,1)=REAL(psi(i,1))
       ENDDO
       mlen=fpar%nnr1
       ALLOCATE(isel(mlen),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ippc=0
       CALL selectp(isel,tau0,ippc)
       mlen=(ions1%nat+1)*ippc
       ALLOCATE(efield(mlen),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,ippc
          efield(i)=rhoe(isel(i),1)
       ENDDO
       CALL atfield(efield,psi(:,1),vtemp,reirop,qphi,isel,ippc)
       CALL espsolv(efield,echrg,ippc)
       IF (paral%parent) CALL prmem('    PROPPT')
       IF (paral%parent) CALL printcg(tau0,achrg,echrg)
       DEALLOCATE(vtemp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(qphi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       DEALLOCATE(achrg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(echrg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(isel,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(efield,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(reirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  LOCALIZE WAVEFUNCTION                                       ==
    ! ==--------------------------------------------------------------==
    IF (prop1%glocl) THEN
       glocal%tgloc=.TRUE.
       CALL phfac(tau0)

       IF (glocal%tg_real) THEN
          DEALLOCATE(cs,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          nc = 2*2*ncpw%ngw*nstate
          lenbk=nxxfun(prop2%numorb)
          nc=MAX(2*2*lenbk*parai%nproc,nc)
          ALLOCATE(cs(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL g_loc_dr_real(c0,cs,nstate,eigv,tau0,taup,iteropt%nfi)
       ELSE
          nc = 2*2*ncpw%ngw*nstate
          lenbk=nxxfun(prop2%numorb)
          ! write(6,*) MEPOS,LENBK,2*2*LENBK*NPROC,NC
          nc=MAX(2*2*lenbk*parai%nproc,nc)
          DEALLOCATE(cs,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cs(nkpt%ngwk,nc/nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! STOP 'NC'
          CALL g_loc_dr_comp(c0,cs,nstate,eigv,tau0,taup,iteropt%nfi)
       ENDIF

    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  PLOT WAVEFUNCTION                                           ==
    ! ==--------------------------------------------------------------==
    IF (prop4%tcubefile_dens.OR.prop4%tcubefile_orb.OR.prop4%tcubefile_pot) THEN
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       CALL give_scr_rhoofr(lrho,tag)
       CALL give_scr_forces(lscr,tag,nstate,.TRUE.,.FALSE.)
       lscr=MAX(il_rhoe_1d*il_rhoe_2d,2*maxfft,lrho,lscr)
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   scr)!,lscr)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   rhoe)!,il_rhoe_1d*il_rhoe_2d)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(psi)!,SIZE(psi))
       ALLOCATE(z_11(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   z_11)!,nstate*nstate)
       ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ALLOCATE(eirop(nkpt%nhgk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eivps(nkpt%nhgk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! malloc complete.
       CALL phfac(tau0)
       CALL eicalc(eivps,eirop)
       ! density?
       IF (group%nogrp .GT. 1) THEN
          CALL stopgm ('PROPPT/CUBEF',&
               'Please program the orbital split by yourself',& 
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (prop4%tcubecenter) THEN
          ! use the center from the input 
          CALL dcopy(3,prop4%cubecin,1,cubecenter,1)
       ELSE
          ! calculate the 'geometric' center of the system
          ! FIXME: check if cntl%proper PBC correction of the ions is done
          ! for this kind of center determination.
          cubecenter(1) = 0._real_8
          cubecenter(2) = 0._real_8
          cubecenter(3) = 0._real_8
          nnat = 0
          DO isp=1,ions1%nsp
             DO iat=1,ions0%na(isp)
                cubecenter(1) = cubecenter(1) + tau0(1,iat,isp)
                cubecenter(2) = cubecenter(2) + tau0(2,iat,isp)
                cubecenter(3) = cubecenter(3) + tau0(3,iat,isp)
                nnat = nnat + 1
             ENDDO
          ENDDO
          cubecenter(1) =  cubecenter(1) / REAL(nnat,kind=real_8)
          cubecenter(2) =  cubecenter(2) / REAL(nnat,kind=real_8)
          cubecenter(3) =  cubecenter(3) / REAL(nnat,kind=real_8)
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,A,3F10.5,A)') ' PROPPT| CUBE CENTER IS AT: ',&
            (cubecenter(i),i=1,3),' BOHR'
       ! ==--------------------------------------------------------------==
       ! TOTAL DENSITY AND SPIN DENSITY (IF AVAILABLE)
       IF (prop4%tcubefile_dens) THEN
          IF (tkpts%tkpnt) THEN
             ! AK: FIXME: CHECK THIS, COMPARE TO OUTPUT FROM rhopri.F
             ! AK: seems like the calls to RNLSM are missing here.
             IF (paral%io_parent)&
                  WRITE(0,*)&
                  ' PROPPT| WARNING: CUBEFILES NOT TESTED WITH k-points.'
             CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
          ELSE
             CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
          ENDIF
          cubefilename = 'RHO_TOT.cube'
          CALL xstring(cubefilename,ia,ie)
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' PROPPT| WRITING DENSITY TO: ',&
               cubefilename(ia:ie)
          CALL cubefile(cubefilename,rhoe,cubecenter,psi,prop4%thalfmesh)
          IF (cntl%tlsd) THEN
             CALL daxpy(fpar%nnr1,-2._real_8,rhoe(1,2),1,rhoe(1,1),1)
             cubefilename = 'RHO_SPIN.cube'
             CALL xstring(cubefilename,ia,ie)
             IF (paral%io_parent)&
                  WRITE(6,'(A,A)') ' PROPPT| WRITING SPIN DENSITY TO: ',&
                  cubefilename(ia:ie)
             CALL cubefile(cubefilename,rhoe,&
                  cubecenter,psi,prop4%thalfmesh)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! orbitals?
       IF (prop4%tcubefile_orb) THEN
          CALL lag_mult(c0,cs,psi,rhoe,z_11,nstate)
          n_cubefiles = 0
          IF (cntl%tlsd) THEN
             DO i=1,prop4%num_cubeorb
                j = spin_mod%nsup - icubeorb(i) + 1
                IF (j .GE. 1) THEN
                   CALL ffttor(c0(1,j,1),scr(1), psi,ncpw%ngw,.TRUE.)
                   IF (paral%io_parent) WRITE (cubefilename,'(A6,I0,A5)')&
                        'PSI_A.',icubeorb(i),'.cube'
                   CALL cubefile(cubefilename,scr,cubecenter,psi,prop4%thalfmesh)
                   n_cubefiles = n_cubefiles + 1
                ELSE
                   IF (paral%io_parent) WRITE(6,*)'wrong alpha state index'
                   CALL stopgm(procedureN,'wrong alpha state index',&
                        __LINE__,__FILE__)
                ENDIF
                j = spin_mod%nsdown - icubeorb(i) + 1
                IF (j .GE. 1) THEN
                   CALL ffttor(c0(1,j+spin_mod%nsup,1),scr(1), psi,ncpw%ngw,.TRUE.)
                   IF (paral%io_parent) WRITE (cubefilename,'(A6,I0,A5)')&
                        'PSI_B.',icubeorb(i),'.cube'
                   CALL cubefile(cubefilename,scr,cubecenter,psi,prop4%thalfmesh)
                   n_cubefiles = n_cubefiles + 1
                ELSE
                   IF (paral%io_parent) WRITE(6,*)'wrong beta state index'
                   CALL stopgm(procedureN,'wrong beta state index',&
                        __LINE__,__FILE__)
                ENDIF
             ENDDO
          ELSE! LDA
             DO i=1,prop4%num_cubeorb
                j = nstate - icubeorb(i) + 1
                IF (j .LT. 1) THEN
                   IF (paral%io_parent) WRITE(6,*)'wrong state index'
                   CALL stopgm(procedureN,'wrong state index',&
                        __LINE__,__FILE__)
                ELSE
                   CALL ffttor(c0(1,j,1),scr,psi,ncpw%ngw,.TRUE.)
                   IF (paral%io_parent)&
                        WRITE (cubefilename,'(A4,I0,A5)')&
                        'PSI.',icubeorb(i),'.cube'
                   CALL cubefile(cubefilename,scr,cubecenter,psi,prop4%thalfmesh)
                   n_cubefiles = n_cubefiles + 1
                ENDIF
             ENDDO
          ENDIF
          IF (paral%io_parent) WRITE(6,*) ' NUMBER OF CUBEFILES PRINTED = ', n_cubefiles
       ENDIF
       IF (prop4%tcubefile_pot) THEN
          CALL rscpot(c0(:,:,1),tau0,fion,rhoe,psi,&
               .FALSE.,.FALSE.,crge%n,nkpt%nkpnt)
          IF (paral%io_parent)&
               WRITE (cubefilename,'(A)')  'V_Hxc.cube'
          CALL cubefile(cubefilename,rhoe,cubecenter,psi,prop4%thalfmesh)
       ENDIF
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(z_11,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.FALSE.,.FALSE.)
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    RETURN
  END SUBROUTINE proppt
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE propcal(c0,c2,cm,sc0,taup,We,f,amu,&
       vpot,psi,nstate,nkpoint,nfi,infi)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE CALCULATES PROPERTIES DURING THE SIMULATION     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:,:), c2(:,:)
    COMPLEX(real_8), TARGET                  :: cm(:)
    COMPLEX(real_8)                          :: sc0(*)
    REAL(real_8)                             :: taup(:,:,:), amu, &
                                                vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: We(nstate,nkpoint), &
                                                f(nstate,nkpoint)
    INTEGER                                  :: nfi, infi

    CHARACTER(*), PARAMETER                  :: procedureN = 'propcal'

    COMPLEX(real_8), DIMENSION(:, :), &
      POINTER                                :: cm_p
    INTEGER                                  :: i, ia, iat, ierr, is, j
    LOGICAL                                  :: ferror
    REAL(real_8), ALLOCATABLE                :: center(:,:), save_fion(:,:,:)

!(nkpt%ngwk,nstate,nkpoint)
! Variables
! ==--------------------------------------------------------------==
! == DIPOLE MOMENT                                                ==
! ==--------------------------------------------------------------==

    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF ((cntl%tdipd.AND.(MOD(nfi-1,cnti%npdip).EQ.0)).OR.&
         (vdwl%vdwd.AND.vdwwfl%twannup)) THEN
       cm_p(1:nstate,1:nstate)=>cm
       CALL ddipo(taup,c0(:,:,1),cm_p,c2,sc0,nstate,center)
       IF (wannl%twann) THEN
          IF (wannl%twmol.OR.wannl%twdos) THEN
             ALLOCATE(save_fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL forcedr(c0(:,:,1),c2,sc0,vpot,psi,taup,save_fion,We,&
                  nstate,1,.FALSE.,.TRUE.)
             IF (wannl%twmol) THEN
                CALL molorb(c0,c2,taup,nstate,center)
                CALL forcedr(c0(:,:,1),c2,sc0,vpot,psi,taup,save_fion,We,&
                     nstate,1,.FALSE.,.TRUE.)
             ENDIF
             IF (wannl%twdos) CALL wc_dos(c0,c2,nstate,center)
             DEALLOCATE(save_fion,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL wannier_print(nfi,c0(:,:,1),taup,nstate,psi(:,1),center)
       ENDIF
    ENDIF
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == Conductivity                                                 ==
    ! ==--------------------------------------------------------------==
    IF (condpa%tconduct.AND.&
         ((MOD(nfi,condpa%iconduct).EQ.0).OR.(infi.EQ.cnti%nomore))) THEN
       IF (infi.EQ.1) THEN
          condpa%nfconduct=0
          condpa%filecond='CONDUCTIVITIES'
          IF (paral%io_parent)&
               CALL fileopen(iucond,condpa%filecond,fo_new,ferror)
          IF (ferror) THEN
             IF (paral%io_parent)&
                  WRITE(iucond,*) condpa%nconduct
          ELSE
             IF (paral%io_parent)&
                  CALL fileopen(iucond,condpa%filecond,fo_app+fo_verb+fo_mark,&
                  ferror)
          ENDIF
       ENDIF
       condpa%nfconduct=condpa%nfconduct+1
       CALL conductivity(c0,We,f,nstate,nkpoint)
       IF (paral%io_parent)&
            WRITE(iucond,*) 'INFI=',infi,' NFI=',nfi
       DO i=1,condpa%nconduct
          IF (paral%io_parent)&
               WRITE(iucond,100) (REAL(i-1,kind=real_8)+0.5_real_8)*condpa%condstep,conduct(i),&
               conduct2(i)-conduct(i)**2,normcon(i)
100       FORMAT(f7.3,3x,2(1pe14.5),3x,i5)
       ENDDO
       IF ((infi.EQ.cnti%nomore).AND.paral%io_parent)&
            CALL fileclose(iucond)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Polarisability                                               ==
    ! ==--------------------------------------------------------------==
    IF (tpolarb.AND.&
         ((MOD(nfi,ipolarise).EQ.0).OR.(infi.EQ.cnti%nomore))) THEN
       IF (infi.EQ.1) THEN
          nfpolar=0
          condpa%filecond='POLARISABILITY'
          IF (paral%io_parent)&
               CALL fileopen(iupola,condpa%filecond,fo_app+fo_verb+fo_mark,&
               ferror)
       ENDIF
       nfpolar=nfpolar+1
       CALL polarproj(c0,We,f,amu,&
            vpot,psi,nstate,nkpoint)
       IF (paral%io_parent)WRITE(iupola,120) infi,((alphap(i,j),j=1,3),i=1,3)
       IF (tzeff) THEN
          IF (paral%io_parent)&
               WRITE(iupola,'(A)') 'EFFECTIVE CHARGE:'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                IF (paral%io_parent)&
                     WRITE(iupola,'(A2,I3/3(3(F11.4,1X)/))')&
                     elem%el(ions0%iatyp(is)),iat,&
                     (zeff(1,j,ia,is),j=1,3),&
                     (zeff(2,j,ia,is),j=1,3),&
                     (zeff(3,j,ia,is),j=1,3)
             ENDDO
          ENDDO
       ENDIF
       IF ((infi.EQ.cnti%nomore).AND.paral%io_parent)&
            CALL fileclose(iupola)
    ENDIF
120 FORMAT(i5/'XX:',f14.5,'   XY:',f14.5,'   XZ:',f14.5/&
         'YX:',f14.5,'   YY:',f14.5,'   YZ:',f14.5/&
         'ZX:',f14.5,'   ZY:',f14.5,'   ZZ:',f14.5)
    RETURN
  END SUBROUTINE propcal
  ! ==================================================================
  SUBROUTINE prdip
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, iex, iexj, ii
    REAL(real_8)                             :: d1, d2, d3, dd

! ==--------------------------------------------------------------==

    IF (prop1%ldip .OR. prop1%dberry) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' DIPOLE MOMENT '
       d1=moment%dmom(1)
       d2=moment%dmom(2)
       d3=moment%dmom(3)
       dd=SQRT(d1*d1+d2*d2+d3*d3)
       IF (paral%io_parent)&
            WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
       IF (paral%io_parent)&
            WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
       IF (paral%io_parent)&
            WRITE(6,'(4F12.5,A)') au_deb*d1, au_deb*d2, au_deb*d3,&
            au_deb*dd,'   Debye'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    IF (prop1%locd) THEN
       DO i=1,numbld
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' LOCAL DIPOLE MOMENT '
          IF (paral%io_parent)&
               WRITE(6,'(A,F8.4,A,F8.4,A)')&
               ' X = ',xminld(i),' - ',xmaxld(i),' a.u.'
          IF (paral%io_parent)&
               WRITE(6,'(A,F8.4,A,F8.4,A)')&
               ' Y = ',yminld(i),' - ',ymaxld(i),' a.u.'
          IF (paral%io_parent)&
               WRITE(6,'(A,F8.4,A,F8.4,A)')&
               ' Z = ',zminld(i),' - ',zmaxld(i),' a.u.'
          IF (paral%io_parent)&
               WRITE(6,'(3(A,F8.4))')&
               ' CENTER OF CHARGE: ',rcc(1,i),' /',&
               rcc(2,i),' /',rcc(3,i)
          d1=dmomlo(1,i)
          d2=dmomlo(2,i)
          d3=dmomlo(3,i)
          dd=SQRT(d1*d1+d2*d2+d3*d3)
          IF (paral%io_parent)&
               WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
          IF (paral%io_parent)&
               WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
          IF (paral%io_parent)&
               WRITE(6,'(4F12.5,A)') au_deb*d1, au_deb*d2, au_deb*d3,&
               au_deb*dd,'   Debye'
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDDO
    ENDIF
    IF (prop1%lext) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A)')&
            ' EXCITED STATE DIPOLE MOMENT DIFFERENCES '
       DO iex=1,extd%nexdip
          iexj=(iex-1)*prop2%numorb
          IF (paral%io_parent)&
               WRITE(6,'(/,A,I4)') ' OCCUPATION OF STATE ',iex
          IF (paral%io_parent)&
               WRITE(6,'(20F4.1)') (focc(ii+iexj),ii=1,prop2%numorb)
          IF (prop1%locd) THEN
             DO i=1,numbld
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' LOCAL DIPOLE MOMENT '
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' X = ',xminld(i),' - ',xmaxld(i),' a.u.'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' Y = ',yminld(i),' - ',ymaxld(i),' a.u.'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' Z = ',zminld(i),' - ',zmaxld(i),' a.u.'
                d1=exd(1,i,iex)
                d2=exd(2,i,iex)
                d3=exd(3,i,iex)
                dd=SQRT(d1*d1+d2*d2+d3*d3)
                IF (paral%io_parent)&
                     WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
                IF (paral%io_parent)&
                     WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
                IF (paral%io_parent)&
                     WRITE(6,'(4F12.5,A)') au_deb*d1,au_deb*d2,au_deb*d3,&
                     au_deb*dd,'   Debye'
                IF (paral%io_parent)&
                     WRITE(6,*)
             ENDDO
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' DIPOLE MOMENT '
             d1=exd(1,1,iex)
             d2=exd(2,1,iex)
             d3=exd(3,1,iex)
             dd=SQRT(d1*d1+d2*d2+d3*d3)
             IF (paral%io_parent)&
                  WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
             IF (paral%io_parent)&
                  WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
             IF (paral%io_parent)&
                  WRITE(6,'(4F12.5,A)') au_deb*d1,au_deb*d2,au_deb*d3,&
                  au_deb*dd,'   Debye'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF
       ENDDO
    ENDIF
    IF (prop1%ltrm) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') ' TRANSITION MOMENTS'
       DO iex=1,extd%ntrmom
          IF (paral%io_parent)&
               WRITE(6,'(/,A,I4,A)',advance="no")&
               ' EXCITATION ',iex,' BETWEEN STATES: '
          IF (paral%io_parent)&
               WRITE(6,'(2I4)') (nsdip(ii,iex),ii=1,2)
          IF (prop1%locd) THEN
             DO i=1,numbld
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' LOCAL TRANSITION MOMENT '
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' X = ',xminld(i),' - ',xmaxld(i),' a.u.'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' Y = ',yminld(i),' - ',ymaxld(i),' a.u.'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F8.4,A,F8.4,A)')&
                     ' Z = ',zminld(i),' - ',zmaxld(i),' a.u.'
                d1=trmom(1,i,iex)
                d2=trmom(2,i,iex)
                d3=trmom(3,i,iex)
                dd=SQRT(d1*d1+d2*d2+d3*d3)
                IF (paral%io_parent)&
                     WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
                IF (paral%io_parent)&
                     WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
                IF (paral%io_parent)&
                     WRITE(6,'(4F12.5,A)') au_deb*d1,au_deb*d2,au_deb*d3,&
                     au_deb*dd,'   Debye'
                IF (paral%io_parent)&
                     WRITE(6,*)
             ENDDO
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' TRANSITION MOMENT '
             d1=trmom(1,1,iex)
             d2=trmom(2,1,iex)
             d3=trmom(3,1,iex)
             dd=SQRT(d1*d1+d2*d2+d3*d3)
             IF (paral%io_parent)&
                  WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
             IF (paral%io_parent)&
                  WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
             IF (paral%io_parent)&
                  WRITE(6,'(4F12.5,A)') au_deb*d1,au_deb*d2,au_deb*d3,&
                  au_deb*dd,'   Debye'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prdip
  ! ==================================================================
  SUBROUTINE give_scr_propcal(lpropcal,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lpropcal
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lddipo

! ==--------------------------------------------------------------==

    lddipo = 0
    IF (cntl%tdipd) CALL give_scr_ddipo(lddipo,tag)
    IF (tpolarb) THEN
       CALL give_scr_polarproj(lpropcal,tag,nstate)
    ELSE
       lpropcal=0
    ENDIF
    lpropcal=MAX(lpropcal,lddipo)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_propcal
  ! ==================================================================
  SUBROUTINE give_scr_espc(lespc,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lespc
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_qphi, lrnlsm

    CALL give_qphi(il_qphi)
    IF (pslo_com%tivan) THEN
       CALL give_scr_rnlsm(lrnlsm,tag,crge%n,.FALSE.)
    ELSE
       lrnlsm=0
    ENDIF
    lespc=MAX(2*ncpw%nhg+il_qphi,lrnlsm)
    tag='2*NHG+IL_QPHI'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_espc
  ! ==================================================================
  SUBROUTINE give_scr_avgp(lavgp,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lavgp
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_qphi, lrnlsm

    CALL give_qphi(il_qphi)
    IF (pslo_com%tivan) THEN
       CALL give_scr_rnlsm(lrnlsm,tag,crge%n,.FALSE.)
    ELSE
       lrnlsm=0
    ENDIF
    lavgp=MAX(2*ncpw%nhg+il_qphi,lrnlsm)
    tag='2*NHG+IL_QPHI'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_avgp


END MODULE proppt_utils
