MODULE secder_driver
  USE adat,                            ONLY: elem
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE bsym,                            ONLY: bsclcs,&
                                             cnstwgt
  USE bsympnt,                         ONLY: fnbs
  USE calc_alm_utils,                  ONLY: calc_alm
  USE cnst,                            ONLY: au_deb
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE copot_utils,                     ONLY: copot
  USE cotr,                            ONLY: cotc0,&
                                             hess,&
                                             lskcor
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: ddipo
  USE detdof_utils,                    ONLY: detdof
  USE dipo_utils,                      ONLY: rsdipo
  USE dipomod,                         ONLY: moment
  USE dynit_utils,                     ONLY: dynit
  USE ehpsi_utils,                     ONLY: set_b2l
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def
  USE fint,                            ONLY: fint1
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE gsize_utils,                     ONLY: gsize
  USE hessin_utils,                    ONLY: hessin
  USE hessout_utils,                   ONLY: hessout
  USE initrun_driver,                  ONLY: initrun
  USE ions,                            ONLY: coord_fdiff,&
                                             ions0,&
                                             ions1,&
                                             r_fdiff,&
                                             tref_fdiff
  USE k_updwf_utils,                   ONLY: k_updwf
  USE kddipo_utils,                    ONLY: kddipo
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01,&
                                             tshl
  USE lr_tddft_utils,                  ONLY: lr_tddft
  USE lsforce_utils,                   ONLY: lsforce
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mp_interface,                    ONLY: mp_bcast
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE ortho_utils,                     ONLY: ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: infi,&
                                             infw,&
                                             iteropt,&
                                             ropt_mod
  USE rrane_utils,                     ONLY: rrane
  USE secder_utils,                    ONLY: give_scr_secder,&
                                             mldgau,&
                                             molvib,&
                                             purged,&
                                             vibeig,&
                                             writeaclimax
  USE setbsstate_utils,                ONLY: setbsstate
  USE setirec_utils,                   ONLY: write_irec
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: restart1,&
                                             rout1
  USE symm,                            ONLY: isymu,&
                                             symmi
  USE symtrz_utils,                    ONLY: symmat
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, parm
  USE testex_utils,                    ONLY: testex
  USE updrho_utils,                    ONLY: updrho
  USE updwf_utils,                     ONLY: updwf
  USE utils,                           ONLY: nxxfun,&
                                             symma
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: secder

CONTAINS

  ! ==================================================================
  SUBROUTINE secder(c0,c1,c2,sc0,cg0,pme,gde,vpp,eigv)
    ! ==--------------------------------------------------------------==
    ! ==  VIBRATIONAL ANALYSIS                                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(:,:,:), c1(:,:,:), c2(:,:,:), &
      sc0(nkpt%ngwk,crge%n,*), cg0(nkpt%ngwk,crge%n,*), pme(*), gde(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), eigv(crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'secder'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/'X','Y','Z'/)

    CHARACTER(len=100)                       :: filen, filen1
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cm(:,:), eirop(:), psi(:,:)
    INTEGER :: i, ia, ia2, iat, iat2, icol, id, idelt, ierr, il_eirop, &
      il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, info, irec(100), is, is2, &
      isymuat, ix, j, k, k2, lscr, mode, naux, nhpsi, nmm, nnx, nsym
    LOGICAL                                  :: cvbswf, cvhswf, ferror, &
                                                readhe, tinfo
    REAL(real_8) :: d1, d2, d3, dd, detot, dip(3,ions1%nat,2), dist, ekin1, &
      ekin2, ekincp, ekinh1, ekinh2, etotbs, etotg, etoths, etoto, &
      ir(3,3*ions1%nat), tcpu, temp1, temp2, thl(2), time1, time2, xlm, &
      xmass, ylm, zlm
    REAL(real_8), ALLOCATABLE :: dipder(:,:), dummy(:,:), fdis(:,:,:,:), &
      rhoe(:,:), rinb(:,:), scr(:), sder(:,:), tscr(:,:,:), vibe(:), xma(:)

    time1 =m_walltime()
    tinfo=.FALSE.
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dipder(3,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dipder)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    IF (paral%parent) THEN
       ALLOCATE(xma(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vibe(3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fdis(3,maxsys%nax,maxsys%nsx,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sder(3*ions1%nat,3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xma)
       CALL zeroing(vibe)
       CALL zeroing(fdis)
       CALL zeroing(sder)
    ENDIF
    IF (cntl%tdiag) THEN
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rinb(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rout0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       rhoo => rin0
    ELSEIF (tkpts%tkpnt) THEN
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%bsymm) THEN
       ALLOCATE(fnbs(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    nacc = 7
    ener_com%ecnstr = 0.0_real_8
    ener_com%erestr = 0.0_real_8
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d, &
         il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmm=1
    IF (cntl%tddft) THEN
       ! TODO refactor RHOO 
       ALLOCATE(rhoo(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(potr(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td01%ns_tri.GT.0.OR.tshl%isc) nmm=2
    ENDIF
    il_psi_1d=il_psi_1d*nmm
    ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_secder(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! EXTRA MEMORY REQUIRED FOR DIPOLE MOMENTS
    IF (prop1%dberry.AND.tkpts%tkpnt) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,A,23X,A,23X,A,/,1X,A,/,1X,A,/)')&
            '******','WARNING','******',&
            'DIPOLE MOMENTS WITH K-POINTS: '//&
            'PLEASE FINISH WRITING KDDIPO',&
            'DISABLING DIPOLE MOMENTS FOR NOW'
       prop1%dberry=.FALSE.
       prop1%ldip=.FALSE.
    ENDIF
    IF (prop1%dberry) THEN
       lenbk=nxxfun(prop2%numorb)
       ALLOCATE(cm(nkpt%ngwk,crge%n*nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE IF (prop1%ldip) THEN
       il_eirop=ncpw%nhg*2
       ALLOCATE(eirop(il_eirop),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! ==--------------------------------------------------------------==
    bsclcs=1
    IF (cntl%bsymm)CALL setbsstate
    IF (cntl%bsymm.AND.paral%io_parent) THEN
       WRITE(6,*)
       WRITE(6,*) 'BSYMM: BS WAVEFUNCTION INITIALIZATION'
    ENDIF
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    IF (cntl%bsymm) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'BSYMM: HS WAVEFUNCTION INITIALIZATION'
       bsclcs=2
       CALL setbsstate
       CALL initrun(irec,c0(:,:,2:),c2(:,:,2),sc0(1,1,2),rhoe,psi,&
            eigv(1,2))
    ENDIF
    CALL write_irec(irec)
    ! END INITIALIZATION
    ! ..  Do not symmetrize density...  
    cntl%tsymrho=.FALSE.
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull).AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)&
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'//&
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' WARNING! USE KPOINTS FULL  !!!'//&
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'//&
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDIF
    ! ..
    etoto=0._real_8
    ener_com%ebogo=0._real_8
    IF (paral%parent) CALL detdof(tau0,tscr)
    cnti%npara=-1
    readhe=.FALSE.
    IF (paral%parent) CALL hessin(tau0,readhe)
    etotg=0.0_real_8
    CALL write_irec(irec)
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==      OPTIMIZATION OF WAVEFUNCTION AT MINIMUM STRUCTURE       ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/," ",22("*"),A,23("*"),/)') ' MINIMUM STRUCTURE '
    ropt_mod%convwf=.FALSE.
    cvbswf=.FALSE.
    cvhswf=.FALSE.
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    time1=m_walltime()

    DO infw=1,cnti%nomore
       IF (cntl%tdiag) THEN
          IF (cntl%bsymm)THEN
             CALL stopgm('SECDER','CNTL%BSYMM .AND. TDIAG NOT IMPLEMENTED',& 
                  __LINE__,__FILE__)
          ELSE
             CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                  rhoe,psi,&
                  crge%n,.FALSE.,.TRUE.,.FALSE.,infw,thl,nhpsi)
          ENDIF
          ! ..store initial density
          CALL dcopy(nnx,rin0,1,rinb,1)
       ELSE
          IF (.NOT.cntl%bsymm)THEN
             ! UPDATE THE WAVEFUNCTIONS
             IF (tkpts%tkpnt) THEN
                CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                     rhoe,psi,crge%n,.FALSE.,infi)
             ELSE
                CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,tscr,pme,gde,vpp,eigv,&
                     rhoe,psi,crge%n,.FALSE.,.TRUE.)
             ENDIF
          ELSE
             ! Update HS WF for Brokensymmetry
             IF (.NOT.cvhswf) THEN
                IF (infw.EQ.1) THEN
                   bsclcs=2
                   CALL setbsstate
                ENDIF
                CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,&
                     fion,pme,gde,vpp,eigv(1,2),rhoe,&
                     psi,crge%n,.FALSE.,.TRUE.)
                cvhswf=ropt_mod%convwf
                ! Update BS WF for Brokensymmetry
             ELSEIF (.NOT.cvbswf) THEN
                IF (bsclcs.EQ.2) THEN
                   IF (paral%io_parent)&
                        WRITE(6,*) 'BSYMM: HS STATE CONVERGED. ',&
                        'SWITCHING TO BS STATE.'
                   etoto=0.0_real_8
                   bsclcs=1
                   CALL setbsstate
                   CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                        vpp,eigv,rhoe,psi,crge%n,&
                        .FALSE.,.TRUE.)
                   IF (tinfo.AND.paral%parent) THEN
                      CALL wrener
                      IF (paral%io_parent)&
                           WRITE(6,'(2A)')&
                           ' NFI      GEMAX       CNORM',&
                           '           ETOT        DETOT      TCPU'
                   ENDIF
                ELSE
                   CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                        vpp,eigv,rhoe,psi,crge%n,&
                        .FALSE.,.TRUE.)
                ENDIF
                cvbswf=ropt_mod%convwf
             ENDIF
             ropt_mod%convwf=cvhswf.AND.cvbswf
          ENDIF
       ENDIF
       ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%parent) THEN
          detot=ener_com%etot-etoto
          IF (infw.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv(1,bsclcs),crge%f,ener_com%amu,crge%n,ener_com%etot,etoto,&
               tcpu,gemax,cnorm,thl,ropt_mod%engpri,infw,infw)
          etoto=ener_com%etot
       ENDIF
       CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf.OR.soft_com%exsoft) GOTO 100
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' NO CONVERGENCE IN ',infw,' STEPS '
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! NUCLEAR GRADIENT
    IF (cntl%tdiag) THEN
       CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
            rhoe,psi,&
            crge%n,.TRUE.,.TRUE.,cntl%tpres,infi,thl,nhpsi)
    ELSE
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       IF (tkpts%tkpnt) THEN
          CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,&
               eigv,rhoe,psi,crge%n,.TRUE.,infi)
       ELSE
          CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,tau0,fion,eigv,&
               crge%n,1,.TRUE.,.TRUE.)
       ENDIF
       etotbs=ener_com%etot
       IF (cntl%bsymm)THEN
          IF (paral%parent)CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
          bsclcs=2
          CALL setbsstate
          CALL forcedr(c0(:,:,2),c2(:,:,2),sc0(:,:,2),rhoe,psi,tau0,&
               fion,eigv,crge%n,1,.TRUE.,.TRUE.)
          etoths=ener_com%etot
          ener_com%etot = (1.0_real_8+cnstwgt)*etotbs-cnstwgt*etoths
          IF (paral%parent)CALL lsforce(fnbs,fion)
          CALL mp_bcast(fion,SIZE(fion),parai%io_source,parai%cp_grp)
       ENDIF
    ENDIF
    CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)
    time2=m_walltime()
    tcpu=(time2-time1)*0.001_real_8
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.3,A)') ' TIME FOR MINIMUM STRUCTURE :',&
            tcpu,' SECONDS'
       CALL wrener
       CALL wrgeof(tau0,fion)
    ENDIF
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    CALL dcopy(2*nkpt%ngwk*crge%n*nkpt%nkpnt*bsclcs,c0,1,cg0,1)
    CALL dcopy(3*maxsys%nax*ions1%nsp,tau0(1,1,1),1,taup(1,1,1),1)
    ! GET DIPOLE MOMENT FOR MINIMUM STRUCTURE
    IF (prop1%dberry) THEN
       IF (tkpts%tkpnt) THEN
          CALL kddipo(tau0,c0,cm,c2,sc0,crge%n)
       ELSE
          ALLOCATE(dummy(4,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL ddipo(tau0,c0(:,:,1),cm,c2(:,:,1),sc0,crge%n,dummy)
          DEALLOCATE(dummy,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ELSE IF (prop1%ldip) THEN
       CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
       IF (tkpts%tkpnt) THEN
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
       ELSE
          CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       ENDIF
       ! CALL EICALC(SCR,EIROP)
       CALL rsdipo(tau0,eirop,psi,rhoe)
    ENDIF
    IF (paral%parent.AND.(prop1%ldip.OR.prop1%dberry)) THEN
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
            WRITE(6,'(4F12.5,A)')au_deb*d1,au_deb*d2,au_deb*d3,&
            au_deb*dd,'  Debye'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    ! Open the FINDIF and FDDIP file
    filen='FINDIF'
    filen1='FDDIP'
    IF (paral%io_parent) THEN
       CALL fileopen(24,filen,fo_def,ferror)
       IF (restart1%rvib) REWIND(24)
       IF (prop1%ldip.OR.prop1%dberry) THEN
          CALL fileopen(25,filen1,fo_def,ferror)
          IF (restart1%rvib.AND.paral%io_parent) REWIND(25)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==         THE BASIC LOOP OVER PERTURBATIONS (COORDINATES)      ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",22("*"),A,22("*"),/)') ' FINITE DIFFERENCES '
       IF (tref_fdiff) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,T30,3(F12.3))') 'REFERENCE:',coord_fdiff
          IF (paral%io_parent)&
               WRITE(6,'(A,T54,F12.3)') 'RADIUS OF THE SPHERE:',r_fdiff
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               'WE COMPUTE ONLY THE FINITE DIFFERENCES FOR ATOMS'
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               'INSIDE THIS SPHERE.'
       ENDIF
    ENDIF
    ! No symmetrization of gradients in this part!
    nsym=symmi%nrot
    symmi%nrot=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          IF (tref_fdiff) THEN
             xlm=tau0(1,ia,is)-coord_fdiff(1)
             ylm=tau0(2,ia,is)-coord_fdiff(2)
             zlm=tau0(3,ia,is)-coord_fdiff(3)
             CALL pbc(xlm,ylm,zlm,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             dist = SQRT(xlm*xlm+ylm*ylm+zlm*zlm)
             IF (dist.GT.r_fdiff) THEN
                GOTO 31416
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A,I7,4X,A4,4X,A,T54,1PE12.5)')&
                     " **** ATOM=",iat,elem%el(ions0%iatyp(is)),"DISTANCE=",dist
             ENDIF
          ENDIF
          DO k=1,3
             DO idelt=-1,1,2
                time1=m_walltime()
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(" **** ATOM=",I7,4X,A4,4X,A4,4X,A,T54,1PE12.5)')&
                     iat,elem%el(ions0%iatyp(is)),cdir(k),&
                     'DISPLACEMENT=',cntr%fdiff*idelt
                ! NEW GEOMETRY            
                IF (restart1%rvib .AND. paral%parent) THEN
                   IF (paral%io_parent)&
                        READ(24,fmt=*,END=300)
                   DO is2=1,ions1%nsp
                      DO ia2=1,ions0%na(is2)
                         IF (paral%io_parent)&
                              READ(24,*) (fion(k2,ia2,is2),k2=1,3)
                      ENDDO
                   ENDDO
                   id=(idelt+1)/2 + 1
                   IF (prop1%ldip.OR.prop1%dberry) READ(25,*) (dip(k2,iat,id),k2=1,3)
                   GOTO 301
300                CONTINUE
                   restart1%rvib=.FALSE.
                   IF (paral%io_parent)&
                        CALL fileclose(24)
                   IF (paral%io_parent)&
                        CALL fileopen(24,filen,fo_app,ferror)
                   IF (prop1%ldip.OR.prop1%dberry) THEN
                      IF (paral%io_parent)&
                           CALL fileclose(25)
                      IF (paral%io_parent)&
                           CALL fileopen(25,filen,fo_app,ferror)
                   ENDIF
                ENDIF
301             CONTINUE
                CALL mp_bcast(restart1%rvib,parai%io_source,parai%cp_grp)
                IF (.NOT.restart1%rvib) THEN
                   isymuat=0
                   IF (symmi%indpg.NE.0) THEN
                      IF (isymu(iat).EQ.0) THEN
                         CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
                         isymuat=1
                      ENDIF
                   ENDIF
                   IF (isymuat.EQ.0.AND.lskcor(k,iat).NE.0) THEN
                      CALL dcopy(3*maxsys%nax*ions1%nsp,taup(1,1,1),1,tau0(1,1,1),1)
                      tau0(k,ia,is)=tau0(k,ia,is)+idelt*cntr%fdiff
                      CALL dcopy(2*nkpt%ngwk*crge%n*nkpt%nkpnt*bsclcs,cg0,1,c0,1)
                      IF (cntl%trane) CALL rrane(c0,c2(:,:,1),crge%n)
                      IF (cntl%bsymm.AND.cntl%trane)CALL rrane(c0(:,:,2:),c2(:,:,2),crge%n)
                      CALL phfac(tau0)
                      IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)
                      IF (.NOT.cntl%bsymm)THEN
                         IF (pslo_com%tivan)CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
                         CALL ortho(crge%n,c0(:,:,1),c2(:,:,1))
                      ELSE
                         bsclcs=1
                         CALL setbsstate
                         IF (pslo_com%tivan)CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
                         CALL ortho(crge%n,c0(:,:,1),c2(:,:,1))
                         bsclcs=2
                         CALL setbsstate
                         IF (pslo_com%tivan)CALL rnlsm(c0(:,:,2),crge%n,&
                              1,1,.FALSE.)
                         CALL ortho(crge%n,c0(:,:,2),c2(:,:,2))
                      ENDIF
                      ropt_mod%convwf=.FALSE.
                      cvhswf=.FALSE.
                      cvbswf=.FALSE.
                      ropt_mod%sdiis=.TRUE.
                      ropt_mod%spcg=.TRUE.
                      ! NL-Projector Overlap Matrix
                      IF (cntl%tdiag.AND.fint1%ttrot) THEN
                         CALL calc_alm
                      ENDIF
                      IF (cntl%tlanc) THEN
                         CALL set_b2l()
                      ENDIF
                      DO infw=1,cnti%nomore
                         IF (cntl%tdiag) THEN
                            CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                                 rhoe,psi,&
                                 crge%n,.FALSE.,tinfo,.FALSE.,&
                                 infw,thl,nhpsi)
                         ELSE
                            ! UPDATE THE WAVEFUNCTIONS
                            IF (.NOT.cntl%bsymm)THEN
                               IF (tkpts%tkpnt) THEN
                                  CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,&
                                       eigv,rhoe,psi,crge%n,.FALSE.,infi)
                               ELSE
                                  CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,tscr,pme,gde,vpp,&
                                       eigv,rhoe,psi,crge%n,.FALSE.,.TRUE.)
                               ENDIF
                            ELSE
                               ! Update HS WF for Brokensymmetry
                               IF (.NOT.cvhswf) THEN
                                  IF (infw.EQ.1) THEN
                                     bsclcs=2
                                     CALL setbsstate
                                  ENDIF
                                  CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),&
                                       tau0,fion,pme,gde,vpp,eigv(1,2),rhoe,&
                                       psi,crge%n,.FALSE.,.TRUE.)
                                  cvhswf=ropt_mod%convwf
                                  ! Update BS WF for Brokensymmetry
                               ELSEIF (.NOT.cvbswf) THEN
                                  IF (bsclcs.EQ.2) THEN
                                     IF (paral%io_parent)&
                                          WRITE(6,*) 'BSYMM: HS STATE CONVERGED. ',&
                                          'SWITCHING TO BS STATE.'
                                     etoto=0.0_real_8
                                     bsclcs=1
                                     CALL setbsstate
                                     CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                          vpp,eigv,rhoe,psi,&
                                          crge%n,.FALSE.,.TRUE.)
                                     IF (tinfo.AND.paral%parent) THEN
                                        CALL wrener
                                        IF (paral%io_parent)&
                                             WRITE(6,'(2A)')&
                                             ' NFI      GEMAX       CNORM',&
                                             '           ETOT        DETOT      TCPU'
                                     ENDIF
                                  ELSE
                                     CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                          vpp,eigv,rhoe,psi,&
                                          crge%n,.FALSE.,.TRUE.)
                                  ENDIF
                                  cvbswf=ropt_mod%convwf
                               ENDIF
                               ropt_mod%convwf=cvhswf.AND.cvbswf
                            ENDIF
                         ENDIF
                         IF (paral%parent) THEN
                            detot=ener_com%etot-etoto
                            IF (infw.EQ.1) detot=0.0_real_8
                            time2=m_walltime()
                            tcpu=(time2-time1)*0.001_real_8
                            IF (cntl%tdiag.AND.tinfo) THEN
                               CALL wrprint_wfopt(eigv(1,bsclcs),crge%f,ener_com%amu,crge%n,ener_com%etot,&
                                    etoto,tcpu,gemax,cnorm,thl,.FALSE.,infw,&
                                    infw)
                            ELSEIF (cntl%tdiag) THEN
                               IF (paral%io_parent)&
                                    WRITE(6,&
                                    '(T9,I4,T15,1PE12.3,T27,0PF12.6,T39,'//&
                                    '1PE12.3,T53,0PF5.2,T59,0PF7.2) ')&
                                    infw,gemax,ener_com%etot,detot,nhpsi/REAL(crge%n,kind=real_8),tcpu
                            ELSE
                               IF (paral%io_parent)&
                                    WRITE(6,&
                                    '(T9,I4,T15,1PE12.3,T27,0PF15.6,T44,'//&
                                    '1PE12.3,T59,0PF7.2) ')&
                                    infw,gemax,ener_com%etot,detot,tcpu
                            ENDIF
                            etoto=ener_com%etot
                         ENDIF
                         CALL testex(soft_com%exsoft)
                         IF (soft_com%exsoft) GOTO 200
                         IF (ropt_mod%convwf) GOTO 101
                      ENDDO
101                   CONTINUE
                      ! NUCLEAR GRADIENT
                      IF (cntl%tdiag) THEN
                         CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
                         CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                              rhoe,psi,&
                              crge%n,.TRUE.,.TRUE.,cntl%tpres,infi,thl,nhpsi)
                         CALL dcopy(fpar%nnr1*clsd%nlsd,rinb,1,rin0,1)
                      ELSE
                         bsclcs=1
                         IF (cntl%bsymm)CALL setbsstate
                         IF (tkpts%tkpnt) THEN
                            CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,&
                                 eigv,rhoe,psi,crge%n,.TRUE.,infi)
                         ELSE
                            CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,tau0,fion,eigv,&
                                 crge%n,1,.TRUE.,.TRUE.)
                         ENDIF
                         IF (cntl%bsymm)THEN
                            IF (paral%parent)CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
                            etotbs=ener_com%etot
                            bsclcs=2
                            CALL setbsstate
                            CALL forcedr(c0(:,:,2),c2(:,:,2),sc0(:,:,2),rhoe,&
                                 psi,tau0,fion,eigv,crge%n,1,&
                                 .TRUE.,.TRUE.)
                            etoths=ener_com%etot
                            ener_com%etot = (1.0_real_8+cnstwgt)*etotbs-cnstwgt*etoths
                            IF (paral%parent)CALL lsforce(fnbs,fion)
                            CALL mp_bcast(fion,SIZE(fion),parai%io_source,parai%cp_grp)
                         ENDIF
                         IF (cntl%tddft)&
                              CALL lr_tddft(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,tau0,fion,&
                              eigv,crge%n,.TRUE.,td01%ioutput)
                      ENDIF
                      CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)
                   ENDIF
                   ! WAVEFUNCTION AND FORCES DONE, NOW DIPOLE MOMENT
                   IF (prop1%dberry) THEN
                      IF (tkpts%tkpnt) THEN
                         CALL kddipo(tau0,c0,cm,c2,sc0,crge%n)
                      ELSE
                         ALLOCATE(dummy(4,1),STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                              __LINE__,__FILE__)
                         CALL ddipo(tau0,c0(:,:,1),cm,c2(:,:,1),sc0,crge%n,dummy)
                         DEALLOCATE(dummy,STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                              __LINE__,__FILE__)
                      ENDIF
                   ELSE IF (prop1%ldip) THEN
                      CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
                      IF (tkpts%tkpnt) THEN
                         CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
                      ELSE
                         CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
                      ENDIF
                      ! CALL EICALC(SCR,EIROP)
                      CALL rsdipo(tau0,eirop,psi,rhoe)
                   ENDIF
                   IF (paral%parent.AND.(prop1%ldip.OR.prop1%dberry)) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A)') ' DIPOLE MOMENT '
                      id=(idelt+1)/2 + 1
                      d1=moment%dmom(1)
                      d2=moment%dmom(2)
                      d3=moment%dmom(3)
                      dip(1,iat,id) = d1
                      dip(2,iat,id) = d2
                      dip(3,iat,id) = d3
                      dd=SQRT(d1*d1+d2*d2+d3*d3)
                      IF (paral%io_parent)&
                           WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
                      IF (paral%io_parent)&
                           WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
                      IF (paral%io_parent)&
                           WRITE(6,'(4F12.5,A)')au_deb*d1,au_deb*d2,au_deb*d3,&
                           au_deb*dd,'  Debye'
                      IF (paral%io_parent)&
                           WRITE(6,*)
                   ENDIF
                   ! WRITE GRADIENT AND DIPOLE MOMENT TO DISK
                   IF (paral%parent) THEN
                      IF (paral%io_parent)&
                           WRITE(24,'(4I10,1PE20.10)') is,ia,k,idelt,cntr%fdiff
                      DO is2=1,ions1%nsp
                         DO ia2=1,ions0%na(is2)
                            IF (paral%io_parent) THEN
                               WRITE(24,*) (fion(k2,ia2,is2),k2=1,3)
                               CALL m_flush(24)
                            ENDIF
                         ENDDO
                      ENDDO
                      IF (prop1%ldip.OR.prop1%dberry) THEN
                         IF (paral%io_parent) THEN
                            WRITE(25,*) (dip(k2,iat,id),k2=1,3)
                            CALL m_flush(25)
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
                IF (paral%parent) THEN
                   CALL gsize(fion,gnmax,gnorm)
                   id=(idelt+1)/2 + 1
                   CALL dcopy(3*maxsys%nax*maxsys%nsx,fion(1,1,1),1,fdis(1,1,1,id),1)
                   time2=m_walltime()
                   tcpu=(time2-time1)*0.001_real_8
                   IF (paral%io_parent)&
                        WRITE(6,'(" ITER=",I6,T22,"ENERGY=",T46,F20.10)')&
                        infw,ener_com%etot
                   IF (paral%io_parent)&
                        WRITE(6,&
                        '(" TCPU=",0PF9.2,T22,"GRADIENT=",T46,1PE20.3,/)')&
                        tcpu,gnorm
                ENDIF
                CALL testex(soft_com%exsoft)
                IF (soft_com%exsoft) GOTO 200
             ENDDO
             ! UPDATE HESSIAN
             IF (paral%parent) THEN
                icol=3*(iat-1)+k
                iat2=0
                DO is2=1,ions1%nsp
                   DO ia2=1,ions0%na(is2)
                      iat2=iat2+1
                      DO k2=1,3
                         ix=3*(iat2-1)+k2
                         sder(ix,icol)=(fdis(k2,ia2,is2,2)-&
                              fdis(k2,ia2,is2,1))/(2.0_real_8*cntr%fdiff)
                      ENDDO
                   ENDDO
                ENDDO
                ! FINITE DIFFERENCES OF DIPOLE MOMENT
                IF (prop1%ldip.OR.prop1%dberry) THEN
                   DO k2=1,3
                      dipder(k2,icol) = ( dip(k2,iat,2) - dip(k2,iat,1) )&
                           / (2.0_real_8*cntr%fdiff)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
31416     CONTINUE
       ENDDO
    ENDDO
    ! Reset coordinates
    CALL dcopy(3*maxsys%nax*ions1%nsp,taup(1,1,1),1,tau0(1,1,1),1)
    symmi%nrot=nsym
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(24)
       IF ((prop1%ldip.OR.prop1%dberry).AND.paral%io_parent)&
            CALL fileclose(25)
       ! Symmetrize Hessian
       IF ( symmi%indpg.NE.0.AND.(symmi%nrot.NE.1.OR.symmi%ntvec.GT.1) ) THEN
          CALL symmat(sder,2)
       ENDIF
       CALL symma(sder,3*ions1%nat)
       ! Write Hessian to file HESSIAN
       IF (cotc0%nodim.EQ.3*ions1%nat) THEN
          CALL dcopy(cotc0%nodim*cotc0%nodim,sder,1,hess,1)
          CALL hessout
       ELSE
          CALL dcopy(9*ions1%nat*ions1%nat,sder,1,hess,1)
       ENDIF
       ! Write output file for MOLVIB program
       CALL molvib(tau0,sder)
       ! Mass weighted force constants
       ix=0
       DO is=1,ions1%nsp
          xmass=SQRT(1._real_8/rmass%pma0(is))
          DO ia=1,ions0%na(is)
             DO k=1,3
                ix=ix+1
                xma(ix)=xmass
             ENDDO
          ENDDO
       ENDDO
       !$omp parallel do private(I,J)
       DO i=1,3*ions1%nat
          DO j=1,3*ions1%nat
             sder(i,j)=sder(i,j)*xma(i)*xma(j)
          ENDDO
       ENDDO
       ! Diagonalization
       naux=9*ions1%nat*ions1%nat
       info=0
       CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,scr,naux,info)
       IF (info.NE.0) CALL stopgm('SECDER','DSYEV INFO',& 
            __LINE__,__FILE__)
       DO i=1,3*ions1%nat
          vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(/," ",64("*"),/,A,/)')&
            ' HARMONIC FREQUENCIES [cm**-1]:'
       IF (paral%io_parent)&
            WRITE(6,'(4(F16.4))') (vibe(i),i=1,3*ions1%nat)
       ! Write output file for Eigenvectors   
       CALL vibeig(vibe,sder,3*ions1%nat,.TRUE.)
       ! PURIFICATION
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') ' PURIFICATION OF DYNAMICAL MATRIX'
       CALL dcopy(9*ions1%nat*ions1%nat,hess,1,sder,1)
       ! tu     purification only AFTER mass weighting of the Hessian
       ! tu     CALL PURGED(SDER,HESS,TAU0)
       !$omp parallel do private(I,J)
       DO i=1,3*ions1%nat
          DO j=1,3*ions1%nat
             sder(i,j)=sder(i,j)*xma(i)*xma(j)
          ENDDO
       ENDDO
       CALL purged(sder,hess,tau0)
       ! Diagonalization
       naux=9*ions1%nat*ions1%nat
       info=0
       CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,scr,naux,info)
       IF (info.NE.0) CALL stopgm('SECDER','DSYEV INFO',& 
            __LINE__,__FILE__)
       !$omp parallel do private(I)
       DO i=1,3*ions1%nat
          vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(/," ",64("*"),/,A,/)')&
            ' HARMONIC FREQUENCIES [cm**-1]:'
       IF (paral%io_parent)&
            WRITE(6,'(4(F16.4))') (vibe(i),i=1,3*ions1%nat)
       IF (paral%io_parent)&
            WRITE(6,'(A,E16.8)') 'ChkSum(FREQ) = ',SUM(ABS(vibe(1:3*ions1%nat)))
       ! Write output file for Eigenvectors   
       CALL vibeig(vibe,sder,3*ions1%nat,.FALSE.)
       IF (rout1%vout) CALL mldgau(vibe,sder,3*ions1%nat,tau0)
       IF (rout1%acout) CALL writeaclimax(vibe,sder,3*ions1%nat,tau0)
       ! CALCULATE (DERIVATIVES OF) DIPOLE MOMENT IN BASIS OF EIGENMODES
       CALL zeroing(ir)!,3*cotc0%nodim)
       DO mode=1,cotc0%nodim
          DO j=1,cotc0%nodim
             ir(1,mode) = ir(1,mode) + dipder(1,j) * sder(j,mode)
             ir(2,mode) = ir(2,mode) + dipder(2,j) * sder(j,mode)
             ir(3,mode) = ir(3,mode) + dipder(3,j) * sder(j,mode)
          ENDDO
       ENDDO
       ! WRITE BORN EFFECTIVE CHARGES AND IR INTENSITIES
       IF (prop1%ldip.OR.prop1%dberry) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/," ",64("*"),/,A,/)')&
               ' BORN EFFECTIVE CHARGE TENSOR (NOT SYMMETRIZED):'
          IF (3*ions1%nat.EQ.cotc0%nodim) THEN
             iat=0
             j=0
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   iat=iat+1
                   IF (paral%io_parent)&
                        WRITE(6,&
                        '(I5,1X,A,3X,3F10.5,/,6X,A,3X,3F10.5,/,'//&
                        '6X,A,3X,3F10.5)')&
                        iat, 'X', (dipder(k2,1+j), k2=1,3),&
                        'Y', (dipder(k2,2+j), k2=1,3),&
                        'Z', (dipder(k2,3+j), k2=1,3)
                   j=j+3
                ENDDO
             ENDDO
          ELSE
             DO j = 1,cotc0%nodim
                IF (paral%io_parent)&
                     WRITE(6,'(3F10.6)') (dipder(k2,j), k2=1,3)
             ENDDO
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(/," ",64("*"),/,A,/)') ' RELATIVE IR INTENSITIES:'
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,3(10X,A),6X,A)') '      FREQ',' X',' Y',' Z',&
               ' TOTAL'
          DO j=1,cotc0%nodim
             d1 = ir(1,j)
             d2 = ir(2,j)
             d3 = ir(3,j)
             dd=SQRT(d1*d1+d2*d2+d3*d3)
             IF (paral%io_parent)&
                  WRITE(6,'(5F12.5)') vibe(j),d1,d2,d3,dd
          ENDDO
       ENDIF
    ENDIF
200 CONTINUE
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dipder,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tdiag) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (tkpts%tkpnt) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) THEN
       DEALLOCATE(fdis,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(sder,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vibe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xma,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%bsymm)DEALLOCATE(fnbs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE secder

END MODULE secder_driver
