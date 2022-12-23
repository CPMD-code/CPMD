MODULE sdlinres_utils
  USE adat,                            ONLY: elem
  USE canon_utils,                     ONLY: canon,&
                                             give_scr_canon
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE copot_utils,                     ONLY: give_scr_copot
  USE cotr,                            ONLY: cotc0,&
                                             hess,&
                                             lskcor
  USE detdof_utils,                    ONLY: detdof
  USE dynit_utils,                     ONLY: dynit
  USE eicalc_utils,                    ONLY: eicalc,&
                                             eicalc1
  USE eind_ii_utils,                   ONLY: eind_ii
  USE eind_loc_utils,                  ONLY: eind_loc
  USE eind_nl_utils,                   ONLY: eind_nl
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE hessin_utils,                    ONLY: hessin
  USE hessout_utils,                   ONLY: hessout
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03
  USE lr_in_utils,                     ONLY: lr_in
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nl_res_utils,                    ONLY: give_scr_nl_res,&
                                             nl_res
  USE nlcc,                            ONLY: corel
  USE nlps,                            ONLY: ndfnl
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE opt_lr_utils,                    ONLY: give_scr_opt_lr,&
                                             opt_lr
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE rho1ofr_utils,                   ONLY: rho1ofr
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_2d_utils,                  ONLY: give_scr_rnlsm_2d,&
                                             rnlsm_2d
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: infw,&
                                             iteropt,&
                                             ropt_mod
  USE sd_ii_utils,                     ONLY: sd_ii
  USE sd_loc2_utils,                   ONLY: sd_loc2
  USE sd_loc_utils,                    ONLY: sd_loc
  USE sd_nl2_utils,                    ONLY: sd_nl2
  USE sd_nl_utils,                     ONLY: sd_nl
  USE secder_utils,                    ONLY: mldgau,&
                                             molvib,&
                                             purged,&
                                             vibeig,&
                                             writeaclimax
  USE setirec_utils,                   ONLY: write_irec
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             fnl,&
                                             fnl2
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: restart1,&
                                             rout1
  USE symm,                            ONLY: isymu,&
                                             symmi
  USE symtrz_utils,                    ONLY: give_scr_symmat,&
                                             symmat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             nacc,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE utils,                           ONLY: symma
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sdlinres
  PUBLIC :: give_scr_sdlinres

CONTAINS

  ! ==================================================================
  SUBROUTINE sdlinres(c0,c2,sc0,c1,eigv)
    ! ==--------------------------------------------------------------==
    ! ==  VIBRATIONAL ANALYSIS USING LINEAR RESPONSE                  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,1), c2(ncpw%ngw,crge%n), &
      sc0(ncpw%ngw,crge%n), c1(ncpw%ngw,crge%n)
    REAL(real_8)                             :: eigv(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'sdlinres'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/'X','Y','Z'/)

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eirop1(:), &
                                                eivps(:), eivps1(:), gde(:), &
                                                h1nl(:), pme(:), psi(:,:)
    EXTERNAL                                 :: dasum
    INTEGER :: i, ia, ia2, iat, iat2, icol, id, ierr, ig, ij, il_psi_1d, &
      il_psi_2d, il_rhoe_1d, il_rhoe_2d, info, irec(100), is, is2, isymuat, &
      ix, j, k, k2, lddxc_1d, lddxc_2d, ldfnl, lscr, naux, ngde, npme, nsym, &
      nvpp
    INTEGER, ALLOCATABLE                     :: iato(:,:)
    LOGICAL                                  :: ferror, readhe
    REAL(real_8) :: dasum, detot, eind, ekin1, ekin2, ekincp, ekinh1, ekinh2, &
      etotg, etoto, tcpu, temp1, temp2, thl(2), time1, time2, trace, vp, xmass
    REAL(real_8), ALLOCATABLE                :: ddxc(:,:), rhoe(:,:), scr(:), &
                                                sder(:,:), tscr(:,:,:), &
                                                vibe(:), vpp(:), xma(:)
    REAL(real_8), POINTER                    :: dfnl1(:,:,:,:,:,:), &
                                                fnl1(:,:,:,:,:), &
                                                fnl3(:,:,:,:,:)

! ==================================================================

    time1 =m_walltime()
    ! ==--------------------------------------------------------------==
    ! ..input for linear response
    CALL lr_in
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' K-POINTS WITH LINEAR RESPONSE NOT IMPLEMENTED'
       CALL stopgm('SECDER','WRONG OPTIONS',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(fion)!,SIZE(fion))
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    IF (paral%parent) THEN
       ALLOCATE(xma(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vibe(3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(sder(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sder)!,9*ions1%nat*ions1%nat)
    nacc = 7
    ener_com%ecnstr = 0.0_real_8
    ener_com%erestr = 0.0_real_8
    ! ..atom pointer
    id=(maxsys%nax*maxsys%nsx)/2+1
    ALLOCATE(iato(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(iato)!,maxsys%nax*maxsys%nsx)
    ij=0
    DO is=1,ions1%nsp
       DO iat=1,ions0%na(is)
          iato(iat,is)=ij
          ij=ij+1
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_sdlinres(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    ALLOCATE(rhoo(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(potr(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ALLOCATE(h1nl(nkpt%ngwk*crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eivps1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eirop1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (lr03%txc_analytic) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=(2*clsd%nlsd-1)
    ELSEIF (lr01%lopti.EQ.0.OR.lr01%lopti.EQ.2) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=(2*clsd%nlsd-1)
    ELSE
       lddxc_1d=1
       lddxc_2d=1
    ENDIF
    ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (isos3%ps_type.EQ.1) CALL stopgm('SDLINRES','HOCKNEY PS not impl.',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! ==--------------------------------------------------------------==
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    ! END INITIALIZATION
    ! ..  Do not symmtrize density...  
    cntl%tsymrho=.FALSE.
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull)) THEN
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    time1=m_walltime()
    ! ..allocate additional memory for optimisers
    IF (cntl%tsde) THEN
       npme = 1
       ngde = 1
       nvpp = 1
    ELSE IF (cntl%diis) THEN
       npme = (ncpw%ngw*crge%n+8)*cnti%mdiis/2 !vw rm 2*
       ngde = ((ncpw%ngw*crge%n+8)*cnti%mdiis)/2 !vw rm 2*
       nvpp = ncpw%ngw
    ELSE IF (cntl%pcg) THEN
       npme = ncpw%ngw*crge%n !vw rm 2*
       ngde = 1
       nvpp = ncpw%ngw
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
       CALL stopgm('SDLINRES',' ',& 
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
    DO infw=1,cnti%nomore
       ! UPDATE THE WAVEFUNCTIONS
       CALL updwf(c0(:,:,1),c2,sc0,tau0,tscr,pme,gde,vpp,eigv,&
            rhoe,psi,crge%n,.FALSE.,.TRUE.)
       ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%parent) THEN
          detot=ener_com%etot-etoto
          IF (infw.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etoto,&
               tcpu,gemax,cnorm,thl,ropt_mod%engpri,infw,infw)
          etoto=ener_com%etot
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf.OR.soft_com%exsoft) GOTO 100
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' NO CONVERGENCE IN ',infw,' STEPS '
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! ..release memory for ground state optimiser
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! NUCLEAR GRADIENT
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.FALSE.,.TRUE.)
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
    ! Open the FINDIF file
    IF (paral%parent) THEN
       filen='FINDIF'
       IF (paral%io_parent)&
            CALL fileopen(24,filen,fo_def,ferror)
       IF (restart1%rvib.AND.paral%io_parent) REWIND(24)
    ENDIF
    ! ..transform to canonical orbitals and get eigenvalues
    CALL canon(c0,c2,crge%f,crge%n,eigv)
    ! ..recalculate FNL and DFNL for transformed orbitals
    ! ..we could also rotate the old ones, but this is easier
    CALL rnlsm(c0(:,:,1),crge%n,1,1,.TRUE.)
    ! ..store potential and density 
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
    CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)
    ! ..NLCC
    IF (corel%tinlc) THEN
       CALL stopgm("SDLINRES","NLCC not implemented",& 
            __LINE__,__FILE__)
    ENDIF
    ! ..alocate, calculate, and save DDFNL
    ldfnl=6*ions1%nat*maxsys%nhxs*ndfnl
    IF (ldfnl.LE.0) ldfnl=1
    ALLOCATE(ddfnl(ions1%nat,maxsys%nhxs,6,ndfnl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rnlsm_2d(c0,crge%n)
    ! ..save FNL
    CALL fnl_set("SAVE")
    CALL fnlalloc(crge%n,.TRUE.,.FALSE.)
    ! ..calculate and store EIVPS and EIROP
    CALL eicalc(eivps,eirop)
    ! ..calculate and store analytic dmu/dn
    CALL lr_xcpot(ddxc,rhoo,.FALSE.)
    ! ==--------------------------------------------------------------==
    ! ..additional memory for optimisers
    nvpp = ncpw%ngw
    IF (lr01%lopti.EQ.0) THEN
       npme = MAX((ncpw%ngw*crge%n+8)*cnti%mdiis/2,ncpw%ngw*crge%n) !vw rm 2*
       npme = MAX(ncpw%ngw*crge%n*lr01%mldiis,npme) !vw rm 2*
       ngde = MAX(((ncpw%ngw*crge%n+8)*cnti%mdiis)/2,1) !vw rm 2*
       ngde = MAX(ncpw%ngw*crge%n*lr01%mldiis,ngde) !vw rm 2*
    ELSEIF (lr01%lopti.EQ.1) THEN
       npme = 1
       ngde = 1
    ELSE IF (lr01%lopti.EQ.2) THEN
       npme = ncpw%ngw*crge%n !vw rm 2*
       ngde = 1
    ELSE IF (lr01%lopti.EQ.3) THEN
       npme = (ncpw%ngw*crge%n+8)*cnti%mdiis/2 !vw rm 2*
       ngde = ((ncpw%ngw*crge%n+8)*cnti%mdiis)/2 !vw rm 2*
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' WRONG OPTION FOR LINEAR RESPONSE OPTIMIZATION'
       CALL stopgm('SDLINRES',' ',& 
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
    cntl%prec=.TRUE.
    CALL ksdiag(vpp)
    trace=dasum(crge%n,eigv,1)/REAL(crge%n,kind=real_8)
    !$omp parallel do private(IG,VP) shared(TRACE)
    DO ig=1,ncpw%ngw
       vp=vpp(ig)+trace
       vpp(ig)=ABS(vp/(vp**2+lr02%lr_hthrs**2))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==         THE BASIC LOOP OVER PERTURBATIONS (COORDINATES)      ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",23("*"),A,24("*"),/)') ' LINEAR RESPONSE '
    ENDIF
    ! No symmetrization in this part!
    nsym=symmi%nrot
    symmi%nrot=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             time1=m_walltime()
             IF (paral%io_parent)&
                  WRITE(6,&
                  '(" **** ATOM=",I7,4X,A4,T49,A,A4)')&
                  iat,elem%el(ions0%iatyp(is)),'PERTURBATION:',cdir(k)
             ! ..new perturbation
             IF (restart1%rvib .AND. paral%parent) THEN
                IF (paral%io_parent)&
                     READ(24,fmt=*,END=300)
                DO is2=1,ions1%nsp
                   DO ia2=1,ions0%na(is2)
                      IF (paral%io_parent)&
                           READ(24,*) (fion(k2,ia2,is2),k2=1,3)
                   ENDDO
                ENDDO
                GOTO 301
300             CONTINUE
                restart1%rvib=.FALSE.
                IF (paral%io_parent)&
                     CALL fileclose(24)
                IF (paral%io_parent)&
                     CALL fileopen(24,filen,fo_app,ferror)
             ENDIF
301          CONTINUE
             CALL mp_bcast(restart1%rvib,parai%source,parai%allgrp)
             IF (restart1%rvib.AND..NOT.paral%parent) CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
             IF (.NOT.restart1%rvib) THEN
                isymuat=0
                IF (symmi%indpg.NE.0) THEN
                   IF (isymu(iat).EQ.0) THEN
                      CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
                      isymuat=1
                   ENDIF
                ENDIF
                IF (isymuat.EQ.0.AND.lskcor(k,iat).NE.0) THEN
                   ! ---------------------------------------------------------------------==
                   ! ..this is the response part
                   eind=0._real_8
                   CALL eind_ii(tau0,is,ia,k,eind,iteropt%iesr)
                   CALL eind_loc(eind,is,iat,k,rhoo,psi(:,1),eirop)
                   ! ..calculate constant force from NL-PP
                   CALL fnl_set("SWITCH")
                   CALL nl_res(is,iat,k,h1nl,crge%n)
                   CALL eind_nl(eind,is,iat,k,crge%n)
                   CALL fnl_set("SWITCH")
                   ! ..calculate the constant part of the local potential
                   CALL eicalc1(k,is,iat,eivps1,eirop1)
                   ! ..calculate first order wavefunction
                   ropt_mod%sdiis=.TRUE.
                   ropt_mod%spcg=.TRUE.
                   CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                        h1nl,eivps1,eirop1,ddxc,vpp,pme,gde,&
                        psi,crge%n,"PHONON","CANON")
                   ! ..calculate first order density
                   CALL rho1ofr(c0,c1,crge%f(:,1),rhoe,psi(:,1),crge%n)
                   ! ..contribution from local potentials
                   CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
                   CALL sd_loc(fion,rhoe,psi(:,1),eirop1)
                   ! ..calculate the FNL and DFNL for C1
                   CALL rnlsm(c1,crge%n,1,1,.TRUE.)
                   fnl1 => fnl
                   fnl3 => fnl2
                   dfnl1 => dfnl
                   CALL fnl_set("SWITCH")
                   CALL sd_nl(fion,crge%f,fnl,fnl2,dfnl,fnl1,fnl3,dfnl1,crge%n)
                   CALL fnl_set("SWITCH")
                   ! ..end of the response part
                   ! ---------------------------------------------------------------------==
                   ! ..Write response to disk
                   CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,tscr,1)
                   CALL mp_sum(tscr,3*maxsys%nax*maxsys%nsx,parai%allgrp)
                   IF (paral%parent) THEN
                      IF (paral%io_parent)&
                           WRITE(24,'(3I10)') is,ia,k
                      DO is2=1,ions1%nsp
                         DO ia2=1,ions0%na(is2)
                            IF (paral%io_parent)&
                                 WRITE(24,*) (tscr(k2,ia2,is2),k2=1,3)
                         ENDDO
                      ENDDO
                      CALL m_flush(24)
                   ENDIF
                ELSE
                   IF (paral%io_parent)&
                        WRITE(6,'(A)') ' Pertubation related by symmetry'
                ENDIF
                IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
                IF (soft_com%exsoft) GOTO 200
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(T33,A)') 'Result read from file FINDIF'
             ENDIF
             ! Update Hessian
             icol=3*(iat-1)+k
             iat2=0
             DO is2=1,ions1%nsp
                DO ia2=1,ions0%na(is2)
                   iat2=iat2+1
                   DO k2=1,3
                      ix=3*(iat2-1)+k2
                      sder(ix,icol)=fion(k2,ia2,is2)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    symmi%nrot=nsym
    ! ==--------------------------------------------------------------==
    CALL fnldealloc(.FALSE.,.FALSE.)
    CALL fnl_set("RECV")
    ! ..contribution from ESR
    CALL sd_ii(tau0,sder,iteropt%iesr,iato)
    ! ..contribution from <psi_0|V(2)|psi_0>
    ! ..local part(there are diagonal terms only)
    CALL sd_loc2(sder,rhoo,psi(:,1),iato,eirop)
    ! ..nonlocal part(there are diagonal terms only)
    CALL sd_nl2(sder,iato,crge%n)
    DEALLOCATE(ddfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL mp_sum(sder,9*ions1%nat*ions1%nat,parai%allgrp)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(24)
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
       !$omp parallel do private(I)
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
    ENDIF
200 CONTINUE
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(iato,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       DEALLOCATE(vibe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xma,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(h1nl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eivps1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sder,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! TODO refactor RHOO 
    DEALLOCATE(rhoo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(potr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ..release memory for LR optimiser
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sdlinres
  ! ==================================================================
  SUBROUTINE give_scr_sdlinres(lsecder,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsecder
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcanon, lcopot, lfnonloc_p, lforces, linitrun, lnl_res, &
      lopt_lr, lortho, lrhoofr, lrnlsm, lrnlsm_2d, lsymmat, lupdate, nstate

    nstate=crge%n
    linitrun=0
    lforces=0
    lcopot=0
    lortho=0
    lupdate=0
    lcanon=0
    lsymmat=0
    lrhoofr=0
    lrnlsm=0
    lrnlsm_2d=0
    lopt_lr=0
    lfnonloc_p=0
    ! 
    CALL give_scr_initrun(linitrun,tag)
    CALL give_scr_forcedr(lforces,tag,nstate,.FALSE.,.TRUE.)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    CALL give_scr_ortho(lortho,tag,nstate)
    CALL give_scr_updwf(lupdate,tag,nstate,.FALSE.)
    CALL give_scr_canon(lcanon,tag,nstate)
    CALL give_scr_symmat(lsymmat,tag)
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_nl_res(lnl_res,nstate,tag)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    CALL give_scr_rnlsm_2d(lrnlsm_2d,tag,nstate)
    CALL give_scr_opt_lr(lopt_lr,"PHONON",tag)
    ! 
    lsecder=MAX(9*ions1%nat*ions1%nat+9*ions1%nat,linitrun,lforces,lcopot,&
         lortho,lupdate,lcanon,lsymmat,lnl_res,lrhoofr,lrnlsm,&
         lrnlsm_2d,lopt_lr,lfnonloc_p)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_sdlinres
  ! ==================================================================

END MODULE sdlinres_utils
