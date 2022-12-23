MODULE ohlr_utils
  USE atwf,                            ONLY: atwp,&
                                             catom,&
                                             loadc_foc_array_size,&
                                             xsmat
  USE canon_utils,                     ONLY: canon,&
                                             give_scr_canon
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE cppt,                            ONLY: nzh
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE dotp_utils,                      ONLY: dotp
  USE dynit_utils,                     ONLY: dynit
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos3
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03,&
                                             lrhd,&
                                             lrhi
  USE localize_utils,                  ONLY: localize
  USE lr_in_utils,                     ONLY: lr_in
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE opt_lr_utils,                    ONLY: give_scr_opt_lr,&
                                             opt_lr
  USE ovlap_utils,                     ONLY: dmatc,&
                                             ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE puttau_utils,                    ONLY: getatm
  USE rho1ofr_utils,                   ONLY: rho1ofr,&
                                             rhoabofr,&
                                             rhosofr
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE setbasis_utils,                  ONLY: loadc,&
                                             rtog
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE utils,                           ONLY: idamin
  USE v1ofrho1_utils,                  ONLY: give_scr_v1ofrho1,&
                                             v1ofrho1
  USE vhk_utils,                       ONLY: give_scr_vhk,&
                                             vhk
  USE vpsi_utils,                      ONLY: vpsimt
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ohlr
  !public :: give_scr_ohlr
  !public :: rhobas
  !public :: tsumhk

CONTAINS

  ! ==================================================================
  SUBROUTINE ohlr(c0,c1,c2,sc0)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), &
                                                c1(ncpw%ngw,*), &
                                                c2(ncpw%ngw,crge%n), &
                                                sc0(ncpw%ngw,crge%n)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ohlr'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: edummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: gde(:), h1nl(:,:), pme(:), &
                                                psi(:,:)
    EXTERNAL                                 :: dasum
    INTEGER :: i, ierr, ig, ii, iii, ik, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, irec(100), is, isub, j, jj, k, kk, knfw, kst, l, lddxc_1d, &
      lddxc_2d, lscr, ngde, npme, nvpp
    INTEGER, ALLOCATABLE                     :: list(:,:), listst(:)
    LOGICAL                                  :: ferror
    REAL(real_8) :: cd(3), cref(3), dasum, detot, e2(2), eind, ekin1, ekin2, &
      ekincp, ekinh1, ekinh2, etoto, tcpu, temp1, temp2, thl(2), time1, &
      time2, trace, vp
    REAL(real_8), ALLOCATABLE :: center(:), ddxc(:,:), eigv(:), fhk(:), &
      fone(:), hmat(:,:), hmate(:,:,:,:), hnsc(:,:), hnsce(:,:,:,:), &
      rho1(:,:), rhoe(:,:), scr(:), tvec(:), tvece(:,:), vpp(:)

    CALL tiset(procedureN,isub)

    ! ==================================================================
    cntl%tsymrho=.FALSE.
    edummy=0._real_8
    ! ..input for linear response
    CALL lr_in
    ! 
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    ! Memory for densities
    ALLOCATE(rho1(fpar%nnr1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ALLOCATE(listst(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (lrhd%local_orb) THEN
       ALLOCATE(eigv(crge%n*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(eigv(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(fone(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fhk(fpar%nnr1*clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (lrhd%diagonal) THEN
       ALLOCATE(hmat(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hnsc(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tvec(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hmat)!,n*n)
       CALL zeroing(hnsc)!,n*n)
       CALL zeroing(tvec)!,n)
    ELSE
       ALLOCATE(list(2,crge%n**2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hmate(crge%n,crge%n,crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hnsce(crge%n,crge%n,crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tvece(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hmate)!,n*n*n*n)
       CALL zeroing(hnsce)!,n*n*n*n)
       CALL zeroing(tvece)!,n*n)
       ii=0
       DO i=1,crge%n
          DO j=i,crge%n
             ii=ii+1
             list(1,ii)=i
             list(2,ii)=j
          ENDDO
       ENDDO
    ENDIF
    DO i=1,crge%n
       fone(i)=1._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (cntl%tlsd.AND.paral%parent) CALL stopgm('OHLR','LSD not implemented',& 
         __LINE__,__FILE__)
    ALLOCATE(h1nl(ncpw%ngw,crge%n*nkpt%ngwk/ncpw%ngw),STAT=ierr)
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
    IF (isos3%ps_type.EQ.1) CALL stopgm('OHLR','HOCKNEY PS not impl.',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d, &
         il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_ohlr(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    ALLOCATE(rhoo(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(potr(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! Set IREC to restart file.
    CALL read_irec(irec)
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            '   REFERENCE POINT'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
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
    ELSEIF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
       CALL stopgm('OHLR',' ',& 
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
    etoto=0._real_8
    DO knfw=1,cnti%nomore
       ! ..UPDATE THE WAVEFUNCTIONS
       CALL updwf(c0(:,:,1),c2,sc0,tau0,taup,pme,gde,vpp,eigv,&
            rhoe,psi,crge%n,.FALSE.,.TRUE.)
       ! ..PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%parent) THEN
          detot=ener_com%etot-etoto
          IF (knfw.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etoto,&
               tcpu,gemax,cnorm,thl,ropt_mod%engpri,knfw,knfw)
          etoto=ener_com%etot
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf.OR.soft_com%exsoft) GOTO 200
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' NO CONVERGENCE IN ',knfw,' STEPS '
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,knfw)
    CALL stopgm('OHLR',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
200 CONTINUE
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
    ! ..localise orbitals
    IF (lrhd%local_orb) CALL localize(tau0,c0,c2,sc0,crge%n)
    ! ..calculate KS-Matrix 
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.FALSE.,.FALSE.)
    time2=m_walltime()
    tcpu=(time2-time1)*0.001_real_8
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,T46,F12.3,A)') ' TIME FOR MINIMUM STRUCTURE :',&
            tcpu,' SECONDS'
       CALL wrener
    ENDIF
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    IF (lrhd%local_orb) THEN
       CALL ovlap(crge%n,eigv,c0(:,:,1),c2)
       CALL mp_sum(eigv,crge%n*crge%n,parai%allgrp)
    ELSE
       ! ..transform to canonical orbitals and get eigenvalues
       CALL canon(c0,c2,crge%f,crge%n,eigv)
       ! ..recalculate FNL for transformed orbitals
       ! ..we could also rotate the old ones, but this is easier
       CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
    ENDIF
    ! ..store potential and density
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
    CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)
    ! ..NLCC
    IF (corel%tinlc) THEN
       CALL stopgm("OHLR","NLCC not implemented",& 
            __LINE__,__FILE__)
    ENDIF
    ! ..HK Potential estimate
    CALL vhk(rhoe,fhk)
    ! ..save FNL
    CALL fnl_set("SAVE")
    CALL fnlalloc(crge%n,.FALSE.,.FALSE.)
    ! ..calculate and store analytic dmu/dn
    CALL lr_xcpot(ddxc,rhoo,.FALSE.)
    ! ==--------------------------------------------------------------==
    ! ..define reference orbitals
    IF (lrhi%numo.LT.0) THEN
       ! ..all orbitals
       lrhi%numo=crge%n
       DO i=1,lrhi%numo
          listst(i)=i
       ENDDO
    ELSE
       IF (lrhi%refat.LT.0 .OR. .NOT.lrhd%local_orb) THEN
          ! ..last numo orbitals in the list
          DO i=1,lrhi%numo
             listst(i)=crge%n-lrhi%numo+i
          ENDDO
       ELSE
          ! ..local orbitals
          IF (paral%io_parent)&
               CALL fileopen(30,'WANNIER_CENTER',fo_old,ferror)
          IF (ferror) CALL stopgm('OHLR','COULD NOT OPEN WANNIER_CENTER'&
               ,& 
               __LINE__,__FILE__)
          IF (paral%io_parent) REWIND(30)
          ALLOCATE(center(crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL getatm(lrhi%refat,tau0,cref)
          DO i=1,crge%n
             IF (paral%io_parent)&
                  READ(30,*) j,cd(1),cd(2),cd(3)
             cd(1)=cd(1)-cref(1)
             cd(2)=cd(2)-cref(2)
             cd(3)=cd(3)-cref(3)
             CALL pbc(cd(1),cd(2),cd(3),cd(1),cd(2),cd(3),1,parm%apbc,parm%ibrav)
             center(i)=cd(1)**2 + cd(2)**2 + cd(3)**2
          ENDDO
          DO i=1,lrhi%numo
             j=idamin(crge%n,center,1)
             listst(i)=j
             center(j)=1.e30_real_8
          ENDDO
          IF (paral%io_parent) CALL fileclose(30)
          DEALLOCATE(center,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF REFERENCE CALCULATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="),/)')
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    ! ..get basis for density expansions
    CALL rhobas("INIT",rhoo,psi(:,1))
    CALL rhobas("EXPAND",rhoo,psi(:,1))
    ! ..additional memory for optimisers
    nvpp = ncpw%ngw
    IF (lr01%lopti.EQ.0) THEN
       npme = MAX((ncpw%ngw*crge%n+8)*cnti%mdiis/2,ncpw%ngw*crge%n) !vw rm 2*
       npme = MAX(ncpw%ngw*crge%n*lr01%mldiis,npme) !vw rm 2*
       ngde = MAX(((ncpw%ngw*crge%n+8)*cnti%mdiis)/2,1) !vw rm 2*
       ngde = MAX(ncpw%ngw*crge%n*lr01%mldiis,ngde) !vw rm 2*
    ELSE IF (lr01%lopti.EQ.1) THEN
       npme = 1
       ngde = 1
    ELSE IF (lr01%lopti.EQ.2) THEN
       npme = ncpw%ngw*crge%n !vw rm 2*
       ngde = 1
    ELSE IF (lr01%lopti.EQ.3) THEN
       npme = (ncpw%ngw*crge%n+8)*cnti%mdiis/2  !vw rm 2*
       ngde = ((ncpw%ngw*crge%n+8)*cnti%mdiis)/2  !vw rm 2*
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' WRONG OPTION FOR LINEAR RESPONSE OPTIMIZATION'
       CALL stopgm('OHLR',' ',& 
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
    DO ig=1,ncpw%ngw
       vp=vpp(ig)+trace
       vpp(ig)=ABS(vp/(vp**2+lr02%lr_hthrs**2))
    ENDDO
    ! Open the FINDIF file
    filen='FINDIF'
    IF (paral%io_parent)&
         CALL fileopen(24,filen,fo_def,ferror)
    IF (restart1%roh.AND.paral%io_parent) REWIND(24)
    ! ..start linear response calculation
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",23("*"),A,24("*"),/)') ' LINEAR RESPONSE '
    ENDIF
    IF (lrhd%diagonal) THEN
       ! ..diagonal orbital hardness
       DO kst=1,lrhi%numo
          is=listst(kst)
          CALL mp_sync(parai%allgrp)
          IF (paral%io_parent)&
               WRITE(6,'(" **** STATE=",I7)') is
          IF (restart1%roh.AND.paral%parent) THEN
             IF (paral%io_parent)&
                  READ(24,fmt=*,END=300)
             IF (paral%io_parent)&
                  READ(24,*) tvec(is)
             IF (paral%io_parent)&
                  READ(24,*) (hmat(i,is),i=1,crge%n)
             GOTO 301
300          CONTINUE
             restart1%roh=.FALSE.
          ENDIF
301       CONTINUE
          CALL mp_bcast(restart1%roh,parai%source,parai%allgrp)
          IF (.NOT.restart1%roh) THEN
             ! ..linear response
             ! ..calculate rho(i)
             CALL rhosofr(c0(:,is,1),rho1(:,1),psi(:,1))
             ! ..calculate EH(2) + EXC(2) for the state density
             ! ..on output RHO1 is the associated potential
             CALL v1ofrho1(e2,rho1,ddxc,psi)
             eind=e2(1)+e2(2)
             ! ..apply this potential to the ground state functions
             CALL zeroing(h1nl)!,ngw*n)
             CALL vpsimt(c0(:,:,1),h1nl,crge%f(:,1),rho1,psi(:,1),crge%n,clsd%nlsd,.FALSE.)
             ! ..calculate first order wavefunction
             ropt_mod%sdiis=.TRUE.
             ropt_mod%spcg=.TRUE.
             IF (lrhd%local_orb) THEN
                CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                     h1nl,edummy,edummy,ddxc,vpp,pme,gde,&
                     psi,crge%n,"ORBHARD","LOCAL")
             ELSE
                CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                     h1nl,edummy,edummy,ddxc,vpp,pme,gde,&
                     psi,crge%n,"ORBHARD","CANON")
             ENDIF
             ! ..calculate first order density
             CALL rho1ofr(c0,c1,crge%f(:,1),rhoe,psi(:,1),crge%n)
             CALL tsumhk(fhk,rhoe,tvec(is))
             CALL rhobas("EXPAND",rhoe,psi(:,1))
             CALL v1ofrho1(e2,rhoe,ddxc,psi)
             CALL zeroing(c2)!,ngw*n)
             CALL vpsimt(c0(:,:,1),c2,fone,rhoe,psi(:,1),crge%n,clsd%nlsd,.FALSE.)
             DO j=1,lrhi%numo
                i=listst(j)
                hnsc(i,is)=-0.5_real_8*dotp(ncpw%ngw,c0(:,i,1),h1nl(:,i))/crge%f(i,1)
                hmat(i,is)=-0.5_real_8*dotp(ncpw%ngw,c0(:,i,1),c2(:,i))
                hmat(i,is)=hmat(i,is)+hnsc(i,is)
             ENDDO
             CALL mp_sum(hnsc(:,is),crge%n,parai%allgrp)
             CALL mp_sum(hmat(:,is),crge%n,parai%allgrp)
             ! ..end of the response part
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(24,'(I10)') is
                IF (paral%io_parent)&
                     WRITE(24,*) tvec(is)
                IF (paral%io_parent)&
                     WRITE(24,*) (hmat(i,is),i=1,crge%n)
                CALL m_flush(24)
             ENDIF
          ENDIF
          IF (soft_com%exsoft) GOTO 100
       ENDDO
    ELSE
       ! ..extended orbital hardness
       DO kst=1,lrhi%numo
          DO kk=kst,lrhi%numo
             is=listst(kst)
             ik=listst(kk)
             CALL mp_sync(parai%allgrp)
             IF (paral%io_parent)&
                  WRITE(6,'(" **** STATES=",I7,I7)') is,ik
             IF (restart1%roh.AND.paral%parent) THEN
                IF (paral%io_parent)&
                     READ(24,fmt=*,END=500)
                IF (paral%io_parent)&
                     READ(24,*) tvece(is,ik)
                IF (paral%io_parent)&
                     READ(24,*) ((hmate(i,j,is,ik),i=1,crge%n),j=1,crge%n)
                GOTO 501
500             CONTINUE
                restart1%roh=.FALSE.
             ENDIF
501          CONTINUE
             CALL mp_bcast(restart1%roh,parai%source,parai%allgrp)
             IF (.NOT.restart1%roh) THEN
                ! ..linear response
                ! ..calculate rho(i)
                CALL zeroing(rho1)!,nnr1)
                CALL rhoabofr(1,c0(:,is:is,1),c0(:,ik:ik,1),rho1(:,1),psi(:,1))
                ! ..calculate EH(2) + EXC(2) for the state density
                ! ..on output RHO1 is the associated potential
                CALL v1ofrho1(e2,rho1,ddxc,psi)
                eind=e2(1)+e2(2)
                ! ..apply this potential to the ground state functions
                CALL zeroing(h1nl)!,ngw*n)
                CALL vpsimt(c0(:,:,1),h1nl,crge%f(:,1),rho1,psi(:,1),crge%n,clsd%nlsd,.FALSE.)
                ! ..calculate first order wavefunction
                ropt_mod%sdiis=.TRUE.
                ropt_mod%spcg=.TRUE.
                IF (lrhd%local_orb) THEN
                   CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                        h1nl,edummy,edummy,ddxc,vpp,pme,gde,&
                        psi,crge%n,"ORBHARD","LOCAL")
                ELSE
                   CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                        h1nl,edummy,edummy,ddxc,vpp,pme,gde,&
                        psi,crge%n,"ORBHARD","CANON")
                ENDIF
                ! ..calculate first order density
                CALL rho1ofr(c0,c1,crge%f(:,1),rhoe,psi(:,1),crge%n)
                CALL tsumhk(fhk,rhoe,tvece(is,ik))
                CALL rhobas("EXPAND",rhoe,psi(:,1))
                CALL v1ofrho1(e2,rhoe,ddxc,psi)
                CALL zeroing(c2)!,ngw*n)
                CALL vpsimt(c0(:,:,1),c2,fone,rhoe,psi(:,1),crge%n,clsd%nlsd,.FALSE.)
                DO ii=1,lrhi%numo
                   DO jj=ii,lrhi%numo
                      i=listst(ii)
                      j=listst(jj)
                      hnsce(i,j,is,ik)=-0.5_real_8*dotp(ncpw%ngw,c0(:,j,1),h1nl(:,i))/&
                           crge%f(i,1)
                      hmate(i,j,is,ik)=-0.5_real_8*dotp(ncpw%ngw,c0(:,j,1),c2(:,i))
                      hmate(i,j,is,ik)=hmate(i,j,is,ik)+hnsce(i,j,is,ik)
                   ENDDO
                ENDDO
                CALL mp_sum(hnsce(:,:,is,ik),crge%n*crge%n,parai%allgrp)
                CALL mp_sum(hmate(:,:,is,ik),crge%n*crge%n,parai%allgrp)
                ! ..end of the response part
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(24,'(2I10)') is,ik
                   IF (paral%io_parent)&
                        WRITE(24,*) tvece(is,ik)
                   IF (paral%io_parent)&
                        WRITE(24,*) ((hmate(i,j,is,ik),i=1,crge%n),j=1,crge%n)
                   CALL m_flush(24)
                ENDIF
             ENDIF
             IF (soft_com%exsoft) GOTO 100
          ENDDO
       ENDDO
    ENDIF
    CALL rhobas("FREE",rhoo,psi(:,1))
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,A)')&
            ' *********************** ORBITAL HARDNESS',&
            ' ***********************'
       IF (lrhd%diagonal) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          DO i=1,lrhi%numo,8
             ii=MIN(8,lrhi%numo-i+1)
             IF (paral%io_parent)&
                  WRITE(6,'(6X,8I8)') (listst(iii),iii=i,i+ii-1)
             DO j=1,lrhi%numo
                IF (paral%io_parent)&
                     WRITE(6,'(I6,8F8.3)') listst(j),&
                     (hmat(listst(j),listst(iii)),iii=i,i+ii-1)
             ENDDO
          ENDDO
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*)
          DO k=1,(lrhi%numo*(lrhi%numo+1))/2,8
             l=MIN(8,(lrhi%numo*(lrhi%numo+1))/2-k+1)
             IF (paral%io_parent)&
                  WRITE(6,'(9X,8(I4,I3))')&
                  (listst(list(1,ii)),listst(list(2,ii)),ii=k,k+l-1)
             DO i=1,lrhi%numo
                DO j=i,lrhi%numo
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,2I4,8F7.3)')&
                        listst(i),listst(j),&
                        (hmate(listst(i),listst(j),&
                        listst(list(1,ii)),listst(list(2,ii))),ii=k,k+l-1)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            ' ****************************************',&
            '************************'
       IF (lrhd%diagonal) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               CALL fileopen(20,'HARDNESS',fo_def,ferror)
          DO i=1,lrhi%numo
             IF (paral%io_parent)&
                  WRITE(20,*) (hmat(listst(i),listst(j)),j=1,lrhi%numo)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(20)
          IF (paral%io_parent)&
               WRITE(6,'(A)')' Hardness in frozen wavefunction approximation'
          DO i=1,lrhi%numo,8
             ii=MIN(8,lrhi%numo-i+1)
             IF (paral%io_parent)&
                  WRITE(6,'(6X,8I8)') (listst(iii),iii=i,i+ii-1)
             DO j=1,lrhi%numo
                IF (paral%io_parent)&
                     WRITE(6,'(I6,8F8.3)') listst(j),&
                     (hnsc(listst(j),listst(iii)),iii=i,i+ii-1)
             ENDDO
          ENDDO
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               CALL fileopen(20,'HARDNESS',fo_def,ferror)
          DO i=1,lrhi%numo
             DO j=i,lrhi%numo
                IF (paral%io_parent)&
                     WRITE(20,*) (hnsce(listst(i),listst(j),&
                     listst(list(1,k)),listst(list(2,k))),&
                     k=1,(lrhi%numo*(lrhi%numo+1))/2)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(20)
          DO k=1,(lrhi%numo*(lrhi%numo+1))/2,8
             l=MIN(8,(lrhi%numo*(lrhi%numo+1))/2-k+1)
             IF (paral%io_parent)&
                  WRITE(6,'(9X,8(I4,I3))')&
                  (listst(list(1,ii)),listst(list(2,ii)),ii=k,k+l-1)
             DO i=1,lrhi%numo
                DO j=i,lrhi%numo
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,2I4,8F7.3)') listst(i),listst(j),&
                        (hnsce(listst(i),listst(j),&
                        listst(list(1,ii)),listst(list(2,ii))),ii=k,k+l-1)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            ' ****************************************',&
            '************************'
       IF (lrhd%diagonal) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' Hohenberg-Kohn Term '
          DO i=1,lrhi%numo
             IF (paral%io_parent)&
                  WRITE(6,'(I6,10X,F8.3)') listst(i),tvec(listst(i))
          ENDDO
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) ' Hohenberg-Kohn Term '
          DO i=1,lrhi%numo
             DO j=1,lrhi%numo
                IF (paral%io_parent)&
                     WRITE(6,'(2I6,10X,F8.3)') listst(i),listst(j),&
                     tvece(listst(i),listst(j))
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            ' ****************************************',&
            '************************'
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(24)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(listst,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fhk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lrhd%diagonal) THEN
       DEALLOCATE(hmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hnsc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tvec,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(list,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hmate,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hnsce,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tvece,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(h1nl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fone,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ohlr
  ! ==================================================================
  SUBROUTINE give_scr_ohlr(lohlr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lohlr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcanon, lddipo, lforces, &
                                                linitrun, lopt_lr, lrhoofr, &
                                                lrnlsm, lupdwf, lv1ofrho1, &
                                                lvhk, nstate

    nstate=crge%n
    CALL give_scr_initrun(linitrun,tag)
    CALL give_scr_updwf(lupdwf,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_forcedr(lforces,tag,nstate,.FALSE.,.FALSE.)
    CALL give_scr_canon(lcanon,tag,nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_v1ofrho1(lv1ofrho1,tag)
    CALL give_scr_opt_lr(lopt_lr,"ORBHARD",tag)
    CALL give_scr_vhk(lvhk,tag)
    IF (lrhd%local_orb)  THEN
       CALL give_scr_ddipo(lddipo,tag)
       lddipo=MAX(2*ncpw%nhg,lddipo)
    ELSE
       lddipo=0
    ENDIF
    ! 
    lohlr=MAX(linitrun,lupdwf,lrhoofr,lforces,lcanon,lrnlsm,&
         lv1ofrho1,lopt_lr,lddipo,lvhk)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ohlr
  ! ==================================================================
  SUBROUTINE rhobas(tag,rhoe,psi)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*)                         :: tag
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhobas'

    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rho1g(:)
    INTEGER                                  :: i, ia, iaorb, iat, ierr, ig, &
                                                info, ir, is, natst, nnrs
    INTEGER, ALLOCATABLE, SAVE               :: ipiv(:)
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: foc(loadc_foc_array_size), &
                                                rhoval, sdtot
    REAL(real_8), ALLOCATABLE, SAVE          :: fxmat(:,:), sdc(:)

    IF (INDEX(tag,"INIT").NE.0) THEN
       CALL rtog(1)
       ALLOCATE(xsmat(atwp%nattot,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fxmat(atwp%nattot+1,atwp%nattot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(catom(ncpw%nhg,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rho1g(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sdc(atwp%nattot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipiv(atwp%nattot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ..Load AO in PW basis
       iaorb=1
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             CALL loadc(catom(1,iaorb),foc,ncpw%nhg,ncpw%nhg,atwp%nattot,SIZE(foc),&
                  is,iat,natst)
             iaorb=iaorb+natst
          ENDDO
       ENDDO
       ! ..Overlap matrix
       CALL dsyrk('U','T',atwp%nattot,2*ncpw%nhg,2._real_8,catom,2*ncpw%nhg,&
            0._real_8,xsmat,atwp%nattot)
       CALL dmatc('U',atwp%nattot,xsmat,atwp%nattot)
       IF (geq0)&
            CALL dger(atwp%nattot,atwp%nattot,-1.0_real_8,catom,2*ncpw%nhg,catom,2*ncpw%nhg,&
            xsmat,atwp%nattot)
       CALL mp_sum(xsmat,atwp%nattot*atwp%nattot,parai%allgrp)
       IF (paral%io_parent)&
            CALL fileopen(22,'DENRHO',fo_def,ferror)
    ELSEIF (INDEX(tag,"EXPAND").NE.0) THEN
       CALL zeroing(fxmat)!,(atwp%nattot+1)*(atwp%nattot+1))
       CALL zeroing(sdc)!,atwp%nattot+1)
       DO i=1,atwp%nattot
          CALL dcopy(atwp%nattot,xsmat(1,i),1,fxmat(1,i),1)
          IF (geq0) sdc(i)=SQRT(parm%omega)*REAL(catom(1,i))
       ENDDO
       CALL mp_sum(sdc,atwp%nattot,parai%allgrp)
       sdtot=0._real_8
       DO i=1,atwp%nattot
          fxmat(i,atwp%nattot+1)=sdc(i)
          fxmat(atwp%nattot+1,i)=sdc(i)
          sdtot=sdtot+ABS(sdc(i))
       ENDDO
       CALL zeroing(psi)!,maxfft)
       rhoval=0._real_8
       DO ir=1,fpar%nnr1
          psi(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          rhoval=rhoval+rhoe(ir,1)
       ENDDO
       CALL mp_sum(rhoval,parai%allgrp)
       nnrs=spar%nr1s*spar%nr2s*spar%nr3s
       rhoval=rhoval*parm%omega/REAL(nnrs,kind=real_8)
       CALL  fwfftn(psi,.FALSE.,parai%allgrp)
       DO ig=1,ncpw%nhg
          rho1g(ig) = psi(nzh(ig))
       ENDDO
       DO i=1,atwp%nattot
          sdc(i)=SQRT(parm%omega)*dotp(ncpw%nhg,catom(:,i),rho1g)
       ENDDO
       sdc(atwp%nattot+1)=rhoval
       CALL mp_sum(sdc,atwp%nattot,parai%allgrp)
       IF (paral%parent) THEN
          IF (sdtot.LT.1.e-10_real_8) THEN
             IF (rhoval.GT.1.e-6_real_8) THEN
                IF (paral%io_parent)&
                     WRITE(6,*)&
                     ' WARNING| BASIS FN INTEGRAL IS ZERO, CHARGE NOT'
                CALL zeroing(sdc)!,atwp%nattot+1)
             ELSE
                CALL dgesv(atwp%nattot,1,fxmat,atwp%nattot+1,ipiv,sdc,atwp%nattot+1,info)
                IF (info.NE.0) CALL stopgm('RHOBAS','INFO 1 /= DGESV',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSE
             CALL dgesv(atwp%nattot+1,1,fxmat,atwp%nattot+1,ipiv,sdc,atwp%nattot+1,info)
             IF (info.NE.0) CALL stopgm('RHOBAS','INFO 2 /= DGESV',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (paral%io_parent)&
               WRITE(22,*) (sdc(i),i=1,atwp%nattot)
       ENDIF
    ELSEIF (INDEX(tag,"FREE").NE.0) THEN
       DEALLOCATE(xsmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fxmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(catom,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rho1g,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(sdc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ipiv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            CALL fileclose(22)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhobas

END MODULE ohlr_utils

! ==================================================================
SUBROUTINE tsumhk(fhk,rhoe,tvec)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE mp_interface, ONLY: mp_sum
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,fpar,parm
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  REAL(real_8)                               :: fhk(fpar%nnr1,*), &
                                                rhoe(fpar%nnr1,*), tvec

  INTEGER                                    :: ir

! ==--------------------------------------------------------------==

  tvec=0._real_8
  DO ir=1,fpar%nnr1
     tvec=tvec+fhk(ir,1)*rhoe(ir,1)
  ENDDO
  IF (cntl%tlsd) THEN
     DO ir=1,fpar%nnr1
        tvec=tvec+fhk(ir,2)*rhoe(ir,2)
     ENDDO
  ENDIF
  CALL mp_sum(tvec,parai%allgrp)
  tvec=tvec*parm%omega/REAL(fpar%nnr1,kind=real_8)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE tsumhk
! ==================================================================
