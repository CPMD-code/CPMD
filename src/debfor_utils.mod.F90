MODULE debfor_utils
  USE adat,                            ONLY: elem
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE calc_alm_utils,                  ONLY: calc_alm,&
                                             give_scr_calc_alm
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cotr,                            ONLY: lskcor
  USE detdof_utils,                    ONLY: detdof
  USE dynit_utils,                     ONLY: dynit
  USE ehpsi_utils,                     ONLY: set_b2l
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: coord_fdiff,&
                                             ions0,&
                                             ions1,&
                                             r_fdiff,&
                                             tref_fdiff
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01,&
                                             tshl
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft,&
                                             lr_tddft
  USE machine,                         ONLY: m_walltime
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE ortho_utils,                     ONLY: give_scr_ortho,&
                                             ortho
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: infi,&
                                             infw,&
                                             iteropt,&
                                             ropt_mod
  USE rrane_utils,                     ONLY: rrane
  USE setirec_utils,                   ONLY: write_irec
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: restart1
  USE symm,                            ONLY: isymu,&
                                             symmi
  USE symtrz_utils,                    ONLY: give_scr_symmat,&
                                             symvec
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, parm
  USE testex_utils,                    ONLY: testex
  USE updrho_utils,                    ONLY: give_scr_updrho,&
                                             updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: debfor
  PUBLIC :: give_scr_debfor

CONTAINS

  ! ==================================================================
  SUBROUTINE debfor(c0,c1,c2,cg0,sc0,pme,gde,vpp,eigv)
    ! ==--------------------------------------------------------------==
    ! == DEBUG FORCES                                                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,1), c1(*), c2(nkpt%ngwk,crge%n), &
      cg0(*), sc0(nkpt%ngwk,crge%n), pme(*), gde(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'debfor'

    CHARACTER(len=1), DIMENSION(1:3)         :: cdir = (/'X','Y','Z'/)
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: ia, iat, id, idelt, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, ir, irec(100), is, isymuat, k, nhpsi, nmm, nnx, nstate, nsym
    LOGICAL                                  :: tinfo
    REAL(real_8) :: dist, ea, ekin1, ekin2, ekincp, ekinh1, ekinh2, er, &
      etoto, tcpu, temp1, temp2, thl(2), time1, time2, xlm, ylm, zlm
    REAL(real_8), ALLOCATABLE                :: fdis(:,:,:,:), rhoe(:,:), &
                                                rinb(:,:), tscr(:,:,:)

    time1 =m_walltime()
    nstate=crge%n
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
    ALLOCATE(fdis(3,maxsys%nax,maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
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
    ENDIF
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    nacc = 7
    ener_com%ecnstr = 0.0_real_8
    ener_com%erestr = 0.0_real_8
    ! ==--------------------------------------------------------------==
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d, &
         & il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
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
    il_psi_1d = il_psi_1d * nmm
    ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! ==--------------------------------------------------------------==
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL write_irec(irec)
    ! END INITIALIZATION
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull)) THEN
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDIF
    ! ..
    etoto=0._real_8
    ener_com%ebogo=0._real_8
    IF (paral%parent) CALL detdof(tau0,tscr)
    cnti%npara=-1
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==      OPTIMIZATION OF WAVEFUNCTION                            ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/," ",20("*"),A,21("*"),/)') ' OPTIMIZE WAVEFUNCTION '
    ropt_mod%convwf=.FALSE.
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    time1=m_walltime()
    DO infw=1,cnti%nomore
       IF (cntl%tdiag) THEN
          CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,tscr,eigv,&
               rhoe,psi,&
               nstate,.FALSE.,.TRUE.,.FALSE.,infw,thl,nhpsi)
          ! ..store initial density
          CALL dcopy(nnx,rin0,1,rinb,1)
       ELSE
          ! UPDATE THE WAVEFUNCTIONS
          CALL updwf(c0(:,:,1),c2,sc0,tau0,tscr,pme,gde,vpp,eigv,&
               rhoe,psi,nstate,.FALSE.,.TRUE.)
       ENDIF
       ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%parent) THEN
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,nstate,ener_com%etot,etoto,&
               tcpu,gemax,cnorm,thl,ropt_mod%engpri,infw,infw)
          etoto=ener_com%etot
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf.OR.soft_com%exsoft) GOTO 100
    ENDDO
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' NO CONVERGENCE IN ',infw,' STEPS '
    ENDIF
    CALL zhwwf(2,irec,c0,c0,nstate,eigv,tau0,tau0,tau0,iteropt%nfi)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! NUCLEAR GRADIENT
    IF (cntl%tddft) THEN
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.TRUE.,.TRUE.)
       CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,.TRUE.,td01%ioutput)
    ELSEIF (cntl%tdiag) THEN
       CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
            rhoe,psi,&
            nstate,.TRUE.,.TRUE.,cntl%tpres,infi,thl,nhpsi)
    ELSE
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.TRUE.,.TRUE.)
    ENDIF
    CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)
    time2=m_walltime()
    tcpu=(time2-time1)*0.001_real_8
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.3,A)') ' TIME FOR ANALYTIC GRADIENT :',&
            tcpu,' SECONDS'
       CALL wrener
       CALL wrgeof(tau0,fion)
    ENDIF
    CALL zhwwf(2,irec,c0,c0,nstate,eigv,tau0,tau0,tau0,iteropt%nfi)
    CALL dcopy(2*nkpt%ngwk*nstate,c0,1,cg0(1),1)
    CALL dcopy(3*maxsys%nax*ions1%nsp,tau0(1,1,1),1,taup(1,1,1),1)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",22("*"),A,22("*"),/)') ' FINITE DIFFERENCES '
       IF (tref_fdiff.AND.paral%io_parent) THEN
          WRITE(6,'(A,T30,3(F12.3))') 'REFERENCE:',coord_fdifF
          WRITE(6,'(A,T54,F12.3)') 'RADIUS OF THE SPHERE:',r_fdifF
          WRITE(6,'(A)')&
               'WE COMPUTE ONLY THE FINITE DIFFERENCES FOR ATOMS'
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
                     WRITE(6,'(" **** ATOM=",I7,4X,A4,4X,A4,4X,A,T54,1PE12.5)')&
                     iat,elem%el(ions0%iatyp(is)),cdir(k),'DISPLACEMENT=',cntr%fdiff*idelt
                ! NEW GEOMETRY
                isymuat=0
                IF (symmi%indpg.NE.0) THEN
                   IF (isymu(iat).EQ.0) THEN
                      isymuat=1
                   ENDIF
                ENDIF
                IF (isymuat.EQ.0.AND.lskcor(k,iat).NE.0) THEN
                   CALL dcopy(3*maxsys%nax*ions1%nsp,taup(1,1,1),1,tau0(1,1,1),1)
                   tau0(k,ia,is)=tau0(k,ia,is)+idelt*cntr%fdiff
                   CALL dcopy(2*nkpt%ngwk*nstate,cg0(1),1,c0,1)
                   IF (cntl%trane) THEN
                      CALL rrane(c0,c2,nstate)
                   ENDIF
                   CALL phfac(tau0)
                   IF (corel%tinlc) THEN
                      CALL copot(rhoe,psi,.FALSE.)
                   ENDIF
                   IF (pslo_com%tivan) THEN
                      CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
                   ENDIF
                   CALL ortho(nstate,c0(:,:,1),c2)
                   ropt_mod%convwf=.FALSE.
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
                         CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,tscr,eigv,&
                              rhoe,psi,&
                              nstate,.FALSE.,tinfo,.FALSE.,&
                              infw,thl, nhpsi)
                      ELSE
                         ! UPDATE THE WAVEFUNCTIONS
                         CALL updwf(c0(:,:,1),c2,sc0,tau0,tscr,pme,gde,vpp,eigv,&
                              rhoe,psi,nstate,.FALSE.,.TRUE.)
                      ENDIF
                      etoto=ener_com%etot
                      IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
                      IF (soft_com%exsoft) GOTO 200
                      IF (ropt_mod%convwf) GOTO 101
                   ENDDO
101                CONTINUE
                   ! cntl%tddft
                   IF (cntl%tddft) THEN
                      CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,tau0,tscr,eigv,&
                           nstate,.FALSE.,td01%ioutput)
                   ENDIF
                ENDIF
                id=(idelt+1)/2 + 1
                fdis(k,ia,is,id)=ener_com%etot
                IF (paral%parent) THEN
                   time2=m_walltime()
                   tcpu=(time2-time1)*0.001_real_8
                   IF (paral%io_parent)&
                        WRITE(6,'(A,I6,T15,A,T24,F20.10,T51,A,0PF9.2)')&
                        " ITER=",infw,"ENERGY=",ener_com%etot," TCPU=",tcpu
                ENDIF
                IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
                IF (soft_com%exsoft) GOTO 200
             ENDDO
          ENDDO
31416     CONTINUE
       ENDDO
    ENDDO
    symmi%nrot=nsym
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             fdis(k,ia,is,1)=(fdis(k,ia,is,2)-fdis(k,ia,is,1))&
                  /(2._real_8*cntr%fdiff)
          ENDDO
       ENDDO
    ENDDO
    CALL symvec(fdis(:,:,:,1))
    IF (paral%parent) CALL wrgeof(taup,fdis(:,:,:,1))
    CALL dcopy(3*maxsys%nax*ions1%nsp,fdis,1,fdis(1,1,1,2),1)
    CALL daxpy(3*maxsys%nax*ions1%nsp,-1._real_8,fion,1,fdis(1,1,1,2),1)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,T10,A,T18,A,T30,A,T49,A)') "ATOM","COMP",&
            "ABSOLUTE ERROR","RELATIVE ERROR[%]"
       iat=0
       ir=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ea=0._real_8
             DO k=1,3
                fdis(k,ia,is,2)=ABS(fdis(k,ia,is,2))
                IF (ea.LT.fdis(k,ia,is,2)) THEN
                   ir=k
                   ea=fdis(k,ia,is,2)
                   er=ABS(ea/(fion(k,ia,is)+1.e-20_real_8))*100._real_8
                ENDIF
             ENDDO
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(I5,T10,A,T18,I4,T30,F14.8,T52,F14.6)')&
                  iat,elem%el(ions0%iatyp(is)),ir,ea,er
          ENDDO
       ENDDO
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
    DEALLOCATE(fdis,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tdiag) THEN
       DEALLOCATE(rinb,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE debfor
  ! ==================================================================
  SUBROUTINE give_scr_debfor(ldebfor,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldebfor
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcalc_alm, lcopot, lforces, &
                                                linitrun, lortho, lrhoofr, &
                                                lrnlsm, lsymmat, ltddft, &
                                                lupdate, nstate

    nstate=crge%n
    lcopot=0
    lortho=0
    lrnlsm=0
    lrhoofr=0
    lcalc_alm=0
    lforces=0
    lsymmat=0
    CALL give_scr_initrun(linitrun,tag)
    IF (restart1%restart)THEN
       IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
       CALL give_scr_ortho(lortho,tag,nstate)
    ENDIF
    IF (cntl%tddft) CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
    IF (cntl%tdiag) THEN
       IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
       CALL give_scr_rhoofr(lrhoofr,tag)
       IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
       CALL give_scr_updrho(lupdate,tag,nstate,.TRUE.,cntl%tpres)
    ELSE
       CALL give_scr_updwf(lupdate,tag,nstate,.FALSE.)
       CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.TRUE.)
    ENDIF
    IF (.NOT.restart1%rvib) CALL give_scr_ortho(lortho,tag,nstate)
    IF (.NOT.(symmi%indpg.EQ.0.OR.symmi%nrot.EQ.1))&
         CALL give_scr_symmat(lsymmat,tag)
    ldebfor=MAX(9*ions1%nat*ions1%nat+9*ions1%nat,&
         linitrun,lcopot,lortho,lrnlsm,lrhoofr,lcalc_alm,&
         lupdate,lforces,lsymmat,ltddft)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_debfor
  ! ==================================================================

END MODULE debfor_utils
