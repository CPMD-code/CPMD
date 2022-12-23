MODULE sh_tddft_utils
  USE adjmu_utils,                     ONLY: adjmu
  USE cnst,                            ONLY: factem,&
                                             pi,&
                                             ry
  USE conv,                            ONLY: nac,&
                                             nbc
  USE coor,                            ONLY: fion,&
                                             tau0
  USE davidson_utils,                  ONLY: davidson
  USE ddip,                            ONLY: lenbk,&
                                             ngwmax
  USE ddipo_utils,                     ONLY: setdip
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileopen
  USE fileopenmod,                     ONLY: fo_app
  USE fint,                            ONLY: fint1
  USE friesner_utils,                  ONLY: friesner
  USE func,                            ONLY: func1
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gsortho
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: &
       c1k, c_sh, ck, nua, shlct, shsigma, shsigmaf, td01, td03, tshf, tshi, tshl, &
       v_sh, xfmqc
  USE machine,                         ONLY: m_flush
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mm_dimmod,                       ONLY: mmdim,&
                                             naq
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE opeigr_p_utils,                  ONLY: opeigr_p
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: supergroup
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE prmem_utils,                     ONLY: prmem
  USE prng_utils,                      ONLY: repprngu
  USE randtowf_utils,                  ONLY: randtowf
  USE readsr_utils,                    ONLY: xstring
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rmas,                            ONLY: rmass
  USE spin,                            ONLY: clsd,&
                                             spin_mod,&
                                             tdsp1
  USE symm,                            ONLY: symmi,&
                                             symmt
  USE system,                          ONLY: &
       cnti, cntl, fpar, mapgp, maxsys, ncpw, nkpt, parap, parm, spar
  USE tauf,                            ONLY: itaur
  USE td_nacvs_utils,                  ONLY: vcouplings
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE utils,                           ONLY: nxxfun
  USE wrener_utils,                    ONLY: wreigen
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lz_trans
  !public :: chk_extrm
  !public :: chk_sh
  !public :: chk_lz
  PUBLIC :: td_spect
  PUBLIC :: resize_f
  PUBLIC :: tmprd
  PUBLIC :: tmpwr
  PUBLIC :: adjustphase
  PUBLIC :: adjustphase_lr
  PUBLIC :: getsigma
  !public :: getoverlap
  !public :: getoverlap_index
  PUBLIC :: update_sh
  PUBLIC :: rscvp_sh
  PUBLIC :: normalize_c_sh
  PUBLIC :: write_shcoeff
  PUBLIC :: read_shcoeff
  PUBLIC :: sh_extpot
  PUBLIC :: sh_lct
  PUBLIC :: ecoupl
  PUBLIC :: xf_qmom
  PUBLIC :: xf_force_history
  PUBLIC :: xf_force_components
  !  PUBLIC :: xf_initialize_all_buffers
  PUBLIC :: xf_fill_all_buffers ! .. MIN changed the name of subroutine
  !public :: td_berry

CONTAINS

  ! ==================================================================
  SUBROUTINE lz_trans(transition,eigv,nstate,nvir,npa,npb,nlr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: transition(3,4,20)
    REAL(real_8)                             :: eigv(*)
    INTEGER                                  :: nstate, nvir, npa, npb, nlr

    CHARACTER(*), PARAMETER                  :: procedureN = 'lz_trans'
    CHARACTER(len=5), DIMENSION(2), &
      PARAMETER                              :: state = (/"ALPHA","BETA "/)
    CHARACTER(len=7), DIMENSION(2), PARAMETER :: &
      stype = (/"HOMO - ","TOP -  "/)
    REAL(real_8), PARAMETER                  :: mine = 0.15_real_8, &
                                                minw = 0.075_real_8

    INTEGER                                  :: fmem, i, ierr, io, is, istat, &
                                                it, iv, ixst, nz
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: extrm_ns, extrm_rs, ostda, &
                                                test_sh = .FALSE.
    REAL(real_8)                             :: maxw
    REAL(real_8), ALLOCATABLE, SAVE          :: delta(:,:), past(:,:)

    IF (ifirst.EQ.0) THEN
       ALLOCATE(past(20,10),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(delta(20,10),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(past)!,200)
       CALL zeroing(delta)!,200)
    ENDIF
    ifirst=ifirst+1
    ! 
    ostda=(td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon)
    IF (paral%parent) THEN
       IF (td03%treorder.OR.td03%tdlocal) THEN
          ixst=2
       ELSE
          ixst=1
       ENDIF
       IF (ostda) THEN
          nz=nua
       ELSE
          nz=nstate
       ENDIF
       ! 
       maxw=0.0_real_8
       fmem=td01%fstate
       ! 
       CALL chk_extrm(td01%fstate,past,delta,eigv,extrm_rs)
       ! 
       DO is=1,nlr
          test_sh=.FALSE.
          IF (is.NE.td01%fstate) CALL chk_extrm(is,past,delta,eigv,extrm_ns)
          IF (ifirst.LT.5) THEN
             test_sh=.FALSE.
          ELSE
             ! CALL CHK_SH(IS,FSTATE,PAST,DELTA,EIGV,MINE,RY,
             ! &                 EXTRM_RS,EXTRM_NS,TEST_SH)
             CALL chk_lz(is,td01%fstate,past,delta,eigv,&
                  extrm_rs,extrm_ns,test_sh)
          ENDIF
          it=0
          DO i=1,4
             IF (transition(1,i,is).EQ.0) GOTO 100
             it=it+1
          ENDDO
100       CONTINUE
          IF (paral%io_parent)&
               WRITE(6,'(/,A,I3,T42,A,F10.3,A)') " STATE=",is,&
               "EIGENVALUE=",2._real_8*ry*eigv(is)," eV"
          ! 
          IF (td03%sh_diag) THEN
             IF (cntl%tlsd) THEN
                DO i=1,it
                   io=transition(1,i,is)
                   IF (io.GT.spin_mod%nsup) THEN
                      istat=2
                      io=nstate-io
                   ELSE
                      istat=1
                      io=spin_mod%nsup-io
                   ENDIF
                   iv=transition(2,i,is)
                   IF (iv.GT.npa) iv=iv-npa
                   IF (paral%io_parent)&
                        WRITE(6,'(T3,A,2X,A,T24,F7.3,T41,A,I3,A,A,I3)')&
                        " TRANSITION ",state(istat),&
                        (REAL(transition(3,i,is))/10000._real_8)**2,&
                        stype(ixst),io," --> ","LUMO + ",iv-1
                ENDDO
             ELSE
                DO i=1,it
                   IF (paral%io_parent)&
                        WRITE(6,'(T10,A,T24,F7.3,T41,A,I3,A,A,I3)')&
                        " TRANSITION ",&
                        (REAL(transition(3,i,is))/10000._real_8)**2,&
                        stype(ixst),nz-transition(1,i,is)," --> ",&
                        "LUMO + ",transition(2,i,is)-1
                ENDDO
                ! 
             ENDIF
          ENDIF
          ! ----------------------------------------------------------
          IF (cntl%tmdbo) THEN
             IF (test_sh) THEN
                td01%fstate=is
                ifirst=1
             ENDIF
          ENDIF
          ! -----------------------------------------------------------
       ENDDO
    ENDIF
    IF ((fmem.NE.td01%fstate).AND.(cntl%tmdbo)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,I4,A,I4,A)')&
            ' ----------  TRANSITION TO NEW STATE: ',&
            fmem,'  TO ', td01%fstate, "   ----------"
    ENDIF
    CALL mp_bcast(td01%fstate,parai%source,parai%allgrp)
    IF ((paral%parent.AND.cntl%tmdbo).AND.paral%io_parent)&
         WRITE(6,'(/,A,I4)')" RUN ON STATE ",td01%fstate
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lz_trans
  ! ==================================================================
  SUBROUTINE chk_extrm(rs,past,delta,eigv,extrm)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: rs
    REAL(real_8)                             :: past(20,*), delta(20,*), &
                                                eigv(*)
    LOGICAL                                  :: extrm

    INTEGER                                  :: i

! Variable
! ==================================================================

    extrm=.FALSE.
    DO i=9,1,-1
       past(rs,i+1)=past(rs,i)
    ENDDO
    past(rs,1)=eigv(rs)
    DO i=1,9
       delta(rs,i)=past(rs,i+1)-past(rs,i)
    ENDDO
    IF ( ((delta(rs,1)*delta(rs,2)).LT.0._real_8) .OR.&
         ((delta(rs,1)*delta(rs,3)).LT.0._real_8) .OR.&
         ((delta(rs,1)*delta(rs,4)).LT.0._real_8) ) extrm=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chk_extrm
  ! ==================================================================
  SUBROUTINE chk_sh(ns,rs,past,delta,eigv,mine,ry,&
       extrm_rs,extrm_ns,test_sh)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns, rs
    REAL(real_8)                             :: past(20,*), delta(20,*), &
                                                eigv(*), mine, ry
    LOGICAL                                  :: extrm_rs, extrm_ns, test_sh

    LOGICAL                                  :: smallgap

    test_sh=.FALSE.
    smallgap=.FALSE.
    IF (ns.EQ.rs) RETURN
    IF (ABS(2.0_real_8*ry*(eigv(rs)-eigv(ns))).LT.mine) smallgap=.TRUE.
    ! IF (EXTRM_RS.AND.EXTRM_NS.AND.SMALLGAP) TEST_SH=.TRUE.
    IF (extrm_rs.AND.smallgap) test_sh=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chk_sh
  ! ==================================================================
  SUBROUTINE chk_lz(ns,rs,past,delta,eigv,extrm_rs,extrm_ns,test_sh)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns, rs
    REAL(real_8)                             :: past(20,*), delta(20,*), &
                                                eigv(*)
    LOGICAL                                  :: extrm_rs, extrm_ns, test_sh

    CHARACTER(len=50)                        :: filen
    INTEGER                                  :: i, i_min_gap
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: extrm_ed, ferror, fexist
    REAL(real_8)                             :: diff_ene(10), massey, &
                                                max_slope, min_gap, prob, &
                                                slope, test

! ==================================================================
! PROB=1-P is the probability to change (adiab) state
! P is def in eq 3 of JACS,120,5499,1998
! ==================================================================

    test_sh=.FALSE.
    extrm_ed=.FALSE.
    IF (ns.EQ.rs) RETURN
    ! IF (.NOT.EXTRM_RS) RETURN
    min_gap=1000.0_real_8
    i_min_gap=10
    DO i=10,1,-1
       diff_ene(i)=ABS(past(rs,i)-past(ns,i))
       ! IF (DIFF_ENE(I).LT.MIN_GAP) THEN
       IF ((diff_ene(i).LT.min_gap).AND.(diff_ene(i).GT.0._real_8)) THEN
          min_gap=diff_ene(i)
          i_min_gap=i
       ENDIF
    ENDDO
    ! the minimum energy difference is between I=2 and I=3
    IF ((i_min_gap.GT.1).AND.(i_min_gap.LE.3)) extrm_ed=.TRUE.
    max_slope=0._real_8
    DO i=1,10-1
       slope=ABS(diff_ene(i)-diff_ene(i+1))/dt_ions
       max_slope=MAX(max_slope,slope)
    ENDDO
    massey=((min_gap/2._real_8)**2)/max_slope
    prob=EXP(-2*pi*massey)
    test=repprngu()
    IF (ifirst.EQ.0) THEN
       filen='LANDAU_ZENER.dat'
       IF (paral%io_parent)&
            INQUIRE(file=filen,exist=fexist)
       IF (paral%io_parent)&
            CALL fileopen(333,filen,fo_app,ferror)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               WRITE(333,'(A50)')&
               '# RS, NS, PROB, TEST MASSEY, MAX_SLOPE, I_MIN_GAP'
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         WRITE(333,'(2I4,4(F15.6),I4)')&
         rs,ns,prob,test,massey,max_slope,i_min_gap
    IF (extrm_ed.AND.(ABS(ns-rs).EQ.1)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,I3,A5,I3,A,F10.6)')&
            ' LZ TRANSITION PROBABILITY BETWEEN STATES ',&
            rs,' AND ',ns, ': ',prob
       ! WRITE(6,'(A,F10.6,A,F10.6)') ' PROBABILITY: ',
       ! &    PROB,' RANDOM NUMBER: ',TEST
    ENDIF
    IF ((test.LT.prob).AND.(extrm_ed).AND.(ABS(ns-rs).EQ.1))&
         test_sh=.TRUE.
    ifirst=1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chk_lz
  ! ==================================================================
  SUBROUTINE td_spect(c0,c2,rhoe,ntddft,cvirt,sc0,psi,edav)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: ntddft
    COMPLEX(real_8)                          :: cvirt(ncpw%ngw,*), sc0(*), &
                                                psi(:)
    REAL(real_8)                             :: edav(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_spect'

    COMPLEX(real_8), ALLOCATABLE             :: cr(:,:), cscr(:,:), cx(:,:)
    INTEGER                                  :: i, ib, ierr, is, nconv, ncr, &
                                                ncscr, ndiag, nhpsi
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: dodavi, dolanc, tadjmu, &
                                                trefine
    REAL(real_8)                             :: wk_fake(1)
    REAL(real_8), ALLOCATABLE                :: vpp(:)

! 

    dolanc=.TRUE.
    dodavi=.NOT.dolanc
    ! 
    ALLOCATE(potr(fpar%nnr1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cnti%nkry_block.LE.0) cnti%nkry_block=ntddft

    CALL m_flush(6)
    IF (ntddft.GT.crge%n.AND.ifirst.EQ.0) THEN
       DO i=crge%n+1,ntddft
          CALL randtowf(c0(:,i:i),1,1,1)
       ENDDO
       ifirst=ifirst+1
    ENDIF
    IF (cntl%tlsd) THEN
       CALL gsortho(c0(1,1),ncpw%ngw,tdsp1%nupel+1,spin_mod%nsup)
       CALL gsortho(c0(1,spin_mod%nsup+1),ncpw%ngw,tdsp1%ndoel+1,spin_mod%nsdown)
    ELSE
       CALL gsortho(c0,ncpw%ngw,crge%n+1,ntddft)
    ENDIF
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
    CALL rhoofr(c0(:,1:ntddft),rhoe,psi,ntddft)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)

    IF (dolanc) THEN
       fint1%tfral=.TRUE.
       cntl%tlanc=.FALSE.
       CALL zeroing(edav(1:ntddft))!,ntddft)
       ncscr=2*ncpw%ngw*MAX(ntddft,cnti%nkry_max*cnti%nkry_block)
       ALLOCATE(cscr(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (func1%mhfx.NE.0) THEN
          ALLOCATE(cx(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(2*ncpw%ngw*ntddft,c0,1,cx,1)
       ENDIF
       ! IF(PARENT) CALL PRMEM('   SPECTRA')
       IF (cntl%tlsd) THEN
          CALL friesner(spin_mod%nsup,REAL(tdsp1%nupel,kind=real_8),nconv,nhpsi,tdsp1%nupel,&
               c0,c2,sc0,cscr,cx,crge%f(1:spin_mod%nsup,1),potr(1,1),psi,edav,&
               0,trefine,.TRUE.)
          nac=nconv
          ib=spin_mod%nsup+1
          itaur=2
          CALL friesner(spin_mod%nsdown,REAL(tdsp1%ndoel,kind=real_8),nconv,nhpsi,tdsp1%ndoel,&
               c0(1,ib),c2,sc0,cscr,cx(1,ib),crge%f(ib:ib+spin_mod%nsdown-1,1),potr(1,2),psi,&
               edav(ib),0,trefine,.TRUE.)
          itaur=1
          nbc=nconv
          ndiag=ntddft
       ELSE
          ndiag=ntddft
          CALL friesner(ndiag,crge%nel,nconv,nhpsi,ntddft,&
               c0,c2,sc0,cscr,cx,crge%f(1:ndiag,1),potr,psi,edav,&
               0,trefine,.TRUE.)
          nac=nconv
       ENDIF
       DEALLOCATE(cscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (func1%mhfx.NE.0) DEALLOCATE(cx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       cnti%ndavv=2*ntddft
       CALL zeroing(edav(1:2*cnti%ndavv))!,2*cnti%ndavv)
       ncscr=2*ncpw%ngw*cnti%ndavv
       ALLOCATE(cscr(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ncr=2*ncpw%ngw*cnti%ndavv
       ALLOCATE(cr(ncpw%ngw,ncr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cx(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL ksdiag(vpp)
       IF (paral%parent) CALL prmem('   SPECTRA')
       IF (cntl%tlsd) THEN
          CALL dcopy(2*ncpw%ngw*spin_mod%nsup,c0,1,cx,1)
          CALL davidson(spin_mod%nsup,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,1),&
               c0(:,1:spin_mod%nsup),crge%f(1:spin_mod%nsup,1),tdsp1%nupel,psi,edav,vpp)
          CALL dcopy(2*ncpw%ngw*spin_mod%nsup,cx,1,c0,1)
          ib = spin_mod%nsup+1
          itaur=2
          CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,c0(1,ib),1,cx,1)
          CALL davidson(spin_mod%nsdown,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,2),&
               c0(:,ib:ib+spin_mod%nsdown-1),crge%f(ib:ib+spin_mod%nsdown-1,1),tdsp1%ndoel,psi,edav(cnti%ndavv+1),vpp)
          itaur=1
          CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,cx,1,c0(1,ib),1)
          DO is=1,spin_mod%nsdown
             edav(spin_mod%nsup+is)=edav(cnti%ndavv+is)
          ENDDO
       ELSE
          ndiag=ntddft
          CALL dcopy(2*ncpw%ngw*ndiag,c0,1,cx,1)
          CALL davidson(ndiag,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,1),&
               c0(:,1:ndiag),crge%f(:,1),crge%n,psi,edav,vpp)
          CALL dcopy(2*ncpw%ngw*ndiag,cx,1,c0,1)
       ENDIF
       DEALLOCATE(cx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) THEN
       tadjmu=.TRUE.
       CALL adjmu(ntddft,1,crge%nel,fint1%betael,edav,wk_fake,cntl%tlsd,ener_com%amu,tadjmu)
       CALL wreigen(edav,crge%f,ener_com%amu,ntddft)
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    ! 
    RETURN
  END SUBROUTINE td_spect
  ! ==================================================================
  SUBROUTINE resize_f
    ! ==--------------------------------------------------------------==
    ! 
    CHARACTER(*), PARAMETER                  :: procedureN = 'resize_f'

    EXTERNAL                                 :: dasum
    INTEGER                                  :: i, ierr, nroot, ntddft, nx
    REAL(real_8)                             :: dasum
    REAL(real_8), ALLOCATABLE                :: fo(:)

    nroot=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri,2)
    ntddft=crge%n+nroot
    ALLOCATE(fo(ntddft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL dcopy(crge%n,crge%f,1,fo,1)
    DO i=crge%n+1,ntddft
       fo(i)=0.0_real_8
    ENDDO
    DEALLOCATE(crge%f,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(crge%f(ntddft,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(crge%f)!,ntddft)
    IF (cntl%tlsd) THEN
       tdsp1%nupel=NINT(dasum(spin_mod%nsup,fo,1))
       tdsp1%ndoel=NINT(dasum(spin_mod%nsdown,fo(spin_mod%nsup+1),1))
       nx=ntddft-spin_mod%nsup-spin_mod%nsdown
       spin_mod%nsup=spin_mod%nsup+(nx+1)/2
       spin_mod%nsdown=spin_mod%nsdown+nx/2
       DO i=1,tdsp1%nupel
          crge%f(i,1)=1._real_8
       ENDDO
       DO i=spin_mod%nsup+1,spin_mod%nsup+tdsp1%ndoel
          crge%f(i,1)=1._real_8
       ENDDO
    ELSE
       CALL dcopy(crge%n,fo,1,crge%f,1)
    ENDIF
    DEALLOCATE(fo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! CALL DCOPY(NTDDFT,F,1,FO,1)

    RETURN
  END SUBROUTINE resize_f
  ! 
  ! ==================================================================
  SUBROUTINE tmprd(phi,nstate,nroot,itime,tag)
    ! ==--------------------------------------------------------------==

    ! -- Arguments --
    INTEGER                                  :: nstate, nroot
    COMPLEX(real_8)                          :: phi(ncpw%ngw,nstate,nroot)
    INTEGER                                  :: itime
    CHARACTER(len=30)                        :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'tmprd'

    INTEGER                                  :: ierr, ig, ip, ipp, ipw, &
                                                iroot, len, lst, msgid
    COMPLEX(real_8), ALLOCATABLE             :: phi2(:), phisx(:)
    INTEGER, ALLOCATABLE                     :: mapw(:)

! -- Variables --

    IF (paral%parent) THEN
       len=2*spar%ngws
    ELSE
       len=2*ncpw%ngw
    ENDIF
    ALLOCATE(phisx(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       ALLOCATE(phi2(spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(mapw(2*(ncpw%ngw+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            OPEN(unit=120,file=tag,status='OLD',&
            form='UNFORMATTED')
       IF (paral%io_parent)&
            READ(120) itime
    ENDIF
    ! MSGLEN=1 * 8/IRAT
    ! CALL MY_BCAST(ITIME,MSGLEN,SOURCE,ALLGRP)

    DO iroot=1,nroot
       DO lst=1,nstate
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  READ(120) (phisx(ig),ig=1,spar%ngws)
          ENDIF

          DO ipp=1,parai%nproc
             ip=parap%pgroup(ipp)
             msgid=ip
             IF (paral%parent) THEN
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      mapw(ig)=mapgp(ig)
                   ENDDO
                ELSE
                   msgid=2
                   !msglen=parap%sparm(3,ipp-1)*8/irat
                   CALL mp_recv(mapw,parap%sparm(3,ipp-1),ip,msgid,parai%allgrp)
                ENDIF
                DO ipw=1,parap%sparm(3,ipp-1)
                   phi2(ipw)=phisx(mapw(ipw))
                ENDDO
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      phi(ig,lst,iroot)=phi2(ig)
                   ENDDO
                ELSE
                   msgid=1
                   !msglen=2*parap%sparm(3,ipp-1)*8
                   CALL mp_send(phi2,parap%sparm(3,ipp-1),ip,msgid,parai%allgrp)
                ENDIF
             ELSE
                IF (ip.EQ.parai%me) THEN
                   msgid=2
                   !msglen=ngw*8/irat
                   CALL mp_send(mapgp,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   msgid=1
                   !msglen=2*ngw*8
                   CALL mp_recv(phisx,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   DO ig=1,ncpw%ngw
                      phi(ig,lst,iroot)=phisx(ig)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF ((paral%parent).AND.paral%io_parent)&
         CLOSE(120)

    IF (paral%parent) THEN
       DEALLOCATE(mapw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(phi2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(phisx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tmprd
  ! ==================================================================
  SUBROUTINE tmpwr(phi,nstate,nroot,itime,tag)
    ! ==--------------------------------------------------------------==

    ! -- Arguments --
    INTEGER                                  :: nstate, nroot
    COMPLEX(real_8)                          :: phi(ncpw%ngw,nstate,nroot)
    INTEGER                                  :: itime
    CHARACTER(len=30)                        :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'tmpwr'

    COMPLEX(real_8), ALLOCATABLE             :: phi2(:), phisx(:)
    INTEGER                                  :: ierr, ig, ip, ipp, ipw, &
                                                iroot, len_g, lst, msgid
    INTEGER, ALLOCATABLE                     :: mapw(:)

! -- Variables --

    IF (paral%parent) THEN
       len_g=2*spar%ngws
    ELSE
       len_g=2*ncpw%ngw
    ENDIF
    ALLOCATE(phisx(len_g),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       ALLOCATE(phi2(spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(mapw(2*(ncpw%ngw+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            OPEN(unit=120,file=tag,status='UNKNOWN',form='UNFORMATTED')
       IF (paral%io_parent)&
            WRITE(120) itime
    ENDIF

    IF (parai%cp_nproc>1)THEN

       DO iroot=1,nroot
          DO lst=1,nstate
             DO ipp=1,parai%nproc
                ip=parap%pgroup(ipp)
                msgid=ip
                IF (paral%parent) THEN
                   IF (ip.EQ.parai%me) THEN
                      DO ig=1,ncpw%ngw
                         phi2(ig)=phi(ig,lst,iroot)
                         mapw(ig)=mapgp(ig)
                      ENDDO
                   ELSE
                      msgid=1
                      !msglen=2*parap%sparm(3,ip)*8
                      CALL mp_recv(phi2,parap%sparm(3,ip),ip,msgid,parai%allgrp)
                      msgid=2
                      !msglen=parap%sparm(3,ip)*8/irat
                      CALL mp_recv(mapw,parap%sparm(3,ip),ip,msgid,parai%allgrp)
                   ENDIF
                   DO ipw=1,parap%sparm(3,ipp-1)
                      phisx(mapw(ipw))=phi2(ipw)
                   ENDDO
                ELSE
                   IF (ip.EQ.parai%me) THEN
                      DO ig=1,ncpw%ngw
                         phisx(ig)=phi(ig,lst,iroot)
                      ENDDO
                      msgid=1
                      !msglen=2*ngw*8
                      CALL mp_send(phisx,ncpw%ngw,parai%source,msgid,parai%allgrp)
                      msgid=2
                      !msglen=ngw*8/irat
                      CALL mp_send(mapgp,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   ENDIF
                ENDIF
             ENDDO
             CALL mp_sync(parai%allgrp)
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE (120) (phisx(ig),ig=1,spar%ngws)
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO iroot=1,nroot
          DO lst=1,nstate
             WRITE (120) (phi(ig,lst,iroot),ig=1,nkpt%ngwk)
          ENDDO
       ENDDO
    ENDIF
    IF ((paral%parent).AND.paral%io_parent)&
         CLOSE(120)

    IF (paral%parent) THEN
       DEALLOCATE(phi2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(mapw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(phisx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tmpwr
  ! ==================================================================
  SUBROUTINE adjustphase(c0,ck,nn,ovlapm)
    ! ==================================================================
    ! 
    INTEGER                                  :: nn
    COMPLEX(real_8)                          :: ck(ncpw%ngw,nn), &
                                                c0(ncpw%ngw,nn)
    REAL(real_8)                             :: ovlapm(nn,nn)

    INTEGER                                  :: i

! 

    CALL zeroing(ovlapm)!,nn*nn)
    CALL ovlap(nn,ovlapm,c0,ck)
    CALL mp_sum(ovlapm,nn*nn,parai%allgrp)
    DO i=1,nn
       ! if (parent) write(6,*) 'ovlap:',I,OVLAPM(I,I)
       IF (ovlapm(i,i).LT.0.0) THEN
          CALL dscal(2*ncpw%ngw,-1.0_real_8,c0(1,i),1)
       ENDIF
    ENDDO
    CALL ovlap(nn,ovlapm,c0,ck)
    CALL mp_sum(ovlapm,nn*nn,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE adjustphase
  ! ==================================================================
  SUBROUTINE adjustphase_lr(c1,c1k,nn,nroot,ovlapm)
    ! ==================================================================
    ! 
    INTEGER                                  :: nn
    COMPLEX(real_8)                          :: c1k(ncpw%ngw,nn,*), &
                                                c1(ncpw%ngw,nn,*)
    INTEGER                                  :: nroot
    REAL(real_8)                             :: ovlapm(nn,nn)

    INTEGER                                  :: i, iroot

! 

    DO iroot=1,nroot
       CALL zeroing(ovlapm)!,nn*nn)
       CALL ovlap(nn,ovlapm,c1(:,:,iroot),c1k(:,:,iroot))
       CALL mp_sum(ovlapm,nn*nn,parai%allgrp)
       DO i=1,nn
          ! if (parent) write(6,*) 'ovlap_lr:',I,OVLAPM(I,I)
          IF (ovlapm(i,i).LT.0.0) THEN
             CALL dscal(2*ncpw%ngw,-1.0_real_8,c1(1,i,iroot),1)
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE adjustphase_lr
  ! ==================================================================
  SUBROUTINE getsigma(c0,c1,nstate,nroot)
    ! ==================================================================

    ! 
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*)
    INTEGER                                  :: nroot

    CHARACTER(*), PARAMETER                  :: procedureN = 'getsigma'

    COMPLEX(real_8), ALLOCATABLE             :: cx1(:,:), cx2(:,:), cy(:,:)
    INTEGER                                  :: i, id, ierr, ii, j, jd, jj, &
                                                k1, k2
    INTEGER, ALLOCATABLE                     :: x1(:), x2(:), y1(:), y2(:)
    REAL(real_8)                             :: sp1, sp2, spu, spx, spy
    REAL(real_8), ALLOCATABLE                :: scalmat(:,:,:,:)

! Variables
! 
! 

    ALLOCATE(cx1(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cx2(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cy(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(scalmat((nroot+1),(nroot+1),(nstate),(nstate)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(x1(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(x2(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(y1(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(y2(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    ! compute all scalar products
    ! i j  indeces of SCALMAT(i,j, ..)
    ! --------------------------------------------------------
    ! 11          12           13          14          ...
    ! [C0|CK]      [C0|C1K1]     [C0|C1K2]    [C0|C1K3]    ...
    ! 21          22           23          24          ...
    ! [C1|CK]      [C1|C1K1]     [C1|C1K2]    [C1|C1K3]    ...
    ! ...          ...           ...          ...
    ! ...          ...           ...          ...
    ! 
    ! where C0 and CK are the ground state orbitals at the two 
    ! different times (t and t+dt, respectively).
    ! C1  are the linear response orbitals for S1 at time t
    ! CK1 are the linear response orbitals for S1 at time t+dt
    ! 
    ! IF (parent) write(6,*) 'start calculation overlap matrices ...'
    ! call m_flush(6)

    DO ii=1,nstate
       DO jj=1,nstate
          scalmat(1,1,ii,jj)=dotp(ncpw%ngw,c0(:,ii),ck(:,jj))
       ENDDO
    ENDDO
    DO i=1,nroot
       id=i+1
       DO ii=1,nstate
          DO jj=1,nstate
             scalmat(id,1,ii,jj)=dotp(ncpw%ngw,c1(:,ii,i),ck(:,jj))
          ENDDO
       ENDDO
    ENDDO
    DO j=1,nroot
       jd=j+1
       DO ii=1,nstate
          DO jj=1,nstate
             scalmat(1,jd,ii,jj)=dotp(ncpw%ngw,c0(:,ii),c1k(:,jj,j))
          ENDDO
       ENDDO
    ENDDO
    DO i=1,nroot
       id=i+1
       DO j=1,nroot
          jd=j+1
          DO ii=1,nstate
             DO jj=1,nstate
                scalmat(id,jd,ii,jj)=dotp(ncpw%ngw,c1(:,ii,i),c1k(:,jj,j))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! 
    ! IF (parent) write(6,*) '... end calculation matrices'
    ! call m_flush(6)
    ! 
    ! 
    ! Scalar product of the unexcited spin manifold
    IF (cntl%tlsd) THEN
       ! not yet implemented
       CALL stopgm('GETSIGMA','LSD NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSEIF (td01%ns_sin.NE.0) THEN
       IF (MOD(NINT(crge%nel),2).EQ.0) THEN
          CALL getoverlap(spu,c0,ck,nstate)
       ELSE
          CALL getoverlap(spu,c0,ck,nstate)
       ENDIF
    ELSEIF (td01%ns_tri.NE.0) THEN
       CALL getoverlap(spu,c0,ck,nstate)
    ENDIF
    ! Scalar product of the excited spin manifold
    DO i=0,nroot
       id=i+1
       DO j=i+1,nroot
          jd=j+1
          sp1=0._real_8
          sp2=0._real_8
          DO k1=1,nstate
             DO ii=1,nstate
                x1(ii)=1
             ENDDO
             IF (i.NE.0) THEN
                x1(k1)=i+1
             ENDIF
             DO ii=1,nstate
                y1(ii)=1
             ENDDO
             IF (i.NE.0) THEN
                y1(k1)=i+1
             ENDIF
             DO k2=1,nstate
                DO ii=1,nstate
                   y2(ii)=1
                ENDDO
                y2(k2)=j+1
                CALL getoverlap_index(spx,x1,y2,nroot,nstate,scalmat,1)
                sp1=sp1+spx
                ! 
                DO ii=1,nstate
                   x2(ii)=1
                ENDDO
                x2(k2)=j+1
                CALL getoverlap_index(spy,y1,x2,nroot,nstate,scalmat,2)
                sp2=sp2+spy
             ENDDO
          ENDDO
          shsigma(1,id,jd)=(sp1*spu-sp2*spu)/(2._real_8*dt_ions)
          ! if (parent) write(6,'(2I3,3F15.9)') 
          ! &                ID,JD,SHSIGMA(1,ID,JD),SP1,SP2
          shsigma(1,jd,id)=-shsigma(1,id,jd)
       ENDDO
    ENDDO
    ! 
    ! CALL GETSIGMA_OLD(C0,C1,NSTATE,NROOT)

    ! IF (PARENT) THEN
    ! WRITE(6,*) 'SIGAMA VALUES '
    ! DO ID=1,NROOT+1
    ! WRITE(6,'(100F15.7)') (SHSIGMA(1,JD,ID),JD=1,NROOT+1)
    ! ENDDO
    ! ENDIF

    DEALLOCATE(y2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(y1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(x2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(x1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scalmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cx2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cx1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getsigma
  ! ==================================================================
  SUBROUTINE getoverlap(sp,c1,c2,n1)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: sp
    COMPLEX(real_8)                          :: c1(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: n1

    CHARACTER(*), PARAMETER                  :: procedureN = 'getoverlap'

    INTEGER                                  :: i, ierr, info, j, msglen
    INTEGER, ALLOCATABLE                     :: iscr(:)
    REAL(real_8), ALLOCATABLE                :: ovlap(:,:)

    ALLOCATE(ovlap(n1,n1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iscr(20*n1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    IF (cntl%tlsd) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'T-SHTDDFT: LSD NOT YET IMPLEMENTED'
       STOP
    ELSE
       CALL zeroing(ovlap)!,n1*n1)
       DO i=1,n1
          DO j=1,n1
             ovlap(i,j)=dotp(ncpw%ngw,c1(:,i),c2(:,j))
          ENDDO
       ENDDO
       msglen=8
       CALL mp_sum(ovlap,n1*n1,parai%allgrp)
       CALL mp_bcast(ovlap,n1*n1,parai%source,parai%allgrp)
       CALL dgetrf(n1,n1,ovlap,n1,iscr,info)
       IF (info.NE.0) THEN
          ! &  call stopgm('ampli','error in matrix inversion (1)')
          IF (paral%io_parent)&
               WRITE(6,*) 'matrix almost singular'
          ! DO I=1,N1
          ! OVLAP(I,I)=0._real_8
          ! ENDDO
       ENDIF
       IF (paral%parent) THEN
          sp=1._real_8
          DO i=1,n1
             sp=sp*ovlap(i,i)
          ENDDO
       ENDIF
    ENDIF
    ! 
#ifdef __PARALLEL
    CALL mp_bcast(sp,parai%source,parai%allgrp)
#endif
    DEALLOCATE(iscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ovlap,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getoverlap
  ! ==================================================================
  SUBROUTINE getoverlap_index(sp,x,y,nroot,n1,scalmat,transp)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: sp
    INTEGER                                  :: x(*), y(*), nroot, n1
    REAL(real_8)                             :: scalmat(nroot+1,nroot+1,n1,n1)
    INTEGER                                  :: transp

    CHARACTER(*), PARAMETER :: procedureN = 'getoverlap_index'
    REAL(real_8), PARAMETER                  :: tol = 1.e-3_real_8 

    INTEGER                                  :: i, ierr, info, j, msglen
    INTEGER, ALLOCATABLE                     :: iscr(:)
    INTEGER, SAVE                            :: totdet = 0, totzero = 0
    LOGICAL                                  :: setzerodet
    REAL(real_8)                             :: maxcol, maxline
    REAL(real_8), ALLOCATABLE                :: ovlap(:,:)

    ALLOCATE(ovlap(n1,n1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iscr(20*n1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    setzerodet=.FALSE.
    ! IF IT FINDS A COLUM OR LINE WITH ALL ENTRIES < TOL
    ! THE DETERMINANT IS NOT COMPUTED AND IT'S VALUE IS SET 
    ! TO ZERO
    ! 
    IF (cntl%tlsd) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'T-SHTDDFT: LSD NOT YET IMPLEMENTED'
       STOP
    ELSE
       CALL zeroing(ovlap)!,n1*n1)
       IF (transp.EQ.1) THEN
          DO j=1,n1
             maxline=0._real_8
             DO i=1,n1
                ovlap(i,j)=scalmat(x(i),y(j),i,j)
                maxline=MAX(maxline,ABS(ovlap(i,j)))
             ENDDO
             CALL mp_max(maxline,parai%allgrp)
             IF (maxline.LE.tol) setzerodet=.TRUE.
          ENDDO
       ELSE
          DO i=1,n1
             maxline=0._real_8
             DO j=1,n1
                ovlap(j,i)=scalmat(y(i),x(j),i,j)
                maxline=MAX(maxline,ABS(ovlap(j,i)))
             ENDDO
             CALL mp_max(maxline,parai%allgrp)
             IF (maxline.LE.tol) setzerodet=.TRUE.
          ENDDO
       ENDIF
       DO j=1,n1
          maxcol=0._real_8
          DO i=1,n1
             maxcol=MAX(maxcol,ABS(ovlap(i,j)))
          ENDDO
          CALL mp_max(maxcol,parai%allgrp)
          IF (maxcol.LE.tol) setzerodet=.TRUE.
       ENDDO

       ! SETZERODET=.FALSE.

       IF (setzerodet) THEN
          sp=0._real_8
          totzero=totzero+1
       ELSE
          totdet=totdet+1
          msglen=8
          CALL mp_sum(ovlap,n1*n1,parai%allgrp)
          CALL mp_bcast(ovlap,n1*n1,parai%source,parai%allgrp)
          CALL dgetrf(n1,n1,ovlap,n1,iscr,info)
          IF (info.NE.0) THEN
             ! &    call stopgm('ampli','error in matrix inversion (1)')
             IF (paral%io_parent)&
                  WRITE(6,*) 'matrix almost singular'
             ! DO I=1,N1
             ! OVLAP(I,I)=0._real_8
             ! ENDDO
          ENDIF
          IF (paral%parent) THEN
             sp=1._real_8
             DO i=1,n1
                sp=sp*ovlap(i,i)
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    ! 
#ifdef __PARALLEL
    CALL mp_bcast(sp,parai%source,parai%allgrp)
#endif
    DEALLOCATE(iscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ovlap,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getoverlap_index
  ! ==================================================================
  SUBROUTINE update_sh(nroot)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nroot

    INTEGER                                  :: i, j

    DO i=1,nroot+1
       c_sh(2,i)=c_sh(1,i)
       v_sh(2,i)=v_sh(1,i)
       DO j=1,nroot+1
          shsigma(2,i,j)=shsigma(1,i,j)
          shsigmaf(2,i,j)=shsigmaf(1,i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE update_sh
  ! ==================================================================
  SUBROUTINE rscvp_sh(dt_sh,tempp,velp)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: dt_sh, tempp, velp(:,:,:)

    INTEGER                                  :: ia, is, k
    REAL(real_8)                             :: alfap, targ_temp

! Variables
! ==--------------------------------------------------------------==
! ==  Dynamical rescaling factor (targ_temp/tempp), due to        ==
! ==  surface hop, where tempp is calculated every step           ==
! ==--------------------------------------------------------------==

    targ_temp=tempp+dt_sh
    alfap=SQRT(targ_temp/tempp)
    IF (lqmmm%qmmm) THEN
       DO is=1,mmdim%nspq
          DO ia=1,NAq(is)
             velp(1,ia,is)=alfap*velp(1,ia,is)
             velp(2,ia,is)=alfap*velp(2,ia,is)
             velp(3,ia,is)=alfap*velp(3,ia,is)
          ENDDO
       ENDDO
    ELSE
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             velp(1,ia,is)=alfap*velp(1,ia,is)
             velp(2,ia,is)=alfap*velp(2,ia,is)
             velp(3,ia,is)=alfap*velp(3,ia,is)
          ENDDO
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*) 'VELOCITIES ARE RESCALED BY',dt_sh,'K'
    IF (paral%io_parent)&
         WRITE(6,*) dt_sh/factem, 'HARTREE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rscvp_sh
  ! ==================================================================
  SUBROUTINE normalize_c_sh(nroot)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nroot

    INTEGER                                  :: i
    LOGICAL                                  :: debug
    REAL(real_8)                             :: sum2

    debug=.FALSE.
    sum2=0._real_8
    DO i=1,nroot+1
       sum2=sum2+REAL(c_sh(1,i)*CONJG(c_sh(1,i)))
       IF (paral%parent.AND.debug) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'c_sh(1,', i,')=',c_sh(1,i)
          IF (paral%io_parent)&
               WRITE(6,*) '|C_SH|^2', i, REAL(c_sh(1,i)*CONJG(c_sh(1,i)))
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) 'SUM OF ALL (C_SH)^2', sum2
    DO i=1,nroot+1
       c_sh(1,i)=c_sh(1,i)/SQRT(sum2)
    ENDDO

    ! SUM2=0._real_8
    ! DO I=1,NROOT+1
    ! SUM2=SUM2+REAL(C_SH(1,I)*CONJG(C_SH(1,I)))
    ! ENDDO
    ! IF (PARENT) WRITE(6,*) 'SUM OF ALL (C_SH)^2 after normalization',
    ! &  SUM2

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE normalize_c_sh
  ! ==================================================================
  SUBROUTINE write_shcoeff(shprob,etot,nroot)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: etot
    INTEGER                                  :: nroot
    REAL(real_8)                             :: shprob(nroot+1,nroot+1)

    CHARACTER(len=100) :: filen, filen_a_sh, filen_c_sh, filen_e_sh, &
      filen_f_sh, filen_f_xf, filen_fbo_xf, filen_fnacv_xf, filen_fqmom_xf, &
      filen_p_sh, filen_s_sh, filen_v_sh
    CHARACTER(len=12)                        :: cflbod
    INTEGER                                  :: i, i1, i2, i3, icount, ir, &
                                                is, iv, j
    LOGICAL                                  :: fexist
    REAL(real_8)                             :: temp(10000),temp2(10000)

! Variables
! 

    IF (paral%parent) THEN
       ! 
       is=0
       DO ir=1,nroot+1
          DO iv=1,nroot+1
             is=is+1
             temp(is)=shsigma(1,ir,iv)
             temp2(is)=shsigmaf(1,ir,iv)
             IF (ir.EQ.iv) temp(is)=0.0_real_8
          ENDDO
       ENDDO
       !
       IF (tshl%txfmqc.AND.tmw) THEN  ! MIN: change the filenames depending on the walker_id (just added the comment)
          cflbod='SH_COEFC_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_c_sh=filen(i1:i2)//'.dat'

          cflbod='SH_COEFA_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_a_sh=filen(i1:i2)//'.dat'

          cflbod='SH_VPOTS_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_v_sh=filen(i1:i2)//'.dat'

          cflbod='SH_COUPL_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_s_sh=filen(i1:i2)//'.dat'

          cflbod='SH_PROBS_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_p_sh=filen(i1:i2)//'.dat'

          cflbod='SH_ENERG_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_e_sh=filen(i1:i2)//'.dat'

          cflbod='SH_STATE_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_f_sh=filen(i1:i2)//'.dat'

          cflbod='XF_ACFOR_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_f_xf=filen(i1:i2)//'.dat'

          cflbod='XF_BOFOR_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_fbo_xf=filen(i1:i2)//'.dat'

          cflbod='XF_NACVFOR_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_fnacv_xf=filen(i1:i2)//'.dat'

          cflbod='XF_QMOMFOR_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_fqmom_xf=filen(i1:i2)//'.dat'
       ELSE
          filen_c_sh="SH_COEFC.dat"
          filen_a_sh="SH_COEFA.dat"
          filen_v_sh="SH_VPOTS.dat"
          filen_s_sh="SH_COUPL.dat"
          filen_p_sh="SH_PROBS.dat"
          filen_e_sh="SH_ENERG.dat"
          filen_f_sh="SH_STATE.dat"
       ENDIF

       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_c_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_c_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_c_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,40F15.9)') tshi%shstep,&
            (REAL(c_sh(1,ir)),AIMAG(c_sh(1,ir)),ir=1,nroot+1)
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_a_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_a_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_a_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,20F15.9)') tshi%shstep,&
            (SQRT(REAL(c_sh(1,ir)*CONJG(c_sh(1,ir)))),ir=1,nroot+1)
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_v_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_v_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_v_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,20F15.9)') tshi%shstep,&
            (v_sh(1,ir),ir=1,nroot+1)
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_s_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_s_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_s_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,1000F15.9)') tshi%shstep,&
            (temp(ir),ir=1,(nroot+1)*(nroot+1))
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_p_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_p_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_p_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,20F15.9)') tshi%shstep,&
            (shprob(ir,td01%fstate+1),ir=1,nroot+1)
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_e_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_e_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_e_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(I10,20F15.9)') tshi%shstep,&
            (etot+v_sh(1,ir),ir=1,nroot+1),etot+v_sh(1,td01%fstate+1)
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_f_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_f_sh,status='NEW')
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_f_sh,status='UNKNOWN',position='APPEND')
       ENDIF
       IF (paral%io_parent)&
            WRITE(50,'(2I10)') tshi%shstep,td01%fstate
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       ! MIN: write f_nli (accumulated BO forces)
       IF (tshl%txfmqc) THEN ! MIN: added IF for the sake of safety
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file=filen_f_xf,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_f_xf,status='NEW')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_f_xf,status='UNKNOWN',position='APPEND')
          ENDIF
          IF (paral%io_parent)&
               WRITE(50,'(10000(E15.9,1X))') xfmqc%f_nli
          IF (paral%io_parent)&
               CLOSE(50)
       ENDIF
       ! FEDE: open files to write the force terms for the xf
       IF (tshl%txfmqc) THEN
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file=filen_fbo_xf,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fbo_xf,status='NEW')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fbo_xf,status='UNKNOWN',position='APPEND')
          ENDIF
          !
          IF (paral%io_parent) THEN
             DO i1=1,tshi%nroot_sh+1 !nstatse
                DO i2=1,3 !x,y,z
                   ! each line contains 6(natoms) terms
                   WRITE(50,'(50(E15.9,1X))') xfmqc%fion_state(i2,:,:,i1)
                ENDDO
             ENDDO
          ENDIF
          IF (paral%io_parent)&
               CLOSE(50)
       ENDIF
       !
       IF (tshl%txfmqc) THEN 
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file=filen_fnacv_xf,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fnacv_xf,status='NEW')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fnacv_xf,status='UNKNOWN',position='APPEND')
          ENDIF
          IF (paral%io_parent) THEN
             !             DO i1=1,6 !off-diagonal non-adiabatic couplings
             !                   DO i2=1,3 ! x,y,z
             !                          WRITE(50,'(50(E15.9,1X))') xfmqc%nacv(i2,:,:,i1)
             icount=0
             DO i=1,tshi%nroot_sh+1
                DO j=i+1,tshi%nroot_sh+1
                   icount=icount+1
                   DO i2=1,3
                      WRITE(50,'(500(E15.9,1X))')                      &
                           2.0d0*dreal(c_sh(1,j)*dconjg(c_sh(1,i))) *       & ! TODO: MIN: c_sh(1 or c_sh(2 ??
                           (xfmqc%eigv(i)-xfmqc%eigv(j)) *      &
                           xfmqc%nacv(i2,:,:,icount)
                   ENDDO
                ENDDO
             ENDDO
             !                   ENDDO
             !             ENDDO
          ENDIF
          IF (paral%io_parent)&
               CLOSE(50)
       ENDIF
       !
       IF (tshl%txfmqc) THEN 
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file=filen_fqmom_xf,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fqmom_xf,status='NEW')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_fqmom_xf,status='UNKNOWN',position='APPEND')
          ENDIF
          IF (paral%io_parent) THEN
             DO i1=1,3 ! x,y,z
                DO i2=1,maxsys%nax
                   DO i3=1,maxsys%nsx
                      !          WRITE(50,'(100(E15.9,1X))') xfmqc%fa_ni(i2,1,1),xfmqc%fa_ni(i2,1,2), &
                      !             xfmqc%fa_ni(i2,1,3),xfmqc%fa_ni(i2,2,3),xfmqc%fa_ni(i2,3,3),xfmqc%fa_ni(i2,4,3)
                      WRITE(50,'(100(E15.9,1X))') xfmqc%fa_ni(i1,i2,i3)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF (paral%io_parent) CLOSE(50)
       ENDIF
       ! 
       IF (td01%fstate.EQ.0) THEN
          tshl%s0_sh=.TRUE.
       ELSE
          tshl%s0_sh=.FALSE.
       ENDIF
    ENDIF
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_shcoeff
  ! ==================================================================
  SUBROUTINE read_shcoeff
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(len=100) :: filen, filen_a_sh, filen_c_sh, filen_e_sh, &
      filen_f_sh, filen_f_xf, filen_p_sh, filen_s_sh, filen_v_sh
    CHARACTER(len=12)                        :: cflbod
    INTEGER                                  :: i1, i2, ir, is, iv, nroot
    LOGICAL                                  :: fexist
    REAL(real_8)                             :: im(100), re(100), &
                                                temp(10000), temp2(10000)

! 

    nroot=tshi%nroot_sh
    ! 
    IF (paral%parent) THEN
       ! 
       IF (tshl%txfmqc.AND.tmw) THEN
          cflbod='SH_COEFC_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_c_sh=filen(i1:i2)//'.dat'

          cflbod='SH_COEFA_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_a_sh=filen(i1:i2)//'.dat'

          cflbod='SH_VPOTS_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_v_sh=filen(i1:i2)//'.dat'

          cflbod='SH_COUPL_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_s_sh=filen(i1:i2)//'.dat'

          cflbod='SH_PROBS_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_p_sh=filen(i1:i2)//'.dat'

          cflbod='SH_ENERG_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_e_sh=filen(i1:i2)//'.dat'

          cflbod='SH_STATE_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_f_sh=filen(i1:i2)//'.dat'

          cflbod='XF_ACFOR_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
          CALL xstring(filen,i1,i2)
          filen_f_xf=filen(i1:i2)//'.dat'
       ELSE
          filen_c_sh="SH_COEFC.dat"
          filen_a_sh="SH_COEFA.dat"
          filen_v_sh="SH_VPOTS.dat"
          filen_s_sh="SH_COUPL.dat"
          filen_p_sh="SH_PROBS.dat"
          filen_e_sh="SH_ENERG.dat"
          filen_f_sh="SH_STATE.dat"
       ENDIF

       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_c_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'ERROR ON READING SURFACE HOPPING FILES'
          IF (paral%io_parent)&
               WRITE(6,*) 'FILE SHOULD EXIST'
          STOP
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_c_sh,status='UNKNOWN',position='APPEND')
          IF (paral%io_parent)&
               BACKSPACE(50)
       ENDIF
       IF (paral%io_parent)&
            READ(50,'(I10,40F15.9)') tshi%shstep,&
            (re(ir),im(ir),ir=1,nroot+1)
       DO ir=1,nroot+1
          c_sh(2,ir)=re(ir)+(0._real_8,1._real_8)*im(ir)
       ENDDO
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_v_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'ERROR ON READING SURFACE HOPPING FILES'
          IF (paral%io_parent)&
               WRITE(6,*) 'FILE SHOULD EXIST'
          STOP
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_v_sh,status='UNKNOWN',position='APPEND')
          IF (paral%io_parent)&
               BACKSPACE(50)
       ENDIF
       IF (paral%io_parent)&
            READ(50,'(I10,20F15.9)') tshi%shstep,&
            (v_sh(2,ir),ir=1,nroot+1)
       DO ir=1,nroot+1
          v_sh(1,ir)=v_sh(2,ir)
       ENDDO
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_s_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'ERROR ON READING SURFACE HOPPING FILES'
          IF (paral%io_parent)&
               WRITE(6,*) 'FILE SHOULD EXIST'
          STOP
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_s_sh,status='UNKNOWN',position='APPEND')
          IF (paral%io_parent)&
               BACKSPACE(50)
       ENDIF
       IF (paral%io_parent)&
            READ(50,'(I10,1000F15.9)') tshi%shstep,&
            (temp(ir),ir=1,(nroot+1)*(nroot+1))
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       is=0
       DO ir=1,nroot+1
          DO iv=1,nroot+1
             is=is+1
             shsigma(1,ir,iv)=temp(is)
             shsigma(2,ir,iv)=temp(is)
             shsigmaf(1,ir,iv)=temp2(is)
             shsigmaf(2,ir,iv)=temp2(is)
          ENDDO
       ENDDO
       ! 
       fexist=.FALSE.
       IF (paral%io_parent)&
            INQUIRE(file=filen_f_sh,exist=fexist)
       IF (.NOT.fexist) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'ERROR ON READING SURFACE HOPPING FILES'
          IF (paral%io_parent)&
               WRITE(6,*) 'FILE SHOULD EXIST'
          STOP
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=50,file=filen_f_sh,status='UNKNOWN',position='APPEND')
          IF (paral%io_parent)&
               BACKSPACE(50)
       ENDIF
       IF (paral%io_parent)&
            READ(50,'(2I10)') tshi%shstep,td01%fstate
       IF (paral%io_parent)&
            CLOSE(50)
       ! 
       IF( tshl%txfmqc ) THEN ! MIN: for the sake of safety
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file=filen_f_xf,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'ERROR ON READING SURFACE HOPPING FILES'
             IF (paral%io_parent)&
                  WRITE(6,*) 'FILE SHOULD EXIST'
             STOP
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file=filen_f_xf,status='UNKNOWN',position='APPEND')
             IF (paral%io_parent)&
                  BACKSPACE(50)
          ENDIF
          IF (paral%io_parent)&
               READ(50,'(10000(E15.9,1X))') xfmqc%f_nli
          ![FEDE
          !IF (paral%io_parent) THEN
          !  WRITE(6,*) 'Im Here'
          !  WRITE(6,'(100000(E15.9,1X))') xfmqc%f_nli
          !END IF
          !FEDE]
          IF (paral%io_parent)&
               CLOSE(50)
       ENDIF

       IF (td01%fstate.EQ.0) THEN
          tshl%s0_sh=.TRUE.
       ELSE
          tshl%s0_sh=.FALSE.
       ENDIF

    ENDIF
    ! 
    CALL mp_bcast(c_sh,SIZE(c_sh),parai%source,parai%allgrp)
    CALL mp_bcast(v_sh,SIZE(v_sh),parai%source,parai%allgrp)
    CALL mp_bcast(shsigma,SIZE(shsigma),parai%source,parai%allgrp)
    CALL mp_bcast(shsigmaf,SIZE(shsigmaf),parai%source,parai%allgrp)
    CALL mp_bcast(td01%fstate,parai%source,parai%allgrp)
    CALL mp_bcast(tshi%shstep,parai%source,parai%allgrp)
    CALL mp_bcast(tshl%s0_sh,parai%source,parai%allgrp)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE read_shcoeff
  ! ==================================================================
  SUBROUTINE sh_extpot
    ! ==--------------------------------------------------------------==
    ! 
    ! -- Variables --
    REAL(real_8), PARAMETER                  :: cau = 137.0359885_real_8 

    REAL(real_8)                             :: freq, period, time_l

! APOT IS THE VECTOR POTENTIAL
! ETPOT IS THE FIELD
! AAMPL IS THE AMPLITUDE OF THE FIELD (!UNITS!)
! 
! THE PULSE IS HARDCODED IN THE FOLLOWING SECTION.
! You can use the following parameter read from input:
! AAMPL : the amplitude of the field
! ADIR  : the direction of the field (1:X, 2:Y, 3:Z, 4: isotropic)
! AFREQ : the frequency of the pulse
! APARA1: additional input parameter that can be used to 
! describe the pulse 
! The actual current time is:
! TIME  : absolute total time

    time_l=tshi%shstep*dt_ions
    ! The user needs to provide the expressions for:
    ! APOT : value of the vector potential in the direction ADIR
    ! ETPOT: corresponding electric field (=- 1/c * dA/dt)
    ! -------------------------------------------------------------------
    freq=tshf%afreq
    period=tshf%apara1
    ! 
    tshf%apot=-tshf%aampl*SIN(freq*time_l)
    ! 
    tshf%etpot=(tshf%aampl/cau)*freq*COS(freq*time_l)
    ! -------------------------------------------------------------------
    ! 
    RETURN
  END SUBROUTINE sh_extpot

  ! ==================================================================
  SUBROUTINE sh_lct
    ! ==--------------------------------------------------------------==
    ! 
    ! -- Variables --
    REAL(real_8)                             :: time_l

    time_l=tshi%shstep*dt_ions
    IF(paral%io_parent)WRITE(*,*) 'sh_lct_seed_t=', shlct%sh_lct_seed_t
    IF(paral%io_parent)WRITE(*,*) 'sh_lct_seed_m=', shlct%sh_lct_seed_m
    IF(paral%io_parent)WRITE(*,*) 'sh_lct_l=', shlct%sh_lct_l
    IF(paral%io_parent)WRITE(*,*) 'adir=', shlct%adir

    IF (time_l.LE.shlct%sh_lct_seed_t) THEN
       tshf%apot=-shlct%sh_lct_seed_m
       tshf%etpot=shlct%sh_lct_seed_m
    ELSE
       IF (paral%io_parent) WRITE(*,*) 'trans. dipole', shlct%trdipl
       ! BASILE: Corrected pulse
       tshf%apot=-shlct%sh_lct_l*(shlct%trdipl* &
            AIMAG(CONJG(c_sh(1,(shlct%lct_state+1)))* &
            c_sh(1,1)))
    ENDIF
    RETURN
  END SUBROUTINE sh_lct
  ! ======================================================================
  SUBROUTINE ecoupl(nstate,c0,nroot,c1,eigs)
    ! ======================================================================
    ! These from raman_p:

    ! Input:
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: nroot
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,nroot)
    REAL(real_8)                             :: eigs(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ecoupl'
    REAL(real_8), PARAMETER                  :: cau = 137.0359885_real_8 

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, id, ierr, j, jd, length
    REAL(real_8)                             :: trdip(3)
    REAL(real_8), ALLOCATABLE                :: scros(:)

! Local:
! ------------------------------------------------------------------
    IF (.NOT.td03%tda) RETURN
    IF (cntl%tlsd) RETURN
    IF (td01%ns_tri.NE.0) RETURN
    IF (td01%msubta.GT.0) RETURN
    ! ==--------------------------------------------------------------==
    IF (symmt%tpgauto .OR. symmt%tmsym .OR. symmi%indpg.NE.0) THEN
       CALL stopgm('td_ext_pot',&
            'POINT GROUP symmetry not implemented.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Allocate SCROS
    length=0
    CALL give_scr_opeigr(length,tag,nstate)
    ALLOCATE(scros(length),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) CALL prmem('     TD_OS')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"LB:",61("-"))')
    ! DO ILRV = 1, NROOT
    ! CALL TD_BERRY(C0,C1(1,1,ILRV),NSTATE,TRDIP)
    ! enddo

    ! DO I=0,NROOT
    i=0
    DO j=1,3
      trdip(j)=0._real_8
    ENDDO
    DO j=i+1,nroot
       CALL td_berry(c0,c1(1,1,j),nstate,trdip)
       IF (j.eq.shlct%lct_state) THEN
        IF (shlct%adir.eq.1) THEN
         shlct%trdipl = dabs(trdip(1))
        ELSEIF (shlct%adir.eq.2) THEN
         shlct%trdipl = dabs(trdip(2))
        ELSEIF (shlct%adir.eq.3) THEN
         shlct%trdipl = dabs(trdip(3))
        ELSEIF (shlct%adir.eq.4) THEN
         shlct%trdipl = (1.d0/dsqrt(3.d0))*(trdip(1)+trdip(2)+&
                                       trdip(3))
        ENDIF
        if (paral%io_parent) write(*,*) 'trans. dipole', shlct%trdipl,&
                               'nroot', nroot
        if (paral%io_parent) write(*,*) 'trans. dipole', shlct%trdipl, 'J=', j,&
                               'LCT_STATE=', shlct%lct_state
       ENDIF
       id=i+1
       jd=j+1
       IF (tshi%adir.EQ.1.OR.shlct%adir.EQ.1) THEN
          ! if (parent) write(6,*) 'ivano',SHSIGMA(1,ID,JD),ETPOT*TRDIP(1)
          ! if (ID.EQ.1 .AND. JD.EQ.2 .AND. parent) 
          ! &     write(111,*) -APOT/CAU*EIGS(J)*TRDIP(1),eigs(J)
          shsigmaf(1,id,jd)=tshf%apot*dabs(trdip(1))
       ELSEIF (tshi%adir.EQ.2.OR.shlct%adir.EQ.2) THEN
          shsigmaf(1,id,jd)=tshf%apot*dabs(trdip(2))
       ELSEIF (tshi%adir.EQ.3.OR.shlct%adir.EQ.3) THEN
          shsigmaf(1,id,jd)=tshf%apot*dabs(trdip(3))
       ELSEIF (tshi%adir.EQ.4.OR.shlct%adir.EQ.4) THEN
          ! 
          IF (i.EQ.0) THEN
             shsigmaf(1,id,jd)=&
                  -(tshf%apot)*(1._real_8/SQRT(3._real_8))*(trdip(1)+trdip(2)+&
                  trdip(3))
          ELSE
             shsigmaf(1,id,jd)=&
                  -(tshf%apot)*(1._real_8/SQRT(3._real_8))*(trdip(1)&
                  +trdip(2)+&
                  trdip(3))
          ENDIF
          !
       ENDIF
       shsigmaf(1,jd,id)=-shsigmaf(1,id,jd)
    ENDDO
    !
    DO i=2,nroot-1
     DO j=i+1,nroot
       shsigmaf(1,i,j)=0.0_real_8
       shsigmaf(1,j,i)=shsigmaf(1,i,j)
     ENDDO
    ENDDO
    DO id=1,nroot+1
       shsigmaf(1,id,id)=0.0_real_8
       shsigmaf(2,id,id)=0.0_real_8
    ENDDO
    ! ENDDO
    DEALLOCATE(scros,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE ecoupl
  ! ==================================================================
  SUBROUTINE td_berry(c0,c1,nstate,trdip)
    ! ==================================================================
    ! Input:
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: trdip(3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_berry'

    COMPLEX(real_8)                          :: cone, czero, dd
    COMPLEX(real_8), ALLOCATABLE             :: cw(:,:), cwork(:,:,:), &
                                                ddmat_1(:,:), ddmat_p(:,:), &
                                                f_raman(:,:,:), sc0(:), &
                                                Work(:)
    INTEGER                                  :: i, ierr, ig, il_cw, il_cw1, &
                                                info, is1, isub, k, n1, &
                                                nmcol, nop1, nop2, nop3, nxx
    INTEGER, ALLOCATABLE                     :: ipiv(:), mapcol(:), mapful(:)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: a, b, c, d, f3, f4, osx, osy, &
                                                osz, ox, spinfac

! ------------------------------------------------------------------

    CALL tiset('    td_os_raman_p',isub)
    cone=CMPLX(1._real_8,0._real_8,kind=real_8)
    czero=CMPLX(0._real_8,0._real_8,kind=real_8)
    ! ------------------------------------------------------------------
    ! Memory allocations
    ! ------------------------------------------------------------------
    !il_framan=6*ncpw%ngw*nstate+12
    il_cw=2*ncpw%ngw*nstate+8
    il_cw1=4*ncpw%ngw*nstate
    ALLOCATE(f_raman(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(f_raman)!,3*ngw*nstate)
    ALLOCATE(ddmat_p(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddmat_1(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipiv(50*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(Work(nstate*50),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    lenbk=nxxfun(nstate)
    nxx=MAX(2*lenbk*parai%nproc,il_cw)
    il_cw=nxx
    ALLOCATE(sc0(nxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cw(il_cw,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cw)!,SIZE(cw))
    !IF (ifirst.EQ.0) THEN
    n1=0
    DO i=0,parai%nproc-1
       n1=MAX(n1,parap%sparm(3,i))
    ENDDO
    ngwmax=n1
    ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmcol=parai%nproc*ngwmax
    ALLOCATE(mapcol(nmcol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cwork(ncpw%ngw,nstate,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setdip(mapful,mapcol)
    ifirst=1
    !ENDIF
    ! ------------------------------------------------------------------
    ! .....Spin degeneracy factor
    IF (cntl%tlsd) THEN
       spinfac=1._real_8
    ELSE
       spinfac=2._real_8
    ENDIF
    ! ------------------------------------------------------------------
    ! Set up phase operator matrices (k=1,3) and f_raman.
    ! ------------------------------------------------------------------
    osx = 0._real_8
    osy = 0._real_8
    osz = 0._real_8
    DO k=1,3
       CALL zeroing(cw)!,SIZE(cw))
       CALL zeroing(ddmat_p)!,nstate*nstate)
       CALL zeroing(ddmat_1)!,nstate*nstate)
       CALL zeroing(cwork)!,2*ngw*nstate)
       CALL zeroing(sc0)!,SIZE(sc0))
       CALL zeroing(f_raman)!,il_framan/2)
       ! .......Compute cwork = opeigr*c(g) and Q(i,j) = <phi_i | opeigr | phi_j>.
       IF (parm%ibrav.EQ.1 .OR. parm%ibrav.EQ.8) THEN
          nop1=k
          nop2=0
          nop3=0
       ELSE IF (parm%ibrav.EQ.2) THEN
          IF (k.EQ.1) THEN
             nop1=1
             nop2=2
             nop3=0
          ELSE IF (k.EQ.2) THEN
             nop1=1
             nop2=3
             nop3=0
          ELSE
             nop1=2
             nop2=3
             nop3=0
          ENDIF
       ELSE IF (parm%ibrav.EQ.3) THEN
          IF (k.EQ.1) THEN
             nop1=1
             nop2=-2
             nop3=-3
          ELSE IF (k.EQ.2) THEN
             nop1=1
             nop2=2
             nop3=-3
          ELSE
             nop1=1
             nop2=2
             nop3=3
          ENDIF
       ELSE
          IF (paral%io_parent)&
               WRITE (6,*)&
               'TD_OS_BERRY: ONLY SUPERCELLS OF TYPE 1,2,3,8 allowed.'
       ENDIF
       CALL mp_sync(parai%allgrp)
       CALL opeigr_p(c0,cw,sc0,nstate,mapful,mapcol,ddmat_p,&
            nop1,nop2,nop3,dd,cwork)
       CALL mp_sync(parai%allgrp)
       ! ........Make a copy of ddmat_p (ddmat_p need be inverted).
       CALL zcopy(nstate*nstate,ddmat_p,1,ddmat_1,1)
       ! ------------------------------------------------------------------
       ! Compute inverse of ddmat_p (Q^-1) for mu=k in Putrino et al.
       ! Note ddmat_p is overwrittwen.
       ! ------------------------------------------------------------------
       CALL zgetrf(nstate,nstate,ddmat_p,nstate,ipiv,info)
       IF (info.EQ.0) THEN
          CALL zgetri(nstate,ddmat_p,nstate,ipiv,Work,nstate,info)
       ELSE
          CALL stopgm('td_os_berry','error in matrix inversion',& 
               __LINE__,__FILE__)
       ENDIF
       ! ------------------------------------------------------------------
       ! Build f_raman.
       ! ------------------------------------------------------------------
       ! ........cwork(1,*,*) = [opeigr*c(g)] * Q^-1;
       ! ........cwork(2,*,*) = [opeigr*c(-g)] * Q^-1.
       DO i=1,2
          CALL zeroing(cw)!,SIZE(cw))
          CALL zgemm('N','N',ncpw%ngw,nstate,nstate,cone,cwork(1,1,i)&
               ,ncpw%ngw,ddmat_p(1,1),nstate,czero,cw(1,1),ncpw%ngw)
          CALL zcopy(ncpw%ngw*nstate,cw(1,1),1,cwork(1,1,i),1)
       ENDDO
       CALL mp_sync(parai%allgrp)
       IF (geq0) THEN
          is1=2
          DO i=1,nstate
             b=AIMAG(cwork(1,i,1))
             f_raman(1,i,k)=CMPLX(b,0._real_8,kind=real_8)
          ENDDO
       ELSE
          is1=1
       ENDIF
       ! This is to have a sum running on ngw below.
       DO i=1,nstate
          DO ig=is1,ncpw%ngw
             a=REAL(cwork(ig,i,1))
             b=AIMAG(cwork(ig,i,1))
             c=REAL(cwork(ig,i,2))
             d=AIMAG(cwork(ig,i,2))
             f3=a-c
             f4=b+d
             f_raman(ig,i,k)=0.5_real_8*CMPLX(f4,-f3,kind=real_8)
          ENDDO
       ENDDO
       ! ------------------------------------------------------------------
       ! Compute matrix elements
       ! ------------------------------------------------------------------
       ! .......Compute os = sum_k sum_i <chi_i| opeigr(k) | phi_i> / Q + c.c.
       IF (geq0) THEN
          is1=2
          DO i=1,nstate
             ox      =&
                  REAL( CONJG(c1(1,i))*f_raman(1,i,k) )&
                  +REAL( CONJG(f_raman(1,i,k))*c1(1,i) )
             IF (k.EQ.1) osx = osx + ox
             IF (k.EQ.2) osy = osy + ox
             IF (k.EQ.3) osz = osz + ox
          ENDDO
       ELSE
          is1=1
       ENDIF

       DO i=1,nstate
          DO ig=is1,ncpw%ngw
             ox      =&
                  2._real_8*REAL( CONJG(c1(ig,i))*f_raman(ig,i,k)  )&
                  +2._real_8*REAL( CONJG(f_raman(ig,i,k))*c1(ig,i)  )
             IF (k.EQ.1) osx = osx + ox
             IF (k.EQ.2) osy = osy + ox
             IF (k.EQ.3) osz = osz + ox
          ENDDO
       ENDDO
       ! ------------------------------------------------------------------
    ENDDO                                    ! On cartesian components
    CALL mp_sync(parai%allgrp)
    CALL mp_sum(osx,parai%allgrp)
    CALL mp_sum(osy,parai%allgrp)
    CALL mp_sum(osz,parai%allgrp)
    ! 
    trdip(1)=spinfac*osx
    trdip(2)=spinfac*osy
    trdip(3)=spinfac*osz
    ! 
    DEALLOCATE(f_raman,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(Work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ipiv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddmat_1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddmat_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapful,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapcol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    td_os_raman_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_berry
  ! ==================================================================

  SUBROUTINE xf_qmom
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'xf_qmom'

    INTEGER                                  :: fileunit, i, i_traj, ia, &
                                                ierr, index_ij, is, istate, &
                                                j, j_traj, jstate, k, power
    INTEGER, ALLOCATABLE                     :: ntr_j(:,:,:)
    REAL(real_8)                             :: dist2, dist_cutoff, n_j, &
                                                threshold
    REAL(real_8), ALLOCATABLE                :: avR(:,:,:), avR2(:,:,:), &
                                                d_ij(:,:), exp_g_i(:,:), &
                                                g_i(:), slope_i(:,:,:,:)

!, slope_i

    IF (.NOT.tshl%txfmqc) RETURN

    ALLOCATE(g_i(xfmqc%n_xftraj),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(g_i)
    ALLOCATE(d_ij(xfmqc%n_xftraj,xfmqc%n_xftraj),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(d_ij)
    ![FEDE
    ALLOCATE(exp_g_i(xfmqc%n_xftraj,xfmqc%n_xftraj),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(exp_g_i)
    ALLOCATE(slope_i(3,maxsys%nax,maxsys%nsx,xfmqc%n_xftraj),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(slope_i)
    ALLOCATE(avR(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(avR)
    ALLOCATE(avR2(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(avR2)
    ALLOCATE(ntr_j(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ntr_j)
    dist_cutoff=0.50_real_8 !to be changed
    threshold=0.25_real_8!xfmqc%threshold
    !FEDE]

    IF (tmw) THEN
       i_traj=mwi%walker_id
    ELSE
       i_traj=1
    ENDIF

    ![FEDE ... calculation of variances
    DO i=1,xfmqc%n_xftraj
       avR=0._real_8
       avR2=0._real_8
       ntr_j=0
       ! in this loop I check the distance among the trajectories from the
       ! selected trajectry (i) and I estimate the variance associated to the
       ! gaussian on the trajectory (i)
       ! "I have 3N_n variances for each trajectories"
       dist2=0._real_8
       !avR=0._real_8
       !avR2=0._real_8
       !ntr_j=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                DO j=1,xfmqc%n_xftraj
                   dist2=(xfmqc%tauall(k,ia,is,i)-xfmqc%tauall(k,ia,is,j))**2
                   IF(SQRT(dist2) .LT. dist_cutoff) THEN
                      avR(k,ia,is)=avR(k,ia,is)+xfmqc%tauall(k,ia,is,j)
                      ! average position of trajectories within a sphere
                      avR2(k,ia,is)=avR2(k,ia,is)+xfmqc%tauall(k,ia,is,j)**2 
                      ! average position squared of trajectories within a sphere
                      ntr_j(k,ia,is)=ntr_j(k,ia,is)+1
                      ! number of trajectories within the sphere
                   ENDIF
                ENDDO ! j loop
                xfmqc%sigma(k,ia,is,i)=(avR2(k,ia,is)/DBLE(ntr_j(k,ia,is))- &
                     (avR(k,ia,is)/DBLE(ntr_j(k,ia,is)))**2)/SQRT(SQRT(DBLE(ntr_j(k,ia,is))))
                IF(xfmqc%sigma(k,ia,is,i) .LT. threshold .OR. ntr_j(k,ia,is)==0) &
                     xfmqc%sigma(k,ia,is,i)=dist_cutoff
                ! variance associated to each trajectory depending on the spreading
                ! of the trajectories close to the trajectory (i)
             ENDDO   ! k loop
          ENDDO     ! ia loop
       ENDDO       ! is loop
    ENDDO         ! i loop
    !IF(paral%io_parent) THEN
    !    WRITE(6,*) 'HERE ARE 3*Natoms*Ntraj values; the distance cutoff is',dist_cutoff
    !    WRITE(6,*) 'SIGMAS for TRAJ', i_traj
    !    WRITE(6,*) xfmqc%sigma(:,1,1,i_traj),xfmqc%sigma(:,1,2,i_traj),                    &
    !        xfmqc%sigma(:,1,3,i_traj),xfmqc%sigma(:,2,3,i_traj),xfmqc%sigma(:,2,3,i_traj), &
    !        xfmqc%sigma(:,4,3,i_traj) 
    !ENDIF

    g_i(i_traj)=0._real_8
    ![FEDE - this part computes the nuclear density as the sum of normalized
    !Gaussians centered at the positions of the trajectories (label j)
    !with variances as computed above
    power=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             power=power+1
          ENDDO
       ENDDO
    ENDDO
    DO j=1,xfmqc%n_xftraj
       n_j=1.0_real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                n_j=n_j*xfmqc%sigma(k,ia,is,j)
             ENDDO
          ENDDO
       ENDDO
       n_j=n_j*(SQRT(2._real_8*pi))**power
       IF(n_j==0.0_real_8) THEN
          WRITE(*,*) "GAUSSIANS NORMALIZATION = 0"
          STOP
       END IF
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                exp_g_i(i_traj,j)=exp_g_i(i_traj,j)- &
                     ((xfmqc%tauall(k,ia,is,i_traj)-xfmqc%tauall(k,ia,is,j))/ &
                     xfmqc%sigma(k,ia,is,j))**2/2._real_8 
             ENDDO
          ENDDO
       ENDDO
       g_i(i_traj)=g_i(i_traj)+DEXP(exp_g_i(i_traj,j))/n_j
       IF(g_i(i_traj) .LT. 0.000001_real_8) g_i(i_traj)=0.0_real_8
    ENDDO
    !FEDE]

    ![FEDE - compute W_ij whose sum over j is the slope of teh quantum momentum
    !when a sum of Gaussians is used to approximate the nuclear density
    DO j=1,xfmqc%n_xftraj
       n_j=1.0_real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                n_j=n_j*xfmqc%sigma(k,ia,is,j)
             ENDDO
          ENDDO
       ENDDO
       n_j=n_j*(SQRT(2._real_8*pi))**power
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                IF(g_i(i_traj) .GT. 0.0000001_real_8) &
                     xfmqc%w_ij(k,ia,is,i_traj,j)=(1.0_real_8/xfmqc%sigma(k,ia,is,j)**2)* & !New - 17.03.16
                     DEXP(exp_g_i(i_traj,j))/(n_j*g_i(i_traj))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !FEDE]

    ![FEDE - compute the slope form W_ij, summing over j
    DO j=1,xfmqc%n_xftraj
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                slope_i(k,ia,is,i_traj)=slope_i(k,ia,is,i_traj)-xfmqc%w_ij(k,ia,is,i_traj,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !IF(paral%IO_parent) &
    !  WRITE (6,'(A6,i4,100f14.4)') "SLOPE",i_traj,slope_i(:,1,1,i_traj),slope_i(:,1,2,i_traj),   &
    !      slope_i(:,1,3,i_traj),slope_i(:,2,3,i_traj),slope_i(:,3,3,i_traj),slope_i(:,4,3,i_traj)
    !FEDE]
    ! Compute the quantity that will be summed up over the trajectories to give
    ! the center of the quantum momentum
    CALL zeroing(xfmqc%bw0_nkli)
    DO is=1,ions1%nsp  ! atomic species
       DO ia=1,ions0%na(is) ! number of atoms per species
          DO k=1,3 ! x,y,z
             index_ij=0
             DO istate=1,tshi%nroot_sh+1
                DO jstate=istate+1,tshi%nroot_sh+1
                   index_ij=index_ij+1
                   ! denomenator
                   xfmqc%d0_ni(k,ia,is,index_ij)=0._real_8
                   DO j_traj=1,xfmqc%n_xftraj
                      IF(xfmqc%w_ij(k,ia,is,i_traj,j_traj).GT.xfmqc%threshold) THEN
                         xfmqc%d0_ni(k,ia,is,index_ij)=xfmqc%d0_ni(k,ia,is,index_ij)+         &
                              xfmqc%brho_l(istate,j_traj) *          &
                              xfmqc%brho_l(jstate,j_traj) *          &
                              (xfmqc%brho_l(istate,j_traj)+xfmqc%brho_l(jstate,j_traj))* &
                              (xfmqc%bf_nli(k,ia,is,istate,j_traj)-  &
                              xfmqc%bf_nli(k,ia,is,jstate,j_traj))*slope_i(k,ia,is,j_traj)
                      ENDIF
                   ENDDO
                   ! numerator
                   xfmqc%r0_ni(k,ia,is,index_ij)=xfmqc%brho_l(istate,i_traj) *          &
                        xfmqc%brho_l(jstate,i_traj) *          &
                        (xfmqc%brho_l(istate,i_traj)+xfmqc%brho_l(jstate,i_traj))* &
                        (xfmqc%bf_nli(k,ia,is,istate,i_traj)-  &
                        xfmqc%bf_nli(k,ia,is,jstate,i_traj))* &
                        xfmqc%tauall(k,ia,is,i_traj)*slope_i(k,ia,is,i_traj)
                   IF(ABS(xfmqc%d0_ni(k,ia,is,index_ij)).LT.0.00001_real_8) THEN !New - 24.03.16
                      xfmqc%bw0_nkli(k,ia,is,index_ij,i_traj)=0._real_8
                   ELSE
                      IF(paral%io_parent) THEN
                         xfmqc%bw0_nkli(k,ia,is,index_ij,i_traj)=xfmqc%r0_ni(k,ia,is,index_ij)/&
                              xfmqc%d0_ni(k,ia,is,index_ij)
                      ELSE
                         xfmqc%bw0_nkli(k,ia,is,index_ij,i_traj)=0._real_8
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO
    ENDDO
    ! the trajectories share the quantity that each has computed in order to sum it up later on
    ! and get the same center for all trajectories
    CALL mp_sum(xfmqc%bw0_nkli, & 
         3*maxsys%nax*maxsys%nsx*((tshi%nroot_sh+1)*(tshi%nroot_sh))/2*xfmqc%n_xftraj, &
         supergroup)
    !xfmqc%bw0_nkli is now summed up over the trajectories to get the center
    DO is=1,ions1%nsp  ! atomic species
       DO ia=1,ions0%na(is) ! number of atoms per species
          DO k=1,3 ! x,y,z
             index_ij=0
             DO istate=1,tshi%nroot_sh+1
                DO jstate=istate+1,tshi%nroot_sh+1
                   index_ij=index_ij+1
                   xfmqc%r0_ni(k,ia,is,index_ij)=0._real_8
                   DO j_traj=1,xfmqc%n_xftraj
                      IF(xfmqc%w_ij(k,ia,is,i_traj,j_traj).GT.xfmqc%threshold) THEN
                         xfmqc%r0_ni(k,ia,is,index_ij)=xfmqc%r0_ni(k,ia,is,index_ij)+         &
                              xfmqc%bw0_nkli(k,ia,is,index_ij,j_traj)
                      ENDIF
                   ENDDO
                   IF(ABS(slope_i(k,ia,is,i_traj)) .LT. 0.0000001_real_8) THEN
                      xfmqc%r0_ni(k,ia,is,index_ij)=xfmqc%tauall(k,ia,is,i_traj)
                      !ELSE
                      !  xfmqc%r0_ni(k,ia,is,index_ij) = xfmqc%r0_ni(k,ia,is,index_ij)!/slope_i(k,ia,is,i_traj)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
    ENDDO

    ![FEDE - output centers
    index_ij=0
    DO istate=1,tshi%nroot_sh+1
       DO jstate=istate+1,tshi%nroot_sh+1
          index_ij=index_ij+1
          IF (paral%io_parent) THEN
             DO is=1,ions1%nsp  
                DO ia=1,ions0%na(is)
                   fileunit=100+10*i_traj+index_ij
                   WRITE(fileunit,'(I10,40F15.9)') tshi%shstep,&
                        (xfmqc%r0_ni(k,ia,is,index_ij),k=1,3)
                ENDDO
             ENDDO
          END IF
       ENDDO
    ENDDO
    !FEDE]

    ! quantum momentum
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             index_ij=0
             DO istate=1,tshi%nroot_sh+1
                DO jstate=istate+1,tshi%nroot_sh+1
                   index_ij=index_ij+1
                   xfmqc%qmoment(k,ia,is,index_ij)= slope_i(k,ia,is,i_traj) *(xfmqc%tauall(k,ia,is,i_traj)-&
                        xfmqc%r0_ni(k,ia,is,index_ij))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! OUTPUT QUANTUM MOMENTUM FOR CHECK on PROTONATED FORMALDIMINE
    index_ij=0
    !DO istate=1,tshi%nroot_sh+1
    !  DO jstate=istate+1,tshi%nroot_sh+1
    !     index_ij=index_ij+1
    !     IF (paral%io_parent) & 
    !        WRITE(200+index_ij,'(I10,40F15.9)') tshi%shstep,&
    !           xfmqc%qmoment(1,1,1,index_ij),xfmqc%qmoment(2,1,1,index_ij),xfmqc%qmoment(3,1,1,index_ij), &
    !           xfmqc%qmoment(1,1,2,index_ij),xfmqc%qmoment(2,1,2,index_ij),xfmqc%qmoment(3,1,2,index_ij), &
    !           xfmqc%qmoment(1,1,3,index_ij),xfmqc%qmoment(2,1,3,index_ij),xfmqc%qmoment(3,1,3,index_ij), &
    !           xfmqc%qmoment(1,2,3,index_ij),xfmqc%qmoment(2,2,3,index_ij),xfmqc%qmoment(3,2,3,index_ij), &
    !           xfmqc%qmoment(1,3,3,index_ij),xfmqc%qmoment(2,3,3,index_ij),xfmqc%qmoment(3,3,3,index_ij), &
    !           xfmqc%qmoment(1,4,3,index_ij),xfmqc%qmoment(2,4,3,index_ij),xfmqc%qmoment(3,4,3,index_ij)
    !  ENDDO
    !ENDDO

    index_ij=0
    DO istate=1,tshi%nroot_sh+1
       DO jstate=istate+1,tshi%nroot_sh+1
          index_ij=index_ij+1
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   xfmqc%k_li(istate,jstate)=xfmqc%k_li(istate,jstate)+2._real_8/rmass%pma(is)*  &
                        xfmqc%qmoment(k,ia,is,index_ij)*xfmqc%f_nli(k,ia,is,istate) !New - 17.03.16 
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             xfmqc%fa_ni(k,ia,is)=0._real_8
             index_ij=0
             DO istate=1,tshi%nroot_sh+1
                DO jstate=istate+1,tshi%nroot_sh+1
                   xfmqc%fa_ni(k,ia,is)=xfmqc%fa_ni(k,ia,is)-                                    &
                        xfmqc%k_li(istate,jstate)*ABS(c_sh(2,istate))**2*ABS(c_sh(2,jstate))**2 * & !New - 17.03.16
                        (ABS(c_sh(2,istate))**2+ABS(c_sh(2,jstate))**2)*                          &
                        (xfmqc%f_nli(k,ia,is,istate)-xfmqc%f_nli(k,ia,is,jstate))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    xfmqc%fa_ni=xfmqc%fa_ni/2.0_real_8/DBLE(tshi%nroot_sh)

    !New - 17.03.16
    DO istate=1,tshi%nroot_sh+1
       xfmqc%cf(istate)=0._real_8
       DO jstate=1,tshi%nroot_sh+1
          IF(jstate/=istate) THEN
             xfmqc%cf(istate)= xfmqc%cf(istate)+                                          &
                  0.50_real_8*(xfmqc%k_li(jstate,istate)-xfmqc%k_li(istate,jstate))*       &
                  (ABS(c_sh(2,istate))**2+ABS(c_sh(2,jstate))**2)*ABS(c_sh(2,jstate))**2*  & 
                  xfmqc%cf(istate)
          ENDIF
       ENDDO
    ENDDO
    xfmqc%cf= xfmqc%cf/2.0_real_8/DBLE(tshi%nroot_sh)

    DEALLOCATE(g_i,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(d_ij,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(exp_g_i,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(slope_i,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(avR,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(avR2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ntr_j,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE xf_qmom
  ! ==================================================================

  SUBROUTINE xf_force_history
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER :: procedureN = 'xf_force_history'

    INTEGER                                  :: ia, is, k, l

    IF (.NOT.tshl%txfmqc) RETURN

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,tshi%nroot_sh+1
                IF( (ABS(c_sh(2,l)))**2 .GT. xfmqc%threshold .AND. &
                     (ABS(c_sh(2,l)))**2 .LT. (1.0_real_8-xfmqc%threshold) ) THEN
                   xfmqc%f_nli(k,ia,is,l)=xfmqc%f_nli(k,ia,is,l) + &
                        xfmqc%fion_state(k,ia,is,l)*dt_ions
                ELSE
                   xfmqc%f_nli(k,ia,is,l)= 0.0_real_8
                END IF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE xf_force_history
  ! ==================================================================

  SUBROUTINE xf_force_components(c0,c1,c2,sc0,rhoe,eigs)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), INTENT(inout)           :: c0(:,:,:), &
                                                c1(nkpt%ngwk,crge%n,*), &
                                                c2(:,:), sc0(:)
    REAL(real_8)                             :: rhoe(:,:), eigs(:)

    CHARACTER(*), PARAMETER :: procedureN = 'xf_force_components'

    INTEGER                                  :: i, i1, ia, icount, is, j, k, &
                                                nact, no, nvir

    IF (.NOT.tshl%txfmqc) RETURN
    ! TODO IVANO
    ! remember to fill array bfion_state, which is the buffer array 
    ! for the local fion_state
    ! the same for bfa_ni and fa_na

    CALL xf_force_history
    CALL xf_qmom
    !
    ! add force in eq. 4 first term
    CALL zeroing(fion)
    DO i1=1,tshi%nroot_sh+1
       DO is = 1,ions1%nsp
          DO ia = 1,ions0%na(is)
             fion(1,ia,is)=fion(1,ia,is) + &
                  (ABS(c_sh(1,i1)))**2 * xfmqc%fion_state(1,ia,is,i1) ! TODO MIN: c_sh(1,i1) or c_sh(2,i1) ???
             fion(2,ia,is)=fion(2,ia,is) + &
                  (ABS(c_sh(1,i1)))**2 * xfmqc%fion_state(2,ia,is,i1)
             fion(3,ia,is)=fion(3,ia,is) + &
                  (ABS(c_sh(1,i1)))**2 * xfmqc%fion_state(3,ia,is,i1)
             !if (paral%io_parent) then
             !write(6,*) ABS(c_sh(1,i1)),xfmqc%fion_state(1,ia,is,i1)
             !endif
          ENDDO
       ENDDO
    ENDDO
    ! add force in eq. 4 second term
    no=0
    DO i=1,crge%n
       IF (crge%f(i,1).GT.0.00001_real_8) no=no+1
    ENDDO
    !ivano
    nvir=crge%n-no
    IF (td03%tda.AND.td01%msubta.GT.0) THEN
       nact=td01%msubta
       !if (.not.td03%tdacanon) tdasubopt=.true.
    ELSE
       IF (cntl%tlsd) THEN
          nact=tdsp1%nupel+tdsp1%ndoel
       ELSE
          nact=no
       ENDIF
    ENDIF
    CALL vcouplings(c0(:,:,1),c1,c2,sc0,rhoe,eigs,xfmqc%eigv,        &
         crge%n,no,nvir,nact,tshi%nroot_sh)
    icount=0
    DO i=tshi%nroot_sh,1,-1
       xfmqc%eigv(i+1)=xfmqc%eigv(i)
    ENDDO
    xfmqc%eigv(1)=0.d0
    DO i=1,tshi%nroot_sh+1
       DO j=i+1,tshi%nroot_sh+1
          icount=icount+1
          DO is = 1,ions1%nsp
             DO ia = 1,ions0%na(is)
                DO k=1,3
                   fion(k,ia,is)=fion(k,ia,is) + &
                        2.0d0*dreal(c_sh(1,j)*dconjg(c_sh(1,i))) * & ! TODO: MIN: c_sh(1 or c_sh(2 ??
                        (xfmqc%eigv(i)-xfmqc%eigv(j)) *      &
                        xfmqc%nacv(k,ia,is,icount)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! add force in eq. 18
    DO is = 1,ions1%nsp
       DO ia = 1,ions0%na(is)
          DO k=1,3
             fion(k,ia,is)=fion(k,ia,is) + &
                  xfmqc%fa_ni(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    !] EXACT FACTORIZATION
  END SUBROUTINE xf_force_components
  ! ==================================================================

  SUBROUTINE xf_fill_all_buffers ! MIN: modified by MIN

    INTEGER                                  :: i_traj, ia, is, ist, k

    IF( tshl%txfmqc .AND. tmw ) THEN
       i_traj = mwi%walker_id
    ELSE
       i_traj = 1
    ENDIF
    ! tauall
    CALL zeroing(xfmqc%tauall)
    DO is = 1,ions1%nsp
       DO ia = 1,ions0%na(is)
          DO k=1,3
             IF(paral%io_parent) THEN
                xfmqc%tauall(k,ia,is,i_traj)=tau0(k,ia,is)
             ELSE
                xfmqc%tauall(k,ia,is,i_traj)=0._real_8
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(xfmqc%tauall,3*maxsys%nax*maxsys%nsx*mwi%nwalk,supergroup)
    ! bf_nli
    CALL zeroing(xfmqc%bf_nli)
    DO is = 1,ions1%nsp
       DO ia = 1,ions0%na(is)
          DO k=1,3
             DO ist = 1, tshi%nroot_sh+1
                IF(paral%io_parent) THEN
                   xfmqc%bf_nli(k,ia,is,ist,i_traj) = xfmqc%f_nli(k,ia,is,ist)
                ELSE
                   xfmqc%bf_nli(k,ia,is,ist,i_traj)=0._real_8
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(xfmqc%bf_nli,3*maxsys%nax*maxsys%nsx*(tshi%nroot_sh+1)*mwi%nwalk,supergroup)
    ! brho_l
    CALL zeroing(xfmqc%brho_l)
    DO ist = 1, tshi%nroot_sh+1
       IF(paral%io_parent) THEN
          xfmqc%brho_l(ist,i_traj) = ABS(c_sh(2,ist))**2
       ELSE
          xfmqc%brho_l(ist,i_traj) = 0._real_8
       END IF
    ENDDO
    CALL mp_sum(xfmqc%brho_l,(tshi%nroot_sh+1)*mwi%nwalk,supergroup)

  END SUBROUTINE xf_fill_all_buffers

END MODULE sh_tddft_utils
