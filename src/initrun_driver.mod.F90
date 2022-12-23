MODULE initrun_driver
  USE ainitwf_utils,                   ONLY: ainitwf
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE bsym,                            ONLY: bsclcs,&
                                             resthswf
  USE calc_alm_utils,                  ONLY: calc_alm
  USE chksym_utils,                    ONLY: updatsym
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE geofile_utils,                   ONLY: geofile
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: lractive
  USE merge_utils,                     ONLY: merge
  USE metr,                            ONLY: metr_com
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: clc,&
                                             lqmmm,&
                                             rtr_l
  USE mp_interface,                    ONLY: mp_bcast
  USE newcell_utils,                   ONLY: newcell
  USE nlcc,                            ONLY: corel,&
                                             vnlcc,&
                                             vnlt
  USE ortho_utils,                     ONLY: ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: ipcurr,&
                                             np_low
  USE pslo,                            ONLY: pslo_com
  USE ranc_utils,                      ONLY: ranc
  USE ranp_utils,                      ONLY: ranp
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rinitwf_driver,                  ONLY: rinitwf
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: iteropt
  USE rrandd_utils,                    ONLY: rrandd
  USE rrane_utils,                     ONLY: rrane
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec
  USE setsc_utils,                     ONLY: ihmat
  USE shop,                            ONLY: sh02
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: irec_co,&
                                             irec_rho,&
                                             irec_wf,&
                                             restart1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpbl,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdw_utils,                       ONLY: rwannier
  USE vdwcmod,                         ONLY: nfragx,&
                                             nwfcx,&
                                             rwann,&
                                             swann,&
                                             tauref,&
                                             trwanncx,&
                                             vdwl,&
                                             vdwwfl
  USE wrgeo_utils,                     ONLY: wrgeof
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initrun

CONTAINS

  ! ==================================================================
  SUBROUTINE initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
    ! ==--------------------------------------------------------------==
    ! == Initialization of a run (read from RESTART file or not)      ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: irec(:)
    COMPLEX(real_8)                          :: c0(:,:,:), cm(*), sc0(*)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'initrun'

    INTEGER                                  :: i, ierr, ikind, ikpt, ipx, &
                                                isub, kbeg, kend, kinc, &
                                                msglen, nkpoint, nstate, nwfc
    LOGICAL                                  :: status, statusdummy

! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    CALL mm_dim(mm_go_mm,status)

    IF (corel%tinlc) THEN
       ALLOCATE(vnlt(ncpw%nhg*2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(vnlcc(ncpw%nhg,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    CALL read_irec(irec)
    IF (cntl%tshop) THEN
       nstate=sh02%nst_s0+sh02%nst_s1
    ELSE
       nstate=crge%n
    ENDIF
    IF (cntl%tddft) THEN
       IF (.NOT.restart1%restart) CALL stopgm('INITRUN','TDDFT NEEDS RESTART',& 
            __LINE__,__FILE__)
       IF (.NOT.restart1%rlr) CALL stopgm('INITRUN','TDDFT NEEDS RESTART LINRES',& 
            __LINE__,__FILE__)
       lractive=.TRUE.
    ENDIF
    IF (.NOT.restart1%restart) THEN
       IF (bsclcs.EQ.1) resthswf = .TRUE.
#if defined (__GROMOS)
       IF (cntl%tqmmm)THEN
          IF (rtr_l%restart_traj) CALL mm_restart_traj(tau0,velp,irec)
       ENDIF
#endif
       IF (cntl%tshop) CALL stopgm('INITRUN','NEED A RESTART FILE',& 
            __LINE__,__FILE__)
       ! Randomize cell parameters
       IF (cntl%tranc) CALL ranc(tau0)
       ! Randomization of the atomic coordinates.
       IF (cntl%tranp.AND.paral%io_parent) CALL ranp(tau0)
       CALL mp_bcast(tau0,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       ! INITIALIZATION OF WAVEFUNCTION
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL rinitwf(c0,cm,sc0,nstate,tau0,taup,rhoe,psi)
    ELSEIF (restart1%restart.AND.restart1%rnon) THEN
       IF (bsclcs.EQ.1) resthswf = .TRUE.
#if defined (__GROMOS)
       IF (cntl%tqmmm)THEN
          IF (rtr_l%restart_traj) CALL mm_restart_traj(tau0,velp,irec)
       ENDIF
#endif
       IF (restart1%rgeo.AND.paral%io_parent) CALL geofile(tau0,velp,'READ')
       ! Randomize cell parameters
       IF (cntl%tranc) CALL ranc(tau0)
       ! Randomization of the atomic coordinates.
       IF (cntl%tranp.AND.paral%parent) CALL ranp(tau0)
       CALL mp_bcast(tau0,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(velp,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       ! INITIALIZATION OF WAVEFUNCTION
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL rinitwf(c0,cm,sc0,nstate,tau0,taup,rhoe,psi)
       IF (restart1%rgeo.AND.paral%io_parent) CALL geofile(tau0,velp,'READ')
       CALL mp_bcast(tau0,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(velp,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    ELSEIF ((restart1%restart.AND. .NOT.restart1%rnon).OR.cntl%tmerge) THEN
       IF (bsclcs.EQ.1) resthswf = .FALSE.
       ! READ RESTART FILE
       IF (cntl%tmerge) THEN
          CALL MERGE(c0,cm,rhoe,psi)
       ELSE
          CALL zhrwf(1,irec,c0,cm,nstate,eigv,tau0,velp,taup,iteropt%nfi)
       ENDIF
       IF (restart1%rgeo.AND.paral%io_parent) CALL geofile(tau0,velp,'READ')
#if defined (__GROMOS)
       IF (cntl%tqmmm)THEN
          IF (rtr_l%restart_traj) CALL mm_restart_traj(tau0,velp,irec)
       ENDIF
#endif
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL ainitwf(c0,nstate,tau0,taup,rhoe,psi)
       CALL mm_dim(mm_go_mm,statusdummy)
       ! Randomize cell parameters
       IF (cntl%tranc) CALL ranc(tau0)
       ! Randomization of the atomic coordinates.
       IF (cntl%tranp.AND.paral%io_parent) CALL ranp(tau0)
       CALL mp_bcast(tau0,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(velp,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       ! Update symmetry operation
       IF (irec(irec_co).EQ.1) CALL updatsym(tau0)
       ! INITIALIZE A NEW WAVEFUNCTION
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (irec(irec_wf).EQ.0) THEN
          IF (cntl%tshop) CALL stopgm('INITRUN','NEED WAVEFUNCTIONS',& 
               __LINE__,__FILE__)
          CALL rinitwf(c0,cm,sc0,nstate,tau0,fion,rhoe,psi)
       ENDIF
       ! Randomize electronic wavefunctions around loaded ones
       IF (cntl%trane) CALL rrane(c0,cm,nstate)
       ! update cell
       IF (cntl%tprcp.OR.cntl%tpres.OR.restart1%rcell) THEN
          IF (isos1%tclust) THEN
             IF (paral%io_parent)&
                  WRITE(0,'(/,1X,A,/,1X,A,/,1X,A,/)')&
                  'INITRUN| WARNING: READ CELL FROM RESTART UNSUPPORTED',&
                  'INITRUN| WARNING: FOR CLUSTER CALCULATIONS',&
                  'INITRUN| WARNING: USING PARAMETERS FROM INPUT FILE'
          ELSE
             msglen = 9 * 8
             CALL mp_bcast_byte(metr_com, size_in_bytes_of(metr_com),parai%io_source,parai%cp_grp)
             DO i=1,3
                parm%a1(i) = metr_com%ht(1,i)
                parm%a2(i) = metr_com%ht(2,i)
                parm%a3(i) = metr_com%ht(3,i)
             ENDDO
             CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
             CALL newcell
          ENDIF
       ENDIF

#if defined (__GROMOS)
       ! (re-)check coordinates, e.g. needed when restarting from GEOMETRY file.
       IF (cntl%tqmmm.AND..NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0,cm,nstate)
       ENDIF
#endif

       CALL phfac(tau0)
       IF (corel%tinlc) CALL copot(rhoe,psi,cntl%tprcp)
       ! Orthogonalize electronic wavefunctions
       IF (irec(irec_wf).EQ.0.OR.cntl%trane.OR.(pslo_com%tivan.AND..NOT.cntl%nonort)) THEN
          CALL inq_swap(kbeg,kend,kinc)
          DO ikpt=kbeg,kend,kinc
             nkpoint=nkpbl(ikpt)
             IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,'C0')
             DO ikind=1,nkpoint
                IF (cntl%tshop) THEN
                   IF (pslo_com%tivan) CALL rnlsm(c0(:,1:sh02%nst_s0,ikind),sh02%nst_s0,&
                        ikpt,ikind,.FALSE.)
                   CALL ortho(sh02%nst_s0,c0(:,1:sh02%nst_s0,ikind),cm)
                   IF (pslo_com%tivan) CALL rnlsm(c0(:,sh02%nst_s0+1:sh02%nst_s0+sh02%nst_s1,ikind),&
                        sh02%nst_s1,&
                        ikpt,ikind,.FALSE.)
                   CALL ortho(sh02%nst_s1,c0(:,sh02%nst_s0+1:sh02%nst_s0+sh02%nst_s1,ikind),cm)
                ELSE
                   IF (pslo_com%tivan) CALL rnlsm(c0(:,1:nstate,ikind),nstate,&
                        ikpt,ikind,.FALSE.)
                   CALL ortho(nstate,c0(:,1:nstate,ikind),cm)
                ENDIF
             ENDDO
             IF (tkpts%tkblock) CALL wkpt_swap(c0,nstate,ikpt,'C0')
          ENDDO
       ENDIF
    ENDIF

    ! Read WF centers & spread from restart file
    IF (vdwl%vdwd) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr-np_low+1
       ELSE
          ipx=1
       ENDIF
       IF (paral%io_parent) THEN
          nwfc=nstate
          vdwwfl%trwannc=trwanncx(ipx)
          CALL rwannier(nwfc,tauref(:,:,:,ipx),rwann(:,:,ipx),swann(:,ipx),vdwwfl%trwannc)
          IF (.NOT.vdwwfl%trwannc) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,tauref(1,1,1,ipx),1)
          ENDIF
       ENDIF
       CALL mp_bcast(tauref(:,:,:,ipx),3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(rwann(:,:,ipx),3*nwfcx*nfragx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(swann(:,ipx),nwfcx*nfragx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(vdwwfl%trwannc,parai%io_source,parai%cp_grp)
       trwanncx(ipx)=vdwwfl%trwannc
    ENDIF
    CALL mm_dim(mm_go_mm,statusdummy)
    IF (paral%io_parent) THEN
       IF (.NOT.lqmmm%qmmm) THEN
          ! we cannot do this since the gromos/cpmd translation
          ! tables have not been set up yet.
          ! SI: To check, if forces are always available

          IF (.NOT.ALLOCATED(fion)) THEN
             ! FIXME SPECT/IR/methan-in-nosymm testcase crashes with 'not
             ! allocated' error here
             ALLOCATE(fion(1,1,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)

          ENDIF

          CALL wrgeof(tau0,fion)

          ! FIXME
          IF (SIZE(fion,1) == 1 .AND.&
               SIZE(fion,2) == 1 .AND. SIZE(fion,3) == 1) THEN
             DEALLOCATE(fion,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF

          CALL geofile(tau0,velp,'WRITE')
       ENDIF
    ENDIF
    ! Initialization of density and potential
    ! for diagonalization schemes
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (cntl%tdiag.AND..NOT.clc%classical) THEN
       IF (cntl%tshop) CALL stopgm('INITRUN','CNTL%TDIAG AND TSHOP NOT SUPPORTED'&
            ,& 
            __LINE__,__FILE__)
       IF (irec(irec_rho).EQ.0) THEN
          IF (pslo_com%tivan) THEN
             CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
          ENDIF
          IF (tkpts%tkpnt) THEN
             CALL rhoofr_c(c0,rhoe,psi(:,1),nstate)
          ELSE
             CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),nstate)
          ENDIF
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rin0,1)
       ENDIF
       IF (andr2%trand) CALL rrandd(rin0,andr2%amprd)
       CALL zeroing(rout0)!,nnr1*clsd%nlsd)
       CALL zeroing(rmix)!,nnr1*clsd%nlsd)
       ! NL-Projector Overlap Matrix
       IF (fint1%ttrot) CALL calc_alm
    ELSEIF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,psi(:,1),nstate)
       CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rin0,1)
    ENDIF
    ! Init some cntl%diis counters
    iteropt%ndisrs=0
    iteropt%ndistp=0
    CALL mm_dim(mm_revert,status)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initrun

END MODULE initrun_driver
