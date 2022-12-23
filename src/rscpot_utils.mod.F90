MODULE rscpot_utils
  USE cdft_utils,                      ONLY: cdft_forces
  USE cnst_dyn,                        ONLY: en_harm_locs,&
                                             lmeta
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com,&
                                             ener_d
  USE epr_efg_utils,                   ONLY: save_rho
  USE hubbardu,                        ONLY: c2u0,hubbu
  USE hubbardu_utils,                  ONLY: hubbardUcorrection
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE kpts,                            ONLY: tkpts
  USE mm_input,                        ONLY: iqmmm,&
                                             lqmmm
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop5
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlfor_utils,                    ONLY: rnlfor
  USE rnlrh_utils,                     ONLY: rnlrh
  USE ropt,                            ONLY: iteropt
  USE rpiiint_utils,                   ONLY: rpiiint
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE stress_utils,                    ONLY: give_scr_stress,&
                                             stress
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdw_utils,                       ONLY: vdw,&
                                             vdw_wf
  USE vdwcmod,                         ONLY: &
       empvdwi, empvdwr, idvdw, ivdw, jvdw, vdwbe, vdwi, vdwl, vdwr, vdwrm, &
       vdwst, vdwwfr
  USE vofrho_utils,                    ONLY: give_scr_vofrho,&
                                             vofrho
  USE vofrhoc_utils,                   ONLY: ctfor,&
                                             vcomp
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rscpot
  PUBLIC :: give_scr_rscpot

CONTAINS

  ! ==================================================================
  SUBROUTINE rscpot(c0,tau0,fion,rhoe,psi,&
       tfor,tstress,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE POT |PSI> AND THE ATOMIC FORCES (POT PART)         ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   C0(NGW,NSTATE) ELECTRONIC WAVEFUNCTIONS                    ==
    ! ==   TAU0           ATOMIC COORDINATES                          ==
    ! ==   TFOR           .TRUE. CALCULATES THE ATOMIC FORCES         ==
    ! ==   TSTRESS        .TRUE. CALCULATES THE STRESS TENSOR         ==
    ! ==   NSTATE         NUMBER OF STATES                            ==
    ! ==   NKPOINT        NUMBER OF K POINTS                          ==
    ! == OUTPUT:                                                      ==
    ! ==   FION      ATOMIC FORCES                                    ==
    ! ==   RHOE      ELECTRONIC POTENTIAL                             ==
    ! ==   PSI       USED FOR FFT                                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    LOGICAL                                  :: tfor, tstress
    INTEGER                                  :: nstate, nkpoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'rscpot'

    INTEGER                                  :: isub, nfto, nwfc
    LOGICAL                                  :: dorho
    REAL(real_8)                             :: e_mm, ebindnovdw, ebindvdw, &
                                                etotnovdw

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==
    ! == Calculate NON-LOCAL PARTS of the energy and forces, which    ==
    ! == depend on the wavefunctions themselves, and not on the       ==
    ! == electronic density.                                          ==
    ! ==--------------------------------------------------------------==
    CALL rnlrh(ener_com%enl,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == Compute the residual part of the ion-ion interaction (ESR)   ==
    ! == due to the overlap of the smeared ionic charge densities     ==
    ! == corresponding to different atoms.                            ==
    ! ==--------------------------------------------------------------==
    nfto=3*maxsys%nax*maxsys%nsx
    CALL zeroing(fion)!,nfto)
    CALL rpiiint(ener_com%esr,tau0,fion,iteropt%iesr,tfor)
    ! ==--------------------------------------------------------------==
    ! == Compute the electronic density and kinetic energy            ==
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,psi(:,1),nstate)
    ELSE
       dorho=.NOT.lqmmm%qmmm .OR. (lqmmm%qmmm.AND.cntl%bsymm)&
            .OR.(lqmmm%qmmm .AND. iqmmm%coupl_model.EQ.0)
       IF (dorho) THEN
          CALL rhoofr(c0,rhoe,psi(:,1),nstate)
       ENDIF
    ENDIF

    IF ((prop5%teprefg.AND.cntl%tlsd).OR.prop5%tefg) THEN
       CALL save_rho(rhoe)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Calculate the Stress Tensor                                  ==
    ! ==--------------------------------------------------------------==
    IF (tstress) THEN
       CALL stress(c0,tau0,crge%f,psi(:,1),nstate)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Compute the LOCAL PARTS of the energy and the forces, which  ==
    ! == depend exclusively on the electronic density.                ==
    ! ==--------------------------------------------------------------==
    IF (cntl%cdft.AND.tfor) CALL cdft_forces(fion,tau0,rhoe(:,1))
    CALL vofrho(tau0,fion,rhoe,psi,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! == Force on the ions from nonlocal pseudopotential              ==
    ! ==--------------------------------------------------------------==
    IF (tfor) THEN
       CALL rnlfor(fion,crge%f,wk,nstate,nkpoint)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions for empirical van der Waals     ==
    ! ==--------------------------------------------------------------==
    vdwr%evdw=0._real_8
    CALL zeroing(vdwr%devdw)!,6)
    IF (vdwl%vdwc) THEN
       CALL vdw(tau0,empvdwi%nvdw,idvdw,ivdw,jvdw,vdwst,vdwrm,vdwbe,&
            empvdwr%vdweps,empvdwr%s6grim,vdwi%nxvdw,vdwi%nyvdw,vdwi%nzvdw,vdwr%evdw,fion,vdwr%devdw)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions for Wannier-based van der Waals ==
    ! ==--------------------------------------------------------------==
    IF (vdwl%vdwd) THEN
       vdwr%evdw=0._real_8
       CALL zeroing(vdwr%devdw)!,6)
       nwfc=nstate
       IF (tfor) CALL vdw_wf(tau0,fion,nwfc)
    ENDIF
    !==--------------------------------------------------------------==
    !== Energy and forces from DFT+U Hamiltonian                     ==
    !== (forces on wavefunctions are stored in C2U0)                 ==
    !==--------------------------------------------------------------==
      hubbu%ehub=0.D0
      IF(cntl%thubb)THEN
        IF(.not.allocated(c2u0))allocate(C2U0(ncpw%NGW,NSTATE))
        CALL hubbardUcorrection(C0,C2U0,TAU0,FION,NSTATE,PSI(:,1),TFOR,0)
      END IF 
    ! ==--------------------------------------------------------------==
    ! == QM/MM of Roethlisberger group                                ==
    ! ==--------------------------------------------------------------==
    e_mm=0._real_8
#if defined (__QMECHCOUPL)
    ! mb      Unknown origin: the routine mm_cpmd_MM_energy does not exist!
    ! mb      call mm_cpmd_MM_energy(E_MM)
    IF (tfor) THEN
       CALL mm_cpmd_add_mm_forces_f77(fion)
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    ! == TOTAL ENERGY                                                 ==
    ! ==--------------------------------------------------------------==
    ener_com%eeig = 0._real_8
    ener_com%etot = ener_com%ekin + ener_com%eht + ener_com%epseu + ener_com%enl + ener_com%exc + e_mm&
         + vdwr%evdw + hubbu%ehub  ! Empirical or WFC-ab initio van der Waals correction
    ! Sum ETOT,EKIN,EPSEU,ENL,EHT,EHEP,EHEE,EHII,EXC,VXC,EGC
    CALL mp_sum(ener_com% etot,parai%allgrp)
    CALL mp_sum(ener_com% ekin,parai%allgrp)
    CALL mp_sum(ener_com% epseu,parai%allgrp)
    CALL mp_sum(ener_com% enl,parai%allgrp)
    CALL mp_sum(ener_com% eht,parai%allgrp)
    CALL mp_sum(ener_com% ehep,parai%allgrp)
    CALL mp_sum(ener_com% ehee,parai%allgrp)
    CALL mp_sum(ener_com% ehii,parai%allgrp)
    CALL mp_sum(ener_com% exc,parai%allgrp)
    CALL mp_sum(ener_com% vxc,parai%allgrp)
    CALL mp_sum(ener_com% egc,parai%allgrp)
    CALL mp_sum(hubbu%ehub   ,parai%allgrp)

    ! _vdwWF
    IF (vdwl%vdwd.AND.tfor.AND.vdwwfr%enmonomer.NE.0.0_real_8) THEN
       etotnovdw=ener_com%etot-vdwr%evdw
       ebindnovdw=(etotnovdw-vdwwfr%enmonomer)*27.211607_real_8
       ebindvdw=(ener_com%etot-vdwwfr%enmonomer)*27.211607_real_8
       IF (paral%io_parent) THEN
          WRITE(6,737) ebindnovdw,ebindvdw,vdwr%evdw*27.211607_real_8
737       FORMAT(1x,'-> ETOT, ETOT+vdW, vdW (eV):',3e12.3)
       ENDIF
    ENDIF
    ! _vdwWF
    ! == Metadynamics : spin den localization                         ==
    IF (lmeta%tlocalizespin) THEN
       ener_com%etot = ener_com%etot + en_harm_locs
    ENDIF
    IF (lspin2%tlse .AND. lspin2%tcas22) THEN
       IF (ABS(e_mm).GT.1.e-10_real_8) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) " ******************************************"
             IF (paral%io_parent)&
                  WRITE(6,*) " WARNING: MM ENERGY CONTRIBUTION NOT TESTED"
             IF (paral%io_parent)&
                  WRITE(6,*) " ******************************************"
          ENDIF
       ENDIF
       CALL mp_sum(ener_c%etot_a,parai%allgrp)
       CALL mp_sum(ener_c%ekin_a,parai%allgrp)
       CALL mp_sum(ener_c%epseu_a,parai%allgrp)
       CALL mp_sum(ener_c%enl_a,parai%allgrp)
       CALL mp_sum(ener_c%eht_a,parai%allgrp)
       CALL mp_sum(ener_c%exc_a,parai%allgrp)
       CALL mp_sum(ener_c%eext_a,parai%allgrp)
       CALL mp_sum(ener_c%etot_2,parai%allgrp)
       CALL mp_sum(ener_c%ekin_2,parai%allgrp)
       CALL mp_sum(ener_c%epseu_2,parai%allgrp)
       CALL mp_sum(ener_c%enl_2,parai%allgrp)
       CALL mp_sum(ener_c%eht_2,parai%allgrp)
       CALL mp_sum(ener_c%exc_2,parai%allgrp)
       CALL mp_sum(ener_c%eext_2,parai%allgrp)
       CALL mp_sum(ener_c%etot_ab,parai%allgrp)
       CALL mp_sum(ener_c%ekin_ab,parai%allgrp)
       CALL mp_sum(ener_c%epseu_ab,parai%allgrp)
       CALL mp_sum(ener_c%enl_ab,parai%allgrp)
       CALL mp_sum(ener_c%eht_ab,parai%allgrp)
       CALL mp_sum(ener_c%exc_ab,parai%allgrp)
       CALL mp_sum(ener_c%eext_ab,parai%allgrp)

       CALL mp_sum(ener_d%etot_b,parai%allgrp)
       CALL mp_sum(ener_d% ecas,parai%allgrp)
       CALL mp_sum(ener_d%etot_t,parai%allgrp)

       ener_c%etot_a = ener_c%ekin_a+ener_c%epseu_a+ener_c%enl_a+ener_c%eht_a+ener_com%esr-ener_com%eself+ener_c%exc_a+&
            (e_mm-ener_com%eext+ener_c%eext_a)
       ener_c%etot_2 = ener_c%ekin_2+ener_c%epseu_2+ener_c%enl_2+ener_c%eht_2+ener_com%esr-ener_com%eself+ener_c%exc_2+&
            (e_mm-ener_com%eext+ener_c%eext_2)
       ener_d%etot_b = ener_com%etot
       ener_c%etot_ab = ener_c%ekin_ab+ener_c%epseu_ab+ener_c%enl_ab+ener_c%eht_ab+ener_c%exc_ab+ener_c%eext_ab
       IF (ener_d%etot_b .LT. ener_c%etot_a) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) " ******************************************"
             IF (paral%io_parent)&
                  WRITE(6,*) " WRONG ORDER OF STATES: NOW SWITCH A AND B"
             IF (paral%io_parent)&
                  WRITE(6,*) " ******************************************"
             CALL dswap(2*ncpw%ngw,c0(1,clsd%ialpha),1,c0(1,clsd%ibeta),1)
          ENDIF
       ENDIF
       ener_d%casang = 0.5_real_8*ATAN(2._real_8*ener_c%etot_ab/(ener_c%etot_a-ener_d%etot_b))
       ener_d%ecas   = SIN(ener_d%casang)**2 * ener_c%etot_a + COS(ener_d%casang)**2 * ener_d%etot_b -&
            2._real_8*SIN(ener_d%casang)*COS(ener_d%casang) * ener_c%etot_ab
       ener_d%ecas   = ener_d%ecas - ener_d%etot_b
       ener_com%etot   = ener_d%etot_b + ener_d%ecas
       ener_d%etot_t = ener_d%etot_b - ener_com%exc + ener_d%etot_t
       ! Now get the final potentials and forces
       CALL vcomp(rhoe(:,1),rhoe(:,4))
       IF (tfor) THEN
          CALL ctfor(fion)
       ENDIF
    ELSEIF (lspin2%tlse .AND. lspin2%tpenal) THEN
       CALL mp_sum(ener_c%ekin_ab,parai%allgrp)
       CALL mp_sum(ener_c%eht_ab,parai%allgrp)
       ener_c%etot_ab = ener_c%ekin_ab + ener_c%eht_ab
       ener_com%etot = ener_com%etot + ener_c%etot_ab*ener_c%etot_ab
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE rscpot
  ! ==================================================================
  SUBROUTINE give_scr_rscpot(lrscpot,tag,tstress)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrscpot
    CHARACTER(len=30)                        :: tag
    LOGICAL                                  :: tstress

    INTEGER                                  :: lrhoofr, lstress, lvofrho

! ==--------------------------------------------------------------==

    CALL give_scr_rhoofr(lrhoofr,tag)
    IF (tstress) THEN
       CALL give_scr_stress(lstress,tag)
    ELSE
       lstress=0
    ENDIF
    CALL give_scr_vofrho(lvofrho,tag)
    lrscpot=MAX(lrhoofr,lstress,lvofrho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rscpot
  ! ==================================================================

END MODULE rscpot_utils
