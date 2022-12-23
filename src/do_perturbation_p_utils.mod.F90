#include "cpmd_global.h"

MODULE do_perturbation_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cotr,                            ONLY: cotc0
  USE detdof_utils,                    ONLY: detdof
  USE eicalc_utils,                    ONLY: eicalc
  USE eigensystem_p_utils,             ONLY: eigensystem_p
  USE epr_p_utils,                     ONLY: epr_p,&
                                             give_scr_epr
  USE error_handling,                  ONLY: stopgm
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_p_utils,                  ONLY: give_scr_forces_p
  USE fukui_p_utils,                   ONLY: fukui_p
  USE hardness_p_utils,                ONLY: hardness_p
  USE interaction_manno_p_utils,       ONLY: interaction_manno_p
  USE interaction_p_utils,             ONLY: give_scr_interaction,&
                                             interaction_p
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lanc_phon_p_utils,               ONLY: lanc_phon_p
  USE machine,                         ONLY: m_flush
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast
  USE nmr_p_utils,                     ONLY: give_scr_nmr,&
                                             nmr_p
  USE nuclear_p_utils,                 ONLY: nuclear_p
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE parac,                           ONLY: parai,&
                                             paral
  USE pert_kpoint_p_utils,             ONLY: give_scr_pert_kpoint_p,&
                                             pert_kpoint_p
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phfac_utils,                     ONLY: phfac
  USE phonons_p_utils,                 ONLY: phonons_p
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE raman_p_utils,                   ONLY: raman_p
  USE response_pmod,                   ONLY: dmbi,&
                                             response1,&
                                             response2,&
                                             rho0,&
                                             vofrho0,&
                                             vxc0
  USE rhoofr_p_utils,                  ONLY: give_scr_rhoofr_p
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rhopri_utils,                    ONLY: rhopri
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE store_types,                     ONLY: rout1
  USE symm,                            ONLY: symmi,&
                                             symmt
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE updwf_p_utils,                   ONLY: give_scr_updwf_p
  USE voa_p_utils,                     ONLY: voa_p
  USE vofrho_utils,                    ONLY: vofrho
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: do_perturbation
  PUBLIC :: give_scr_perturbation

CONTAINS

  ! ==================================================================
  SUBROUTINE do_perturbation(c0,c1,nstate)
    ! ==--------------------------------------------------------------==
    ! The main routine that controls the linear response calculations.
    ! INPUT: C0: Converged ground state wavefunctions
    ! NSTATE
    ! OUTPUT: C1: LINEAR RESPONSE wavefunctions
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate)
    COMPLEX(real_8), TARGET                  :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'do_perturbation'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), psi(:,:)
    COMPLEX(real_8), POINTER                 :: c0_3d_p(:,:,:)
    INTEGER                                  :: ierr, il_psi_1d, il_psi_2d, &
                                                il_rhoe_1d, il_rhoe_2d, lscr
    REAL(real_8), ALLOCATABLE                :: scr(:), tscr(:,:,:), z11(:,:)
    REAL(real_8), ALLOCATABLE, TARGET        :: rhoe(:,:)
    REAL(real_8), POINTER                    :: drhoe(:,:)

! ==--------------------------------------------------------------==
! internal variables
! ==--------------------------------------------------------------==
! Wavefunctions, densities and related:
! potentials and related:
! others & scratch:
! (nstate,nstate)
! ==--------------------------------------------------------------==

    IF (symmt%tpgauto .OR. symmt%tmsym .OR. symmi%indpg.NE.0) THEN
       CALL stopgm('do_perturb',&
            'POINT GROUP symmetry not implemented.',& 
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tlsd) THEN
       IF (response1%traman .OR. response1%tlanphon .OR.&
            response1%tnmr .OR. response1%tkpert .OR. response1%tinteraction .OR.&
            response1%teigensystem .OR. response1%thardness) THEN
          IF (paral%parent.AND.paral%io_parent) THEN
             WRITE (6,*) 'LSD calculation'&
                  //' not yet tested. Please complain to'
             WRITE (6,*) 'D.S. (sebastia@mpip-mainz.mpg.de)'
          ENDIF
       ENDIF
       IF (response1%tlanphon .OR.&
            response1%tkpert .OR. response1%tinteraction .OR.&
            response1%teigensystem .OR. response1%thardness) THEN
          CALL stopgm('do_perturb','LSD: drugs are illegal'&
               //'in the U.S. while guns are not.',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.lqmmm%qmmm) THEN
       WRITE (6,*) '!! QM/MM + LINEAR RESPONSE: This combination is'&
            // ' not yet fully tested.  !!'
       WRITE (6,*) '!! QM/MM + LINEAR RESPONSE: Please use with'&
            // 'special care.              !!'
    ENDIF
    IF (tkpts%tkpnt .AND. .NOT. response1%tkpert)&
         CALL stopgm('do_perturb','K POINTS NOT IMPLEMENTED.',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm('do_perturb','VDB USPP NOT IMPLEMENTED.',& 
         __LINE__,__FILE__)
    IF (ropt_mod%tesr .OR. iteropt%iesr .NE.0)&
         CALL stopgm('do_perturb','TESR/IESR NOT IMPLEMENTED.',& 
         __LINE__,__FILE__)
    IF (isos3%ps_type.EQ.1)&
         CALL stopgm('do_perturb','HOCKNEY NOT IMPLEMENTED.',& 
         __LINE__,__FILE__)
    IF (paral%io_parent.AND.isos1%tclust .AND. .NOT.lqmmm%qmmm) THEN
       WRITE (6,*) '!! BE CAREFUL WITH ISOLATED MOLECULES        !!'
       WRITE (6,*) '!! NOT FULLY TESTED WITH PERTURBATION THEORY !!'
    ENDIF
    IF (lspin2%tlse) THEN
       CALL stopgm('do_perturb','LSE NOT IMPLEMENTED.',& 
            __LINE__,__FILE__)
    ENDIF
    response1%response_running = .TRUE. ! This is for WRENER
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,'(" *",19X,A,20X,"*")') ' RESPONSE CALCULATIONS'
       WRITE(6,'(1X,64("*"),/)')
       ! Format: 65 Spalten
81     FORMAT ("* RESPONSE TYPE:        ",a35,5x,"*")
82     FORMAT ("* PRECONDITIONER TYPE:  ",a35,5x,"*")
83     FORMAT ("* PRECONDITIONER CUTOFF:",f10.5,30x,"*")
84     FORMAT ("* WFN CONVERGENCE:      ",e12.3,28x,"*")
85     FORMAT ("* ALGORITHM:            ",a35,5x,"*")
86     FORMAT ("* CG: ANALYTIC STEPS:   ",i4,36x,"*" /,&
            "* CG: LENGTH PREFACTOR: ",f10.5,30x,"*")
87     FORMAT ("* ORTHO-SCHEME:         ",a35,5x,"*")
88     FORMAT ("* GND STATE WFN:        ",a35,5x,"*")
89     FORMAT ("* ",22x,a35,5x,"*")
90     FORMAT (56("*"),"RESPONSE*")
       WRITE (6,*)
       WRITE (6,90)
       IF (response1%traman)   WRITE (6,81)'RAMAN spectra                   '
       IF (response1%phonon)   WRITE (6,81)'PHONON frequencies              '
       IF (response1%tlanphon) WRITE (6,81)'PHONON frequencies (Lanczos)    '
       IF (response1%tnmr)     WRITE (6,81)'NMR chemical shifts             '
       IF (response1%tepr)     WRITE (6,81)'EPR G tensor / Hyperfine tensor '
       IF (response1%tvoa)     WRITE (6,81)'Vibrational Optical Activity    '
       IF (response1%tkpert)   WRITE (6,81)'K point effect as linear response'
       IF (response1%tinteraction) WRITE (6,81)'INTERACTION of Wannier orbitals '
       IF (response1%teigensystem)&
            WRITE (6,81)          'EIGENSYSTEM or softness states  '
       IF (response1%thardness)&
            WRITE (6,81)          'HARDNESS kernel diagonalization '
       IF (response1%preconditioner_p.EQ.1) THEN
          WRITE (6,82)            'smooth (ds)                     '
       ELSE IF (response1%preconditioner_p.EQ.2) THEN
          WRITE (6,82)            'tight (ap)                      '
       ELSE IF (response1%preconditioner_p.EQ.3) THEN
          WRITE  (6,82) 'State dependent + DS        '
       ELSE IF (response1%preconditioner_p.EQ.4) THEN
          WRITE  (6,82) 'State dependent + AP        '
       ELSE IF (response1%preconditioner_p.EQ.5) THEN
          WRITE (6,82)  'precondition with cutoff    '
       ENDIF
       WRITE (6,83) response2%hthrs_p
       WRITE (6,84) response2%tolog_p
       IF (response1%tkeeprealspacewfn) THEN
          WRITE (6,88) 'KEPT in REAL SPACE'
       ELSE
          WRITE (6,88) 'reciprocal space only'
       ENDIF
       IF (response1%pcg_p) THEN
          IF (response1%prec_p) THEN
             WRITE (6,85) 'preconditioned conj gradient'
          ELSE
             WRITE (6,85) 'bare conj gradient        '
          ENDIF
          IF (response1%tpolak_p) THEN
             WRITE (6,89) 'Polak-Ribiere variant     '
          ELSE
             WRITE (6,89) 'Fletcher-Reeves variant   '
          ENDIF
          IF (response1%pcgmin_p) THEN
             WRITE (6,86) response1%cg_analytic,response2%cg_factor
          ELSE
             WRITE (6,89) 'using FIXED STEP LENGTH   '
          ENDIF
       ELSEIF (response1%tsde_p) THEN
          IF (response1%prec_p) THEN
             WRITE (6,85) 'preconditioned steepest descent'
          ELSE
             WRITE (6,85) 'bare steepest descent     '
          ENDIF
       ELSEIF (response1%diis_p) THEN
          IF (response1%prec_p) THEN
             WRITE (6,85) 'PRECONDITIONED DIIS     '
          ELSE
             WRITE (6,85) 'BARE DIIS OPTIMIZATION    '
          ENDIF
       ENDIF
       IF (response1%tinteraction) THEN
          IF (dmbi%tlinearscaling) THEN
             WRITE (6,87) 'MOLECULAR'
          ELSE
             WRITE (6,87) 'complete (std)'
          ENDIF
          IF (dmbi%fastdexc) WRITE (6,89) 'Exc second derivative using LDA '
          IF (dmbi%tatomicwavefunction)&
               WRITE (6,89) 'Atomic wavefunction for W0      '
          IF (dmbi%tatomicrho)&
               WRITE (6,89) 'Atomic density for Rho0         '
          IF (.NOT.dmbi%torthog_wannier)&
               WRITE (6,89) 'Using NON ORTHOGONAL W0         '
          IF (dmbi%tsimple_model)&
               WRITE (6,89) 'Using SIMPLE interaction model  '
          IF (dmbi%tmanno)&
               WRITE (6,89) 'Using OLD MANNOPT code (serial) '
       ENDIF
       WRITE (6,90)
    ENDIF                     ! Parent.
    ! ==--------------------------------------------------------------==
    ! First, ALLOCATE MEMORY which is needed for any kind of perturbation.
    ! Then, the specific perturbation subroutine is called.
    ! ==--------------------------------------------------------------==
    ! Wavefunctions, densities and related:
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d, &
         il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d)
    ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(psi)!,SIZE(psi))

    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! 2D size is not known here
    CALL zeroing(rhoe)!,size(rhoe))

    IF ((response1%tnmr .OR. response1%tkpert).OR.response1%tepr) THEN
       rho0 => rhoe
       drhoe => rhoe
       vxc0 => rhoe

       response1%trho1 = .FALSE.
    ELSE
       ALLOCATE(rho0(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! 2D size is not known here
       CALL zeroing(rho0)!,size(rho0))

       ALLOCATE(drhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! 2D size is not known here
       CALL zeroing(drhoe)!,size(drhoe))

       ALLOCATE(vxc0(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(vxc0)!,size(vxc0))
       response1%trho1 = .TRUE.
    ENDIF

    ! potentials and related:
    ALLOCATE(vofrho0(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vofrho0)!,size(vofrho0))

    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop)!,SIZE(eirop))

    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eivps)!,SIZE(eivps))


    ! others & scratch:
    ALLOCATE(z11(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(z11)!,nstate*nstate)

    IF (.NOT.cntl%tmdfile) THEN
       ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL give_scr_perturbation(lscr,tag,nstate)

    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(scr)!,lscr)
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    CALL eicalc(eivps,eirop)

    ! the potential v[n0]:
    IF (.NOT. response1%tinteraction) & ! in the INTERACTION case, the wfns
         THEN                  ! are restarted lateron
       CALL rhoofr(c0,rho0,psi(:,1),nstate)
       IF (cntl%ttau) CALL stopgm("inr_dr","no tau functionals",& 
            __LINE__,__FILE__)
       CALL dcopy(fpar%nnr1*clsd%nlsd, rho0, 1, vofrho0, 1)
       CALL vofrho(tau0,fion,vofrho0,psi,.FALSE.,.FALSE.)
    ENDIF

    IF (paral%parent) CALL prmem('  DO_PERTURBATION')
    ! calculate the lagrange parameters z11
    IF ((.NOT. response1%tinteraction) .AND. (.NOT. (response1%tnmr.OR.response1%tepr))) THEN
       CALL lag_mult(c0,c1,psi,rhoe,z11,nstate)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (response1%tdummyatom) THEN
       CALL nuclear_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tdummyatom_ref) THEN
       CALL nuclear_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
       response1%tdummyatom_ref=.FALSE.
       IF (rout1%rhoout) THEN
#if defined( _HASNT_F08_POINTER_REMAPPING )
          CALL stopgm(procedureN,'compiler needs to support pointer remapping!',&
               __LINE__,__FILE__)
#else
          c0_3d_p(1:SIZE(c0,1),1:SIZE(c0,2),1:1) => c0
#endif
          CALL rhopri(c0_3d_p,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       ENDIF
    ELSEIF (response1%tdummyatom_tec) THEN
       CALL nuclear_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%toacp4f) THEN
       CALL nuclear_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%traman) THEN
       IF (paral%io_parent)&
            WRITE(6,*) '*** RAMAN calculation'
       ! PAGLIAI/CINECA ADDED TAU0 TO THE CALL 
       CALL raman_p(tau0,c0,c1,psi(:,1),rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tinteraction) THEN
       IF (dmbi%tmanno) THEN
          CALL interaction_manno_p(c0,c1,psi,rhoe,drhoe,&
               eirop,eivps,z11,nstate)
       ELSE
          CALL interaction_p(c0,c1,psi,rhoe,drhoe,&
               eirop,eivps,z11,nstate)
       ENDIF
    ELSEIF (response1%phonon) THEN
       IF (paral%io_parent) THEN
          WRITE(6,'(2x,a)') '*** PHONON calculation '//&
               '(classical diagonalization). '
          IF (response1%projout.AND.response1%rotout) THEN
             WRITE(6,'(2x,a)') 'TRANSLATIONS AND ROTATIONS HAVE'//&
                  ' BEEN DISCARDED'
          ELSEIF (response1%projout.AND..NOT.response1%rotout) THEN
             WRITE(6,'(2x,a)') 'TRANSLATIONS ONLY HAVE BEEN'//&
                  ' DISCARDED'
          ELSE
             WRITE(6,'(2x,a)') 'TRIVIAL EIGENMODES ARE INCLUDED'
          ENDIF
       ENDIF
       CALL phonons_p(c0,c1,psi(:,1),rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tlanphon) THEN
       IF (paral%parent) CALL detdof(tau0,tscr)
       CALL mp_bcast(cotc0%nodim,parai%source,parai%allgrp)

       IF (paral%io_parent) THEN
          WRITE(6,'(2x,a)') '*** PHONON calculation  '//&
               '(LANCZOS iterative method). '
          IF (response1%projout.AND.response1%rotout) THEN
             WRITE(6,'(2x,a)') 'TRANSLATIONS AND ROTATIONS HAVE'//&
                  ' BEEN DISCARDED'
          ELSEIF (response1%projout.AND..NOT.response1%rotout) THEN
             WRITE(6,'(2x,a)') 'TRANSLATIONS ONLY HAVE BEEN'//&
                  ' DISCARDED'
          ELSE
             WRITE(6,'(2x,a)') 'TRIVIAL EIGENMODES ARE INCLUDED'
          ENDIF
       ENDIF
       CALL lanc_phon_p(c0,c1,psi(:,1),rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tkpert) THEN
       IF (paral%io_parent) WRITE(6,*)' *** KPERT CALCULATION ***'
       CALL pert_kpoint_p(c0,c1,psi(:,1),rhoe,eirop,eivps,z11,nstate)
    ELSEIF (response1%teigensystem) THEN
       IF (paral%io_parent) WRITE(6,*) ' *** EIGENSYSTEM CALCULATION ***'
       CALL eigensystem_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%thardness) THEN
       IF (paral%io_parent) WRITE(6,*) ' *** HARDNESS CALCULATION ***'
       CALL hardness_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tfukui) THEN
       IF (paral%io_parent) WRITE(6,*)' *** FUKUI FUNCTIONS CALCULATION ***'
       CALL m_flush(6)
       CALL fukui_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tnmr) THEN
       CALL nmr_p(c0,c1,psi,rhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tepr) THEN
       CALL epr_p(c0,c1,psi,rhoe,&
            eirop,eivps,z11,nstate)
    ELSEIF (response1%tvoa) THEN
       CALL voa_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! DE-ALLOCATE MEMORY
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr,tscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)  ! last allocated = first de-allocated
    IF (.NOT. (response1%tinteraction.OR.cntl%tmdfile)) THEN
       DEALLOCATE(fion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(z11,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vofrho0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (.NOT. (response1%tepr .OR. response1%tnmr .OR. response1%tkpert)) THEN
       DEALLOCATE(vxc0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(drhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rho0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    response1%response_running = .FALSE.
    RETURN
  END SUBROUTINE do_perturbation
  ! ==================================================================
  SUBROUTINE give_scr_perturbation(l_max,tag,nstate)
    INTEGER                                  :: l_max
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate

    INTEGER :: l_forcedr, l_opeigr, l_updwf_p, lepr, lforce1, linteraction, &
      lnmr, lpert_kpoint_p, lrho, lrhoofr, lrnlsm, lrscpot

! ==--------------------------------------------------------------==

    l_max=nstate
    ropt_mod%calste=.FALSE.
    CALL give_scr_nmr(lnmr,tag)
    CALL give_scr_epr(lepr,tag)
    CALL give_scr_interaction(linteraction,tag,nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    CALL give_scr_rhoofr(lrho,tag)
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    CALL give_scr_rhoofr_p(lrhoofr,tag)
    CALL give_scr_forces_p(lforce1,tag,nstate)
    CALL give_scr_forcedr(l_forcedr,tag,nstate,.FALSE.,.FALSE.)
    CALL give_scr_opeigr(l_opeigr,tag,nstate)
    CALL give_scr_updwf_p(l_updwf_p,tag,nstate)
    CALL give_scr_pert_kpoint_p(lpert_kpoint_p,tag,nstate)

    l_max = MAX(linteraction,lrho,lrscpot,lrnlsm,lrhoofr,l_max,&
         lforce1,lnmr,l_forcedr,l_opeigr,l_updwf_p,lepr,&
         lpert_kpoint_p)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_perturbation
  ! ==================================================================

END MODULE do_perturbation_p_utils
