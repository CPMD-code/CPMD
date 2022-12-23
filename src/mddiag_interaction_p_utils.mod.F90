MODULE mddiag_interaction_p_utils
  USE cppt,                            ONLY: gk
  USE ddip,                            ONLY: lenbk
  USE do_perturbation_p_utils,         ONLY: give_scr_perturbation
  USE eicalc_utils,                    ONLY: eicalc
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_driver,                   ONLY: forces
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: localize
  USE lowdin_utils,                    ONLY: lowdin
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phfac_utils,                     ONLY: phfac
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE response_pmod,                   ONLY: &
       dmbi, nmr_options, response1, response2, rho0, vofrho0, wc_array1, &
       wc_array2, wc_velocity
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE utils,                           ONLY: nxxfun
  USE vofrho_utils,                    ONLY: vofrho
  USE wann,                            ONLY: wannl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mddiag_interaction_p
  !public :: store_wc
  !public :: calculate_wc_velocity
  !public :: move_all_wannier

CONTAINS

  SUBROUTINE mddiag_interaction_p(nstate,c0,c1,cr,sc0,cscr,vpp,eigv,&
       rhoe,psi,&
       tau0,velp,taui,fion,ifcalc,&
       irec,tfor,tinfo)
    ! ==================================================================
    ! ** THE FUNCTIONALITY OF THIS ROUTINE IS IDENTICAL TO THAT OF    **
    ! ** FORCES_DIAG, except that the wfn optimization is done        **
    ! ** in the perturbation theory framework.                        **
    ! **                                                              **
    ! ** This routine is in fact an interface between MD_DIAG         **
    ! ** and INTERACTION_P. Most of the arguments are NOT USED        **
    ! ** but kept for the sake of compatibility to FORCES_DIAG.       **
    ! ==--------------------------------------------------------------==
    ! == ENTER WITH THE IONIC POSITIONS IN TAU0 AND AN INPUT GUESS    ==
    ! == FOR THE DENSITY IN RIN0, THIS ROUTINE RETURNS THE CONVERGED  ==
    ! == DENSITY, WAVEFUNCTIONS AND IONIC FORCES.                     ==
    ! ==--------------------------------------------------------------==
    ! == TAU0: ATOMIC POSITION                                        ==
    ! == VELP and TAUI are necessary for RESTART FILE                 ==
    ! == FION: IONIC FORCES                                           ==
    ! == IFCALC: total number of iterations                           ==
    ! ** C0:   Initial guess for the wavefunction.                    **
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate,1), &
                                                c1(ncpw%ngw,nstate), cr(*), &
                                                sc0(*), cscr(*)
    REAL(real_8)                             :: vpp(*), eigv(*), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(*), &
                                                taui(*), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER :: procedureN = 'mddiag_interaction_p'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: dummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), eirop(:), eivps(:), &
                                                h1psi0(:), sct(:,:)
    INTEGER                                  :: ierr, il_rhoe_1d, il_rhoe_2d, &
                                                istate, lscr, nc
    REAL(real_8)                             :: e0, e2
    REAL(real_8), ALLOCATABLE                :: drhoe(:,:), scr(:), z11(:,:)

! Variables
! PT variables:
! potentials and related:
! others & scratch:
! Scratch arrays for the localisation
! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
! DMB E2 calc
! DMB E2 calc
! ==--------------------------------------------------------------==
! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST

    IF (dmbi%inter_pt_firstcall) THEN
       ! testing WCpredict and fastdex cflag
       IF (paral%io_parent)&
            WRITE(6,*) 'PROC',parai%me,'WCpredict =',dmbi%wcpredict
       IF (paral%io_parent)&
            WRITE(6,*) 'PROC',parai%me,'FASTDEXC =',dmbi%fastdexc
    ENDIF
    ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",19X,A,20X,"*")') ' RESPONSE CALCULATIONS '
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"),/)')
       ! Format: 65 Spalten
81     FORMAT ("* RESPONSE TYPE:        ",a35,5x,"*")
82     FORMAT ("* PRECONDITIONER TYPE:  ",a35,5x,"*")
83     FORMAT ("* PRECONDITIONER CUTOFF:",f10.5,30x,"*")
84     FORMAT ("* WFN CONVERGENCE:      ",e12.3,28x,"*")
85     FORMAT ("* ALGORITHM:            ",a35,5x,"*")
86     FORMAT ("* CG: ANALYTIC STEPS:   ",i4,36x,"*" /, &
            "* CG: LENGTH PREFACTOR: ",f10.5,30x,"*")
89     FORMAT ("* ",22x,a35,5x,"*")
90     FORMAT (56("*"),"RESPONSE*")

       response1%tinteraction = .TRUE.

       ! print out information on first call
       IF (dmbi%inter_pt_firstcall) THEN
          IF (paral%io_parent)&
               WRITE (6,*)
          IF (paral%io_parent)&
               WRITE (6,90)
          IF (paral%io_parent)&
               WRITE (6,81)'INTERACTION of Wannier orbitals '

          IF (response1%preconditioner_p.EQ.1) THEN
             IF (paral%io_parent)&
                  WRITE (6,82) 'smooth (ds)                     '
          ELSE IF (response1%preconditioner_p.EQ.2) THEN
             IF (paral%io_parent)&
                  WRITE (6,82) 'tight (ap)                      '
          ELSE IF (response1%preconditioner_p.EQ.3) THEN
             IF (paral%io_parent)&
                  WRITE  (6,82) 'State dependent + DS            '
          ELSE IF (response1%preconditioner_p.EQ.4) THEN
             IF (paral%io_parent)&
                  WRITE  (6,82) 'State dependent + AP            '
          ENDIF
          IF (paral%io_parent)&
               WRITE (6,83) response2%hthrs_p
          IF (paral%io_parent)&
               WRITE (6,84) response2%tolog_p
          IF (response1%pcg_p) THEN
             IF (response1%prec_p) THEN
                IF (paral%io_parent)&
                     WRITE (6,85) 'preconditioned conj gradient    '
             ELSE
                IF (paral%io_parent)&
                     WRITE (6,85) 'bare conj gradient              '
             ENDIF
             IF (response1%tpolak_p) THEN
                IF (paral%io_parent)&
                     WRITE (6,89) 'Polak-Ribiere variant           '
             ELSE
                IF (paral%io_parent)&
                     WRITE (6,89) 'Fletcher-Reeves variant         '
             ENDIF
             IF (response1%pcgmin_p) THEN
                IF (paral%io_parent)&
                     WRITE (6,86) response1%cg_analytic,response2%cg_factor
             ELSE
                IF (paral%io_parent)&
                     WRITE (6,89) 'using FIXED STEP LENGTH         '
             ENDIF
          ELSEIF (response1%tsde_p) THEN
             IF (response1%prec_p) THEN
                IF (paral%io_parent)&
                     WRITE (6,85) 'preconditioned steepest descent '
             ELSE
                IF (paral%io_parent)&
                     WRITE (6,85) 'bare steepest descent           '
             ENDIF
          ENDIF
          IF (dmbi%wcpredict) THEN
             IF (paral%io_parent)&
                  WRITE (6,*) 'WC Pos prediction (reset every',&
                  dmbi%wc_pred_step,'steps)'
          ENDIF

          IF (paral%io_parent)&
               WRITE (6,90)

          ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
          IF (dmbi%wcpredict) THEN
             ! Allocate the memory necessary for the WCs prediction
             ALLOCATE(WC_array1(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(WC_array2(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF           ! WCpredict
          ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
       ENDIF               ! firstcall
    ENDIF                     ! Parent.

    ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
    IF (dmbi%wcpredict.AND.dmbi%inter_pt_firstcall) THEN
       ! Allocate memory for WC velocities (parent and others)
       ALLOCATE(WC_velocity(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
    ! ==--------------------------------------------------------------==
    ! First, ALLOCATE MEMORY which is needed for any kind of perturbation.
    ! Then, the specific perturbation subroutine is called.
    ! might need to be separated and variable passed to the routine..
    ! ==--------------------------------------------------------------==
    ! Wavefunctions, densities and related:
    ALLOCATE(h1psi0(ncpw%ngw*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(h1psi0)!,ngw*nstate)
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    CALL zeroing(psi)!,il_psi)
    CALL zeroing(rhoe)!,il_rhoe)

    ! NB: why compiler does not issue error about mismatching dimensions?..
    ! rho0 is real(8), pointer :: rho0(:) declared in response_p.mod.F
    ALLOCATE(rho0(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rho0)!,il_rhoe)
    ALLOCATE(drhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(drhoe)!,il_rhoe)

    ! potentials and related:
    ALLOCATE(VofRho0(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(VofRho0)!,il_rhoe)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop)!,nhg)
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eivps)!,nhg)

    ! others & scratch:
    ALLOCATE(z11(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(z11)!,nstate*nstate)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)

    CALL give_scr_perturbation(lscr,tag,nstate)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(scr)!,lscr)
    ! ==--------------------------------------------------------------==

    CALL phfac(tau0)
    CALL eicalc(eivps,eirop)

    ! Localising the orbitals
    IF (nmr_options%tlocalize) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'LOCALISING THE KS ORBITALS...'
       wannl%twann=.TRUE.
       prop1%locl=.TRUE.
       prop2%numorb=MAX(prop2%numorb,nstate)
       nc=2*nkpt%ngwk*prop2%numorb
       lenbk=nxxfun(prop2%numorb)
       nc=MAX(2*lenbk*parai%nproc,nc)
       ALLOCATE(sct(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL localize(tau0,c0,cs,sct,nstate)
       CALL eicalc(eivps,eirop)
       DEALLOCATE(sct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,*) 'Done'

       ! -TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
       ! only the parent knows the position of the WCs
       IF (paral%parent.AND.dmbi%wcpredict) THEN
          ! storing the first WCs (old)
          IF (dmbi%wc_sample.EQ.1) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Storing 1st WCs'
             ! SCR contains the coordinates of the wannier centres...
             ! copying the WCs from the scratch to WC_array1
             CALL store_wc(nstate,scr,WC_array1)

             ! little printout just for fun.... works OK
             ! write(6,*) 'WC_array1'
             ! DO istate=1,NSTATE
             ! write(6,*) WC_array1(1,istate),
             !                      WC_array1(2,istate),
             !                      WC_array1(3,istate)
             ! ENDDO

             ! storage of array1 completed, next call will store array2
             dmbi%wc_sample=2
             ! storing the second WCs (new)
          ELSE IF (dmbi%wc_sample.EQ.2) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Storing 2nd WCs'
             ! SCR contains the coordinates of the wannier centres...
             ! copying the WCs from the scratch to WC_array2
             CALL store_wc(nstate,scr,WC_array2)

             ! little printout just for fun.... works OK
             ! write(6,*) 'WC_array2'
             ! DO istate=1,NSTATE
             ! write(6,*) WC_array2(1,istate),
             !                      WC_array2(2,istate),
             !                      WC_array2(3,istate)
             ! ENDDO

             ! Now we can calculate the velocity of the WCs.
             CALL calculate_wc_velocity(nstate,WC_array1,&
                  WC_array2,WC_velocity)

             ! little printout just for fun....
             IF (paral%io_parent)&
                  WRITE(6,*) 'velocities'
             DO istate=1,nstate
                IF (paral%io_parent)&
                     WRITE(6,*) WC_velocity(1,istate),&
                     WC_velocity(2,istate),&
                     WC_velocity(3,istate)
             ENDDO

             ! Storage + velocity calc completed, next call will predict
             dmbi%wc_sample=0
             ! here one should broadcast WC_velocity to all processors...
             ! to be done by parent??
          ENDIF
       ENDIF

       IF (dmbi%wcpredict) THEN
          ! Broadcast of switch?
          CALL mp_bcast(dmbi%wc_sample,parai%source,parai%allgrp)
          ! here one should broadcast WC_velocity to all processors...
          CALL mp_bcast(WC_velocity,SIZE(WC_velocity),parai%source,parai%allgrp)
          IF (paral%io_parent)&
               WRITE(6,*) 'PROC:',parai%me,'WC_sample =',dmbi%wc_sample
       ENDIF

       ! De-allocating the memory for WC prediction
       ! CALL FREEM(IP_WC_array1)
       ! CALL FREEM(IP_WC_array2)
       ! CALL FREEM(IP_WC_velocity)
       ! -TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
    ENDIF
    ! end of localisation

    ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST
    ! Now (if requested) predict the position of WCs
    ! and orthogonalise

    ! Predicting the position of the WCs 
    IF (dmbi%wcpredict) THEN
       ! check for prediction(0) and no localise
       IF (dmbi%wc_sample.EQ.0) THEN
          IF (.NOT.nmr_options%tlocalize) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'PROC',parai%me,' moving WCs to predicted pos'
             ! Shift each state in C0 by the constant vector ammount WC_velocity
             CALL move_all_wannier(nstate,c0,WC_velocity)
          ELSE
             ! if we end up here, it means that we just calculated the 2nd point...
             ! so disable the localisation, and next call will predict WC position
             nmr_options%tlocalize=.FALSE.
          ENDIF
          ! orthogonalising c0 (C1 here is a scratch)
          CALL lowdin(c0,c1,nstate)
       ENDIF
    ENDIF

    ! Synchronise the processors before the PT start
    CALL mp_sync(parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,*) 'BPT:process',parai%me,' ready'
    ! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTTEST

    ! PERTURBATION OPTIMISATION OF W1

    IF (cntl%ttau) CALL stopgm("mddiag-inter","no tau functionals",& 
         __LINE__,__FILE__)
    ! initialisation of the PT vectors (forces,densities)
    CALL rhoofr(c0(:,:,1),rho0,psi(:,1),nstate)
    CALL dcopy(fpar%nnr1, rho0, 1, VofRho0, 1)
    CALL vofrho(tau0,fion,VofRho0,psi,.FALSE.,.FALSE.)

    IF (paral%parent) CALL prmem('  DO_PERTURBATION')

    ! calculate the lagrange parameters z11 and the forces on the wfns:
    CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)

    ! DMB E2 calc
    IF (paral%io_parent)&
         WRITE(6,*) 'ETOT(0) = ',ener_com%etot
    ! storing the starting energy
    e0=ener_com%etot
    ! DMB E2 calc

    ! ==--------------------------------------------------------------==
    ! optimise the response wavefunction
    CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
         eirop,eivps,dummy,h1psi0,&
         z11,nstate,dummy)

    ! END OF PERTURBATION OPTIMISATION OF W1

    ! DMB E2 calc
    IF (paral%io_parent)&
         WRITE(6,*) 'ETOT(2) = ',ener_com%etot
    ! storing the correction energy
    e2=ener_com%etot
    ! DMB E2 calc

    ! ==--------------------------------------------------------------==
    ! DE-allocate memory
    ! ==--------------------------------------------------------------==
    DEALLOCATE(h1psi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho0, STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(drhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(VofRho0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(z11,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    response1%response_running = .FALSE.

    ! add C1 to C0
    CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c1,1,c0,1)

    ! orthogonalising c0+c1 (C1 here is a scratch)
    CALL lowdin(c0,c1,nstate)

    ! calculate the forces on the ions with the new C0
    CALL forces(c0,c1,tau0,fion,rhoe,psi,&
         nstate,1,.FALSE.,.TRUE.)

    ! DMB E2 calc
    ! if(parent) write(6,*) 'ETOT= ',ETOT
    ! if(parent) write(6,*) 'EPT = ',E0+E2
    ! replacing etot by the second order functional energy
    ! ETOT=E0+E2
    ! DMB E2 calc
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mddiag_interaction_p
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE store_wc(nstate,center,WC_array)
    ! ==--------------------------------------------------------------==
    ! This routine stores the coordinates of the Wannier centres
    ! into WC_array, just to get rid of using SCR as storage space..
    INTEGER                                  :: nstate
    REAL(real_8)                             :: center(4,nstate), &
                                                WC_array(3,nstate)

    INTEGER                                  :: istate

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,*) 'Storing wannier centres'
    DO istate=1,nstate
       WC_array(1,istate)=center(1,istate)
       WC_array(2,istate)=center(2,istate)
       WC_array(3,istate)=center(3,istate)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE store_wc
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE calculate_wc_velocity(nstate,WC_array1,WC_array2,&
       WC_velocity)
    ! ==--------------------------------------------------------------==
    ! This routine calculates the velocity of the Wannier centres
    ! using the positions contained in WC_array1(old) and WC_array2(new)
    ! The method used is a simple crude finite difference and the result
    ! is put in WC_velocity
    INTEGER                                  :: nstate
    REAL(real_8)                             :: WC_array1(3,nstate), &
                                                WC_array2(3,nstate), &
                                                WC_velocity(3,nstate)

    INTEGER                                  :: istate

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,*) 'Calculating velocity of wannier centres'
    DO istate=1,nstate
       WC_velocity(1,istate)=WC_array2(1,istate)-WC_array1(1,istate)
       WC_velocity(2,istate)=WC_array2(2,istate)-WC_array1(2,istate)
       WC_velocity(3,istate)=WC_array2(3,istate)-WC_array1(3,istate)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calculate_wc_velocity
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE move_all_wannier(nstate,c0,WC_velocity)
    ! ==--------------------------------------------------------------==
    ! This routine moves all the Wannier functions state-by-state
    ! by a vector displacement given in WC_velocity for each state.
    ! This constitutes a crude propagation scheme for localised functions


    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: WC_velocity(3,nstate)

    COMPLEX(real_8)                          :: EXPiGR
    INTEGER                                  :: Gindx, istate
    REAL(real_8)                             :: GbyR

! Main loop on the states

    DO istate=1,nstate
       ! Loop on the G vectors
       DO Gindx=1,ncpw%ngw
          ! Calculating the product G.r
          GbyR = gk(1,Gindx)*WC_velocity(1,istate)&
               +GK(2,Gindx)*WC_velocity(2,istate)&
               +GK(3,Gindx)*WC_velocity(3,istate)

          ! Adjusting the units
          GbyR=-gbyr*parm%tpiba

          ! Calculating E^i(G.r) using cos(x)+i*sin(x)
          EXPiGR=CMPLX(COS(GbyR),SIN(gbyr),kind=real_8)

          ! multiplying the Wannier coeffs
          c0(Gindx,istate)=c0(gindx,istate)*EXPiGR
       ENDDO
    ENDDO
  END SUBROUTINE move_all_wannier
  ! ==================================================================

END MODULE mddiag_interaction_p_utils
