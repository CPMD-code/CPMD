MODULE interaction_p_utils
  USE atomwf_utils,                    ONLY: atomwf,&
                                             give_scr_atomwf
  USE atwf,                            ONLY: atwp,&
                                             catom
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE csize_utils,                     ONLY: csize
  USE ddip,                            ONLY: lenbk
  USE eicalc_utils,                    ONLY: eicalc
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE forces_driver,                   ONLY: forces
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE localize_utils,                  ONLY: localize
  USE lowdin_utils,                    ONLY: lowdin
  USE machine,                         ONLY: m_walltime
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: dmbi,&
                                             dmbr,&
                                             emono,&
                                             nmr_options,&
                                             response1,&
                                             rho0,&
                                             timetag,&
                                             vofrho0
  USE rhoofr_utils,                    ONLY: rhoofr
  USE ropt,                            ONLY: iteropt
  USE rotate_my_wannier_para_p_utils,  ONLY: rotatemywannier_para
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sdlinres_utils,                  ONLY: give_scr_sdlinres
  USE setirec_utils,                   ONLY: write_irec
  USE sort_utils,                      ONLY: sort2
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: icopy,&
                                             nxxfun
  USE vofrho_utils,                    ONLY: vofrho
  USE wann,                            ONLY: wannl
  USE wv30_utils,                      ONLY: zhwwf
!!use densto_utils, only : densto
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interaction_p
  !public :: print_orbital_resp_densities
  PUBLIC :: give_scr_interaction
  PUBLIC :: sort_array

CONTAINS

  ! ==================================================================
  SUBROUTINE interaction_p(c0,c1,psi,rhoe,drhoe,&
       eirop,eivps,z11,nstate)
    ! ==================================================================
    ! PARALLEL
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(maxfftn,1)
    REAL(real_8)                             :: rhoe(fpar%nnr1,1), drhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate,1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'interaction_p'

    CHARACTER(len=20)                        :: name
    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE :: C0plusC1(:,:,:), cattmp(:,:,:), cs(:,:), &
      h1psi0(:,:), originalC0(:,:), scr(:), sct(:,:), sumC1(:,:), tempC1(:,:)
    INTEGER                                  :: ierr, irec(100), is, &
                                                ISCF_step, isub, isub2, nc
    REAL(real_8)                             :: dummies(1), e0, Efunct, &
                                                Last_Etot, scf_cnorm, &
                                                scf_gemax, sum_mono, tcpu, &
                                                time1, time2

! Common variables
! ==--------------------------------------------------------------==
! vector containing the final result (C0+C1)
! tmp C1 vector for cutoff changes...
! Copy of the original C0 vector for SCF calculations
! Accumulated C1 vector for SCF calc.
! Scratch arrays for the localisation
! Convergence variables for SCF calculations
! temporary: atomic wf for orthogonalisation

99  FORMAT ("*  Calculation of responses DONE.",31x,"*")
98  FORMAT (47("*"),"PT-INTER*RESPONSE*")
    ! ==--------------------------------------------------------------==
    CALL tiset('   PTINTER',isub)
    ALLOCATE(scr(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    timetag=' '
    response1%pr_energy = .TRUE.

    ! Allocating memory
    ALLOCATE(h1psi0(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(C0plusC1(ncpw%ngw,nstate,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(emono(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (dmbi%bptscfiter.GT.1) THEN
       ! SCF arrays...
       ALLOCATE(originalC0(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sumC1(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('  PT-INTER')
    ! ==--------------------------------------------------------------==
    ! If we use atomic wavefunction then ititialise it..
    IF (dmbi%tatomicwavefunction) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Initialising ATOMIC WAVEFUNCTION for W0'

       CALL atomwf(c0,nstate,tau0,fion,rhoe,psi)
       IF (paral%io_parent)&
            WRITE(6,*) 'DONE...'
    ENDIF

    ! Localisation of the wavefunction if we re not loading from disk
    IF (nmr_options%tlocalize.AND..NOT.dmbi%wann_load) THEN
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

       ! freeing memory
       DEALLOCATE(sct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,*) 'Done'
    ENDIF

    ! End of localisation

    IF (dmbi%bptscfiter.GT.1 .AND. paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) '---------- SELF-CONSISTENT PT ----------'
       IF (paral%io_parent)&
            WRITE(6,*) '----- PERFORMING',dmbi%bptscfiter,'SCF iterations -----'
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,E8.3,A)') '----- SCF-PT CONVERGENCE: ',&
            dmbr%scf_tol,' -----'
       IF (paral%io_parent)&
            WRITE(6,*) '----------------------------------------'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! If we are not using an atomic wavefunction then use wannier functions...
    IF (.NOT.dmbi%tatomicwavefunction) THEN
       ! LOAD THE ORBITALS FROM DATABASE into c0
       IF (paral%io_parent)&
            WRITE(6,*) 'Initialising the Wannier orbitals (W0)'

       ! Use the parallel version instead... has no rotations though
       CALL rotatemywannier_para(nstate,c0,psi)

       IF (paral%io_parent)&
            WRITE(6,*) 'Wannier orbitals are now initialised'
    ENDIF

    ! Were we saving a WF? if yes then EXIT
    IF (dmbi%wann_save) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'WANNIER FUNCTION SAVING OPTION ENABLED'
          IF (paral%io_parent)&
               WRITE(6,*) 'PROGRAM TERMINATES HERE....'
       ENDIF
       CALL tihalt('   PTINTER',isub)

       ! DE-allocating the arrays....
       DEALLOCATE(h1psi0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(C0plusC1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (dmbi%cutoff_trick) THEN
          DEALLOCATE(tempC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (dmbi%bptscfiter.GT.1) THEN
          ! SCF arrays...
          DEALLOCATE(originalC0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(sumC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       RETURN
    ENDIF

    ! ---saving the starting guess made of localised wf
    IF (dmbi%save_localised_w0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'SAVING RESTART FOR ORTHOGONALISED W0'
       store1%swf     = .TRUE.
       CALL write_irec(irec)
       CALL zhwwf(2,irec,c0,h1psi0,nstate,dummies,&
            TAU0,VELP,TAUP,iteropt%nfi)
       IF (paral%io_parent)&
            WRITE(6,*) 'DONE'

       ! OK, now we are done....
       IF (paral%io_parent)&
            WRITE(6,*) 'PROGRAM TERMINATES HERE....'
       CALL tihalt('   PTINTER',isub)

       ! DE-allocating the arrays....
       DEALLOCATE(h1psi0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(C0plusC1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (dmbi%cutoff_trick) THEN
          DEALLOCATE(tempC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (dmbi%bptscfiter.GT.1) THEN
          ! SCF arrays...
          DEALLOCATE(originalC0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(sumC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       RETURN
    ENDIF
    ! ---end of saving the starting guess made of localised wf

    ! Calculate the sum of monomer energies on parent
    IF (paral%parent) THEN
       sum_mono=0._real_8
       DO is=1,dmbi%nmol
          sum_mono=sum_mono+emono(is)
       ENDDO
    ENDIF
    ! ==-START OF PT CALCULATION--------------------------------------==
    ! initialise timer
    time1 =m_walltime()

    ! and don t forget to update the constant gnd-state density and
    ! Coulomb-potential arrays:

    IF (dmbi%tatomicrho) THEN
       ! Here we use the atomic density as constant gnd state density
       IF (paral%io_parent)&
            WRITE(6,*) 'Using atomic RHO0 instead of W0^2'

       ! just allocate memory for the atomic WF...
       ALLOCATE(catom(nkpt%ngwk,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cattmp(nkpt%ngwk,atwp%nattot,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ! Use the full atomic WF routine (to fix later...)
       CALL atomwf(cattmp,nstate,tau0,fion,rhoe,psi)
       ! generate the atomic density from the KS-atomic wf... 
       CALL rhoofr(cattmp(:,:,1),rho0,psi(:,1),nstate)
       IF (cntl%ttau) CALL stopgm("interaction_p","no tau functionals",& 
            __LINE__,__FILE__)

       CALL dcopy(fpar%nnr1, rho0, 1, vofrho0, 1)
       CALL vofrho(tau0,fion,vofrho0,psi,.FALSE.,.FALSE.)

       ! Free the atomic WF..
       DEALLOCATE(cattmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       IF (paral%io_parent)&
            WRITE(6,*) 'Done...'
    ELSE
       ! Here the constant gnd state density is computed exactly
       CALL rhoofr(c0(:,:,1),rho0,psi(:,1),nstate)
       CALL dcopy(fpar%nnr1, rho0, 1, vofrho0, 1)
       CALL vofrho(tau0,fion,vofrho0,psi,.FALSE.,.FALSE.)
    ENDIF

    ! saving the total density to disk.... if we don t do SCF
    IF (dmbi%bptscfiter.LE.1) THEN
       CALL ffttog(rho0,scr,psi,ncpw%nhg,.FALSE.)
       name='total.dens'
       CALL densto(scr,tau0,name)
    ENDIF

    ! Total energy of the starting Wannier orbitals (W0)
    IF (dmbi%torthog_wannier) THEN
       CALL forces(c0,h1psi0,tau0,fion,rhoe,psi,&
            nstate,1,.FALSE.,.FALSE.)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*)'***** ENERGY DECOMPOSITION *****'
          IF (paral%io_parent)&
               WRITE(6,*)'EKIN  (       W0)= ',ener_com%ekin,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EHT   (       W0)= ',ener_com%eht,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EPSEU (       W0)= ',ener_com%epseu,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'ENL   (       W0)= ',ener_com%enl,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EXC   (       W0)= ',ener_com%exc,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'------------------'
          IF (paral%io_parent)&
               WRITE(6,*)'ETOT  (       W0)= ',ener_com%etot,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF

       ! Calculate convergence criteria on W0
       CALL csize(h1psi0,nstate,scf_gemax,scf_cnorm)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(a,E10.4,a,E10.4)')&
               'INITIAL GRADIENT OF W0: Gemax = ',scf_gemax,&
               ' Cnorm = ',scf_cnorm
       ENDIF
    ENDIF

    ! Storing the value of E(0)
    e0=ener_com%etot

    ! ==--------------------------------------------------------------==
    CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE (6,98)
    CALL tiset('PSI1-WF-OPT',isub2)
    CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
         eirop,eivps,cdummy,h1psi0,z11,nstate,cdummy)
    CALL tihalt('PSI1-WF-OPT',isub2)
    timetag='calculation of psi^1.'

    ! ==-END OF PT CALCULATION----------------------------------------==
    ! ETOT contains the perturbed energy from rwfopt_p

    ! Calculating the total energy using the functional
    Efunct=e0+ener_com%etot
    IF (dmbi%bptscfiter.GT.1) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'TOTAL ENERGY at step 1 (functional) =',Efunct
    ELSE
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'TOTAL ENERGY (functional) =',Efunct
          IF (paral%io_parent)&
               WRITE(6,*) 'Functional INTERACTION Energy =',&
               (Efunct-sum_mono)*1000,'mH'
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF
    ENDIF

    ! NOW, c1 contains the response orbitals. ANALYSE THEM:

    ! Saving the perturbed density, W0, and rho0 if we don t do SCF
    IF (dmbi%bptscfiter.LE.1) THEN
       CALL print_orbital_resp_densities(c0,c1,nstate,psi,scr)
    ENDIF

    CALL dcopy(2*ncpw%ngw*nstate,c0,1,C0plusC1,1)

    ! Now adding the perturbation (C1) to C0plusC1
    CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c1,1,C0plusC1,1)

    IF (dmbi%bptscfiter.LE.1) THEN
       ! Total energy of the non-orthogonalised W0+W1 if we don t do SCF
       CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
            nstate,1,.FALSE.,.FALSE.)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*)'TOTAL   ENERGY of W0+W1: ',ener_com%etot,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*)'***** ENERGY DECOMPOSITION *****'
          IF (paral%io_parent)&
               WRITE(6,*)'EKIN  (    W0+W1)= ',ener_com%ekin,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EHT   (    W0+W1)= ',ener_com%eht,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EPSEU (    W0+W1)= ',ener_com%epseu,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'ENL   (    W0+W1)= ',ener_com%enl,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'EXC   (    W0+W1)= ',ener_com%exc,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)'------------------'
          IF (paral%io_parent)&
               WRITE(6,*)'ETOT  (    W0+W1)= ',ener_com%etot,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF
    ENDIF

    ! orthogonalising the (w0+w1) states
    ! Allocating scratch array (CS)
    ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Parallel Lowdin
    CALL lowdin(C0plusC1,cs,nstate)
    ! de-alloc CS
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! Now calculating the total energy of the orthogonalised w0+w1

    CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
         nstate,1,.TRUE.,.FALSE.)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (dmbi%bptscfiter.GT.1) THEN
          IF (paral%io_parent)&
               WRITE(6,*) '*-=> SCF-PT STEP 1 <=-*'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*)'TOTAL   ENERGY of orthog. W0+W1: ',&
            ener_com%etot,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERACTION Energy =',&
            (ener_com%etot-sum_mono)*1000,'mH ',&
            '(from orthog. W0+W1)'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'Efunct-EKS(orth.) =',(efunct-ener_com%etot)*1000,'mH'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*)'***** ENERGY DECOMPOSITION *****'
       IF (paral%io_parent)&
            WRITE(6,*)'EKIN  (ort.W0+W1)= ',ener_com%ekin,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)'EHT   (ort.W0+W1)= ',ener_com%eht,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)'EPSEU (ort.W0+W1)= ',ener_com%epseu,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)'ENL   (ort.W0+W1)= ',ener_com%enl,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)'EXC   (ort.W0+W1)= ',ener_com%exc,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)'------------------'
       IF (paral%io_parent)&
            WRITE(6,*)'ETOT  (ort.W0+W1)= ',ener_com%etot,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    ! Calculate convergence criteria
    CALL csize(h1psi0,nstate,scf_gemax,scf_cnorm)
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(a,E10.4,a,E10.4,a,f8.2,a)')&
            '  1: Gemax =',scf_gemax,&
            ' Cnorm =',scf_cnorm,' Time: ',tcpu,' seconds'
    ENDIF

    ! Do we need SCF iterations?
    IF (dmbi%bptscfiter.GT.1) THEN
       ! saving the original C0
       CALL dcopy(2*ncpw%ngw*nstate,c0,1,originalC0,1)
       ! --------- START OF SCF cycles ------------------
       DO ISCF_step=2,dmbi%bptscfiter
          ! Are we converged ?
          IF (scf_gemax.LT.dmbr%scf_tol) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Tolerance reached...',scf_gemax
             GOTO 100
          ENDIF

          ! initialise timer
          time1 =m_walltime()

          ! storing the last Etot
          Last_Etot=ener_com%etot
          ! Storing the value of E(0)
          e0=ener_com%etot

          ! copying C0+C1 to C0...
          CALL dcopy(2*ncpw%ngw*nstate,C0plusC1,1,c0,1)

          ! ---------- START OF PT CALC ------------------------------
          ! initialisation
          CALL rhoofr(c0(:,:,1),rho0,psi(:,1),nstate)
          CALL dcopy(fpar%nnr1, rho0, 1, vofrho0, 1)
          CALL vofrho(tau0,fion,vofrho0,psi,&
               .FALSE.,.FALSE.)

          ! Lagrange multipliers
          CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)

          ! Optimisation of W1
          IF (paral%io_parent)&
               WRITE (6,98)
          CALL tiset('PSI1-WF-OPT',isub2)
          CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,&
               cdummy,h1psi0,z11,nstate,cdummy)
          CALL tihalt('PSI1-WF-OPT',isub2)
          timetag='calculation of psi^1.'
          ! ---------- END OF PT CALC ------------------------------

          ! Is there no interaction force?
          IF (dmbi%tnoresponse) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'NEGLIGIBLE RESPONSE...'
             GOTO 100
          ENDIF

          ! Calculating the total energy using the functional
          Efunct=e0+ener_com%etot
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)
             IF (paral%io_parent)&
                  WRITE(6,*) '*-=> SCF-PT STEP',ISCF_step,'<=-*'
             IF (paral%io_parent)&
                  WRITE(6,*) 'TOTAL ENERGY at step',ISCF_step,&
                  '(functional) =',Efunct
             IF (paral%io_parent)&
                  WRITE(6,*) 'Functional INTERACTION Energy =',&
                  (Efunct-sum_mono)*1000,'mH'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF

          ! copying C0 to C0plusC1
          CALL dcopy(2*ncpw%ngw*nstate,c0,1,C0plusC1,1)

          ! Now adding the perturbation (C1) to C0plusC1
          CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c1,1,C0plusC1,1)

          ! Adding the perturbation to the previous C1
          ! call DAXPY(2*ngw*nstate,1._real_8,C1,1,sumC1,1)

          ! orthogonalising c0+c1 (C1 here is a scratch)
          CALL lowdin(C0plusC1,c1,nstate)

          ! Now calculating the total energy of the orthogonalised w0+w1
          CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
               nstate,1,.TRUE.,.FALSE.)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'TOTAL ENERGY of orthog. W0+W1 at step',&
                  ISCF_step,' : ',ener_com%etot,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*) 'Ediff E(',ISCF_step,')-E(',iscf_step-1,&
                  '): ',ener_com%etot-Last_Etot,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*) 'Efunct-EKS =',(efunct-ener_com%etot)*1000,'mH'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF

          ! Calculate convergence criteria
          CALL csize(h1psi0,nstate,scf_gemax,scf_cnorm)
          IF (paral%parent) THEN
             time2 =m_walltime()
             tcpu = (time2 - time1)*0.001_real_8
             IF (paral%io_parent)&
                  WRITE(6,'(i3,a,E10.4,a,E10.4,a,f8.2,a)') ISCF_step,&
                  ': Gemax =',scf_gemax,&
                  ' Cnorm =',scf_cnorm,' Time: ',tcpu,' seconds'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF
       ENDDO
       ! --------- END OF SCF cycles ------------------
       ! We have exhausted all the SCF cycles without convergence...
       IF (paral%parent.AND.scf_gemax.GT.dmbr%scf_tol) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1x,64("!"))')
          IF (paral%io_parent)&
               WRITE(6,'(" !!",a,t64,"!!")')&
               '   SCF-PT| the maximum number of SCF step is reached'
          IF (paral%io_parent)&
               WRITE(6,'(" !!",a,f10.6,a,t64,"!!")')&
               '         but no convergence (d E / d psi =',scf_gemax,')'
          IF (paral%io_parent)&
               WRITE(6,'(1x,64("!"))')
       ENDIF

       ! We are converged: end of SCF...
100    CONTINUE
    ENDIF

    ! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,98)
       IF (paral%io_parent)&
            WRITE (6,99)
       IF (paral%io_parent)&
            WRITE (6,98)
    ENDIF

    ! ==--------------------------------------------------------------==
    CALL tihalt('   PTINTER',isub)

    ! DE-allocating the arrays....
    DEALLOCATE(h1psi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(C0plusC1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (dmbi%bptscfiter.GT.1) THEN
       ! SCF arrays...
       DEALLOCATE(originalC0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(sumC1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! TODO deallocate at every return statement
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE interaction_p
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE give_scr_interaction(l_max,tag,nstate)
    INTEGER                                  :: l_max
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: latomwf, lsecder

! ==--------------------------------------------------------------==

    l_max=2*fpar%nnr1               ! for the gaussian fitting routine

    ! Atomic wf scratch...
    CALL give_scr_atomwf(latomwf,tag,nstate)
    ! simple model scratch
    CALL give_scr_sdlinres(lsecder,tag)

    l_max=MAX(l_max,latomwf)
    l_max=MAX(l_max,lsecder)

    l_max=MAX(l_max,4*ncpw%ngw*nstate)  ! for fun. repeat if necessary.
    tag = '2*nnr1'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_interaction
  ! ================================================================== 

  ! ==================================================================
  SUBROUTINE sort_array(a,b,n)
    ! ==================================================================
    ! Sorts Array A using an index and arranges A,B according to it
    ! Including the "essential" variables...

    INTEGER                                  :: n, b(n)
    REAL(real_8)                             :: a(n)

    INTEGER                                  :: i, INDEX(n), tempB(n)
    REAL(real_8)                             :: tempA(n)

! sorting array A

    CALL sort2(a,n,index)

    ! sort B accordingly using index
    DO i = 1 , n
       tempB(i)=b(INDEX(i))
    ENDDO
    CALL icopy(n,tempB,1,b,1)

    ! Arrange A and B in descending order
    DO i = 1 , n
       tempA(i)=a(n-i+1)
       tempB(i)=b(n-i+1)
    ENDDO
    CALL dcopy(n,tempA,1,a,1)
    CALL icopy(n,tempB,1,b,1)
    RETURN

  END SUBROUTINE sort_array
  ! ==================================================================

END MODULE interaction_p_utils


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! Subroutine to save the perturbed DENSITY      C
! (2 W0W1 -> rho1), the starting                C 
! orthogonalised wannier functions (W0 -> psi0) C
! and its corresponding density (rho0 = W0*W0)  C
! C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE print_orbital_resp_densities(c0,c1,nstate,psi,scr)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum, mp_bcast
  USE system , ONLY: ncpw,fpar
  USE parac, ONLY : paral,parai
  USE coor , ONLY:tau0
  USE geq0mod , ONLY:geq0
  USE response_pmod , ONLY:dmbi
  USE fft_maxfft,                      ONLY: maxfftn, maxfft
  USE readsr_utils, ONLY : readsr, xstring, readsi
  USE interaction_p_utils, ONLY : sort_array
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
  REAL(real_8)                               :: psi(fpar%nnr1,2), &
                                                scr(fpar%nnr1,2)

  CHARACTER(*), PARAMETER :: procedureN = 'print_orbital_resp_densities'

  CHARACTER(len=8)                           :: numb
  CHARACTER(len=80)                          :: filename
  INTEGER                                    :: curr_pos, ia, ie, ierr, ir, is
  INTEGER, ALLOCATABLE                       :: charge_disp_pos(:)
  REAL(real_8)                               :: charge_disp, sum_C0, sum_C0C1
  REAL(real_8), ALLOCATABLE                  :: charge_disp_val(:), &
                                                rho1_tot(:)

!!use densto_utils, only : densto
! Storage for the response density      
! Storage for the most important charge displacements     
! and their position
! ===== Allocate memory for rho1_tot and initialise ====================   

  ALLOCATE(rho1_tot(fpar%nnr1),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(rho1_tot)!,nnr1)

  ! ===== Allocate memory for charge disp arrays =========================
  ALLOCATE(charge_disp_val(dmbi%max_num_disp),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ALLOCATE(charge_disp_pos(dmbi%max_num_disp),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(charge_disp_val)!,dmbi%max_num_disp)
  DO is=1,dmbi%max_num_disp
     charge_disp_pos(is)=-1
     charge_disp_val(is)=-dmbi%max_num_disp+is
  ENDDO

  ! ===== First loop on all states and constructing rho1 ================= 
  DO is=1,nstate
     ! saving the perturbed density (2W0W1)
     CALL fft2tor(c0(1,is),psi(1,1),&
          c1(1,is),psi(1,2),scr,ncpw%ngw,.FALSE.)

     sum_C0=0._real_8
     sum_C0C1=0._real_8

     DO ir=1,fpar%nnr1
        scr(ir,1)  =  psi(ir,1)*psi(ir,2) *2._real_8
        sum_C0=sum_c0+psi(ir,1)*psi(ir,1)
        IF (scr(ir,1).GT.0) sum_C0C1=sum_c0c1+scr(ir,1)
        rho1_tot(ir)=rho1_tot(ir)+scr(ir,1)
     ENDDO

     CALL mp_sum(sum_C0,parai%allgrp)
     CALL mp_sum(sum_C0C1,parai%allgrp)
     ! call mp_sum(rho1_tot,nnr1,allgrp)
     ! NOT NECESSARY FOR SOME REASONS...
     IF (paral%parent) THEN
        ! calculate charge displacement on parent
        charge_disp=100*(sum_C0C1/REAL(nstate*sum_C0,kind=real_8))
        ! if charge displacmeent is bigger than the last element
        ! of the charge disp array, include it
        IF (charge_disp.GT.charge_disp_val(dmbi%max_num_disp)) THEN
           charge_disp_val(dmbi%max_num_disp)=charge_disp
           ! store position
           charge_disp_pos(dmbi%max_num_disp)=is
           ! re-sort both arrays
           CALL sort_array(charge_disp_val,charge_disp_pos,&
                dmbi%max_num_disp)
        ENDIF
     ENDIF
  ENDDO

  ! ===== All states have been done, saving rho1 =========================
  ! rho1_tot back in G space
  CALL ffttog(rho1_tot,scr(1,1),psi,ncpw%nhg,.TRUE.)

  ! Now save the total rho1 density
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,*) 'Saving total perturbation density (2W0W1)'
     filename = 'rho1.total'
  ENDIF
  CALL densto(scr(1,1),tau0,filename)

  ! Finally clearing rho1 from the memory
  DEALLOCATE(rho1_tot,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)

  ! ===== Now make sure all processor know the disp values and order =====
  CALL mp_bcast(charge_disp_pos,dmbi%max_num_disp,parai%source,parai%allgrp)
  CALL mp_bcast(charge_disp_val,dmbi%max_num_disp,parai%source,parai%allgrp)

  ! ===== Now saving on the states with the largest charge disp ==========
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,*) 'Saving the',dmbi%max_num_disp,'states contributing the',&
          ' largest charge displacements'
     IF (paral%io_parent)&
          WRITE(6,*) '- (in decreasing order of importance) -'
  ENDIF
  DO curr_pos=1,dmbi%max_num_disp
     ! looking up state number
     is=charge_disp_pos(curr_pos)

     IF (is.GT.0) THEN
        ! saving the perturbed density (2W0W1)
        CALL fft2tor(c0(1,is),psi(1,1),&
             c1(1,is),psi(1,2),scr,ncpw%ngw,.FALSE.)

        DO ir=1,fpar%nnr1
           scr(ir,1)  =  psi(ir,1)*psi(ir,2) *2._real_8
        ENDDO

        IF (paral%io_parent)&
             WRITE(6,fmt='(1X,A,I3,A,F8.5,A)') 'STATE',is,&
             ' contributes',charge_disp_val(curr_pos),&
             '% of charge displacement'

        ! rho1 (W0W1) back in G space
        CALL ffttog(scr(1,1),scr(1,1),psi,ncpw%nhg,.TRUE.)
        ! Is the norm of the perturbed density equal to zero?
        IF ((ABS(REAL(scr(1,1),kind=real_8)).GT.1.e-13_real_8).AND.(geq0)) THEN
           IF (paral%io_parent)&
                WRITE(6,*) 'WARNING: R-Norm of state',is,&
                '[n1(g=0)] =',scr(1,1),'! = 0'
        ENDIF
        ! Convert the number of the state to chars (from wannier_print)
        IF (paral%parent) THEN
           IF (paral%io_parent)&
                WRITE(numb,'(I8)') is
           CALL xstring(numb,ia,ie)
        ENDIF

        IF (paral%parent)&
             filename = 'rho1-'//NUMB(IA:IE)//'.dens'

        CALL densto(scr(1,1),tau0,filename)

        ! saving the starting wannier orbitals (W0)
        CALL zeroing(psi)!,2*maxfft)
        CALL dcopy(2*ncpw%ngw,c0(1,is),1,psi,1)
        IF (paral%parent)&
             filename = 'psi0-'//NUMB(IA:IE)//'.wf'
        CALL densto(psi,tau0,filename)

        ! saving the starting wannier density (rho0)
        ! G to R space
        CALL ffttor(c0(1,is),psi(1,1),scr,ncpw%ngw,.FALSE.)

        ! calculate the square
        DO ir=1,fpar%nnr1
           scr(ir,1)  =  psi(ir,1)*psi(ir,1)
        ENDDO

        ! R to G space
        CALL ffttog(scr(1,1),scr(1,1),psi,ncpw%nhg,.TRUE.)

        ! Save rho0
        IF (paral%parent)&
             filename = 'rho0-'//NUMB(IA:IE)//'.dens'
        CALL densto(scr(1,1),tau0,filename)
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE print_orbital_resp_densities
! ================================================================== 
