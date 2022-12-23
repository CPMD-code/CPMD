MODULE interaction_manno_p_utils
  USE atomwf_utils,                    ONLY: atomwf
  USE atwf,                            ONLY: atwp,&
                                             catom
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE cppt,                            ONLY: hg
  USE csize_utils,                     ONLY: csize
  USE ddip,                            ONLY: lenbk,&
                                             ngwmax
  USE ddipo_utils,                     ONLY: give_scr_ddipo,&
                                             set_operator,&
                                             setdip
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE forces_driver,                   ONLY: forces
  USE forces_utils,                    ONLY: give_scr_forces
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE localize_utils,                  ONLY: localize
  USE lowdin_utils,                    ONLY: lowdin
  USE machine,                         ONLY: m_walltime
  USE opeigr_utils,                    ONLY: opeigr
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: dmbi,&
                                             dmbr,&
                                             nmr_options,&
                                             response1,&
                                             rho0,&
                                             timetag,&
                                             vofrho0
  USE rhoofr_utils,                    ONLY: rhoofr
  USE ropt,                            ONLY: iteropt
  USE rotate_my_wannier_manno_p_utils, ONLY: idelta,&
                                             ortho_my_wannier,&
                                             rotatemywannier_manno
  USE rwfopt_p_utils,                  ONLY: print_matrix,&
                                             rwfopt_p,&
                                             simple_model_p
  USE setirec_utils,                   ONLY: write_irec
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy,&
                                             nxxfun
  USE vofrho_utils,                    ONLY: vofrho
  USE wann,                            ONLY: wannc,&
                                             wannl
  USE wannier_center_utils,            ONLY: wannier_center
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interaction_manno_p
  !public :: print_orbital_resp_dens_manno
  !public :: save_aldo_wannier
  !public :: save_norm
  !public :: calc_my_wannier_centre


CONTAINS

  ! ==================================================================
  SUBROUTINE interaction_manno_p(c0,c1,psi,rhoe,drhoe,&
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

    CHARACTER(*), PARAMETER :: procedureN = 'interaction_manno_p'

    CHARACTER(len=20)                        :: name
    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE :: C0plusC1(:,:,:), cattmp(:,:,:), cs(:,:), &
      h1psi0(:,:), originalC0(:,:), scr(:), sct(:,:), sumC1(:,:), tempC1(:,:)
    INTEGER                                  :: i, ierr, irec(100), &
                                                ISCF_step, istate, isub, &
                                                isub2, j, k, nc, New_NGW
    REAL(real_8) :: centre(4,nstate), dummies(1), e0, E_total, Efunct, &
      eigenval(nstate), eigenvecs(nstate,nstate), Last_Etot, my_ecut, &
      My_New_Gcut, scf_cnorm, scf_gemax, tcpu, time1, time2, work(3*nstate), &
      z(nstate*(nstate+1)/2)
    REAL(real_8), ALLOCATABLE                :: c0h0c0_ij(:,:), Emat(:,:), &
                                                inverse_of_S(:,:), W0_Sij(:,:)

! ==--------------------------------------------------------------==
! vector containing the final result (C0+C1)
! tmp C1 vector for cutoff changes...
! Copy of the original C0 vector for SCF calculations
! Accumulated C1 vector for SCF calc.
! Scratch arrays for the localisation
! variables for the articial cutoff reduction
! Convergence variables for SCF calculations
! temporary: atomic wf for orthogonalisation
! temporary: non-orthogonal w0, overlap matrix
! Storage array for Wannier centres

99  FORMAT ("*  Calculation of responses DONE.",31x,"*")
98  FORMAT (47("*"),"PT-INTER*RESPONSE*")
    ! ==--------------------------------------------------------------==
    CALL tiset('   PTINTER',isub)
    timetag=' '
    response1%pr_energy = .TRUE.

    ALLOCATE(scr(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO check scr length
    ! Allocating memory
    ALLOCATE(h1psi0(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(C0plusC1(ncpw%ngw,nstate,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (dmbi%cutoff_trick) THEN
       ALLOCATE(tempC1(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (dmbi%bptscfiter.GT.1) THEN
       ! SCF arrays...
       ALLOCATE(originalC0(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sumC1(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! overlap matrix for non-orthog W0 and so on..
    IF (.NOT.dmbi%torthog_wannier) THEN
       ALLOCATE(W0_Sij(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0h0c0_ij(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(inverse_of_S(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(Emat(nstate,nstate),STAT=ierr)
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

       ! save a RESTART file (disabled)
       ! SWF     = .TRUE.
       ! CALL WRITE_IREC(IREC)
       ! CALL ZHWWF(2,IREC,C0,CS,NSTATE,scr,TAU0,VELP,TAUP,NFI)

       ! freeing memory
       DEALLOCATE(cs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(sct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,*) 'Done'
    ENDIF

    ! End of localisation

    IF (dmbi%cutoff_trick) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'CUTOFF TRICK enabled...'
    ENDIF
    ! Are we using a lower cutoff, if yes assign the variables...
    IF (dmbi%cutoff_restr) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'using a lower cutoff for the perturbation'
       ENDIF
       ! MY_ECUT=ECUT/2._real_8
       my_ecut=40._real_8
       dmbi%ngw_zero=get_new_NGW(my_ecut)
    ENDIF
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
    ! If we are not using an atomic wavefuncttion then used wannier functions...
    IF (.NOT.dmbi%tatomicwavefunction) THEN
       ! LOAD THE ORBITALS FROM DATABASE into c0
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Initialising the Wannier orbitals (W0)'
       ENDIF
       CALL rotatemywannier_manno(nstate,c0,psi)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Wannier orbitals are now initialised'
       ENDIF
    ENDIF

    ! Were we saving a WF? if yes then EXIT
    IF (dmbi%wann_save) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'WANNIER FUNCTION SAVING OPTION ENABLED'
       IF (paral%io_parent)&
            WRITE(6,*) 'PROGRAM TERMINATES HERE....'
       CALL tihalt('   PTINTER',isub)

       ! DE-allocating the arrays....
       IF (dmbi%bptscfiter.GT.1) THEN
          DEALLOCATE(sumC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(originalC0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (dmbi%cutoff_trick) THEN
          DEALLOCATE(tempC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       DEALLOCATE(C0plusC1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(h1psi0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
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
       IF (dmbi%bptscfiter.GT.1) THEN
          DEALLOCATE(sumC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(originalC0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (dmbi%cutoff_trick) THEN
          DEALLOCATE(tempC1,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       DEALLOCATE(C0plusC1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(h1psi0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       RETURN
    ENDIF
    ! ---end of saving the starting guess made of localised wf


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
    IF (dmbi%torthog_wannier.AND..NOT.dmbi%tsimple_model) THEN
       ! IF (TORTHOG_WANNIER.OR.TSIMPLE_MODEL) THEN
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

    ! simple interaction model 
    IF (dmbi%tsimple_model) THEN
       CALL simple_model_p(c0,c1,psi,rhoe,nstate,eirop,eivps)
       RETURN
    ENDIF

    ! Storing the value of E(0)
    e0=ener_com%etot

    ! non orthogonal energy <c0|H0|c0>
    IF (.NOT.dmbi%torthog_wannier) THEN
       ! copying F... (state occupation)

       !vw>>>
       !vw this is buggus     call icopy(nstate*nkpt%nkpnt, f, 1, F_tmp, 1)
       CALL stopgm(procedureN,'fix me',&
            __LINE__,__FILE__)
       !vw<<<
       ! change F to single occupation for all states
       DO i=1,nstate
          crge%f(i,1)=1
       ENDDO

       ! Total energy of the starting Wannier orbitals (W0)
       CALL forces(c0,h1psi0,tau0,fion,rhoe,psi,&
            nstate,1,.FALSE.,.FALSE.)

       ! copying original F back to its normal location
       !vw>>>
       !vw this is buggus          call icopy(nstate*nkpt%nkpnt, F_tmp, 1, f, 1)
       CALL stopgm(procedureN,'fix me',&
            __LINE__,__FILE__)
       !vw<<<


       ! CALL MEMORY(IP_CS,NC,'CS/tmp   ')
       ! CALL hpsi(c0,h1psi0,CS,vofrho0,PSI,NSTATE,1,1)
       ! CALL FREEM(IP_CS)

       ! Use OVLAP to build hamiltonian matrix
       CALL ovlap(nstate,c0h0c0_ij,c0(:,:,1),h1psi0)

       ! now printing...
       CALL print_matrix(c0h0c0_ij,nstate,nstate,'c0h0c0_ij ')

       ! overlap matrix
       CALL ovlap(nstate,W0_Sij,c0(:,:,1),c0(:,:,1))

       ! now printing...
       CALL print_matrix(W0_Sij,nstate,nstate,'W0_Sij    ')

       ! calculating inverse_of_S...
       !$omp parallel do private(i,j)
       DO i=1,nstate
          DO j=1,nstate
             inverse_of_S(i,j)=REAL(2*idelta(i-j),kind=real_8)-W0_Sij(i,j)
          ENDDO
       ENDDO

       ! now printing...
       CALL print_matrix(inverse_of_S,nstate,nstate,'inv of S  ')

       ! calculating energy...
       CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,c0h0c0_ij,&
            NSTATE,inverse_of_S,NSTATE,0._real_8,Emat,NSTATE)

       ! now printing...
       CALL print_matrix(Emat,nstate,nstate,'Emat      ')

       ! filling up Z packed matrix for diagonalisation
       k=1
       DO i=1,nstate
          DO j=i,nstate
             z(k)=Emat(i,j)
             k=k+1
          ENDDO
       ENDDO

       ! Diagonalising Z matrix 
       CALL dspevy(1,z,eigenval,eigenvecs,nstate,nstate,work,3*nstate)

       E_total=0
       ! printing eigenval and eigenvecs
       IF (paral%io_parent)&
            WRITE(6,*) 'Eigenvalues + eigenvecs'
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,A10,10(1X,I12))') 'ISTATE',' EIGENVAL',&
            (j,j=1,NSTATE)
       DO i=1,nstate
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,F10.6,10(1X,F12.6))') i,&
               eigenval(i),&
               (eigenvecs(i,j),j=1,NSTATE)
          E_total=e_total-eigenval(i)
       ENDDO

       IF (paral%io_parent)&
            WRITE(6,*) 'Total energy (band energy) =',E_total
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==

    CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)
    ! if (parent) 
    ! &    WRITE(6,*)'interact: <W0_i| H0 |W0_j> = ',
    ! &    (z11(i,i),i=1,nstate)
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
       IF (paral%io_parent)&
            WRITE(6,*) 'TOTAL ENERGY (functional) =',Efunct
    ENDIF

    ! NOW, c1 contains the response orbitals. ANALYSE THEM:

    ! Saving the perturbed density, W0, and rho0 if we don t do SCF
    IF (dmbi%bptscfiter.LE.1) THEN
       CALL print_orbital_resp_dens_manno(c0,c1,nstate,psi,scr)
    ENDIF

    ! if there is no further work to do, add PT + calc energy
    IF (.NOT.dmbi%cutoff_trick) THEN
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
       CALL ortho_my_wannier(nstate,C0plusC1)

       ! Now calculating the total energy of the orthogonalised w0+w1

       CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
            nstate,1,.TRUE.,.FALSE.)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (dmbi%cutoff_restr) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'PT-CUTOFF = ',my_ecut
          ENDIF
          IF (dmbi%bptscfiter.GT.1) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) '*-=> SCF-PT STEP 1 <=-*'
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,*)'TOTAL   ENERGY of orthog. W0+W1: ',&
               ener_com%etot,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*) 'Efunct-EKS =',efunct-ener_com%etot
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

       ! ------temporary.... for charge calc..
       ! write(6,*) 'SAVING restart for orthog W0+W1'
       ! SWF     = .TRUE.
       ! CALL WRITE_IREC(IREC)
       ! CALL ZHWWF(2,IREC,C0plusC1,h1psi0,NSTATE,scr,
       !              TAU0,VELP,TAUP,NFI)
       ! write(6,*) 'DONE'
       ! ------temporary.... for charge calc..

       ! Do we need SCF iterations?
       IF (dmbi%bptscfiter.GT.1) THEN
          ! saving the original C0
          CALL dcopy(2*ncpw%ngw*nstate,c0,1,originalC0,1)
          ! copying C1 from the first cycle
          ! call DCOPY(2*ngw*nstate,C1,1,sumC1,1)
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
                     WRITE(6,*) 'Efunct-EKS =',efunct-ener_com%etot
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

          ! We are converged: save and finish
100       CONTINUE

          ! Saving the final perturbed density and the W0
          ! call print_orbital_resp_dens_manno(originalC0,sumC1,NSTATE,
          !                                  PSI,SCR)
          ! TESTTESTTEST
          ! call DCOPY(2*ngw*nstate,originalC0,1,C0plusC1,1)
          ! call DAXPY(2*ngw*nstate,1._real_8,sumC1,1,C0plusC1,1)
          ! CALL Ortho_My_Wannier(NSTATE,C0plusC1)

          ! call forces(c0plusC1,h1psi0,tau0,FION,RHOE,PSI,
          !                     SCR,LSCR,NSTATE,1,.false.,.false.)
          ! if (parent) then
          ! write(6,*)
          ! write(6,*) 'TEST....'
          ! WRITE(6,*)'TOTAL   ENERGY of orthog. W0+W1: ',
          !               ETOT,' A.U.'
          ! write(6,*) 'Ediff: ',ETOT-Last_Etot,' A.U.'
          ! endif
          ! TESTTESTTEST
       ENDIF
    ELSE
       ! ===== CUTOFF LOWERING EXPERIMENT !!! =====
       ! copying C1 in a temporary location 
       CALL dcopy(2*ncpw%ngw*nstate,c1,1,tempC1,1)
       DO j=1,6
          ! copying C0 to C0plusC1
          CALL dcopy(2*ncpw%ngw*nstate,c0,1,C0plusC1,1)

          ! defining the new cutoff and finding the corresponding G vector
          My_New_Gcut = (cntr%ecut-10._real_8*REAL(j-1,kind=real_8))/parm%tpiba2
          DO i=1,ncpw%ngw
             IF (hg(i).LE.My_New_Gcut) THEN
                New_NGW=i
             ENDIF
          ENDDO
          ! checking if G=0 is selected + correcting
          IF (New_NGW.EQ.1) THEN
             New_NGW=0
          ENDIF

          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'original ECUT = ',cntr%ecut
             IF (paral%io_parent)&
                  WRITE(6,*) 'new_ECUT = ',cntr%ecut-10._real_8*REAL(j-1,kind=real_8)
             IF (paral%io_parent)&
                  WRITE(6,*) 'My_New_Cut = ',My_New_Gcut
             IF (paral%io_parent)&
                  WRITE(6,*) 'new NGW = ',New_NGW
          ENDIF

          ! Copying from tempC1 to C1....
          CALL dcopy(2*ncpw%ngw*nstate,tempC1,1,c1,1)
          ! Zeroing the coefficients of C1 above the cutoff energy for each state
          DO istate=1,nstate
             DO i=New_NGW+1,ncpw%ngw
                c1(i,istate)=0._real_8
             ENDDO
          ENDDO

          ! Now adding the perturbation (C1) to C0plusC1
          CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c1,1,C0plusC1,1)

          ! Total energy of the non-orthogonalised W0+W1
          CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
               nstate,1,.FALSE.,.FALSE.)
          IF (paral%parent)  THEN
             IF (paral%io_parent)&
                  WRITE(6,*)
             IF (paral%io_parent)&
                  WRITE(6,*)'TOTAL   ENERGY of W0+W1: ',ener_com%etot,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)'KINETIC ENERGY of W0+W1: ',ener_com%ekin,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)'HARTREE ENERGY of W0+W1: ',ener_com%eht,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF

          ! orthogonalising the (w0+w1) states
          CALL ortho_my_wannier(nstate,C0plusC1)

          ! Now calculating the total energy of w0+w1

          CALL forces(C0plusC1,h1psi0,tau0,fion,rhoe,psi,&
               nstate,1,.FALSE.,.FALSE.)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)
             IF (paral%io_parent)&
                  WRITE(6,*)'TOTAL   ENERGY of orthog. W0+W1: ',&
                  ener_com%etot,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)'KINETIC ENERGY of orthog. W0+W1: ',&
                  ener_com%ekin,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)'HARTREE ENERGY of orthog. W0+W1: ',&
                  ener_com%eht,' A.U.'
             IF (paral%io_parent)&
                  WRITE(6,*)

             ! Printing a summary...
             IF (paral%io_parent)&
                  WRITE(6,*) 'SUM Ecut',cntr%ecut-10._real_8*REAL(j-1,kind=real_8),&
                  '-> Etot =',ener_com%etot 
          ENDIF

          ! end of Cutoff loop...
       ENDDO
       ! ===== END OF CUTOFF LOWERING EXPERIMENT !!! =====
    ENDIF

    ! Saving the necessary data for Aldo s program..
    ! write(6,*) 'Saving the resulting orbitals ',
    ! '(W0+W1) to disk'

    ! CALL Save_Aldo_Wannier(C0plusC1,NSTATE)

    ! write(6,*) 'Done!'
    ! END of file saving part...

    ! Calculating the WF profile
    ! write(6,*) 'calculating+ saving total wf profile'
    ! CALL PROFILE_WAN(C0plusC1,..,nstate)
    ! write(6,*) 'done'

    ! calculating the dipole moment+WC''s
    CALL calc_my_wannier_centre(C0plusC1,nstate,centre)

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
    ! overlap matrix..
    IF (.NOT.dmbi%torthog_wannier) THEN
       DEALLOCATE(W0_Sij,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c0h0c0_ij,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(Emat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! TODO deallocate at every return statement in this subroutine
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE interaction_manno_p
  ! ==================================================================
  SUBROUTINE save_aldo_wannier(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! Subroutine to save the resulting WF to use
    ! with Aldo s program (Wannier profiles)
    ! FORMAT: one file per state
    ! NGW
    ! C0(Gvec,state)
    ! ==--------------------------------------------------------------==

    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(len=20)                        :: name
    CHARACTER(len=8)                         :: number
    INTEGER                                  :: Gindx, ia, ie, istate
    LOGICAL                                  :: ferror

    CALL stopgm('INTERACTION_P','SAVE-ALDO-WANNIER not parallel.',& 
         __LINE__,__FILE__)

    ! main loop on the states...
    DO istate=1,nstate

       ! First creating the file name for the state istate
       IF (paral%io_parent)&
            WRITE(number,'(I8)')istate
       ! removing the spaces before and after the number (CPMD routine)
       CALL xstring(number,ia,ie)

       name='WANNIER-'//number(ia:ie)
       IF (paral%io_parent)&
            WRITE(6,*) 'opening file: ',name

       ! opening the file to save to
       IF (paral%io_parent)&
            CALL fileopen(500,name,fo_def,ferror)
       ! number of G vectors
       IF (paral%io_parent)&
            WRITE(500,*) ncpw%ngw

       ! now writing the coefficients (c0) (g-space,ascii)
       DO Gindx=1,ncpw%ngw
          IF (paral%io_parent)&
               WRITE(500,*) c0(Gindx,istate)
       ENDDO

       ! Closing the file...
       IF (paral%io_parent)&
            CALL fileclose(500)
    ENDDO

    RETURN
  END SUBROUTINE save_aldo_wannier
  ! ==================================================================
  SUBROUTINE save_norm(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...


    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(len=20)                        :: name
    CHARACTER(len=8)                         :: number
    INTEGER                                  :: Gindx, ia, ie, istate
    LOGICAL                                  :: ferror

    CALL stopgm('INTERACTION_P','SAVE-NORM not parallel.',& 
         __LINE__,__FILE__)

    ! main loop on the states...
    DO istate=1,nstate

       ! First creating the file name for the state istate
       IF (paral%io_parent)&
            WRITE(number,'(I8)')istate
       ! removing the spaces before and after the number (CPMD routine)
       CALL xstring(number,ia,ie)

       name='norm-'//number(ia:ie)
       IF (paral%io_parent)&
            WRITE(6,*) 'opening file: ',name

       ! opening the file to save to
       IF (paral%io_parent)&
            CALL fileopen(500,name,fo_def,ferror)
       IF (paral%io_parent)&
            WRITE(6,*) 'Saving data...'

       ! now writing the coefficients (c0) (g-space,ascii)
       DO Gindx=1,ncpw%ngw
          IF (paral%io_parent)&
               WRITE(500,*) REAL(c0(Gindx,istate)*CONJG(c0(gindx,istate)))
       ENDDO

       ! Closing the file...
       IF (paral%io_parent)&
            CALL fileclose(500)
    ENDDO

    RETURN
  END SUBROUTINE save_norm
  ! ==================================================================
  INTEGER FUNCTION get_new_NGW(my_ecut)
    ! ==--------------------------------------------------------------==
    ! This function returns the higest possible G vector index (<=NGW)
    ! for a given cutoff energy in Rydberg (<=ECUT)
    ! ==--------------------------------------------------------------==

    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)
    IMPLICIT NONE
    REAL(real_8) :: my_ecut

    INTEGER :: Gindx,New_NGW
    REAL(real_8) :: My_New_Gcut

    ! Checking if the cutoff energy is valid
    IF (my_ecut.GT.cntr%ecut) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'INVALID CUTOFF (>ECUT), returning NGW'
       get_new_NGW=ncpw%ngw
       RETURN
    ENDIF

    ! defining the new cutoff 
    My_New_Gcut = my_ecut/parm%tpiba2

    ! finding the corresponding G vector
    DO Gindx=1,ncpw%ngw
       IF (hg(Gindx).LE.My_New_Gcut) THEN
          New_NGW=Gindx
       ENDIF
    ENDDO

    ! checking if G=0 is selected + correcting
    IF (New_NGW.EQ.1) THEN
       New_NGW=0
    ENDIF

    ! Assigning the result + exiting..
    get_new_NGW=New_NGW

    RETURN
  END FUNCTION get_new_NGW
  ! ==================================================================
  SUBROUTINE calc_my_wannier_centre(c0,nstate,centre)
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)


    ! input
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: centre(4,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_my_wannier_centre'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: dd1, dd2, dd3, ddx
    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), sc0(:,:), scr(:), &
                                                xyzmat(:,:,:)
    INTEGER                                  :: i, ierr, ip1, ip2, lforces, &
                                                lscr, lscr2, nmcol
    INTEGER, ALLOCATABLE                     :: mapcol(:), mapful(:)

! variables
! Wannier centre calculation from Aldo (orth_wannier.F)

    ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! first the XYZMAT needs to be calculated
    CALL give_scr_ddipo(lscr2,tag)
    CALL give_scr_forces(lforces,tag,nstate,.FALSE.,.FALSE.)
    IF (paral%io_parent)&
         WRITE(6,*) 'Dimensions = ',lscr2,lforces
    lscr=MAX(lscr2,lforces)
    ALLOCATE(scr(lscr/2+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmcol=parai%nproc*ngwmax
    ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(spar%ngws,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mapcol(nmcol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setdip(mapful,mapcol)
    CALL set_operator(.TRUE.)
    ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ip1=1
    ip2=1+nstate*nstate!vw removed 2*
    !  lscr2=lscr-ip2

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO OPEIGR'
    CALL opeigr(c0,cs,sc0,nstate,mapful,mapcol,scr(ip1),&
         1,0,dd1)

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO COPY'
    CALL dcopy(2*nstate*nstate,scr(ip1),1,xyzmat(1,1,1),1)

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO OPEIGR'
    CALL opeigr(c0,cs,sc0,nstate,mapful,mapcol,scr(ip1),&
         2,0,dd2)

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO COPY'
    CALL dcopy(2*nstate*nstate,scr(ip1),1,xyzmat(1,1,2),1)

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO OPEIGR'
    CALL opeigr(c0,cs,sc0,nstate,mapful,mapcol,scr(ip1),&
         3,0,dd3)

    IF (paral%io_parent)&
         WRITE(6,*) 'ENTERING TO COPY'
    CALL dcopy(2*nstate*nstate,scr(ip1),1,xyzmat(1,1,3),1)

    DO i=4,wannc%nwanopt
       CALL opeigr(c0,cs,sc0,nstate,mapful,mapcol,scr(ip1),&
            wannc%iow(1,i-3),wannc%iow(2,i-3),ddx)
       CALL dcopy(2*nstate*nstate,scr(ip1),1,xyzmat(1,1,i),1)
    ENDDO

    ! then WC''s are calculated
    ! are you sure that TAU0 are the correct current coordinates?

    CALL wannier_center(xyzmat,nstate,nstate,centre,tau0)

    IF (paral%io_parent)&
         WRITE(6,*)'Wannier centers properly generated'

    RETURN
  END SUBROUTINE calc_my_wannier_centre
  ! ==================================================================

END MODULE interaction_manno_p_utils

SUBROUTINE print_orbital_resp_dens_manno(c0,c1,nstate,psi,scr)
  ! ==--------------------------------------------------------------==`
  ! Subroutine to save the perturbed DENSITY
  ! (2 W0W1 -> rho1), the starting
  ! orthogonalised wannier functions (W0 -> psi0)
  ! and its corresponding density (rho0 = W0*W0)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw,fpar
  USE parac, ONLY : paral,parai
  USE coor , ONLY:tau0
  USE fft_maxfft,                      ONLY: maxfftn, maxfft
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
  REAL(real_8)                               :: psi(maxfftn,2), &
                                                scr(fpar%nnr1,2)

  CHARACTER(len=2), DIMENSION(99), PARAMETER :: zahlenchars = (/'01','02','03'&
      ,'04','05','06','07','08','09','10','11','12','13','14','15','16','17',&
      '18','19','20','21','22','23','24','25','26','27','28','29','30','31',&
      '32','33','34','35','36','37','38','39','40','41','42','43','44','45',&
      '46','47','48','49','50','51','52','53','54','55','56','57','58','59',&
      '60','61','62','63','64','65','66','67','68','69','70','71','72','73',&
      '74','75','76','77','78','79','80','81','82','83','84','85','86','87',&
      '88','89','90','91','92','93','94','95','96','97','98','99' /)

  CHARACTER(len=80)                          :: filename
  INTEGER                                    :: ir, is, numstates
  REAL(real_8)                               :: sum_C0, sum_C0C1

  numstates = nstate
  IF (nstate .GT. 50) THEN
     IF (paral%io_parent)&
          WRITE (6,*)&
          ' PRINT RHO0: TOO MANY STATES. PRINTING ONLY 1..50.'
     numstates = 50
  ENDIF

  DO is=1,numstates
     ! saving the perturbed density (2W0W1)
     CALL fft2tor(c0(1,is),psi(1,1),&
          c1(1,is),psi(1,2),scr,ncpw%ngw,.FALSE.)

     sum_C0=0._real_8
     sum_C0C1=0._real_8

     DO ir=1,fpar%nnr1
        scr(ir,1)  =  psi(ir,1)*psi(ir,2)*2._real_8
        sum_C0=sum_c0+psi(ir,1)*psi(ir,1)
        IF (scr(ir,1).GT.0) sum_C0C1=sum_c0c1+scr(ir,1)
     ENDDO

     ! write(6,*) 'sum_C0 =',sum_C0
     ! write(6,*) 'sum_C0C1 =',sum_C0C1

     IF (paral%io_parent)&
          WRITE(6,fmt='(1X,A,I3,A,F8.5,A)') 'STATE',is,' contributes',&
          100*(sum_C0C1/REAL(NSTATE*sum_C0,kind=real_8)),&
          '% of charge displacement'

     CALL ffttog(scr(1,1),scr(1,1),psi,ncpw%nhg,.TRUE.)
     ! Is the norm of the perturbed density equal to zero?
     IF (ABS(REAL(scr(1,1),kind=real_8)).GT.1.e-13_real_8) THEN
        IF (paral%io_parent)&
             WRITE(6,*) 'WARNING: R-Norm of state',is,&
             '[n1(g=0)] =',scr(1,1),'! = 0'
     ENDIF
     filename = 'rho1-'//zahlenchars(is)//'.dens'
     CALL densto(scr(1,1),tau0,filename)

     ! saving the starting wannier orbitals (W0)
     CALL zeroing(psi)!,2*maxfft)
     CALL dcopy(2*ncpw%ngw,c0(1,is),1,psi,1)
     filename = 'psi0-'//zahlenchars(is)//'.wf'
     CALL densto(psi,tau0,filename)

     ! saving the starting wannier density (rho0)
     ! G to R space
     CALL ffttor(c0(1,is),psi(1,1),scr,ncpw%ngw,.FALSE.)

     ! calculate the square
     !$omp parallel do private(ir)
     DO ir=1,fpar%nnr1
        scr(ir,1)=psi(ir,1)*psi(ir,1)
     ENDDO

     ! R to G space
     CALL ffttog(scr(1,1),scr(1,1),psi,ncpw%nhg,.TRUE.)

     ! Save rho0
     filename = 'rho0-'//zahlenchars(is)//'.dens'
     CALL densto(scr(1,1),tau0,filename)
  ENDDO

  RETURN
END SUBROUTINE print_orbital_resp_dens_manno
! ==================================================================
